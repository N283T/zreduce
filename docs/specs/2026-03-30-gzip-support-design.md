# Gzip Input Support Design

## Context

Many PDB/AlphaFold distributions provide structures as `.cif.gz`. The CCD dictionary is also distributed as `components.cif.gz` (~44MB compressed, ~1GB decompressed). Currently zreduce requires pre-decompressed files.

Zig 0.15's native `std.compress.flate.Decompress` has a known panic bug (ziglang/zig#25035) on certain valid gzip files. The workaround is C zlib via `@cImport`.

## Design

### New file: `src/gzip.zig`

Ported from zsasa `src/gzip.zig` (proven in production with 250K+ file batch processing).

- `readGzip(allocator, path) ![]u8` — decompress a gzip file, caller owns returned slice
- `readGzipLimited(allocator, path, max_size) ![]u8` — with size guard
- Uses `gzopen`/`gzread`/`gzclose` via `@cImport(@cInclude("zlib.h"))`
- CRC validation on `gzclose` (detects corrupt data)
- 64KB chunk reading
- Default max size: 4GB

Error set: `GzipError = error{ GzipOpenFailed, GzipReadFailed, FileTooLarge, OutOfMemory }`

### Modified: `src/run.zig`

`readFile` gains `.gz` extension detection:

```zig
pub fn readFile(allocator: Allocator, path: []const u8) ![]const u8 {
    if (std.mem.endsWith(u8, path, ".gz")) {
        return gzip.readGzip(allocator, path);
    }
    const file = try std.fs.cwd().openFile(path, .{});
    defer file.close();
    return try file.readToEndAlloc(allocator, 1024 * 1024 * 1024);
}
```

This automatically handles all callers:
- `run.processFile` — input CIF
- `main.runSubcommand` — CCD dictionary
- `batch.run` — CCD dictionary

### Modified: `src/root.zig`

Add: `pub const gzip = @import("gzip.zig");`

### Modified: `build.zig`

The library module (`mod`) needs libc + libz for `@cImport` to resolve:

```zig
const mod = b.addModule("zreduce", .{
    .root_source_file = b.path("src/root.zig"),
    .target = target,
});
mod.link_libc = true;
mod.linkSystemLibrary("z", .{});
```

The exe module already links both, but since gzip.zig is part of the library module (imported via root.zig), the library must also link them. The exe inherits these via module dependency.

### Tests (in gzip.zig)

1. **readGzip decompresses gzip store block** — minimal in-memory gzip containing "Hello world\n"
2. **readGzip returns GzipOpenFailed for nonexistent file**
3. **readGzipLimited returns FileTooLarge when limit exceeded**

### No CLI changes

No new flags. Detection is automatic via file extension. Both `run` and `batch` subcommands gain gzip support transparently.
