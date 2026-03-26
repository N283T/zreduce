const std = @import("std");
const build_options = @import("build_options");

pub fn main() !void {
    const args = try std.process.argsAlloc(std.heap.page_allocator);
    defer std.process.argsFree(std.heap.page_allocator, args);

    if (args.len < 2) {
        printUsage(args[0]);
        std.process.exit(1);
    }

    const subcmd = args[1];
    if (std.mem.eql(u8, subcmd, "--help") or std.mem.eql(u8, subcmd, "-h")) {
        printUsage(args[0]);
        return;
    }
    if (std.mem.eql(u8, subcmd, "--version") or std.mem.eql(u8, subcmd, "-V")) {
        std.debug.print("zreduce {s}\n", .{build_options.version});
        return;
    }

    std.debug.print("Unknown command: {s}\n", .{subcmd});
    std.process.exit(1);
}

fn printUsage(program_name: []const u8) void {
    std.debug.print(
        \\zreduce {s} - Hydrogen placement for mmCIF structures
        \\
        \\USAGE:
        \\    {s} [OPTIONS] <input.cif> [-o output.cif]
        \\
        \\OPTIONS:
        \\    -h, --help       Show this help message
        \\    -V, --version    Show version
        \\    -d, --dict PATH  Path to components.cif[.gz]
        \\    -o, --output PATH  Output file (default: stdout)
        \\    --no-opt         Skip optimization (placement only)
        \\    --no-flip        Disable Asn/Gln/His flips
        \\
    , .{ build_options.version, program_name });
}
