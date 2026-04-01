# Batch Processing Design

## Context

zreduce currently processes one mmCIF file per invocation. Users with directories of structures (e.g., AlphaFold proteome dumps, PDB mirrors) must script sequential invocations, losing opportunities for:
- Shared CCD dictionary loading (hundreds of MB, read once)
- File-level parallelism across CPU cores
- Aggregated logging/reporting

## Design Decisions

| Decision | Choice | Rationale |
|----------|--------|-----------|
| CLI structure | Subcommand (`run` / `batch`) | Clear separation, extensible |
| Parallelism | File-level only | Simple, avoids nested thread pools; optimizer internal threading disabled in batch |
| CCD dictionary | Loaded once, shared read-only | Dictionary is immutable after parse; avoids re-reading hundreds of MB per file |
| JSONL log | Opt-in (`--jsonl PATH`) | Not everyone needs it; omit = no log overhead |
| Output | Directory-based (`-o out_dir/`) | Each input file produces a corresponding output file |

## CLI Interface

### Subcommand: `run` (single file — current behavior)

```
zreduce run [OPTIONS] <input.cif>

OPTIONS:
    -o, --output PATH    Output file (default: stdout)
    -d, --dict PATH      CCD dictionary
    --json PATH          JSON log file
    --no-opt             Skip optimization
    --no-flip            Disable Asn/Gln/His flips
    --validate           Print validation diagnostics
```

### Subcommand: `batch` (directory)

```
zreduce batch [OPTIONS] <input_dir>

OPTIONS:
    -o, --output PATH    Output directory (default: <input_dir>_reduced/)
    -d, --dict PATH      CCD dictionary (loaded once, shared across all files)
    -j, --threads N      Thread count (default: 0 = auto-detect CPU count)
    --jsonl PATH         Aggregated JSONL log file
    --no-opt             Skip optimization
    --no-flip            Disable flips
    --quiet              Suppress progress output
```

### Global Flags

```
zreduce --help           Show top-level help with subcommand list
zreduce --version        Show version
zreduce <cmd> --help     Show subcommand-specific help
```

## Architecture

### File Processing Pipeline (per file)

Each file is processed independently through the existing pipeline:
1. Read mmCIF source
2. Parse CIF document + model
3. Apply chemistry annotations
4. Place hydrogens (with shared CCD dict)
5. Optimize (with `n_threads=1` to disable internal parallelism)
6. Mark absent H, validate
7. Write output mmCIF
8. (Optional) Write JSONL log line

### Parallel Execution Model

```
Main thread:
  1. Parse args
  2. Load CCD dictionary (once)
  3. Scan input directory for *.cif files
  4. Allocate FileResult[] array
  5. Spawn N worker threads
  6. Progress monitor loop (atomic counter)
  7. Join all threads
  8. Print summary

Worker thread (each):
  1. Arena allocator (backed by smp_allocator)
  2. Loop: atomic fetch-add next_file index
  3. Process file through full pipeline
  4. Write result to FileResult[idx] (disjoint index)
  5. (Optional) Mutex-lock → write JSONL line → unlock
  6. Atomic increment processed_count
```

### Thread Safety

| Resource | Strategy |
|----------|----------|
| CCD dictionary | Read-only after parse, safe to share |
| Input files | Each thread reads its own file |
| Output files | Each thread writes to a unique output path |
| FileResult array | Disjoint index access (no overlap) |
| JSONL output | Mutex-protected streaming writer |
| Progress counter | `std.atomic.Value(usize)` |
| Memory | Per-thread `ArenaAllocator` over `smp_allocator` |

### JSONL Log Format

One JSON object per line, written as each file completes (order is non-deterministic):

```json
{"file":"input.cif","status":"ok","hydrogens":1234,"movers":567,"residues":890,"time_ms":45.2}
{"file":"bad.cif","status":"error","error":"failed to parse mmCIF: InvalidSyntax"}
```

Fields:
- `file`: input filename (basename, not full path)
- `status`: `"ok"` or `"error"`
- `hydrogens`: number of H atoms placed (ok only)
- `movers`: number of movers generated (ok only)
- `residues`: number of residues (ok only)
- `time_ms`: processing time in milliseconds (ok only)
- `error`: error message string (error only)

### Batch Summary (stderr)

```
zreduce batch: 100/100 files processed (2 failed)
  Total time: 12.3s
  Throughput: 8.1 files/s
```

## File Structure

```
src/
  main.zig          Subcommand dispatch (run/batch) + backward compat
  batch.zig         NEW: batch processing (scan, parallel exec, JSONL, summary)
  run.zig           NEW: extracted single-file pipeline from current main.zig
```

### Refactoring main.zig

Current `main.zig` (~280 lines) contains both CLI parsing and the processing pipeline. Refactor:

1. Extract the processing pipeline (steps 2-10) into `run.zig` as a reusable `processFile()` function
2. `main.zig` becomes a thin dispatcher: parse top-level args → call `run.run()` or `batch.run()`
3. `batch.zig` calls `run.processFile()` per file in its worker threads

### `run.processFile()` signature

```zig
pub const ProcessResult = struct {
    n_placed: u32,
    n_residues: u32,
    n_skipped: u32,
    n_movers: u32,
    n_singletons: u32,
    n_brute_force: u32,
    n_vertex_cut: u32,
};

pub const ProcessConfig = struct {
    input_path: []const u8,
    output_path: ?[]const u8 = null,
    dict: ?*ccd.ComponentDict = null,  // shared, read-only
    json_path: ?[]const u8 = null,
    no_opt: bool = false,
    no_flip: bool = false,
    validate_flag: bool = false,
    opt_threads: u32 = 0,  // 0 = auto; batch sets to 1
};

pub fn processFile(allocator: Allocator, config: ProcessConfig) !ProcessResult
```

## Error Handling

- Per-file errors are caught and recorded in `FileResult` (status=error, error_msg set)
- The batch continues processing remaining files
- Exit code: 0 if all succeed, 1 if any file fails
- Fatal errors (can't read input dir, can't create output dir) exit immediately

## Directory Scanning

- Scan `input_dir` for files matching `*.cif` (case-insensitive)
- Non-recursive (flat directory only)
- Sort filenames alphabetically for deterministic progress display
- Skip directories, symlinks to directories, hidden files (`.` prefix)

## Output Directory

- Default: `<input_dir>_reduced/` (append `_reduced` suffix)
- Created automatically if it doesn't exist
- Output filename: same as input filename (e.g., `input_dir/foo.cif` → `out_dir/foo.cif`)
- Overwrite existing files without warning (user controls output dir)
