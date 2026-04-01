# Batch Processing Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add `batch` subcommand with file-level parallel processing and optional JSONL aggregated logging.

**Architecture:** Refactor main.zig into subcommand dispatch + extracted pipeline in run.zig. New batch.zig implements directory scanning, thread pool, JSONL writer, and summary. Optimizer gains configurable thread count to disable internal parallelism in batch mode.

**Tech Stack:** Zig 0.15, std.Thread, std.atomic.Value, std.Thread.Mutex

---

### Task 1: Add configurable thread count to OptConfig

**Files:**
- Modify: `src/optimize/optimizer.zig:22-26` (OptConfig struct)
- Modify: `src/optimize/optimizer.zig:79` (hardcoded getCpuCount)

- [ ] **Step 1: Add n_threads field to OptConfig**

In `src/optimize/optimizer.zig`, add `n_threads` to OptConfig:

```zig
pub const OptConfig = struct {
    n_threads: u32 = 0, // 0 = auto-detect CPU count
    brute_force_limit: u64 = 100_000,
    interaction_cutoff: f32 = 6.0,
    scoring_params: scorer_mod.ScoringParams = .{},
};
```

- [ ] **Step 2: Use config.n_threads in optimize()**

Replace the hardcoded line:

```zig
const n_threads: u32 = @intCast(@min(std.Thread.getCpuCount() catch 1, 8));
```

with:

```zig
const n_threads: u32 = if (config.n_threads > 0)
    config.n_threads
else
    @intCast(@min(std.Thread.getCpuCount() catch 1, 8));
```

- [ ] **Step 3: Run tests**

Run: `zig build test --summary all`
Expected: All 246 tests pass (no behavior change, default is 0 = auto).

- [ ] **Step 4: Commit**

```bash
git add src/optimize/optimizer.zig
git commit -m "feat: add configurable n_threads to OptConfig"
```

---

### Task 2: Extract run.zig from main.zig

**Files:**
- Create: `src/run.zig`
- Modify: `src/main.zig`
- Modify: `src/root.zig` (add run module export)

This task extracts the file processing pipeline (main.zig steps 2-10) into a reusable function. The main.zig `main()` function calls `run.processFile()` after parsing args.

- [ ] **Step 1: Create `src/run.zig` with ProcessConfig and ProcessResult**

```zig
//! Single-file processing pipeline: parse, place, optimize, validate, write.

const std = @import("std");
const Allocator = std.mem.Allocator;
const zreduce = @import("root.zig");
const build_options = @import("build_options");

pub const ProcessConfig = struct {
    input_path: []const u8,
    output_path: ?[]const u8 = null,
    dict: ?*zreduce.ccd.ComponentDict = null,
    json_path: ?[]const u8 = null,
    no_opt: bool = false,
    no_flip: bool = false,
    validate_flag: bool = false,
    opt_threads: u32 = 0, // 0 = auto; batch sets to 1
};

pub const ProcessResult = struct {
    n_placed: u32,
    n_residues: u32,
    n_skipped: u32,
    n_movers: u32 = 0,
    n_singletons: u32 = 0,
    n_brute_force: u32 = 0,
    n_vertex_cut: u32 = 0,
};
```

- [ ] **Step 2: Implement processFile()**

Extract the pipeline from main.zig (steps 2-10). Key differences from the current main.zig:
- Takes `allocator` parameter (batch uses per-thread arena)
- CCD dictionary is passed in via `config.dict` (not loaded here)
- Optimizer uses `config.opt_threads` for thread count
- Returns `ProcessResult` instead of printing directly
- Does NOT call `std.process.exit()` on errors — returns errors

```zig
fn readFile(allocator: Allocator, path: []const u8) ![]const u8 {
    const file = try std.fs.cwd().openFile(path, .{});
    defer file.close();
    return try file.readToEndAlloc(allocator, 1024 * 1024 * 1024);
}

fn markAbsentHydrogens(mdl: *zreduce.model.Model) void {
    for (mdl.atoms.items) |*atom| {
        if (zreduce.optimize.mover.isAbsentH(atom.*)) {
            atom.is_added = false;
        }
    }
}

fn countAddedHydrogens(mdl: *const zreduce.model.Model) u32 {
    var count: u32 = 0;
    for (mdl.atoms.items) |atom| {
        if (atom.is_added and atom.is_hydrogen) count += 1;
    }
    return count;
}

pub fn processFile(allocator: Allocator, config: ProcessConfig) !ProcessResult {
    // 1. Read input mmCIF
    const source = try readFile(allocator, config.input_path);
    defer allocator.free(source);

    // 2. Parse CIF document
    var doc = try zreduce.cif.readString(allocator, source);
    defer doc.deinit();

    // 3. Extract model
    var mdl = try zreduce.mmcif.parseModel(allocator, source);
    defer mdl.deinit();

    // 4. Apply chemistry
    zreduce.place.applyChemistry(&mdl);

    // 5. Place hydrogens
    const place_result = try zreduce.place.addHydrogens(
        &mdl,
        if (config.dict) |d| d else null,
    );

    var result = ProcessResult{
        .n_placed = place_result.n_placed,
        .n_residues = place_result.n_residues,
        .n_skipped = place_result.n_skipped,
    };

    // 6. Optimize
    var movers: []zreduce.optimize.Mover = &.{};
    var movers_owned = false;
    defer {
        for (0..movers.len) |i| movers[i].deinit();
        if (movers_owned) allocator.free(movers);
    }

    if (!config.no_opt) {
        const gen_result = try zreduce.optimize.generateMovers(
            allocator,
            &mdl,
            config.no_flip,
            if (config.dict) |d| d else null,
        );
        movers = gen_result.movers;
        movers_owned = true;
        result.n_movers = @intCast(movers.len);

        if (movers.len > 0) {
            const opt_config = zreduce.optimize.OptConfig{
                .n_threads = config.opt_threads,
            };
            const opt_result = try zreduce.optimize.optimizer.optimize(
                allocator,
                movers,
                &mdl,
                opt_config,
            );
            result.n_singletons = opt_result.n_singletons;
            result.n_brute_force = opt_result.n_brute_force;
            result.n_vertex_cut = opt_result.n_vertex_cut;
        }
    }

    // 7. Mark absent H
    markAbsentHydrogens(&mdl);

    // 8. Validate
    {
        var validation = try zreduce.validate.validateModel(allocator, &mdl);
        defer validation.deinit();
        if (!validation.ok() and config.validate_flag) {
            zreduce.validate.reportIssues(validation.issues, &mdl);
        }
    }

    // 9. Write output mmCIF
    var out_buf: [4096]u8 = undefined;
    if (config.output_path) |out_path| {
        const file = try std.fs.cwd().createFile(out_path, .{});
        defer file.close();
        var fw = file.writer(&out_buf);
        try zreduce.writer.mmcif_writer.writeWithDocument(&fw.interface, &mdl, &doc);
        try fw.interface.flush();
    } else {
        const stdout = std.fs.File.stdout();
        var sw = stdout.writer(&out_buf);
        try zreduce.writer.mmcif_writer.writeWithDocument(&sw.interface, &mdl, &doc);
        try sw.interface.flush();
    }

    // 10. Write JSON log
    if (config.json_path) |json_path| {
        var json_buf: [4096]u8 = undefined;
        const file = try std.fs.cwd().createFile(json_path, .{});
        defer file.close();
        var jw = file.writer(&json_buf);
        try zreduce.writer.json_writer.writeLog(
            &jw.interface,
            build_options.version,
            config.input_path,
            countAddedHydrogens(&mdl),
            movers,
            mdl.residues.items,
            mdl.chains.items,
        );
        try jw.interface.flush();
    }

    return result;
}
```

- [ ] **Step 3: Update root.zig to export run module**

Add to `src/root.zig`:

```zig
pub const run = @import("run.zig");
```

- [ ] **Step 4: Rewrite main.zig to use run.processFile()**

Replace the pipeline code in `main()` with a call to `run.processFile()`. Keep arg parsing and error printing with `std.process.exit()` in main.zig. The `Config` struct stays in main.zig and is mapped to `ProcessConfig`.

```zig
const std = @import("std");
const zreduce = @import("zreduce");
const build_options = @import("build_options");

const Allocator = std.mem.Allocator;

const Config = struct {
    input_path: []const u8,
    output_path: ?[]const u8 = null,
    dict_path: ?[]const u8 = null,
    json_path: ?[]const u8 = null,
    no_opt: bool = false,
    no_flip: bool = false,
    validate: bool = false,
};

fn readFile(allocator: Allocator, path: []const u8) ![]const u8 {
    const file = try std.fs.cwd().openFile(path, .{});
    defer file.close();
    return try file.readToEndAlloc(allocator, 1024 * 1024 * 1024);
}

fn parseRunArgs(args: []const []const u8) ?Config {
    var config = Config{ .input_path = undefined };
    var input_set = false;
    var i: usize = 0;

    while (i < args.len) : (i += 1) {
        const arg = args[i];

        if (std.mem.eql(u8, arg, "-h") or std.mem.eql(u8, arg, "--help")) {
            printRunUsage(args);
            return null;
        } else if (std.mem.eql(u8, arg, "-o") or std.mem.eql(u8, arg, "--output")) {
            i += 1;
            if (i >= args.len) {
                std.debug.print("Error: {s} requires a path argument\n", .{arg});
                std.process.exit(1);
            }
            config.output_path = args[i];
        } else if (std.mem.eql(u8, arg, "-d") or std.mem.eql(u8, arg, "--dict")) {
            i += 1;
            if (i >= args.len) {
                std.debug.print("Error: {s} requires a path argument\n", .{arg});
                std.process.exit(1);
            }
            config.dict_path = args[i];
        } else if (std.mem.eql(u8, arg, "--json")) {
            i += 1;
            if (i >= args.len) {
                std.debug.print("Error: --json requires a path argument\n", .{});
                std.process.exit(1);
            }
            config.json_path = args[i];
        } else if (std.mem.eql(u8, arg, "--no-opt")) {
            config.no_opt = true;
        } else if (std.mem.eql(u8, arg, "--no-flip")) {
            config.no_flip = true;
        } else if (std.mem.eql(u8, arg, "--validate")) {
            config.validate = true;
        } else if (arg.len > 0 and arg[0] == '-') {
            std.debug.print("Error: unknown option '{s}'\n", .{arg});
            std.process.exit(1);
        } else {
            if (input_set) {
                std.debug.print("Error: unexpected positional argument '{s}'\n", .{arg});
                std.process.exit(1);
            }
            config.input_path = arg;
            input_set = true;
        }
    }

    if (!input_set) {
        std.debug.print("Error: missing input file\n", .{});
        std.process.exit(1);
    }

    return config;
}

fn printUsage(program_name: []const u8) void {
    std.debug.print(
        \\zreduce {s} - Hydrogen placement for mmCIF structures
        \\
        \\USAGE:
        \\    {s} <command> [OPTIONS] <args>
        \\
        \\COMMANDS:
        \\    run      Process a single mmCIF file
        \\    batch    Process all mmCIF files in a directory
        \\
        \\GLOBAL OPTIONS:
        \\    -h, --help       Show this help message
        \\    -V, --version    Show version
        \\
        \\Use '{s} <command> --help' for more information.
        \\
    , .{ build_options.version, program_name, program_name });
}

fn printRunUsage(all_args: []const []const u8) void {
    _ = all_args;
    std.debug.print(
        \\USAGE:
        \\    zreduce run [OPTIONS] <input.cif>
        \\
        \\OPTIONS:
        \\    -h, --help         Show this help message
        \\    -d, --dict PATH    CCD dictionary
        \\    -o, --output PATH  Output file (default: stdout)
        \\    --json PATH        Write JSON log to file
        \\    --no-opt           Skip optimization
        \\    --no-flip          Disable Asn/Gln/His flips
        \\    --validate         Print validation diagnostics
        \\
    , .{});
}

pub fn main() !void {
    var gpa = std.heap.GeneralPurposeAllocator(.{}){};
    defer _ = gpa.deinit();
    const allocator = gpa.allocator();

    const args = std.process.argsAlloc(std.heap.page_allocator) catch {
        std.debug.print("Fatal: cannot read process arguments\n", .{});
        std.process.exit(2);
    };

    if (args.len < 2) {
        printUsage(args[0]);
        std.process.exit(1);
    }

    const subcmd = args[1];

    if (std.mem.eql(u8, subcmd, "-h") or std.mem.eql(u8, subcmd, "--help")) {
        printUsage(args[0]);
        return;
    }
    if (std.mem.eql(u8, subcmd, "-V") or std.mem.eql(u8, subcmd, "--version")) {
        std.debug.print("zreduce {s}\n", .{build_options.version});
        return;
    }

    if (std.mem.eql(u8, subcmd, "run")) {
        runSubcommand(allocator, args[2..]);
    } else if (std.mem.eql(u8, subcmd, "batch")) {
        // TODO: Task 4
        std.debug.print("Error: batch subcommand not yet implemented\n", .{});
        std.process.exit(1);
    } else {
        std.debug.print("Error: unknown subcommand '{s}'\n", .{subcmd});
        printUsage(args[0]);
        std.process.exit(1);
    }
}

fn runSubcommand(allocator: Allocator, args: []const []const u8) void {
    const config = parseRunArgs(args) orelse return;

    // Load CCD dictionary
    var ccd_dict: ?zreduce.ccd.ComponentDict = null;
    if (config.dict_path) |dict_path| {
        const dict_source = readFile(allocator, dict_path) catch |err| {
            std.debug.print("Error: cannot read dictionary '{s}': {s}\n", .{ dict_path, @errorName(err) });
            std.process.exit(1);
        };
        defer allocator.free(dict_source);
        ccd_dict = zreduce.ccd.parseComponentDict(allocator, dict_source) catch |err| {
            std.debug.print("Error: failed to parse CCD dictionary: {s}\n", .{@errorName(err)});
            std.process.exit(1);
        };
    }
    defer if (ccd_dict) |*d| d.deinit();

    const proc_config = zreduce.run.ProcessConfig{
        .input_path = config.input_path,
        .output_path = config.output_path,
        .dict = if (ccd_dict) |*d| d else null,
        .json_path = config.json_path,
        .no_opt = config.no_opt,
        .no_flip = config.no_flip,
        .validate_flag = config.validate,
    };

    const result = zreduce.run.processFile(allocator, proc_config) catch |err| {
        std.debug.print("Error: processing failed: {s}\n", .{@errorName(err)});
        std.process.exit(1);
    };

    std.debug.print("zreduce: placed {d} H atoms on {d} residues ({d} skipped)\n", .{
        result.n_placed,
        result.n_residues,
        result.n_skipped,
    });
}
```

- [ ] **Step 5: Run tests and verify `run` subcommand works**

Run: `zig build test --summary all`
Expected: All tests pass.

Manual check:
```bash
zig build -Doptimize=ReleaseFast
./zig-out/bin/zreduce run examples/data/AF-C1P619-F1-model_v4.cif -o /tmp/test_run.cif
./zig-out/bin/zreduce --help
./zig-out/bin/zreduce run --help
```

- [ ] **Step 6: Commit**

```bash
git add src/run.zig src/main.zig src/root.zig
git commit -m "refactor: extract run.zig pipeline, add subcommand dispatch to main.zig"
```

---

### Task 3: Implement batch.zig — directory scanning and sequential processing

**Files:**
- Create: `src/batch.zig`
- Modify: `src/root.zig` (add batch export)
- Modify: `src/main.zig` (wire batch subcommand)

This task implements the batch subcommand with sequential processing first. Parallel execution is added in Task 4.

- [ ] **Step 1: Create `src/batch.zig` with BatchConfig and arg parsing**

```zig
//! Batch processing: process all mmCIF files in a directory.

const std = @import("std");
const Allocator = std.mem.Allocator;
const zreduce = @import("root.zig");
const run_mod = @import("run.zig");
const build_options = @import("build_options");

pub const BatchConfig = struct {
    input_dir: []const u8,
    output_dir: ?[]const u8 = null,
    dict_path: ?[]const u8 = null,
    jsonl_path: ?[]const u8 = null,
    n_threads: u32 = 0, // 0 = auto-detect
    no_opt: bool = false,
    no_flip: bool = false,
    quiet: bool = false,
};

pub const FileResult = struct {
    filename: []const u8, // owned copy
    status: Status,
    result: ?run_mod.ProcessResult = null,
    error_msg: ?[]const u8 = null, // owned copy
    time_ns: u64 = 0,

    pub const Status = enum { ok, err };
};

pub const BatchResult = struct {
    total_files: usize,
    successful: usize,
    failed: usize,
    total_time_ns: u64,
    file_results: []FileResult,
    allocator: Allocator,

    pub fn deinit(self: *BatchResult) void {
        for (self.file_results) |r| {
            self.allocator.free(r.filename);
            if (r.error_msg) |msg| self.allocator.free(msg);
        }
        self.allocator.free(self.file_results);
    }

    pub fn printSummary(self: BatchResult) void {
        const total_ms = @as(f64, @floatFromInt(self.total_time_ns)) / 1_000_000.0;
        const throughput = if (total_ms > 0)
            @as(f64, @floatFromInt(self.successful)) / (total_ms / 1000.0)
        else
            0.0;
        std.debug.print("\nzreduce batch: {d}/{d} files processed", .{ self.successful, self.total_files });
        if (self.failed > 0) std.debug.print(" ({d} failed)", .{self.failed});
        std.debug.print("\n  Total time: {d:.1}s\n  Throughput: {d:.1} files/s\n", .{
            total_ms / 1000.0, throughput,
        });
    }
};
```

- [ ] **Step 2: Implement scanDirectory()**

Scan input directory for `*.cif` files, sorted alphabetically.

```zig
fn endsWithCif(name: []const u8) bool {
    if (name.len < 4) return false;
    const ext = name[name.len - 4 ..];
    return std.mem.eql(u8, ext, ".cif");
}

fn scanDirectory(allocator: Allocator, dir_path: []const u8) ![][]const u8 {
    var dir = try std.fs.cwd().openDir(dir_path, .{ .iterate = true });
    defer dir.close();

    var filenames = std.ArrayListUnmanaged([]const u8).empty;
    errdefer {
        for (filenames.items) |f| allocator.free(f);
        filenames.deinit(allocator);
    }

    var it = dir.iterate();
    while (try it.next()) |entry| {
        if (entry.kind != .file) continue;
        if (entry.name.len > 0 and entry.name[0] == '.') continue;
        if (!endsWithCif(entry.name)) continue;
        const copy = try allocator.dupe(u8, entry.name);
        try filenames.append(allocator, copy);
    }

    const items = try filenames.toOwnedSlice(allocator);
    std.mem.sort([]const u8, items, {}, struct {
        fn lessThan(_: void, a: []const u8, b: []const u8) bool {
            return std.mem.order(u8, a, b) == .lt;
        }
    }.lessThan);

    return items;
}
```

- [ ] **Step 3: Implement processFileInBatch() helper**

Wraps `run_mod.processFile()` with error capture:

```zig
fn processFileInBatch(
    allocator: Allocator,
    input_dir: []const u8,
    output_dir: []const u8,
    filename: []const u8,
    dict: ?*zreduce.ccd.ComponentDict,
    config: BatchConfig,
) FileResult {
    const input_path = std.fs.path.join(allocator, &.{ input_dir, filename }) catch
        return makeErrorResult(allocator, filename, "failed to construct input path");
    defer allocator.free(input_path);

    const output_path = std.fs.path.join(allocator, &.{ output_dir, filename }) catch
        return makeErrorResult(allocator, filename, "failed to construct output path");
    defer allocator.free(output_path);

    var timer = std.time.Timer.start() catch
        return makeErrorResult(allocator, filename, "failed to start timer");

    const proc_config = run_mod.ProcessConfig{
        .input_path = input_path,
        .output_path = output_path,
        .dict = dict,
        .no_opt = config.no_opt,
        .no_flip = config.no_flip,
        .opt_threads = 1, // disable internal parallelism
    };

    const proc_result = run_mod.processFile(allocator, proc_config) catch |err| {
        return makeErrorResult(allocator, filename, @errorName(err));
    };

    return FileResult{
        .filename = allocator.dupe(u8, filename) catch filename,
        .status = .ok,
        .result = proc_result,
        .time_ns = timer.read(),
    };
}

fn makeErrorResult(allocator: Allocator, filename: []const u8, msg: []const u8) FileResult {
    return FileResult{
        .filename = allocator.dupe(u8, filename) catch filename,
        .status = .err,
        .error_msg = allocator.dupe(u8, msg) catch null,
    };
}
```

- [ ] **Step 4: Implement runBatchSequential()**

```zig
fn runBatchSequential(
    allocator: Allocator,
    files: []const []const u8,
    input_dir: []const u8,
    output_dir: []const u8,
    dict: ?*zreduce.ccd.ComponentDict,
    config: BatchConfig,
) !BatchResult {
    var total_timer = try std.time.Timer.start();

    const file_results = try allocator.alloc(FileResult, files.len);

    for (files, 0..) |filename, i| {
        if (!config.quiet) {
            std.debug.print("\rProcessing: {d}/{d}", .{ i, files.len });
        }
        file_results[i] = processFileInBatch(
            allocator, input_dir, output_dir, filename, dict, config,
        );
    }
    if (!config.quiet and files.len > 0) {
        std.debug.print("\rProcessing: {d}/{d}\n", .{ files.len, files.len });
    }

    var successful: usize = 0;
    var failed: usize = 0;
    for (file_results) |r| {
        if (r.status == .ok) successful += 1 else failed += 1;
    }

    return BatchResult{
        .total_files = files.len,
        .successful = successful,
        .failed = failed,
        .total_time_ns = total_timer.read(),
        .file_results = file_results,
        .allocator = allocator,
    };
}
```

- [ ] **Step 5: Implement run() entry point**

```zig
pub fn run(allocator: Allocator, config: BatchConfig) !void {
    // Load CCD dictionary once
    var ccd_dict: ?zreduce.ccd.ComponentDict = null;
    if (config.dict_path) |dict_path| {
        const dict_source = readFile(allocator, dict_path) catch |err| {
            std.debug.print("Error: cannot read dictionary '{s}': {s}\n", .{ dict_path, @errorName(err) });
            return error.DictLoadFailed;
        };
        defer allocator.free(dict_source);
        ccd_dict = zreduce.ccd.parseComponentDict(allocator, dict_source) catch |err| {
            std.debug.print("Error: failed to parse CCD dictionary: {s}\n", .{@errorName(err)});
            return error.DictParseFailed;
        };
    }
    defer if (ccd_dict) |*d| d.deinit();

    // Scan directory
    const files = scanDirectory(allocator, config.input_dir) catch |err| {
        std.debug.print("Error: cannot scan directory '{s}': {s}\n", .{ config.input_dir, @errorName(err) });
        return error.ScanFailed;
    };
    defer {
        for (files) |f| allocator.free(f);
        allocator.free(files);
    }

    if (files.len == 0) {
        std.debug.print("No .cif files found in '{s}'\n", .{config.input_dir});
        return;
    }

    if (!config.quiet) {
        std.debug.print("zreduce batch: {d} files in '{s}'\n", .{ files.len, config.input_dir });
    }

    // Determine output directory
    const output_dir = config.output_dir orelse blk: {
        break :blk try std.fmt.allocPrint(allocator, "{s}_reduced", .{config.input_dir});
    };
    defer if (config.output_dir == null) allocator.free(output_dir);

    try std.fs.cwd().makePath(output_dir);

    // Process files (sequential for now; Task 4 adds parallel)
    var batch_result = try runBatchSequential(
        allocator, files, config.input_dir, output_dir,
        if (ccd_dict) |*d| d else null, config,
    );
    defer batch_result.deinit();

    // Write JSONL log
    if (config.jsonl_path) |jsonl_path| {
        try writeJsonlLog(allocator, batch_result.file_results, jsonl_path);
    }

    // Print summary
    batch_result.printSummary();

    if (batch_result.failed > 0) return error.SomeFilesFailed;
}

fn readFile(allocator: Allocator, path: []const u8) ![]const u8 {
    const file = try std.fs.cwd().openFile(path, .{});
    defer file.close();
    return try file.readToEndAlloc(allocator, 1024 * 1024 * 1024);
}
```

- [ ] **Step 6: Implement writeJsonlLog()**

Sequential version (writes after all files complete):

```zig
fn writeJsonlLog(
    allocator: Allocator,
    file_results: []const FileResult,
    path: []const u8,
) !void {
    const file = try std.fs.cwd().createFile(path, .{});
    defer file.close();
    var buf: [4096]u8 = undefined;
    var w = file.writer(&buf);

    for (file_results) |r| {
        try writeJsonlLine(&w.interface, r);
    }
    try w.interface.flush();
}

fn writeJsonlLine(writer: anytype, r: FileResult) !void {
    try writer.print("{{\"file\":\"{s}\",\"status\":\"{s}\"", .{
        r.filename,
        if (r.status == .ok) "ok" else "error",
    });
    if (r.status == .ok) {
        if (r.result) |res| {
            const time_ms = @as(f64, @floatFromInt(r.time_ns)) / 1_000_000.0;
            try writer.print(",\"hydrogens\":{d},\"movers\":{d},\"residues\":{d},\"time_ms\":{d:.1}", .{
                res.n_placed, res.n_movers, res.n_residues, time_ms,
            });
        }
    } else {
        if (r.error_msg) |msg| {
            try writer.print(",\"error\":\"{s}\"", .{msg});
        }
    }
    try writer.writeAll("}\n");
}
```

- [ ] **Step 7: Wire batch into main.zig and root.zig**

In `src/root.zig`, add:
```zig
pub const batch = @import("batch.zig");
```

In `src/main.zig`, replace the batch TODO:
```zig
} else if (std.mem.eql(u8, subcmd, "batch")) {
    batchSubcommand(allocator, args[2..]);
}
```

With `parseBatchArgs()` and `batchSubcommand()`:

```zig
fn parseBatchArgs(args: []const []const u8) ?zreduce.batch.BatchConfig {
    var config = zreduce.batch.BatchConfig{ .input_dir = undefined };
    var input_set = false;
    var i: usize = 0;

    while (i < args.len) : (i += 1) {
        const arg = args[i];
        if (std.mem.eql(u8, arg, "-h") or std.mem.eql(u8, arg, "--help")) {
            printBatchUsage();
            return null;
        } else if (std.mem.eql(u8, arg, "-o") or std.mem.eql(u8, arg, "--output")) {
            i += 1;
            if (i >= args.len) { std.debug.print("Error: {s} requires a path\n", .{arg}); std.process.exit(1); }
            config.output_dir = args[i];
        } else if (std.mem.eql(u8, arg, "-d") or std.mem.eql(u8, arg, "--dict")) {
            i += 1;
            if (i >= args.len) { std.debug.print("Error: {s} requires a path\n", .{arg}); std.process.exit(1); }
            config.dict_path = args[i];
        } else if (std.mem.eql(u8, arg, "--jsonl")) {
            i += 1;
            if (i >= args.len) { std.debug.print("Error: --jsonl requires a path\n", .{arg}); std.process.exit(1); }
            config.jsonl_path = args[i];
        } else if (std.mem.eql(u8, arg, "-j") or std.mem.eql(u8, arg, "--threads")) {
            i += 1;
            if (i >= args.len) { std.debug.print("Error: {s} requires a number\n", .{arg}); std.process.exit(1); }
            config.n_threads = std.fmt.parseInt(u32, args[i], 10) catch {
                std.debug.print("Error: invalid thread count '{s}'\n", .{args[i]});
                std.process.exit(1);
            };
        } else if (std.mem.eql(u8, arg, "--no-opt")) {
            config.no_opt = true;
        } else if (std.mem.eql(u8, arg, "--no-flip")) {
            config.no_flip = true;
        } else if (std.mem.eql(u8, arg, "--quiet")) {
            config.quiet = true;
        } else if (arg.len > 0 and arg[0] == '-') {
            std.debug.print("Error: unknown option '{s}'\n", .{arg});
            std.process.exit(1);
        } else {
            if (input_set) { std.debug.print("Error: unexpected argument '{s}'\n", .{arg}); std.process.exit(1); }
            config.input_dir = arg;
            input_set = true;
        }
    }

    if (!input_set) {
        std.debug.print("Error: missing input directory\n", .{});
        std.process.exit(1);
    }
    return config;
}

fn printBatchUsage() void {
    std.debug.print(
        \\USAGE:
        \\    zreduce batch [OPTIONS] <input_dir>
        \\
        \\OPTIONS:
        \\    -h, --help         Show this help message
        \\    -d, --dict PATH    CCD dictionary (loaded once)
        \\    -o, --output PATH  Output directory (default: <input>_reduced/)
        \\    -j, --threads N    Thread count (default: auto-detect)
        \\    --jsonl PATH       Aggregated JSONL log file
        \\    --no-opt           Skip optimization
        \\    --no-flip          Disable flips
        \\    --quiet            Suppress progress output
        \\
    , .{});
}

fn batchSubcommand(allocator: Allocator, args: []const []const u8) void {
    const config = parseBatchArgs(args) orelse return;
    zreduce.batch.run(allocator, config) catch |err| {
        if (err != error.SomeFilesFailed) {
            std.debug.print("Error: batch processing failed: {s}\n", .{@errorName(err)});
        }
        std.process.exit(1);
    };
}
```

- [ ] **Step 8: Test sequential batch**

Run: `zig build test --summary all`

Manual test:
```bash
zig build -Doptimize=ReleaseFast
mkdir -p /tmp/batch_test
cp examples/data/AF-C1P619-F1-model_v4.cif /tmp/batch_test/
cp examples/data/AF-P0A9J6-F1-model_v4.cif /tmp/batch_test/
./zig-out/bin/zreduce batch /tmp/batch_test -o /tmp/batch_out --jsonl /tmp/batch.jsonl
cat /tmp/batch.jsonl
```

- [ ] **Step 9: Commit**

```bash
git add src/batch.zig src/main.zig src/root.zig
git commit -m "feat: add batch subcommand with sequential file processing"
```

---

### Task 4: Add parallel execution to batch.zig

**Files:**
- Modify: `src/batch.zig`

- [ ] **Step 1: Add JsonlStreamWriter with mutex**

```zig
const JsonlStreamWriter = struct {
    mutex: std.Thread.Mutex = .{},
    file: std.fs.File,

    fn writeResult(self: *JsonlStreamWriter, allocator: Allocator, r: FileResult) void {
        // Serialize to buffer first (outside mutex)
        var line_buf = std.ArrayListUnmanaged(u8).empty;
        defer line_buf.deinit(allocator);
        writeJsonlLine(line_buf.writer(allocator), r) catch return;

        // Write under mutex
        self.mutex.lock();
        defer self.mutex.unlock();

        var buf: [4096]u8 = undefined;
        var w = self.file.writer(&buf);
        w.interface.writeAll(line_buf.items) catch return;
        w.interface.flush() catch {};
    }
};
```

- [ ] **Step 2: Add ParallelContext and parallelWorker**

```zig
const ParallelContext = struct {
    files: []const []const u8,
    input_dir: []const u8,
    output_dir: []const u8,
    config: BatchConfig,
    dict: ?*zreduce.ccd.ComponentDict,
    results: []FileResult,
    result_allocator: Allocator,
    next_file: std.atomic.Value(usize),
    processed_count: std.atomic.Value(usize),
    jsonl_stream: ?*JsonlStreamWriter,
};

fn parallelWorker(ctx: *ParallelContext) void {
    var arena = std.heap.ArenaAllocator.init(std.heap.smp_allocator);
    defer arena.deinit();

    while (true) {
        const file_idx = ctx.next_file.fetchAdd(1, .monotonic);
        if (file_idx >= ctx.files.len) break;

        const filename = ctx.files[file_idx];

        // Reset arena for each file
        _ = arena.reset(.retain_capacity);
        const alloc = arena.allocator();

        var result = processFileInBatch(
            alloc,
            ctx.input_dir,
            ctx.output_dir,
            filename,
            ctx.dict,
            ctx.config,
        );

        // Copy filename to result_allocator (thread-safe: each index is unique)
        result.filename = ctx.result_allocator.dupe(u8, filename) catch filename;
        if (result.error_msg) |msg| {
            result.error_msg = ctx.result_allocator.dupe(u8, msg) catch null;
        }

        ctx.results[file_idx] = result;

        // Stream JSONL if enabled
        if (ctx.jsonl_stream) |stream| {
            stream.writeResult(alloc, result);
        }

        _ = ctx.processed_count.fetchAdd(1, .monotonic);
    }
}
```

- [ ] **Step 3: Implement runBatchParallel()**

```zig
fn runBatchParallel(
    allocator: Allocator,
    files: []const []const u8,
    input_dir: []const u8,
    output_dir: []const u8,
    dict: ?*zreduce.ccd.ComponentDict,
    config: BatchConfig,
    jsonl_stream: ?*JsonlStreamWriter,
) !BatchResult {
    var total_timer = try std.time.Timer.start();

    const file_results = try allocator.alloc(FileResult, files.len);
    errdefer allocator.free(file_results);

    const cpu_count: u32 = @intCast(std.Thread.getCpuCount() catch 1);
    const n_threads = if (config.n_threads == 0)
        cpu_count
    else
        @min(config.n_threads, cpu_count);

    if (files.len == 1 or n_threads <= 1) {
        // Fall back to sequential
        allocator.free(file_results);
        var batch_result = try runBatchSequential(
            allocator, files, input_dir, output_dir, dict, config,
        );
        batch_result.total_time_ns = total_timer.read();
        return batch_result;
    }

    var ctx = ParallelContext{
        .files = files,
        .input_dir = input_dir,
        .output_dir = output_dir,
        .config = config,
        .dict = dict,
        .results = file_results,
        .result_allocator = allocator,
        .next_file = std.atomic.Value(usize).init(0),
        .processed_count = std.atomic.Value(usize).init(0),
        .jsonl_stream = jsonl_stream,
    };

    const actual_threads: u32 = @min(n_threads, @as(u32, @intCast(files.len)));
    const threads = try allocator.alloc(std.Thread, actual_threads);
    defer allocator.free(threads);

    for (threads) |*thread| {
        thread.* = try std.Thread.spawn(.{}, parallelWorker, .{&ctx});
    }

    // Progress monitor
    if (!config.quiet) {
        while (ctx.processed_count.load(.acquire) < files.len) {
            const processed = ctx.processed_count.load(.acquire);
            std.debug.print("\rProcessing: {d}/{d}", .{ processed, files.len });
            std.Thread.sleep(50 * std.time.ns_per_ms);
        }
        std.debug.print("\rProcessing: {d}/{d}\n", .{ files.len, files.len });
    }

    for (threads) |thread| {
        thread.join();
    }

    var successful: usize = 0;
    var failed: usize = 0;
    for (file_results) |r| {
        if (r.status == .ok) successful += 1 else failed += 1;
    }

    return BatchResult{
        .total_files = files.len,
        .successful = successful,
        .failed = failed,
        .total_time_ns = total_timer.read(),
        .file_results = file_results,
        .allocator = allocator,
    };
}
```

- [ ] **Step 4: Update run() to use parallel execution**

Replace the call to `runBatchSequential` in `run()` with:

```zig
    // Set up JSONL streaming writer
    var jsonl_file: ?std.fs.File = null;
    var jsonl_file_needs_close = false;
    if (config.jsonl_path) |path| {
        jsonl_file = try std.fs.cwd().createFile(path, .{});
        jsonl_file_needs_close = true;
    }
    defer if (jsonl_file_needs_close) {
        if (jsonl_file) |f| f.close();
    };

    var jsonl_stream_storage: JsonlStreamWriter = if (jsonl_file) |jf|
        JsonlStreamWriter{ .file = jf }
    else
        undefined;
    const jsonl_stream_ptr: ?*JsonlStreamWriter = if (jsonl_file != null) &jsonl_stream_storage else null;

    var batch_result = try runBatchParallel(
        allocator, files, config.input_dir, output_dir,
        if (ccd_dict) |*d| d else null, config, jsonl_stream_ptr,
    );
    defer batch_result.deinit();

    batch_result.printSummary();

    if (batch_result.failed > 0) return error.SomeFilesFailed;
```

And remove the separate `writeJsonlLog` call (JSONL is now streamed during processing).

- [ ] **Step 5: Test parallel batch**

Run: `zig build test --summary all`

Manual test with multiple files:
```bash
zig build -Doptimize=ReleaseFast
./zig-out/bin/zreduce batch examples/data -o /tmp/batch_par --jsonl /tmp/batch_par.jsonl
cat /tmp/batch_par.jsonl
./zig-out/bin/zreduce batch examples/data -o /tmp/batch_seq -j 1
```

- [ ] **Step 6: Commit**

```bash
git add src/batch.zig
git commit -m "feat: add parallel file processing to batch subcommand"
```

---

### Task 5: Tests and polish

**Files:**
- Modify: `src/batch.zig` (add unit tests)
- Modify: `src/run.zig` (add unit test)

- [ ] **Step 1: Add scanDirectory test**

```zig
test "scanDirectory finds .cif files" {
    // Use test_data directory which has known .cif files
    const files = try scanDirectory(std.testing.allocator, "src/test_data");
    defer {
        for (files) |f| std.testing.allocator.free(f);
        std.testing.allocator.free(files);
    }
    try std.testing.expect(files.len > 0);
    // Verify sorted
    for (files[0 .. files.len - 1], files[1..]) |a, b| {
        try std.testing.expect(std.mem.order(u8, a, b) == .lt);
    }
    // Verify all end with .cif
    for (files) |f| {
        try std.testing.expect(endsWithCif(f));
    }
}
```

- [ ] **Step 2: Add endsWithCif test**

```zig
test "endsWithCif" {
    try std.testing.expect(endsWithCif("test.cif"));
    try std.testing.expect(endsWithCif("path/to/file.cif"));
    try std.testing.expect(!endsWithCif("test.pdb"));
    try std.testing.expect(!endsWithCif("cif"));
    try std.testing.expect(!endsWithCif(".ci"));
    try std.testing.expect(!endsWithCif(""));
}
```

- [ ] **Step 3: Add writeJsonlLine test**

```zig
test "writeJsonlLine ok result" {
    var buf = std.ArrayListUnmanaged(u8).empty;
    defer buf.deinit(std.testing.allocator);
    const r = FileResult{
        .filename = "test.cif",
        .status = .ok,
        .result = .{ .n_placed = 100, .n_residues = 50, .n_skipped = 2, .n_movers = 30 },
        .time_ns = 1_500_000_000, // 1.5s
    };
    try writeJsonlLine(buf.writer(std.testing.allocator), r);
    const output = buf.items;
    try std.testing.expect(std.mem.indexOf(u8, output, "\"status\":\"ok\"") != null);
    try std.testing.expect(std.mem.indexOf(u8, output, "\"hydrogens\":100") != null);
}

test "writeJsonlLine error result" {
    var buf = std.ArrayListUnmanaged(u8).empty;
    defer buf.deinit(std.testing.allocator);
    const r = FileResult{
        .filename = "bad.cif",
        .status = .err,
        .error_msg = "InvalidSyntax",
    };
    try writeJsonlLine(buf.writer(std.testing.allocator), r);
    const output = buf.items;
    try std.testing.expect(std.mem.indexOf(u8, output, "\"status\":\"error\"") != null);
    try std.testing.expect(std.mem.indexOf(u8, output, "\"error\":\"InvalidSyntax\"") != null);
}
```

- [ ] **Step 4: Run full test suite**

Run: `zig build test --summary all`
Expected: All tests pass (existing 246 + new batch/run tests).

- [ ] **Step 5: End-to-end manual verification**

```bash
zig build -Doptimize=ReleaseFast

# Single file via run subcommand
./zig-out/bin/zreduce run examples/data/AF-C1P619-F1-model_v4.cif -o /tmp/e2e_run.cif
diff <(wc -l < /tmp/e2e_run.cif) <(echo "should have content")

# Batch with parallelism
./zig-out/bin/zreduce batch examples/data -o /tmp/e2e_batch --jsonl /tmp/e2e.jsonl
cat /tmp/e2e.jsonl | wc -l  # should match number of .cif files in examples/data

# Batch single-threaded
./zig-out/bin/zreduce batch examples/data -o /tmp/e2e_seq -j 1

# Help
./zig-out/bin/zreduce --help
./zig-out/bin/zreduce run --help
./zig-out/bin/zreduce batch --help
```

- [ ] **Step 6: Commit**

```bash
git add src/batch.zig src/run.zig
git commit -m "test: add unit tests for batch processing"
```
