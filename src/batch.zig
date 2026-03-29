//! Batch processing: scan a directory of mmCIF files, process in parallel,
//! and produce an aggregated JSONL log.

const std = @import("std");
const Allocator = std.mem.Allocator;
const zreduce = @import("root.zig");
const run_mod = zreduce.run;

// ---------------------------------------------------------------------------
// Configuration
// ---------------------------------------------------------------------------

pub const BatchConfig = struct {
    input_dir: []const u8,
    output_dir: ?[]const u8 = null, // default: <input_dir>_reduced/
    dict_path: ?[]const u8 = null,
    jsonl_path: ?[]const u8 = null,
    n_threads: u32 = 0, // 0 = auto-detect (for Task 4 parallelism)
    no_opt: bool = false,
    no_flip: bool = false,
    quiet: bool = false,
    json_version: []const u8 = "", // passed from main
};

// ---------------------------------------------------------------------------
// Result types
// ---------------------------------------------------------------------------

pub const FileResult = struct {
    filename: []const u8, // owned copy
    status: Status,
    result: ?run_mod.ProcessResult = null,
    error_msg: ?[]const u8 = null, // owned copy
    time_ns: u64 = 0,

    pub const Status = enum { ok, err };
};

pub const BatchResult = struct {
    total_files: u32,
    successful: u32,
    failed: u32,
    total_time_ns: u64,
    file_results: []FileResult,
    allocator: Allocator,

    pub fn printSummary(self: *const BatchResult) void {
        std.debug.print("\n--- Batch Summary ---\n", .{});
        std.debug.print("  Total files : {d}\n", .{self.total_files});
        std.debug.print("  Successful  : {d}\n", .{self.successful});
        std.debug.print("  Failed      : {d}\n", .{self.failed});
        const total_ms = @as(f64, @floatFromInt(self.total_time_ns)) / 1_000_000.0;
        std.debug.print("  Total time  : {d:.1} ms\n", .{total_ms});
    }

    pub fn deinit(self: *BatchResult) void {
        for (self.file_results) |r| {
            self.allocator.free(r.filename);
            if (r.error_msg) |msg| self.allocator.free(msg);
        }
        self.allocator.free(self.file_results);
    }
};

// ---------------------------------------------------------------------------
// Directory scanning
// ---------------------------------------------------------------------------

fn endsWithCif(name: []const u8) bool {
    return std.mem.endsWith(u8, name, ".cif");
}

pub fn scanDirectory(allocator: Allocator, dir_path: []const u8) ![][]const u8 {
    var dir = try std.fs.cwd().openDir(dir_path, .{ .iterate = true });
    defer dir.close();

    var names: std.ArrayListUnmanaged([]const u8) = .empty;
    errdefer {
        for (names.items) |n| allocator.free(n);
        names.deinit(allocator);
    }

    var iter = dir.iterate();
    while (try iter.next()) |entry| {
        if (entry.kind != .file) continue;
        if (entry.name.len == 0 or entry.name[0] == '.') continue;
        if (!endsWithCif(entry.name)) continue;
        const owned = try allocator.dupe(u8, entry.name);
        try names.append(allocator, owned);
    }

    const items = try names.toOwnedSlice(allocator);
    std.mem.sort([]const u8, items, {}, struct {
        fn lessThan(_: void, a: []const u8, b: []const u8) bool {
            return std.mem.order(u8, a, b) == .lt;
        }
    }.lessThan);
    return items;
}

// ---------------------------------------------------------------------------
// Per-file processing
// ---------------------------------------------------------------------------

fn makeErrorResult(allocator: Allocator, filename: []const u8, err_name: []const u8, elapsed_ns: u64) !FileResult {
    const owned_name = try allocator.dupe(u8, filename);
    errdefer allocator.free(owned_name);
    const owned_err = try allocator.dupe(u8, err_name);
    return FileResult{
        .filename = owned_name,
        .status = .err,
        .error_msg = owned_err,
        .time_ns = elapsed_ns,
    };
}

fn processFileInBatch(
    allocator: Allocator,
    filename: []const u8,
    input_dir: []const u8,
    output_dir: []const u8,
    dict: ?*const zreduce.ccd.ComponentDict,
    config: *const BatchConfig,
) !FileResult {
    const input_path = try std.fs.path.join(allocator, &.{ input_dir, filename });
    defer allocator.free(input_path);

    const output_path = try std.fs.path.join(allocator, &.{ output_dir, filename });
    defer allocator.free(output_path);

    const proc_config = run_mod.ProcessConfig{
        .input_path = input_path,
        .output_path = output_path,
        .dict = dict,
        .json_version = config.json_version,
        .no_opt = config.no_opt,
        .no_flip = config.no_flip,
        .opt_threads = 1,
    };

    var timer = try std.time.Timer.start();

    const proc_result = run_mod.processFile(allocator, proc_config) catch |err| {
        const elapsed = timer.read();
        return try makeErrorResult(allocator, filename, @errorName(err), elapsed);
    };

    const elapsed = timer.read();

    return FileResult{
        .filename = try allocator.dupe(u8, filename),
        .status = .ok,
        .result = proc_result,
        .time_ns = elapsed,
    };
}

// ---------------------------------------------------------------------------
// Sequential batch runner
// ---------------------------------------------------------------------------

fn runBatchSequential(
    allocator: Allocator,
    files: []const []const u8,
    input_dir: []const u8,
    output_dir: []const u8,
    dict: ?*const zreduce.ccd.ComponentDict,
    config: *const BatchConfig,
) !BatchResult {
    const results = try allocator.alloc(FileResult, files.len);
    var completed: usize = 0;
    errdefer {
        for (results[0..completed]) |r| {
            allocator.free(r.filename);
            if (r.error_msg) |msg| allocator.free(msg);
        }
        allocator.free(results);
    }

    var successful: u32 = 0;
    var failed: u32 = 0;
    var total_time_ns: u64 = 0;

    for (files, 0..) |filename, i| {
        if (!config.quiet) {
            std.debug.print("\rProcessing: {d}/{d}", .{ i + 1, files.len });
        }

        results[i] = try processFileInBatch(
            allocator,
            filename,
            input_dir,
            output_dir,
            dict,
            config,
        );
        completed = i + 1;

        total_time_ns += results[i].time_ns;
        switch (results[i].status) {
            .ok => successful += 1,
            .err => failed += 1,
        }
    }

    if (!config.quiet and files.len > 0) {
        std.debug.print("\n", .{});
    }

    return BatchResult{
        .total_files = @intCast(files.len),
        .successful = successful,
        .failed = failed,
        .total_time_ns = total_time_ns,
        .file_results = results,
        .allocator = allocator,
    };
}

// ---------------------------------------------------------------------------
// JSONL logging
// ---------------------------------------------------------------------------

const JsonlStreamWriter = struct {
    mutex: std.Thread.Mutex = .{},
    file: std.fs.File,

    fn writeResult(self: *JsonlStreamWriter, allocator: Allocator, file_result: FileResult) void {
        // Serialize to a buffer first (outside mutex)
        var line_buf = std.ArrayListUnmanaged(u8).empty;
        defer line_buf.deinit(allocator);
        writeJsonlLine(line_buf.writer(allocator), file_result) catch return;

        // Write under mutex (unbuffered to avoid stale buffered-writer state)
        self.mutex.lock();
        defer self.mutex.unlock();

        self.file.writeAll(line_buf.items) catch {};
    }
};

// ---------------------------------------------------------------------------
// Parallel execution
// ---------------------------------------------------------------------------

const ParallelContext = struct {
    files: []const []const u8,
    input_dir: []const u8,
    output_dir: []const u8,
    config: *const BatchConfig,
    dict: ?*const zreduce.ccd.ComponentDict,
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

        // Reset arena per file to bound memory usage
        _ = arena.reset(.retain_capacity);
        const alloc = arena.allocator();

        var result = processFileInBatch(
            alloc,
            filename,
            ctx.input_dir,
            ctx.output_dir,
            ctx.dict,
            ctx.config,
        ) catch |err| blk: {
            // processFileInBatch can only fail on OOM (path join or dupe)
            break :blk FileResult{
                .filename = filename, // will be duped below
                .status = .err,
                .error_msg = @errorName(err),
            };
        };

        // Copy owned strings to result_allocator (thread-safe: disjoint index)
        result.filename = ctx.result_allocator.dupe(u8, result.filename) catch filename;
        if (result.error_msg) |msg| {
            result.error_msg = ctx.result_allocator.dupe(u8, msg) catch null;
        }

        ctx.results[file_idx] = result;

        // Stream JSONL line if enabled
        if (ctx.jsonl_stream) |stream| {
            stream.writeResult(alloc, result);
        }

        _ = ctx.processed_count.fetchAdd(1, .monotonic);
    }
}

fn runBatchParallel(
    allocator: Allocator,
    files: []const []const u8,
    input_dir: []const u8,
    output_dir: []const u8,
    dict: ?*const zreduce.ccd.ComponentDict,
    config: *const BatchConfig,
    jsonl_stream: ?*JsonlStreamWriter,
) !BatchResult {
    const file_results = try allocator.alloc(FileResult, files.len);
    errdefer allocator.free(file_results);

    const cpu_count: u32 = @intCast(std.Thread.getCpuCount() catch 1);
    const n_threads = if (config.n_threads == 0)
        cpu_count
    else
        @min(config.n_threads, cpu_count);

    // Fall back to sequential for 1 thread or 1 file
    if (files.len <= 1 or n_threads <= 1) {
        allocator.free(file_results);
        return runBatchSequential(allocator, files, input_dir, output_dir, dict, config);
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

    // Progress monitor on main thread
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

    var successful: u32 = 0;
    var failed: u32 = 0;
    var total_time_ns: u64 = 0;
    for (file_results) |r| {
        switch (r.status) {
            .ok => successful += 1,
            .err => failed += 1,
        }
        total_time_ns += r.time_ns;
    }

    return BatchResult{
        .total_files = @intCast(files.len),
        .successful = successful,
        .failed = failed,
        .total_time_ns = total_time_ns,
        .file_results = file_results,
        .allocator = allocator,
    };
}

fn writeJsonString(writer: anytype, s: []const u8) !void {
    try writer.writeByte('"');
    for (s) |c| switch (c) {
        '"' => try writer.writeAll("\\\""),
        '\\' => try writer.writeAll("\\\\"),
        '\n' => try writer.writeAll("\\n"),
        '\r' => try writer.writeAll("\\r"),
        '\t' => try writer.writeAll("\\t"),
        else => try writer.writeByte(c),
    };
    try writer.writeByte('"');
}

fn writeJsonlLine(writer: anytype, file_result: FileResult) !void {
    try writer.writeAll("{\"file\":");
    try writeJsonString(writer, file_result.filename);
    try writer.writeAll(",\"status\":");
    switch (file_result.status) {
        .ok => {
            const r = file_result.result.?;
            const time_ms = @as(f64, @floatFromInt(file_result.time_ns)) / 1_000_000.0;
            try writer.print("\"ok\",\"hydrogens\":{d},\"movers\":{d},\"residues\":{d},\"time_ms\":{d:.1}}}\n", .{
                r.n_placed,
                r.n_movers,
                r.n_residues,
                time_ms,
            });
        },
        .err => {
            try writer.writeAll("\"error\",\"error\":");
            try writeJsonString(writer, file_result.error_msg orelse "unknown");
            try writer.writeAll("}\n");
        },
    }
}

pub fn writeJsonlLog(allocator: Allocator, file_results: []const FileResult, path: []const u8) !void {
    _ = allocator;
    const file = try std.fs.cwd().createFile(path, .{});
    defer file.close();
    var buf: [4096]u8 = undefined;
    var fw = file.writer(&buf);
    for (file_results) |r| {
        try writeJsonlLine(&fw.interface, r);
    }
    try fw.interface.flush();
}

// ---------------------------------------------------------------------------
// Entry point
// ---------------------------------------------------------------------------

pub fn run(allocator: Allocator, config: BatchConfig) !void {
    // 1. Load CCD dictionary (once)
    var ccd_dict: ?zreduce.ccd.ComponentDict = null;
    if (config.dict_path) |dict_path| {
        const dict_source = try run_mod.readFile(allocator, dict_path);
        defer allocator.free(dict_source);
        ccd_dict = try zreduce.ccd.parseComponentDict(allocator, dict_source);
    }
    defer if (ccd_dict) |*d| d.deinit();

    // 2. Scan directory
    const files = try scanDirectory(allocator, config.input_dir);
    defer {
        for (files) |f| allocator.free(f);
        allocator.free(files);
    }

    if (files.len == 0) {
        std.debug.print("No .cif files found in '{s}'\n", .{config.input_dir});
        return;
    }

    std.debug.print("Found {d} .cif files in '{s}'\n", .{ files.len, config.input_dir });

    // 3. Determine output directory
    const output_dir: []const u8 = if (config.output_dir) |od|
        od
    else blk: {
        const default = try std.fmt.allocPrint(allocator, "{s}_reduced", .{config.input_dir});
        break :blk default;
    };
    defer if (config.output_dir == null) allocator.free(output_dir);

    // 4. Create output directory
    std.fs.cwd().makePath(output_dir) catch |err| {
        std.debug.print("Error: cannot create output directory '{s}': {s}\n", .{ output_dir, @errorName(err) });
        return err;
    };

    // 5. Set up JSONL streaming writer (if requested)
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

    // 6. Process files (parallel or sequential)
    var batch_result = try runBatchParallel(
        allocator,
        files,
        config.input_dir,
        output_dir,
        if (ccd_dict) |*d| d else null,
        &config,
        jsonl_stream_ptr,
    );
    defer batch_result.deinit();

    if (config.jsonl_path) |jsonl_path| {
        std.debug.print("JSONL log written to '{s}'\n", .{jsonl_path});
    }

    // 7. Print summary
    batch_result.printSummary();

    // 8. Return error if any files failed
    if (batch_result.failed > 0) {
        return error.SomeFilesFailed;
    }
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

test "endsWithCif" {
    try std.testing.expect(endsWithCif("test.cif"));
    try std.testing.expect(endsWithCif("path/to/file.cif"));
    try std.testing.expect(!endsWithCif("test.pdb"));
    try std.testing.expect(!endsWithCif("cif"));
    try std.testing.expect(!endsWithCif(".ci"));
    try std.testing.expect(!endsWithCif(""));
}

test "writeJsonlLine ok result" {
    var buf = std.ArrayListUnmanaged(u8).empty;
    defer buf.deinit(std.testing.allocator);
    const r = FileResult{
        .filename = "test.cif",
        .status = .ok,
        .result = .{ .n_placed = 100, .n_residues = 50, .n_skipped = 2, .n_movers = 30 },
        .time_ns = 1_500_000_000,
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

test "scanDirectory finds cif files" {
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
    for (files) |f| {
        try std.testing.expect(endsWithCif(f));
    }
}
