//! Single-file processing pipeline: parse, place, optimize, validate, write.

const std = @import("std");
const Allocator = std.mem.Allocator;
const zreduce = @import("root.zig");

pub const ProcessConfig = struct {
    input_path: []const u8,
    output_path: ?[]const u8 = null, // must be non-null when called from batch workers (stdout is not thread-safe)
    dict: ?*const zreduce.ccd.ComponentDict = null, // shared read-only; caller owns lifetime
    json_path: ?[]const u8 = null,
    json_version: []const u8 = "", // version string for JSON log (passed from main)
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

pub fn readFile(allocator: Allocator, path: []const u8) ![]const u8 {
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

/// Process a single mmCIF file through the full pipeline.
/// Errors are returned (not exit'd) -- the caller decides how to handle them.
pub fn processFile(allocator: Allocator, config: ProcessConfig) !ProcessResult {
    // 1. Read input mmCIF
    const source = try readFile(allocator, config.input_path);
    defer allocator.free(source);

    // 2. Parse CIF document (for preserving non-atom_site categories in output)
    var doc = try zreduce.cif.readString(allocator, source);
    defer doc.deinit();

    // 3. Extract model from CIF
    var mdl = try zreduce.mmcif.parseModel(allocator, source);
    defer mdl.deinit();

    // 4. Apply chemistry annotations
    zreduce.place.applyChemistry(&mdl);

    // 5. Place hydrogens
    const place_result = try zreduce.place.addHydrogens(
        &mdl,
        config.dict,
    );

    var result = ProcessResult{
        .n_placed = place_result.n_placed,
        .n_residues = place_result.n_residues,
        .n_skipped = place_result.n_skipped,
    };

    // 6. Optimize (unless --no-opt)
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
            config.dict,
        );
        movers = gen_result.movers;
        movers_owned = true;
        result.n_movers = @intCast(movers.len);

        if (gen_result.n_skipped > 0) {
            std.debug.print("  Mover generation: {d} skipped (missing atoms or incomplete groups)\n", .{gen_result.n_skipped});
        }

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

    // 7. Mark absent H atoms
    markAbsentHydrogens(&mdl);

    // 8. Validate
    {
        var validation = try zreduce.validate.validateModel(allocator, &mdl);
        defer validation.deinit();

        if (!validation.ok()) {
            std.debug.print("  Validation: {d} issue(s) found\n", .{validation.issues.len});
            if (config.validate_flag) {
                zreduce.validate.reportIssues(validation.issues, &mdl);
            }
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

    // 10. Write JSON log (optional)
    if (config.json_path) |json_path| {
        var json_buf: [4096]u8 = undefined;
        const file = try std.fs.cwd().createFile(json_path, .{});
        defer file.close();
        var jw = file.writer(&json_buf);
        try zreduce.writer.json_writer.writeLog(
            &jw.interface,
            config.json_version,
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
