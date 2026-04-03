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
    quiet: bool = false, // suppress diagnostic prints (batch mode)
    water: zreduce.place.WaterConfig = .{},
    bond_policy: zreduce.place.BondPolicy = .{},
    protonation_path: ?[]const u8 = null,
    fix_path: ?[]const u8 = null,
    dump_movers_path: ?[]const u8 = null,
    strip_h: bool = false,
    format: InputFormat = .mmcif,
    model_filter: ModelFilter = .all,
};

pub const ProcessResult = struct {
    n_placed: u32,
    n_residues: u32,
    n_skipped_existing: u32,
    n_skipped_inter_residue: u32,
    n_skipped_missing_ref: u32,
    n_skipped_quality_filter: u32 = 0,
    n_movers: u32 = 0,
    n_singletons: u32 = 0,
    n_brute_force: u32 = 0,
    n_vertex_cut: u32 = 0,

    pub fn totalSkipped(self: ProcessResult) u32 {
        return self.n_skipped_existing + self.n_skipped_inter_residue + self.n_skipped_missing_ref + self.n_skipped_quality_filter;
    }
};

pub const InputFormat = enum {
    mmcif,
    pdb,
};

pub const ModelFilter = zreduce.mmcif.ModelFilter;

/// Detect input format from file extension.
pub fn detectFormat(path: []const u8) InputFormat {
    if (std.mem.endsWith(u8, path, ".pdb") or
        std.mem.endsWith(u8, path, ".pdb.gz") or
        std.mem.endsWith(u8, path, ".ent") or
        std.mem.endsWith(u8, path, ".ent.gz"))
    {
        return .pdb;
    }
    return .mmcif;
}

pub fn readFile(allocator: Allocator, path: []const u8) ![]const u8 {
    if (std.mem.endsWith(u8, path, ".gz")) {
        return zreduce.gzip.readGzip(allocator, path);
    }
    const file = try std.fs.cwd().openFile(path, .{});
    defer file.close();
    return try file.readToEndAlloc(allocator, 1024 * 1024 * 1024);
}

/// Remove hydrogen rows from the CIF document's _atom_site loop.
/// This keeps the document in sync with the model after stripHydrogens().
fn stripDocumentHydrogens(block: *zreduce.cif.Block) void {
    const loop = block.findLoopMut("_atom_site.type_symbol") orelse return;
    const type_col = loop.findTag("_atom_site.type_symbol") orelse return;
    const w = loop.width();
    const n_rows = loop.length();

    var write: usize = 0;
    var read: usize = 0;
    while (read < n_rows) : (read += 1) {
        const type_sym = loop.val(read, type_col) orelse {
            // Missing type_symbol — not hydrogen, keep the row
            if (write != read) {
                for (0..w) |col| {
                    loop.values.items[write * w + col] = loop.values.items[read * w + col];
                }
            }
            write += 1;
            continue;
        };
        const is_h = std.ascii.eqlIgnoreCase(type_sym, "H") or
            std.ascii.eqlIgnoreCase(type_sym, "D");
        if (is_h) continue;

        if (write != read) {
            for (0..w) |col| {
                loop.values.items[write * w + col] = loop.values.items[read * w + col];
            }
        }
        write += 1;
    }
    loop.values.items.len = write * w;
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

/// Process a single structure file (mmCIF or PDB) through the full pipeline.
/// Errors are returned (not exit'd) -- the caller decides how to handle them.
pub fn processFile(allocator: Allocator, config: ProcessConfig) !ProcessResult {
    // 1. Read input
    const source = try readFile(allocator, config.input_path);
    defer allocator.free(source);

    return switch (config.format) {
        .mmcif => processFileMmcif(allocator, config, source),
        .pdb => processFilePdb(allocator, config, source),
    };
}

/// mmCIF pipeline: supports multi-model processing.
fn processFileMmcif(allocator: Allocator, config: ProcessConfig, source: []const u8) !ProcessResult {
    // 2. Parse CIF document (for preserving non-atom_site categories in output)
    var doc = try zreduce.cif.readString(allocator, source);
    defer doc.deinit();

    // 2a. Strip existing hydrogens from the document BEFORE parsing models.
    // This ensures cif_row_start/end in ModelEntry are consistent with the
    // (possibly stripped) document loop used by the writer and atom lookup.
    if (config.strip_h) {
        stripDocumentHydrogens(&doc.blocks.items[0]);
    }

    // 3. Extract models from CIF block (shares the same doc — no double parse).
    // When strip_h is active, models are parsed from the already-stripped loop,
    // so they contain no H atoms and row indices match the stripped document.
    var entries = try zreduce.mmcif.parseModelsFromBlock(allocator, &doc.blocks.items[0], config.model_filter);
    defer {
        for (entries.items) |*e| e.model.deinit();
        entries.deinit(allocator);
    }

    if (entries.items.len == 0) {
        if (!config.quiet) std.debug.print("  No models found matching filter\n", .{});
        return ProcessResult{
            .n_placed = 0,
            .n_residues = 0,
            .n_skipped_existing = 0,
            .n_skipped_inter_residue = 0,
            .n_skipped_missing_ref = 0,
        };
    }

    if (config.strip_h and !config.quiet) {
        // Count how many H atoms were in the original source (before strip)
        // by comparing original_atom_count with a hypothetical full parse.
        // Since we stripped the doc before parsing, the models already lack H.
        std.debug.print("  Stripped existing H atoms from document\n", .{});
    }

    // 3b. Parse inline components (once — model-independent)
    var inline_dict: ?zreduce.ccd.ComponentDict = null;
    defer if (inline_dict) |*d| d.deinit();
    const block = &doc.blocks.items[0];
    inline_dict = try zreduce.mmcif.parseInlineComponents(allocator, block);

    // Parse overrides (once — shared across models)
    var protonation_overrides: ?zreduce.place.ProtonationOverrides = null;
    defer if (protonation_overrides) |*ov| ov.deinit();
    if (config.protonation_path) |path| {
        protonation_overrides = zreduce.place.protonation.parseFile(allocator, path) catch |err| {
            std.debug.print("Error: failed to load protonation override file '{s}': {s}\n", .{ path, @errorName(err) });
            return err;
        };
    }

    var fix_overrides: ?zreduce.optimize.fix.FixOverrides = null;
    defer if (fix_overrides) |*ov| ov.deinit();
    if (config.fix_path) |path| {
        fix_overrides = zreduce.optimize.fix.parseFile(allocator, path) catch |err| {
            std.debug.print("Error: failed to load fix override file '{s}': {s}\n", .{ path, @errorName(err) });
            return err;
        };
    }

    // Process each model independently
    var result = ProcessResult{
        .n_placed = 0,
        .n_residues = 0,
        .n_skipped_existing = 0,
        .n_skipped_inter_residue = 0,
        .n_skipped_missing_ref = 0,
    };

    for (entries.items) |*entry| {
        const mdl = &entry.model;

        if (!config.quiet and entries.items.len > 1) {
            std.debug.print("  Processing model {d} ({d} atoms)\n", .{ entry.model_num, mdl.atoms.items.len });
        }

        // Build per-model atom lookup for bond parsing
        var atom_lookup = try zreduce.mmcif.buildAtomLookupForRange(allocator, block, entry.cif_row_start, entry.cif_row_end);
        defer atom_lookup.deinit();

        // 4. Apply chemistry annotations
        zreduce.place.applyChemistryWithConfig(mdl, .{
            .protonation = if (protonation_overrides) |*ov| ov else null,
        });

        // 4a-4c. mmCIF-specific bond parsing
        try zreduce.mmcif.parseStructConn(mdl, block, &atom_lookup);
        try zreduce.mmcif.parseBranchLinks(allocator, mdl, block, &atom_lookup);
        zreduce.mmcif.flagLeavingAtoms(mdl, if (inline_dict) |*d| d else null, config.dict);

        // 5. Place hydrogens
        const place_result = try zreduce.place.addHydrogensWithConfig(
            mdl,
            config.dict,
            if (inline_dict) |*d| d else null,
            .{
                .water = config.water,
                .bond_policy = config.bond_policy,
                .protonation = if (protonation_overrides) |*ov| ov else null,
            },
        );

        result.n_placed += place_result.n_placed;
        result.n_residues += place_result.n_residues;
        result.n_skipped_existing += place_result.n_skipped_existing;
        result.n_skipped_inter_residue += place_result.n_skipped_inter_residue;
        result.n_skipped_missing_ref += place_result.n_skipped_missing_ref;
        result.n_skipped_quality_filter += place_result.n_skipped_quality_filter;

        // 6. Generate movers and optimize
        const needs_movers = !config.no_opt or config.fix_path != null or config.dump_movers_path != null;
        if (needs_movers) {
            const gen_result = try zreduce.optimize.generateMovers(
                allocator,
                mdl,
                config.no_flip,
                config.dict,
                if (inline_dict) |*d| d else null,
                if (protonation_overrides) |*ov| ov else null,
                config.bond_policy.mode,
            );
            var movers = gen_result.movers;
            defer {
                for (0..movers.len) |i| movers[i].deinit();
                allocator.free(movers);
            }
            result.n_movers += @intCast(movers.len);

            if (fix_overrides) |*ov| {
                try zreduce.optimize.fix.applyFixes(ov, mdl, movers);
                for (movers) |*m| {
                    if (m.is_fixed) m.applyOrientation(mdl.atoms.items, m.best_orientation);
                }
            }

            if (!config.no_opt and movers.len > 0) {
                const opt_result = try zreduce.optimize.optimizer.optimize(
                    allocator,
                    movers,
                    mdl,
                    .{ .n_threads = config.opt_threads },
                );
                result.n_singletons += opt_result.n_singletons;
                result.n_brute_force += opt_result.n_brute_force;
                result.n_vertex_cut += opt_result.n_vertex_cut;
            }
        }

        // 7. Mark absent H atoms
        markAbsentHydrogens(mdl);

        // 8. Validate
        {
            var validation = try zreduce.validate.validateModel(allocator, mdl);
            defer validation.deinit();

            if (!validation.ok()) {
                if (!config.quiet) std.debug.print("  Model {d}: {d} validation issue(s)\n", .{ entry.model_num, validation.issues.len });
                if (config.validate_flag) {
                    zreduce.validate.reportIssues(validation.issues, mdl);
                }
            }
        }
    }

    // Warn about unmatched overrides (after all models processed)
    if (protonation_overrides) |*ov| {
        if (!config.quiet) {
            // Use first model for name resolution (overrides are file-level, not per-model)
            ov.warnUnmatched(&entries.items[0].model);
        }
    }

    // 9. Write output (all models into one file)
    var out_buf: [4096]u8 = undefined;
    if (config.output_path) |out_path| {
        if (std.mem.endsWith(u8, out_path, ".gz")) {
            var gw = try zreduce.gzip.GzipWriter.init(allocator, out_path);
            errdefer gw.close() catch {};
            const aw = gw.anyWriter();
            try zreduce.writer.mmcif_writer.writeMultiModelWithDocumentWithPolicy(&aw, entries.items, &doc, config.bond_policy);
            try gw.close();
        } else {
            const file = try std.fs.cwd().createFile(out_path, .{});
            defer file.close();
            var fw = file.writer(&out_buf);
            try zreduce.writer.mmcif_writer.writeMultiModelWithDocumentWithPolicy(&fw.interface, entries.items, &doc, config.bond_policy);
            try fw.interface.flush();
        }
    } else {
        const stdout = std.fs.File.stdout();
        var sw = stdout.writer(&out_buf);
        try zreduce.writer.mmcif_writer.writeMultiModelWithDocumentWithPolicy(&sw.interface, entries.items, &doc, config.bond_policy);
        try sw.interface.flush();
    }

    // 10. Write JSON log (optional, uses first model for summary)
    if (config.json_path) |json_path| {
        var json_buf: [4096]u8 = undefined;
        const file = try std.fs.cwd().createFile(json_path, .{});
        defer file.close();
        var jw = file.writer(&json_buf);
        var total_added: u32 = 0;
        for (entries.items) |entry| total_added += countAddedHydrogens(&entry.model);
        try zreduce.writer.json_writer.writeLog(
            &jw.interface,
            config.json_version,
            config.input_path,
            total_added,
            config.bond_policy,
            &.{}, // movers not available here (freed per model)
            entries.items[0].model.residues.items,
            entries.items[0].model.chains.items,
        );
        try jw.interface.flush();
    }

    return result;
}

/// PDB pipeline: single-model processing (multi-model PDB support is a future task).
fn processFilePdb(allocator: Allocator, config: ProcessConfig, source: []const u8) !ProcessResult {
    var pdb_result = try zreduce.pdb.parse(allocator, source);
    var mdl = pdb_result.model;
    defer mdl.deinit();
    defer pdb_result.records.deinit(allocator);

    if (config.strip_h) {
        const n_stripped = mdl.stripHydrogens();
        if (!config.quiet and n_stripped > 0) {
            std.debug.print("  Stripped {d} existing H atoms\n", .{n_stripped});
        }
    }

    var protonation_overrides: ?zreduce.place.ProtonationOverrides = null;
    defer if (protonation_overrides) |*ov| ov.deinit();
    if (config.protonation_path) |path| {
        protonation_overrides = zreduce.place.protonation.parseFile(allocator, path) catch |err| {
            std.debug.print("Error: failed to load protonation override file '{s}': {s}\n", .{ path, @errorName(err) });
            return err;
        };
    }

    var fix_overrides: ?zreduce.optimize.fix.FixOverrides = null;
    defer if (fix_overrides) |*ov| ov.deinit();
    if (config.fix_path) |path| {
        fix_overrides = zreduce.optimize.fix.parseFile(allocator, path) catch |err| {
            std.debug.print("Error: failed to load fix override file '{s}': {s}\n", .{ path, @errorName(err) });
            return err;
        };
    }

    zreduce.place.applyChemistryWithConfig(&mdl, .{
        .protonation = if (protonation_overrides) |*ov| ov else null,
    });

    const place_result = try zreduce.place.addHydrogensWithConfig(
        &mdl,
        config.dict,
        null, // no inline dict for PDB
        .{
            .water = config.water,
            .bond_policy = config.bond_policy,
            .protonation = if (protonation_overrides) |*ov| ov else null,
        },
    );

    if (protonation_overrides) |*ov| {
        if (!config.quiet) ov.warnUnmatched(&mdl);
    }

    var result = ProcessResult{
        .n_placed = place_result.n_placed,
        .n_residues = place_result.n_residues,
        .n_skipped_existing = place_result.n_skipped_existing,
        .n_skipped_inter_residue = place_result.n_skipped_inter_residue,
        .n_skipped_missing_ref = place_result.n_skipped_missing_ref,
        .n_skipped_quality_filter = place_result.n_skipped_quality_filter,
    };

    var movers: []zreduce.optimize.Mover = &.{};
    var movers_owned = false;
    defer {
        for (0..movers.len) |i| movers[i].deinit();
        if (movers_owned) allocator.free(movers);
    }

    const needs_movers = !config.no_opt or config.fix_path != null or config.dump_movers_path != null;
    if (needs_movers) {
        const gen_result = try zreduce.optimize.generateMovers(
            allocator,
            &mdl,
            config.no_flip,
            config.dict,
            null,
            if (protonation_overrides) |*ov| ov else null,
            config.bond_policy.mode,
        );
        movers = gen_result.movers;
        movers_owned = true;
        result.n_movers = @intCast(movers.len);

        if (!config.quiet and gen_result.n_skipped > 0) {
            std.debug.print("  Mover generation: {d} skipped (missing atoms or incomplete groups)\n", .{gen_result.n_skipped});
        }

        if (fix_overrides) |*ov| {
            try zreduce.optimize.fix.applyFixes(ov, &mdl, movers);
            for (movers) |*m| {
                if (m.is_fixed) m.applyOrientation(mdl.atoms.items, m.best_orientation);
            }
            if (!config.quiet) ov.warnUnmatched(&mdl, movers);
        }

        if (config.dump_movers_path) |dump_path| {
            var dump_buf: [4096]u8 = undefined;
            const file = try std.fs.cwd().createFile(dump_path, .{});
            defer file.close();
            var fw = file.writer(&dump_buf);
            try zreduce.optimize.fix.dumpMovers(&fw.interface, &mdl, movers);
            try fw.interface.flush();
        }

        if (!config.no_opt and movers.len > 0) {
            const opt_result = try zreduce.optimize.optimizer.optimize(
                allocator,
                movers,
                &mdl,
                .{ .n_threads = config.opt_threads },
            );
            result.n_singletons = opt_result.n_singletons;
            result.n_brute_force = opt_result.n_brute_force;
            result.n_vertex_cut = opt_result.n_vertex_cut;
        }
    }

    markAbsentHydrogens(&mdl);

    {
        var validation = try zreduce.validate.validateModel(allocator, &mdl);
        defer validation.deinit();

        if (!validation.ok()) {
            if (!config.quiet) std.debug.print("  Validation: {d} issue(s) found\n", .{validation.issues.len});
            if (config.validate_flag) {
                zreduce.validate.reportIssues(validation.issues, &mdl);
            }
        }
    }

    var out_buf: [4096]u8 = undefined;
    if (config.output_path) |out_path| {
        if (std.mem.endsWith(u8, out_path, ".gz")) {
            var gw = try zreduce.gzip.GzipWriter.init(allocator, out_path);
            errdefer gw.close() catch {};
            const aw = gw.anyWriter();
            try zreduce.writer.pdb_writer.writeModel(&aw, &mdl, pdb_result.records.items, config.bond_policy.output_isotope);
            try gw.close();
        } else {
            const file = try std.fs.cwd().createFile(out_path, .{});
            defer file.close();
            var fw = file.writer(&out_buf);
            try zreduce.writer.pdb_writer.writeModel(&fw.interface, &mdl, pdb_result.records.items, config.bond_policy.output_isotope);
            try fw.interface.flush();
        }
    } else {
        const stdout = std.fs.File.stdout();
        var sw = stdout.writer(&out_buf);
        try zreduce.writer.pdb_writer.writeModel(&sw.interface, &mdl, pdb_result.records.items, config.bond_policy.output_isotope);
        try sw.interface.flush();
    }

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
            config.bond_policy,
            movers,
            mdl.residues.items,
            mdl.chains.items,
        );
        try jw.interface.flush();
    }

    return result;
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

test "detectFormat identifies PDB extensions" {
    try std.testing.expectEqual(InputFormat.pdb, detectFormat("input.pdb"));
    try std.testing.expectEqual(InputFormat.pdb, detectFormat("input.pdb.gz"));
    try std.testing.expectEqual(InputFormat.pdb, detectFormat("input.ent"));
    try std.testing.expectEqual(InputFormat.pdb, detectFormat("input.ent.gz"));
    try std.testing.expectEqual(InputFormat.mmcif, detectFormat("input.cif"));
    try std.testing.expectEqual(InputFormat.mmcif, detectFormat("input.cif.gz"));
    try std.testing.expectEqual(InputFormat.mmcif, detectFormat("input.mmcif"));
}
