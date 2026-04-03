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

    // 2-3. Parse based on format
    // mmCIF-specific objects (only used when format == .mmcif)
    var doc: ?zreduce.cif.types.Document = null;
    defer if (doc) |*d| d.deinit();

    var inline_dict: ?zreduce.ccd.ComponentDict = null;
    defer if (inline_dict) |*d| d.deinit();

    var atom_lookup: ?zreduce.mmcif.AtomLookup = null;
    defer if (atom_lookup) |*al| al.deinit();

    // PDB-specific: records for passthrough output.
    // Ownership note: when format == .pdb, pdb_result.model is moved into `mdl`
    // (value copy). `mdl.deinit()` frees the model arrays. This defer only
    // frees the records list. Do NOT call pdb_result.?.deinit() as that would
    // double-free the model.
    var pdb_result: ?zreduce.pdb.PdbParseResult = null;
    defer if (pdb_result) |*r| {
        r.records.deinit(allocator);
    };

    var mdl: zreduce.model.Model = switch (config.format) {
        .mmcif => blk: {
            // 2. Parse CIF document (for preserving non-atom_site categories in output)
            doc = try zreduce.cif.readString(allocator, source);

            // 3. Extract model from CIF
            break :blk try zreduce.mmcif.parseModel(allocator, source);
        },
        .pdb => blk: {
            // 2-3. Parse PDB file
            pdb_result = try zreduce.pdb.parse(allocator, source);
            break :blk pdb_result.?.model;
        },
    };
    defer mdl.deinit();

    // 3a. Strip existing hydrogens if requested
    if (config.strip_h) {
        const n_stripped = mdl.stripHydrogens();
        if (!config.quiet and n_stripped > 0) {
            std.debug.print("  Stripped {d} existing H atoms\n", .{n_stripped});
        }
        // Also strip H rows from the CIF document's _atom_site loop so that
        // the writer's orig_loop row indices stay in sync with the model.
        if (doc) |*d| {
            stripDocumentHydrogens(&d.blocks.items[0]);
        }
    }

    // 3b-4a. mmCIF-specific: inline dict, atom lookup, struct_conn, branch links
    if (config.format == .mmcif) {
        const block = &doc.?.blocks.items[0];

        inline_dict = try zreduce.mmcif.parseInlineComponents(allocator, block);

        // Build atom lookup once for bond parsing
        atom_lookup = try zreduce.mmcif.buildAtomLookup(allocator, block);
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

    // 4. Apply chemistry annotations
    zreduce.place.applyChemistryWithConfig(&mdl, .{
        .protonation = if (protonation_overrides) |*ov| ov else null,
    });

    // 4a-4c. mmCIF-specific bond parsing (after applyChemistry which replaces atom.flags)
    if (config.format == .mmcif) {
        const block = &doc.?.blocks.items[0];

        try zreduce.mmcif.parseStructConn(&mdl, block, &atom_lookup.?);
        try zreduce.mmcif.parseBranchLinks(allocator, &mdl, block, &atom_lookup.?);
        zreduce.mmcif.flagLeavingAtoms(&mdl, if (inline_dict) |*d| d else null, config.dict);
    }

    // 5. Place hydrogens (per-component fallback: inline dict first, then external CCD)
    const place_result = try zreduce.place.addHydrogensWithConfig(
        &mdl,
        config.dict,
        if (inline_dict) |*d| d else null,
        .{
            .water = config.water,
            .bond_policy = config.bond_policy,
            .protonation = if (protonation_overrides) |*ov| ov else null,
        },
    );

    // Warn about unmatched protonation overrides
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

    // 6. Generate movers, apply fixes, and optimize unless --no-opt
    var movers: []zreduce.optimize.Mover = &.{};
    var movers_owned = false;
    defer {
        for (0..movers.len) |i| movers[i].deinit();
        if (movers_owned) allocator.free(movers);
    }

    const needs_movers = !config.no_opt or config.fix_path != null or config.dump_movers_path != null;
    if (needs_movers) {
        // Movers use per-component fallback: inline_dict first, then external CCD dict.
        const gen_result = try zreduce.optimize.generateMovers(
            allocator,
            &mdl,
            config.no_flip,
            config.dict,
            if (inline_dict) |*d| d else null,
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
            if (!config.quiet) std.debug.print("  Validation: {d} issue(s) found\n", .{validation.issues.len});
            if (config.validate_flag) {
                zreduce.validate.reportIssues(validation.issues, &mdl);
            }
        }
    }

    // 9. Write output (format-aware, plain or gzip-compressed based on extension)
    var out_buf: [4096]u8 = undefined;
    if (config.output_path) |out_path| {
        if (std.mem.endsWith(u8, out_path, ".gz")) {
            var gw = try zreduce.gzip.GzipWriter.init(allocator, out_path);
            errdefer gw.close() catch {};
            const aw = gw.anyWriter();
            try dispatchWriter(&aw, &mdl, config, &doc, &pdb_result);
            try gw.close();
        } else {
            const file = try std.fs.cwd().createFile(out_path, .{});
            defer file.close();
            var fw = file.writer(&out_buf);
            try dispatchWriter(&fw.interface, &mdl, config, &doc, &pdb_result);
            try fw.interface.flush();
        }
    } else {
        const stdout = std.fs.File.stdout();
        var sw = stdout.writer(&out_buf);
        try dispatchWriter(&sw.interface, &mdl, config, &doc, &pdb_result);
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
            config.bond_policy,
            movers,
            mdl.residues.items,
            mdl.chains.items,
        );
        try jw.interface.flush();
    }

    return result;
}

/// Dispatch to the correct writer based on input format.
fn dispatchWriter(
    writer: anytype,
    mdl: *const zreduce.model.Model,
    config: ProcessConfig,
    doc: *const ?zreduce.cif.types.Document,
    pdb_result: *const ?zreduce.pdb.PdbParseResult,
) !void {
    switch (config.format) {
        .mmcif => {
            try zreduce.writer.mmcif_writer.writeWithDocumentWithPolicy(writer, mdl, if (doc.*) |*d| d else null, config.bond_policy);
        },
        .pdb => {
            // pdb_result is always set when format == .pdb (set in the parse switch above)
            const pr = pdb_result.*.?;
            try zreduce.writer.pdb_writer.writeModel(writer, mdl, pr.records.items, config.bond_policy.output_isotope);
        },
    }
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
