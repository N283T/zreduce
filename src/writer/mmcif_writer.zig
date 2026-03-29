//! mmCIF writer: outputs model atoms plus a custom _zreduce_log category.

const std = @import("std");
const model_mod = @import("../model.zig");
const Model = model_mod.Model;
const Atom = model_mod.Atom;
const Residue = model_mod.Residue;
const Chain = model_mod.Chain;
const mover_mod = @import("../optimize/mover.zig");
const element = @import("../element.zig");
const AtomType = element.AtomType;
const cif = @import("../cif.zig");
const mmcif = @import("../mmcif.zig");
const place = @import("../place.zig");

/// Convert AtomType back to a 1-2 char element symbol string.
fn elementSymbol(atom_type: AtomType) []const u8 {
    return switch (atom_type) {
        .H, .Har, .Hpol, .Ha_p, .HOd => "H",
        .C, .Car, .C_eq_O => "C",
        .N, .Nacc => "N",
        .O => "O",
        .P => "P",
        .S => "S",
        .Se => "Se",
        .F => "F",
        .Cl => "Cl",
        .Br => "Br",
        .I => "I",
        .Li => "Li",
        .Na => "Na",
        .Mg => "Mg",
        .K => "K",
        .Ca => "Ca",
        .Mn => "Mn",
        .Fe => "Fe",
        .Co => "Co",
        .Ni => "Ni",
        .Cu => "Cu",
        .Zn => "Zn",
        .As => "As",
        .Rb => "Rb",
        .Sr => "Sr",
        .Mo => "Mo",
        .Ag => "Ag",
        .Cd => "Cd",
        .Sn => "Sn",
        .Cs => "Cs",
        .Ba => "Ba",
        .W => "W",
        .Pt => "Pt",
        .Au => "Au",
        .Hg => "Hg",
        .Pb => "Pb",
        .U => "U",
        .unknown => "X",
    };
}

const cif_types = cif.types;
const Document = cif_types.Document;
const Block = cif_types.Block;
const Loop = cif_types.Loop;
const Pair = cif_types.Pair;
const Item = cif_types.Item;

/// Write a model as mmCIF, preserving all original categories from the source document.
/// The _atom_site loop is replaced with the model's atoms (including added H).
/// If doc is null, writes atom_site only (atom-site-only mode).
pub fn writeWithDocument(writer: anytype, model: *const Model, doc: ?*const Document) !void {
    if (doc) |d| {
        if (d.blocks.items.len > 0) {
            const block = &d.blocks.items[0];
            try writer.print("data_{s}\n", .{block.name});

            // Write all items, replacing _atom_site loop
            var atom_site_written = false;
            for (block.items.items) |item| {
                switch (item) {
                    .pair => |p| {
                        try writePairValue(writer, p.tag, p.value);
                    },
                    .loop => |l| {
                        if (isAtomSiteLoop(&l)) {
                            try writer.writeAll("#\n");
                            try writeAtomSitePreserving(writer, model, &l);
                            try writer.writeAll("#\n");
                            atom_site_written = true;
                        } else {
                            try writeLoop(writer, &l);
                        }
                    },
                }
            }

            // If no atom_site loop was found in original, append it
            if (!atom_site_written) {
                try writer.writeAll("#\n");
                try writeAtomSite(writer, model);
                try writer.writeAll("#\n");
            }
            return;
        }
    }

    // Fallback: atom-site-only mode
    try writer.writeAll("data_ZREDUCE\n#\n");
    try writeAtomSite(writer, model);
    try writer.writeAll("#\n");
}

/// Write atom-site-only output (no original document preservation).
pub fn write(writer: anytype, model: *const Model, block_name: []const u8) !void {
    try writer.print("data_{s}\n#\n", .{block_name});
    try writeAtomSite(writer, model);
    try writer.writeAll("#\n");
}

fn isAtomSiteLoop(loop: *const Loop) bool {
    for (loop.tags.items) |tag| {
        if (std.ascii.startsWithIgnoreCase(tag, "_atom_site.")) return true;
    }
    return false;
}

/// Write _atom_site loop preserving the original column structure.
/// Heavy atom rows are written from the original loop data.
/// Added H atom rows fill known columns and use '.' for the rest.
fn writeAtomSitePreserving(writer: anytype, model: *const Model, orig_loop: *const Loop) !void {
    const w = orig_loop.width();
    if (w == 0) return;

    // Write loop header with original tags
    try writer.writeAll("loop_\n");
    for (orig_loop.tags.items) |tag| {
        try writer.print("{s}\n", .{tag});
    }

    // Map column indices for H atom generation
    const ColMap = struct {
        group_pdb: ?usize = null,
        id: ?usize = null,
        type_symbol: ?usize = null,
        label_atom_id: ?usize = null,
        label_alt_id: ?usize = null,
        label_comp_id: ?usize = null,
        label_asym_id: ?usize = null,
        label_entity_id: ?usize = null,
        label_seq_id: ?usize = null,
        cartn_x: ?usize = null,
        cartn_y: ?usize = null,
        cartn_z: ?usize = null,
        occupancy: ?usize = null,
        b_factor: ?usize = null,
        auth_seq_id: ?usize = null,
        auth_comp_id: ?usize = null,
        auth_asym_id: ?usize = null,
        auth_atom_id: ?usize = null,
        pdb_model_num: ?usize = null,
    };
    var cm = ColMap{};
    for (orig_loop.tags.items, 0..) |tag, i| {
        if (std.ascii.eqlIgnoreCase(tag, "_atom_site.group_PDB")) cm.group_pdb = i;
        if (std.ascii.eqlIgnoreCase(tag, "_atom_site.id")) cm.id = i;
        if (std.ascii.eqlIgnoreCase(tag, "_atom_site.type_symbol")) cm.type_symbol = i;
        if (std.ascii.eqlIgnoreCase(tag, "_atom_site.label_atom_id")) cm.label_atom_id = i;
        if (std.ascii.eqlIgnoreCase(tag, "_atom_site.label_alt_id")) cm.label_alt_id = i;
        if (std.ascii.eqlIgnoreCase(tag, "_atom_site.label_comp_id")) cm.label_comp_id = i;
        if (std.ascii.eqlIgnoreCase(tag, "_atom_site.label_asym_id")) cm.label_asym_id = i;
        if (std.ascii.eqlIgnoreCase(tag, "_atom_site.label_entity_id")) cm.label_entity_id = i;
        if (std.ascii.eqlIgnoreCase(tag, "_atom_site.label_seq_id")) cm.label_seq_id = i;
        if (std.ascii.eqlIgnoreCase(tag, "_atom_site.Cartn_x")) cm.cartn_x = i;
        if (std.ascii.eqlIgnoreCase(tag, "_atom_site.Cartn_y")) cm.cartn_y = i;
        if (std.ascii.eqlIgnoreCase(tag, "_atom_site.Cartn_z")) cm.cartn_z = i;
        if (std.ascii.eqlIgnoreCase(tag, "_atom_site.occupancy")) cm.occupancy = i;
        if (std.ascii.eqlIgnoreCase(tag, "_atom_site.B_iso_or_equiv")) cm.b_factor = i;
        if (std.ascii.eqlIgnoreCase(tag, "_atom_site.auth_seq_id")) cm.auth_seq_id = i;
        if (std.ascii.eqlIgnoreCase(tag, "_atom_site.auth_comp_id")) cm.auth_comp_id = i;
        if (std.ascii.eqlIgnoreCase(tag, "_atom_site.auth_asym_id")) cm.auth_asym_id = i;
        if (std.ascii.eqlIgnoreCase(tag, "_atom_site.auth_atom_id")) cm.auth_atom_id = i;
        if (std.ascii.eqlIgnoreCase(tag, "_atom_site.pdbx_PDB_model_num")) cm.pdb_model_num = i;
    }

    // Track serial for renumbering
    var serial: u32 = 1;

    // Write rows grouped by residue: original heavy atoms then added H
    for (model.residues.items, 0..) |res, res_idx| {
        const chain = model.chains.items[res.chain_idx];

        // Original heavy atoms: write rows from original loop data
        for (res.atom_start..res.atom_end) |orig_row_idx| {
            if (orig_row_idx < orig_loop.length()) {
                for (0..w) |col| {
                    if (col > 0) try writer.writeByte(' ');
                    if (cm.id != null and col == cm.id.?) {
                        try writer.print("{d}", .{serial});
                    } else {
                        const val = orig_loop.val(orig_row_idx, col) orelse ".";
                        try writeCifValue(writer, val);
                    }
                }
                try writer.writeByte('\n');
            }
            serial += 1;
        }

        // Added H atoms for this residue
        for (model.atoms.items) |atom| {
            if (!atom.is_added) continue;
            if (atom.residue_idx != res_idx) continue;
            // Skip absent H atoms (flipper sentinel)
            if (mover_mod.isAbsentH(atom)) continue;

            for (0..w) |col| {
                if (col > 0) try writer.writeByte(' ');
                if (cm.group_pdb != null and col == cm.group_pdb.?) {
                    try writer.writeAll("ATOM");
                } else if (cm.id != null and col == cm.id.?) {
                    try writer.print("{d}", .{serial});
                } else if (cm.type_symbol != null and col == cm.type_symbol.?) {
                    try writer.writeAll(elementSymbol(atom.element_type));
                } else if (cm.label_atom_id != null and col == cm.label_atom_id.?) {
                    try writeAtomName(writer, atom.nameSlice());
                } else if (cm.label_alt_id != null and col == cm.label_alt_id.?) {
                    try writeAltId(writer, atom.altloc);
                } else if (cm.label_comp_id != null and col == cm.label_comp_id.?) {
                    try writer.writeAll(res.compIdSlice());
                } else if (cm.label_asym_id != null and col == cm.label_asym_id.?) {
                    try writer.writeAll(chain.labelSlice());
                } else if (cm.label_entity_id != null and col == cm.label_entity_id.?) {
                    // Copy from first original row of this residue if available
                    const first_orig = res.atom_start;
                    if (first_orig < orig_loop.length()) {
                        const val = orig_loop.val(first_orig, col) orelse ".";
                        try writeCifValue(writer, val);
                    } else {
                        try writer.writeAll(".");
                    }
                } else if (cm.label_seq_id != null and col == cm.label_seq_id.?) {
                    try writer.print("{d}", .{res.seq_id});
                } else if (cm.cartn_x != null and col == cm.cartn_x.?) {
                    try writer.print("{d:.3}", .{atom.pos.x});
                } else if (cm.cartn_y != null and col == cm.cartn_y.?) {
                    try writer.print("{d:.3}", .{atom.pos.y});
                } else if (cm.cartn_z != null and col == cm.cartn_z.?) {
                    try writer.print("{d:.3}", .{atom.pos.z});
                } else if (cm.occupancy != null and col == cm.occupancy.?) {
                    try writer.print("{d:.2}", .{atom.occupancy});
                } else if (cm.b_factor != null and col == cm.b_factor.?) {
                    try writer.print("{d:.2}", .{atom.b_factor});
                } else if (cm.auth_seq_id != null and col == cm.auth_seq_id.?) {
                    // Copy from first original row of this residue
                    const first_orig = res.atom_start;
                    if (first_orig < orig_loop.length()) {
                        const val = orig_loop.val(first_orig, col) orelse ".";
                        try writeCifValue(writer, val);
                    } else {
                        try writer.print("{d}", .{res.seq_id});
                    }
                } else if (cm.auth_comp_id != null and col == cm.auth_comp_id.?) {
                    try writer.writeAll(res.compIdSlice());
                } else if (cm.auth_asym_id != null and col == cm.auth_asym_id.?) {
                    try writer.writeAll(chain.authSlice());
                } else if (cm.auth_atom_id != null and col == cm.auth_atom_id.?) {
                    try writeAtomName(writer, atom.nameSlice());
                } else if (cm.pdb_model_num != null and col == cm.pdb_model_num.?) {
                    try writer.writeAll("1");
                } else {
                    // Unknown column — copy from first atom of residue or use '.'
                    const first_orig = res.atom_start;
                    if (first_orig < orig_loop.length()) {
                        const val = orig_loop.val(first_orig, col) orelse ".";
                        try writeCifValue(writer, val);
                    } else {
                        try writer.writeAll(".");
                    }
                }
            }
            try writer.writeByte('\n');
            serial += 1;
        }
    }
}

/// Write a CIF loop as-is (non-atom_site categories).
fn writeLoop(writer: anytype, loop: *const Loop) !void {
    try writer.writeAll("#\nloop_\n");
    for (loop.tags.items) |tag| {
        try writer.print("{s}\n", .{tag});
    }
    const w = loop.width();
    if (w > 0) {
        for (0..loop.length()) |row| {
            for (0..w) |col| {
                if (col > 0) try writer.writeByte(' ');
                const val = loop.val(row, col) orelse ".";
                try writeCifValue(writer, val);
            }
            try writer.writeByte('\n');
        }
    }
    try writer.writeAll("#\n");
}

/// Write an atom name, trimming trailing spaces.
/// Atom names in CIF are typically unquoted even when they contain leading spaces.
fn writeAtomName(writer: anytype, name: []const u8) !void {
    // Trim trailing spaces
    var end = name.len;
    while (end > 0 and name[end - 1] == ' ') end -= 1;
    // Trim leading spaces
    var start: usize = 0;
    while (start < end and name[start] == ' ') start += 1;
    if (start >= end) {
        try writer.writeByte('.');
    } else {
        try writer.writeAll(name[start..end]);
    }
}

/// Write a CIF value, quoting if it contains spaces, quotes, or special characters.
/// Note: '.' and '?' are written unquoted as CIF null/unknown markers.
/// This is correct for round-tripping parsed CIF values where the parser
/// already stripped quotes from actual data values.
fn writeCifValue(writer: anytype, val: []const u8) !void {
    if (val.len == 0) {
        try writer.writeByte('.');
        return;
    }
    // Check if quoting is needed
    var needs_quote = false;
    var has_single = false;
    var has_double = false;
    var has_newline = false;
    for (val) |c| {
        if (c == ' ' or c == '\t') needs_quote = true;
        if (c == '\'') has_single = true;
        if (c == '"') has_double = true;
        if (c == '\n' or c == '\r') has_newline = true;
    }
    // Starts with special char?
    if (val[0] == '_' or val[0] == '#' or val[0] == '$' or val[0] == ';' or
        val[0] == '[' or val[0] == ']' or val[0] == '{' or val[0] == '}') needs_quote = true;
    // Could be confused with CIF keyword?
    if (std.ascii.startsWithIgnoreCase(val, "data_") or
        std.ascii.startsWithIgnoreCase(val, "save_") or
        std.ascii.eqlIgnoreCase(val, "loop_") or
        std.ascii.eqlIgnoreCase(val, "stop_") or
        std.ascii.eqlIgnoreCase(val, "global_"))
    {
        needs_quote = true;
    }
    // Bare '?' and '.' are CIF missing/inapplicable markers — quote if literal
    if (std.mem.eql(u8, val, "?") or std.mem.eql(u8, val, ".")) needs_quote = true;

    if (has_newline) needs_quote = true;

    if (!needs_quote) {
        try writer.writeAll(val);
    } else if (!has_single) {
        try writer.writeByte('\'');
        for (val) |c| {
            if (c == '\n' or c == '\r') {
                try writer.writeByte(' ');
            } else {
                try writer.writeByte(c);
            }
        }
        try writer.writeByte('\'');
    } else if (!has_double) {
        try writer.writeByte('"');
        for (val) |c| {
            if (c == '\n' or c == '\r') {
                try writer.writeByte(' ');
            } else {
                try writer.writeByte(c);
            }
        }
        try writer.writeByte('"');
    } else {
        // Both quote types — single-quote and escape the inner single quotes
        // by splitting: 'can'"'"'t' → valid but ugly. Simpler: use double quote
        // and accept that the inner double-quotes may break. In practice this
        // almost never happens in CIF data. Fall back to replacing quotes.
        try writer.writeByte('"');
        for (val) |c| {
            if (c == '\n' or c == '\r') {
                try writer.writeByte(' ');
            } else if (c == '"') {
                try writer.writeByte('\'');
            } else {
                try writer.writeByte(c);
            }
        }
        try writer.writeByte('"');
    }
}

/// Write a CIF pair tag-value. For pairs (not in loops), semicolon text
/// fields are allowed since they appear on their own lines.
fn writePairCifValue(writer: anytype, val: []const u8) !void {
    if (val.len == 0) {
        try writer.writeByte('.');
        return;
    }
    var has_newline = false;
    for (val) |c| {
        if (c == '\n' or c == '\r') {
            has_newline = true;
            break;
        }
    }
    if (has_newline) {
        try writer.writeAll("\n;");
        try writer.writeAll(val);
        if (val[val.len - 1] != '\n') try writer.writeByte('\n');
        try writer.writeAll(";\n");
    } else {
        try writeCifValue(writer, val);
    }
}

/// Write a CIF pair tag-value, quoting the value if needed.
fn writePairValue(writer: anytype, tag: []const u8, val: []const u8) !void {
    try writer.writeAll(tag);
    var has_newline = false;
    for (val) |c| {
        if (c == '\n' or c == '\r') {
            has_newline = true;
            break;
        }
    }

    if (has_newline) {
        try writePairCifValue(writer, val);
        return;
    }

    try writer.writeByte(' ');
    try writePairCifValue(writer, val);
    try writer.writeByte('\n');
}

/// Write the _atom_site loop with all atoms.
fn writeAtomSite(writer: anytype, model: *const Model) !void {
    try writer.writeAll(
        \\loop_
        \\_atom_site.group_PDB
        \\_atom_site.id
        \\_atom_site.type_symbol
        \\_atom_site.label_atom_id
        \\_atom_site.label_comp_id
        \\_atom_site.label_asym_id
        \\_atom_site.label_seq_id
        \\_atom_site.Cartn_x
        \\_atom_site.Cartn_y
        \\_atom_site.Cartn_z
        \\_atom_site.occupancy
        \\_atom_site.B_iso_or_equiv
        \\_atom_site.label_alt_id
        \\
    );

    // Write atoms grouped by residue: heavy atoms first, then added H atoms.
    // The placer appends H atoms at the end of model.atoms, so we need to
    // reorder output so each residue's H atoms follow its heavy atoms.
    var serial: u32 = 1;
    for (model.residues.items, 0..) |res, res_idx| {
        // First pass: original heavy atoms in this residue's range
        for (model.atoms.items[res.atom_start..res.atom_end]) |atom| {
            try writeAtomRow(writer, model, atom, res, serial);
            serial += 1;
        }
        // Second pass: added H atoms belonging to this residue (appended at end)
        for (model.atoms.items) |atom| {
            if (!atom.is_added) continue;
            if (atom.residue_idx != res_idx) continue;
            // Skip absent H atoms (flipper sentinel position)
            if (mover_mod.isAbsentH(atom)) continue;
            try writeAtomRow(writer, model, atom, res, serial);
            serial += 1;
        }
    }
}

fn writeAtomRow(writer: anytype, model: *const Model, atom: Atom, res: Residue, serial: u32) !void {
    const chain = model.chains.items[res.chain_idx];
    const group = switch (res.entity_type) {
        .water, .non_polymer => "HETATM",
        else => "ATOM",
    };
    const elem_str = elementSymbol(atom.element_type);
    try writer.print(
        "{s} {d} {s} {s} {s} {s} {d} {d:.3} {d:.3} {d:.3} {d:.2} {d:.2} ",
        .{
            group,
            serial,
            elem_str,
            atom.nameSlice(),
            res.compIdSlice(),
            chain.labelSlice(),
            res.seq_id,
            atom.pos.x,
            atom.pos.y,
            atom.pos.z,
            atom.occupancy,
            atom.b_factor,
        },
    );
    try writeAltId(writer, atom.altloc);
    try writer.writeByte('\n');
}

fn writeAltId(writer: anytype, altloc: u8) !void {
    if (altloc == ' ') {
        try writer.writeByte('.');
    } else {
        try writer.writeByte(altloc);
    }
}

/// Write optimization log as a custom mmCIF category.
pub fn writeZreduceLog(
    writer: anytype,
    movers: []const mover_mod.Mover,
    residues: []const Residue,
    chains: []const Chain,
) !void {
    try writer.writeAll(
        \\loop_
        \\_zreduce_log.residue_id
        \\_zreduce_log.action
        \\_zreduce_log.orientation
        \\_zreduce_log.score
        \\
    );
    for (movers) |m| {
        const res = residues[m.residue_idx];
        const chain = chains[res.chain_idx];
        const action = switch (m.kind) {
            .single_h_rotator, .nh3_rotator, .methyl_rotator, .aromatic_methyl => "rotate",
            .amide_flip => "flip_amide",
            .his_flip => "flip_his",
        };
        try writer.print("{s}.{s}.{d} {s} {d} .\n", .{
            chain.labelSlice(),
            res.compIdSlice(),
            res.seq_id,
            action,
            m.best_orientation,
        });
    }
}

// ── Tests ─────────────────────────────────────────────────────────────────────

const testing = std.testing;

test "write mmCIF round-trip" {
    var mdl = Model.init(testing.allocator);
    defer mdl.deinit();

    // Add one chain
    var chain = Chain{};
    chain.setLabelAsymId("A");
    try mdl.chains.append(mdl.allocator, chain);

    // Add one residue
    var res = Residue{};
    res.setCompId("ALA");
    res.seq_id = 1;
    res.chain_idx = 0;
    res.atom_start = 0;
    res.atom_end = 2;
    try mdl.residues.append(mdl.allocator, res);

    // Add two atoms: N and CA
    var n_atom = Atom{ .pos = .{ .x = 1.0, .y = 2.0, .z = 3.0 } };
    n_atom.setName("N");
    n_atom.element_type = .N;
    n_atom.residue_idx = 0;
    try mdl.atoms.append(mdl.allocator, n_atom);

    var ca_atom = Atom{ .pos = .{ .x = 2.5, .y = 3.5, .z = 4.5 } };
    ca_atom.setName("CA");
    ca_atom.element_type = .C;
    ca_atom.residue_idx = 0;
    try mdl.atoms.append(mdl.allocator, ca_atom);

    var buf = std.ArrayList(u8).empty;
    defer buf.deinit(testing.allocator);

    try write(buf.writer(testing.allocator), &mdl, "TEST");

    const output = buf.items;
    try testing.expect(std.mem.indexOf(u8, output, "data_TEST") != null);
    try testing.expect(std.mem.indexOf(u8, output, "_atom_site.Cartn_x") != null);
    try testing.expect(std.mem.indexOf(u8, output, "ALA") != null);
    try testing.expect(std.mem.indexOf(u8, output, "ATOM") != null);
}

test "writeWithDocument preserves multiline pair and loop values" {
    const source =
        \\data_SAMPLE
        \\_struct.entry_id SAMPLE
        \\_struct.title
        \\;first line
        \\second line
        \\;
        \\#
        \\loop_
        \\_pdbx_data_usage.details
        \\_pdbx_data_usage.id
        \\;loop line 1
        \\loop line 2
        \\;
        \\1
        \\#
        \\loop_
        \\_atom_site.group_PDB
        \\_atom_site.id
        \\_atom_site.type_symbol
        \\_atom_site.label_atom_id
        \\_atom_site.label_comp_id
        \\_atom_site.label_asym_id
        \\_atom_site.label_seq_id
        \\_atom_site.Cartn_x
        \\_atom_site.Cartn_y
        \\_atom_site.Cartn_z
        \\_atom_site.occupancy
        \\_atom_site.B_iso_or_equiv
        \\_atom_site.label_alt_id
        \\ATOM 1 N N ALA A 1 0.0 0.0 0.0 1.0 10.0 .
        \\ATOM 2 C CA ALA A 1 1.5 0.0 0.0 1.0 10.0 .
        \\ATOM 3 C C ALA A 1 2.0 1.4 0.0 1.0 10.0 .
        \\ATOM 4 O O ALA A 1 3.2 1.5 0.0 1.0 10.0 .
        \\ATOM 5 C CB ALA A 1 1.9 -0.8 1.2 1.0 10.0 .
        \\#
    ;

    var doc = try cif.readString(testing.allocator, source);
    defer doc.deinit();

    var mdl = try mmcif.parseModel(testing.allocator, source);
    defer mdl.deinit();
    _ = try place.addHydrogens(&mdl, null);

    var buf: std.ArrayListUnmanaged(u8) = .empty;
    defer buf.deinit(testing.allocator);
    try writeWithDocument(buf.writer(testing.allocator), &mdl, &doc);

    var reparsed = try cif.readString(testing.allocator, buf.items);
    defer reparsed.deinit();

    const block = &reparsed.blocks.items[0];
    try testing.expectEqualStrings("SAMPLE", block.name);
    try testing.expectEqualStrings("first line\nsecond line\n", block.findValue("_struct.title").?);

    const data_usage = block.findLoop("_pdbx_data_usage.details").?;
    try testing.expectEqual(@as(usize, 1), data_usage.length());
    try testing.expectEqualStrings("loop line 1 loop line 2 ", data_usage.val(0, 0).?);

    const atom_site = block.findLoop("_atom_site.Cartn_x").?;
    try testing.expect(atom_site.length() > 5);
}

test "writeWithDocument preserves added hydrogen altloc" {
    const source = @embedFile("../test_data/ala_altloc.cif");

    var doc = try cif.readString(testing.allocator, source);
    defer doc.deinit();

    var mdl = try mmcif.parseModel(testing.allocator, source);
    defer mdl.deinit();
    place.applyChemistry(&mdl);
    _ = try place.addHydrogens(&mdl, null);

    var buf: std.ArrayListUnmanaged(u8) = .empty;
    defer buf.deinit(testing.allocator);
    try writeWithDocument(buf.writer(testing.allocator), &mdl, &doc);

    var reparsed = try cif.readString(testing.allocator, buf.items);
    defer reparsed.deinit();

    const block = &reparsed.blocks.items[0];
    const atom_site = block.findLoop("_atom_site.label_atom_id").?;
    const atom_name_col = atom_site.findTag("_atom_site.label_atom_id").?;
    const alt_id_col = atom_site.findTag("_atom_site.label_alt_id").?;

    var ha_a_count: u32 = 0;
    var ha_b_count: u32 = 0;
    for (0..atom_site.length()) |row| {
        const atom_name = atom_site.val(row, atom_name_col) orelse continue;
        if (!std.mem.eql(u8, atom_name, "HA")) continue;
        const alt_id = atom_site.val(row, alt_id_col) orelse continue;
        if (std.mem.eql(u8, alt_id, "A")) ha_a_count += 1;
        if (std.mem.eql(u8, alt_id, "B")) ha_b_count += 1;
    }

    try testing.expectEqual(@as(u32, 1), ha_a_count);
    try testing.expectEqual(@as(u32, 1), ha_b_count);
}

test "elementSymbol returns correct symbols" {
    try testing.expectEqualStrings("H", elementSymbol(.H));
    try testing.expectEqualStrings("H", elementSymbol(.Hpol));
    try testing.expectEqualStrings("C", elementSymbol(.Car));
    try testing.expectEqualStrings("N", elementSymbol(.Nacc));
    try testing.expectEqualStrings("Fe", elementSymbol(.Fe));
    try testing.expectEqualStrings("Cl", elementSymbol(.Cl));
}

test "writeCifValue quotes special-char-prefixed values" {
    // Test that values starting with [, ], {, } get quoted
    const cases = [_][]const u8{ "[bracket", "]close", "{brace", "}close", "?", "." };
    for (cases) |input| {
        var buf = std.ArrayList(u8).empty;
        defer buf.deinit(testing.allocator);
        try writeCifValue(buf.writer(testing.allocator), input);
        const output = buf.items;
        // All should be quoted (start with ' or ")
        try testing.expect(output[0] == '\'' or output[0] == '"');
    }
    // Plain value should NOT be quoted
    {
        var buf = std.ArrayList(u8).empty;
        defer buf.deinit(testing.allocator);
        try writeCifValue(buf.writer(testing.allocator), "hello");
        try testing.expectEqualStrings("hello", buf.items);
    }
}

test "writeZreduceLog produces expected output" {
    const allocator = testing.allocator;

    var chains = [_]Chain{blk: {
        var c = Chain{};
        c.setLabelAsymId("A");
        break :blk c;
    }};
    var residues = [_]Residue{blk: {
        var r = Residue{};
        r.setCompId("ASN");
        r.seq_id = 5;
        r.chain_idx = 0;
        break :blk r;
    }};

    const positions = try allocator.alloc(@import("../math.zig").Vec3(f32), 1);
    defer allocator.free(positions);
    positions[0] = .{ .x = 0, .y = 0, .z = 0 };

    const orientations = try allocator.alloc(mover_mod.Orientation, 1);
    defer allocator.free(orientations);
    orientations[0] = .{ .positions = positions };

    const atom_indices = try allocator.alloc(u32, 1);
    defer allocator.free(atom_indices);
    atom_indices[0] = 0;

    const movers = [_]mover_mod.Mover{.{
        .kind = .amide_flip,
        .residue_idx = 0,
        .atom_indices = atom_indices,
        .orientations = orientations,
        .best_orientation = 1,
        .current_orientation = 0,
        .allocator = allocator,
    }};

    var buf = std.ArrayList(u8).empty;
    defer buf.deinit(allocator);

    try writeZreduceLog(buf.writer(allocator), &movers, &residues, &chains);

    const output = buf.items;
    try testing.expect(std.mem.indexOf(u8, output, "_zreduce_log.residue_id") != null);
    try testing.expect(std.mem.indexOf(u8, output, "flip_amide") != null);
    try testing.expect(std.mem.indexOf(u8, output, "A.ASN.5") != null);
}
