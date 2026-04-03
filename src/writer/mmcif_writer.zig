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
const format = @import("format.zig");

const elementSymbol = format.elementSymbol;
const writeFixedFloat3 = format.writeFixedFloat3;
const writeFixedFloat2 = format.writeFixedFloat2;
const writeAtomName = format.writeAtomName;
const writeAltId = format.writeAltId;
const writeCifValue = format.writeCifValue;
const writeCifValueInLoop = format.writeCifValueInLoop;
const writePairCifValue = format.writePairCifValue;
const writePairValue = format.writePairValue;

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
    return writeWithDocumentWithPolicy(writer, model, doc, .{});
}

/// Write a model as mmCIF with explicit bond/output policy.
pub fn writeWithDocumentWithPolicy(writer: anytype, model: *const Model, doc: ?*const Document, bond_policy: place.BondPolicy) !void {
    if (doc) |d| {
        if (d.blocks.items.len > 0) {
            const block = &d.blocks.items[0];
            try writer.print("data_{s}\n", .{block.name});

            // Write all items, replacing _atom_site loop.
            // Each item is preceded by a single '#' separator line.
            var atom_site_written = false;
            for (block.items.items) |item| {
                try writer.writeAll("#\n");
                switch (item) {
                    .pair => |p| {
                        try writePairValue(writer, p.tag, p.value);
                    },
                    .loop => |l| {
                        if (isAtomSiteLoop(&l)) {
                            try writeAtomSitePreserving(writer, model, &l, bond_policy.output_isotope);
                            atom_site_written = true;
                        } else {
                            try writeLoopBody(writer, &l);
                        }
                    },
                }
            }
            try writer.writeAll("#\n");

            // If no atom_site loop was found in original, append it
            if (!atom_site_written) {
                try writeAtomSite(writer, model, bond_policy.output_isotope);
                try writer.writeAll("#\n");
            }
            return;
        }
    }

    // Fallback: atom-site-only mode
    try writer.writeAll("data_ZREDUCE\n#\n");
    try writeAtomSite(writer, model, bond_policy.output_isotope);
    try writer.writeAll("#\n");
}

/// Write atom-site-only output (no original document preservation).
pub fn write(writer: anytype, model: *const Model, block_name: []const u8) !void {
    return writeWithPolicy(writer, model, block_name, .{});
}

/// Write atom-site-only output with explicit bond/output policy.
pub fn writeWithPolicy(writer: anytype, model: *const Model, block_name: []const u8, bond_policy: place.BondPolicy) !void {
    try writer.print("data_{s}\n#\n", .{block_name});
    try writeAtomSite(writer, model, bond_policy.output_isotope);
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
fn writeAtomSitePreserving(writer: anytype, model: *const Model, orig_loop: *const Loop, output_isotope: place.OutputIsotope) !void {
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

    // Pre-index added H atoms by residue_idx to avoid O(N*R) inner scan.
    // Build a counting-sort index: added_indices[added_offsets[r]..added_offsets[r+1]]
    // gives the atom indices for residue r, in original order.
    const n_residues = model.residues.items.len;
    const allocator = model.allocator;

    // Phase 1: count added atoms per residue.
    const added_counts = try allocator.alloc(u32, n_residues + 1);
    defer allocator.free(added_counts);
    @memset(added_counts, 0);
    for (model.atoms.items) |atom| {
        if (!atom.is_added) continue;
        added_counts[atom.residue_idx] += 1;
    }

    // Phase 2: build prefix-sum offsets (length n_residues + 1).
    const added_offsets = try allocator.alloc(u32, n_residues + 1);
    defer allocator.free(added_offsets);
    added_offsets[0] = 0;
    for (0..n_residues) |r| {
        added_offsets[r + 1] = added_offsets[r] + added_counts[r];
    }
    const total_added = added_offsets[n_residues];

    // Phase 3: fill sorted atom indices (reuse added_counts as cursor).
    const added_indices = try allocator.alloc(u32, total_added);
    defer allocator.free(added_indices);
    @memset(added_counts, 0);
    for (model.atoms.items, 0..) |atom, atom_idx| {
        if (!atom.is_added) continue;
        const r = atom.residue_idx;
        added_indices[added_offsets[r] + added_counts[r]] = @intCast(atom_idx);
        added_counts[r] += 1;
    }

    // Write rows grouped by residue: original heavy atoms then added H
    for (model.residues.items, 0..) |res, res_idx| {
        const chain = model.chains.items[res.chain_idx];

        // Original heavy atoms: write rows from original loop data,
        // but use model atom coordinates (which may have been updated by
        // optimizers, e.g. amide flip swapping O/N positions).
        for (res.atom_start..res.atom_end) |orig_row_idx| {
            if (orig_row_idx < orig_loop.length()) {
                const atom = model.atoms.items[orig_row_idx];
                for (0..w) |col| {
                    if (col > 0) try writer.writeByte(' ');
                    if (cm.id != null and col == cm.id.?) {
                        try writer.print("{d}", .{serial});
                    } else if (cm.cartn_x != null and col == cm.cartn_x.?) {
                        try writer.print("{d:.3}", .{atom.pos.x});
                    } else if (cm.cartn_y != null and col == cm.cartn_y.?) {
                        try writer.print("{d:.3}", .{atom.pos.y});
                    } else if (cm.cartn_z != null and col == cm.cartn_z.?) {
                        try writer.print("{d:.3}", .{atom.pos.z});
                    } else {
                        const val = orig_loop.val(orig_row_idx, col) orelse ".";
                        try writeCifValueInLoop(writer, val);
                    }
                }
                try writer.writeByte('\n');
                serial += 1;
            }
        }

        // Added H atoms for this residue — look up via pre-built index.
        const h_start = added_offsets[res_idx];
        const h_end = added_offsets[res_idx + 1];
        for (added_indices[h_start..h_end]) |atom_idx| {
            const atom = model.atoms.items[atom_idx];
            // Skip absent H atoms (flipper sentinel)
            if (mover_mod.isAbsentH(atom)) continue;

            for (0..w) |col| {
                if (col > 0) try writer.writeByte(' ');
                if (cm.group_pdb != null and col == cm.group_pdb.?) {
                    try writer.writeAll("ATOM");
                } else if (cm.id != null and col == cm.id.?) {
                    try writer.print("{d}", .{serial});
                } else if (cm.type_symbol != null and col == cm.type_symbol.?) {
                    try writer.writeAll(atomTypeSymbol(atom, output_isotope));
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
                        try writeCifValueInLoop(writer, val);
                    } else {
                        try writer.writeAll(".");
                    }
                } else if (cm.label_seq_id != null and col == cm.label_seq_id.?) {
                    try writer.print("{d}", .{res.seq_id});
                } else if (cm.cartn_x != null and col == cm.cartn_x.?) {
                    try writeFixedFloat3(writer, atom.pos.x);
                } else if (cm.cartn_y != null and col == cm.cartn_y.?) {
                    try writeFixedFloat3(writer, atom.pos.y);
                } else if (cm.cartn_z != null and col == cm.cartn_z.?) {
                    try writeFixedFloat3(writer, atom.pos.z);
                } else if (cm.occupancy != null and col == cm.occupancy.?) {
                    try writeFixedFloat2(writer, atom.occupancy);
                } else if (cm.b_factor != null and col == cm.b_factor.?) {
                    try writeFixedFloat2(writer, atom.b_factor);
                } else if (cm.auth_seq_id != null and col == cm.auth_seq_id.?) {
                    // Copy from first original row of this residue
                    const first_orig = res.atom_start;
                    if (first_orig < orig_loop.length()) {
                        const val = orig_loop.val(first_orig, col) orelse ".";
                        try writeCifValueInLoop(writer, val);
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
                        try writeCifValueInLoop(writer, val);
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

/// Write a CIF loop body (tags + data rows) without surrounding '#' separators.
fn writeLoopBody(writer: anytype, loop: *const Loop) !void {
    try writer.writeAll("loop_\n");
    for (loop.tags.items) |tag| {
        try writer.print("{s}\n", .{tag});
    }
    const w = loop.width();
    if (w > 0) {
        for (0..loop.length()) |row| {
            for (0..w) |col| {
                if (col > 0) try writer.writeByte(' ');
                const val = loop.val(row, col) orelse ".";
                try writeCifValueInLoop(writer, val);
            }
            try writer.writeByte('\n');
        }
    }
}

/// Write the _atom_site loop with all atoms.
fn writeAtomSite(writer: anytype, model: *const Model, output_isotope: place.OutputIsotope) !void {
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
    // Pre-index added H atoms by residue_idx to avoid O(N*R) inner scan.
    const n_residues = model.residues.items.len;
    const allocator = model.allocator;
    const added_counts = try allocator.alloc(u32, n_residues + 1);
    defer allocator.free(added_counts);
    @memset(added_counts, 0);
    for (model.atoms.items) |atom| {
        if (!atom.is_added) continue;
        added_counts[atom.residue_idx] += 1;
    }
    const added_offsets = try allocator.alloc(u32, n_residues + 1);
    defer allocator.free(added_offsets);
    added_offsets[0] = 0;
    for (0..n_residues) |r| {
        added_offsets[r + 1] = added_offsets[r] + added_counts[r];
    }
    const total_added = added_offsets[n_residues];
    const added_indices = try allocator.alloc(u32, total_added);
    defer allocator.free(added_indices);
    @memset(added_counts, 0);
    for (model.atoms.items, 0..) |atom, atom_idx| {
        if (!atom.is_added) continue;
        const r = atom.residue_idx;
        added_indices[added_offsets[r] + added_counts[r]] = @intCast(atom_idx);
        added_counts[r] += 1;
    }

    var serial: u32 = 1;
    for (model.residues.items, 0..) |res, res_idx| {
        // Original heavy atoms in this residue's range
        for (model.atoms.items[res.atom_start..res.atom_end]) |atom| {
            try writeAtomRow(writer, model, atom, res, serial, output_isotope);
            serial += 1;
        }
        // Added H atoms for this residue — look up via pre-built index.
        const h_start = added_offsets[res_idx];
        const h_end = added_offsets[res_idx + 1];
        for (added_indices[h_start..h_end]) |atom_idx| {
            const atom = model.atoms.items[atom_idx];
            // Skip absent H atoms (flipper sentinel position)
            if (mover_mod.isAbsentH(atom)) continue;
            try writeAtomRow(writer, model, atom, res, serial, output_isotope);
            serial += 1;
        }
    }
}

fn writeAtomRow(writer: anytype, model: *const Model, atom: Atom, res: Residue, serial: u32, output_isotope: place.OutputIsotope) !void {
    const chain = model.chains.items[res.chain_idx];
    const group = switch (res.entity_type) {
        .water, .non_polymer, .branched => "HETATM",
        .polymer, .unknown => "ATOM",
    };
    const elem_str = atomTypeSymbol(atom, output_isotope);
    try writer.print(
        "{s} {d} {s} {s} {s} {s} {d} ",
        .{
            group,
            serial,
            elem_str,
            atom.nameSlice(),
            res.compIdSlice(),
            chain.labelSlice(),
            res.seq_id,
        },
    );
    try writeFixedFloat3(writer, atom.pos.x);
    try writer.writeByte(' ');
    try writeFixedFloat3(writer, atom.pos.y);
    try writer.writeByte(' ');
    try writeFixedFloat3(writer, atom.pos.z);
    try writer.writeByte(' ');
    try writeFixedFloat2(writer, atom.occupancy);
    try writer.writeByte(' ');
    try writeFixedFloat2(writer, atom.b_factor);
    try writer.writeByte(' ');
    try writeAltId(writer, atom.altloc);
    try writer.writeByte('\n');
}

fn atomTypeSymbol(atom: Atom, output_isotope: place.OutputIsotope) []const u8 {
    if (output_isotope == .deuterium and atom.is_added and atom.is_hydrogen) return "D";
    return elementSymbol(atom.element_type);
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
    _ = try place.addHydrogens(&mdl, null, null);

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
    _ = try place.addHydrogens(&mdl, null, null);

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

test "writeWithPolicy outputs D for added hydrogens when isotope is deuterium" {
    const source = @embedFile("../test_data/tiny.cif");

    var mdl = try mmcif.parseModel(testing.allocator, source);
    defer mdl.deinit();
    _ = try place.addHydrogens(&mdl, null, null);

    var out: std.ArrayListUnmanaged(u8) = .empty;
    defer out.deinit(testing.allocator);
    try writeWithPolicy(out.writer(testing.allocator), &mdl, "TEST", .{
        .output_isotope = .deuterium,
    });

    var doc = try cif.readString(testing.allocator, out.items);
    defer doc.deinit();

    const block = &doc.blocks.items[0];
    const atom_site = block.findLoop("_atom_site.label_atom_id").?;
    const atom_name_col = atom_site.findTag("_atom_site.label_atom_id").?;
    const type_symbol_col = atom_site.findTag("_atom_site.type_symbol").?;

    var found_added_h: bool = false;
    for (0..atom_site.length()) |row| {
        const atom_name = atom_site.val(row, atom_name_col) orelse continue;
        if (!std.mem.startsWith(u8, atom_name, "H")) continue;
        found_added_h = true;
        const type_symbol = atom_site.val(row, type_symbol_col) orelse continue;
        try testing.expectEqualStrings("D", type_symbol);
    }
    try testing.expect(found_added_h);
}

test "writeWithDocumentWithPolicy outputs D in preserving mode" {
    const source = @embedFile("../test_data/tiny.cif");

    var doc = try cif.readString(testing.allocator, source);
    defer doc.deinit();

    var mdl = try mmcif.parseModel(testing.allocator, source);
    defer mdl.deinit();
    _ = try place.addHydrogens(&mdl, null, null);

    var out: std.ArrayListUnmanaged(u8) = .empty;
    defer out.deinit(testing.allocator);
    try writeWithDocumentWithPolicy(out.writer(testing.allocator), &mdl, &doc, .{
        .output_isotope = .deuterium,
    });

    var parsed = try cif.readString(testing.allocator, out.items);
    defer parsed.deinit();

    const block = &parsed.blocks.items[0];
    const atom_site = block.findLoop("_atom_site.label_atom_id").?;
    const atom_name_col = atom_site.findTag("_atom_site.label_atom_id").?;
    const type_symbol_col = atom_site.findTag("_atom_site.type_symbol").?;

    var found_d: bool = false;
    for (0..atom_site.length()) |row| {
        const atom_name = atom_site.val(row, atom_name_col) orelse continue;
        const type_symbol = atom_site.val(row, type_symbol_col) orelse continue;
        if (std.mem.startsWith(u8, atom_name, "H")) {
            // Added H atoms should have type_symbol "D"
            try testing.expectEqualStrings("D", type_symbol);
            found_d = true;
        }
    }
    try testing.expect(found_d);
}

test "deuterium mode preserves existing H atoms as H" {
    const source = @embedFile("../test_data/ala_with_h.cif");

    var mdl = try mmcif.parseModel(testing.allocator, source);
    defer mdl.deinit();
    _ = try place.addHydrogens(&mdl, null, null);

    var out: std.ArrayListUnmanaged(u8) = .empty;
    defer out.deinit(testing.allocator);
    try writeWithPolicy(out.writer(testing.allocator), &mdl, "TEST", .{
        .output_isotope = .deuterium,
    });

    var doc = try cif.readString(testing.allocator, out.items);
    defer doc.deinit();

    const block = &doc.blocks.items[0];
    const atom_site = block.findLoop("_atom_site.label_atom_id").?;
    const type_symbol_col = atom_site.findTag("_atom_site.type_symbol").?;
    const atom_name_col = atom_site.findTag("_atom_site.label_atom_id").?;

    var found_original_h = false;
    var found_added_d = false;
    for (0..atom_site.length()) |row| {
        const atom_name = atom_site.val(row, atom_name_col) orelse continue;
        const type_symbol = atom_site.val(row, type_symbol_col) orelse continue;
        if (!std.mem.startsWith(u8, atom_name, "H")) continue;

        // Check if this is an original atom (in the first residue's range) or added
        // ala_with_h.cif has pre-existing H atoms; added H will have type_symbol D
        if (std.mem.eql(u8, type_symbol, "H")) {
            found_original_h = true;
        } else if (std.mem.eql(u8, type_symbol, "D")) {
            found_added_d = true;
        }
    }
    // ala_with_h.cif has pre-existing H, so we should see both H and D
    try testing.expect(found_original_h);
}

test "writeAtomSitePreserving serial numbers are contiguous even with out-of-range rows" {
    // Construct a document with an atom_site loop that has fewer rows than
    // the model expects (simulating orig_row_idx >= orig_loop.length()).
    // The serial counter must only increment for actually-written rows.
    const source =
        \\data_TEST
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
        \\#
    ;

    var doc = try cif.readString(testing.allocator, source);
    defer doc.deinit();

    var mdl = try mmcif.parseModel(testing.allocator, source);
    defer mdl.deinit();

    // Extend the residue's atom range beyond loop length to simulate gap
    mdl.residues.items[0].atom_end = 4; // only 2 rows exist in loop

    var buf: std.ArrayListUnmanaged(u8) = .empty;
    defer buf.deinit(testing.allocator);
    try writeWithDocument(buf.writer(testing.allocator), &mdl, &doc);

    // Parse output and verify serial numbers are 1 and 2 (no gap)
    var reparsed = try cif.readString(testing.allocator, buf.items);
    defer reparsed.deinit();
    const block = &reparsed.blocks.items[0];
    const atom_site = block.findLoop("_atom_site.id").?;
    const id_col = atom_site.findTag("_atom_site.id").?;
    try testing.expectEqual(@as(usize, 2), atom_site.length());
    try testing.expectEqualStrings("1", atom_site.val(0, id_col).?);
    try testing.expectEqualStrings("2", atom_site.val(1, id_col).?);
}

test "writeWithDocument preserves bare null alt ids on original rows" {
    const source = @embedFile("../test_data/tiny.cif");

    var doc = try cif.readString(testing.allocator, source);
    defer doc.deinit();

    var mdl = try mmcif.parseModel(testing.allocator, source);
    defer mdl.deinit();

    var buf: std.ArrayListUnmanaged(u8) = .empty;
    defer buf.deinit(testing.allocator);
    try writeWithDocument(buf.writer(testing.allocator), &mdl, &doc);

    try testing.expect(std.mem.indexOf(u8, buf.items, "'.'") == null);
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
