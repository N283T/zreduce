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

/// Write a model as mmCIF to the given writer.
/// Includes all original atoms plus added hydrogens.
pub fn write(writer: anytype, model: *const Model, block_name: []const u8) !void {
    try writer.print("data_{s}\n#\n", .{block_name});
    try writeAtomSite(writer, model);
    try writer.writeAll("#\n");
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
    const alt: [1]u8 = if (atom.altloc != ' ') [_]u8{atom.altloc} else [_]u8{'.'};

    try writer.print(
        "{s} {d} {s} {s} {s} {s} {d} {d:.3} {d:.3} {d:.3} {d:.2} {d:.2} {s}\n",
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
            &alt,
        },
    );
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

test "elementSymbol returns correct symbols" {
    try testing.expectEqualStrings("H", elementSymbol(.H));
    try testing.expectEqualStrings("H", elementSymbol(.Hpol));
    try testing.expectEqualStrings("C", elementSymbol(.Car));
    try testing.expectEqualStrings("N", elementSymbol(.Nacc));
    try testing.expectEqualStrings("Fe", elementSymbol(.Fe));
    try testing.expectEqualStrings("Cl", elementSymbol(.Cl));
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
