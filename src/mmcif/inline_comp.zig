//! Inline component dictionary parsing and leaving atom flagging.

const std = @import("std");
const Allocator = std.mem.Allocator;

const cif = @import("../cif.zig");
const model_mod = @import("../model.zig");
const ccd_mod = @import("../ccd.zig");

const Model = model_mod.Model;

/// Helper: find atom index in a slice by name.
fn findAtomIdx(atoms: []const ccd_mod.CompAtom, name: []const u8) ?u16 {
    for (atoms, 0..) |*a, i| {
        if (std.mem.eql(u8, a.nameSlice(), name)) return @intCast(i);
    }
    return null;
}

/// Temporary bond record storing atom names (strings) before index resolution.
const RawBond = struct {
    atom1: []const u8,
    atom2: []const u8,
    order: ccd_mod.BondOrder,
    aromatic: bool,
};

/// Parse inline `_chem_comp_atom` and `_chem_comp_bond` loops from a structure
/// CIF block into a `ccd_mod.ComponentDict`.
///
/// Returns `null` when neither loop is present in the block.
/// The returned ComponentDict owns all its memory.
pub fn parseInlineComponents(allocator: Allocator, block: *const cif.Block) !?ccd_mod.ComponentDict {
    // Check that at least the bond loop is present; if not, return null.
    const bond_loop = block.findLoop("_chem_comp_bond.comp_id") orelse return null;

    // ── Phase 1: group atoms by comp_id ──────────────────────────────────────

    // comp_id -> list of CompAtom (temporary, unowned strings)
    var atom_groups = std.StringHashMap(std.ArrayListUnmanaged(ccd_mod.CompAtom)).init(allocator);
    defer {
        var it = atom_groups.valueIterator();
        while (it.next()) |list| list.deinit(allocator);
        atom_groups.deinit();
    }

    if (block.findLoop("_chem_comp_atom.comp_id")) |atom_loop| {
        const col_comp = atom_loop.findTag("_chem_comp_atom.comp_id") orelse return null;
        const col_atom = atom_loop.findTag("_chem_comp_atom.atom_id");
        const col_sym = atom_loop.findTag("_chem_comp_atom.type_symbol");
        const col_charge = atom_loop.findTag("_chem_comp_atom.charge");
        const col_leaving = atom_loop.findTag("_chem_comp_atom.pdbx_leaving_atom_flag");
        const col_arom = atom_loop.findTag("_chem_comp_atom.pdbx_aromatic_flag");

        const nrows = atom_loop.length();
        for (0..nrows) |row| {
            const comp_id = cif.asString(atom_loop.val(row, col_comp) orelse continue);

            var atom = ccd_mod.CompAtom{};

            if (col_atom) |c| {
                if (atom_loop.val(row, c)) |v| {
                    const s = cif.asString(v);
                    const len = @min(s.len, 4);
                    atom.name_len = @intCast(len);
                    @memcpy(atom.name[0..len], s[0..len]);
                }
            }

            if (col_sym) |c| {
                if (atom_loop.val(row, c)) |v| {
                    const s = cif.asString(v);
                    const len = @min(s.len, 2);
                    @memcpy(atom.element_symbol[0..len], s[0..len]);
                }
            }

            if (col_charge) |c| {
                if (atom_loop.val(row, c)) |v| {
                    atom.charge = cif.value.asIntOr(i8, v, 0);
                }
            }

            if (col_leaving) |c| {
                if (atom_loop.val(row, c)) |v| {
                    atom.leaving = std.ascii.eqlIgnoreCase(cif.asString(v), "Y");
                }
            }

            if (col_arom) |c| {
                if (atom_loop.val(row, c)) |v| {
                    atom.aromatic = std.ascii.eqlIgnoreCase(cif.asString(v), "Y");
                }
            }

            const gop = try atom_groups.getOrPut(comp_id);
            if (!gop.found_existing) gop.value_ptr.* = .empty;
            try gop.value_ptr.append(allocator, atom);
        }
    }

    // ── Phase 2: group bonds by comp_id ──────────────────────────────────────

    // comp_id -> list of RawBond (strings are slices into the CIF source)
    var bond_groups = std.StringHashMap(std.ArrayListUnmanaged(RawBond)).init(allocator);
    defer {
        var it = bond_groups.valueIterator();
        while (it.next()) |list| list.deinit(allocator);
        bond_groups.deinit();
    }

    {
        const col_comp = bond_loop.findTag("_chem_comp_bond.comp_id") orelse return null;
        const col_a1 = bond_loop.findTag("_chem_comp_bond.atom_id_1");
        const col_a2 = bond_loop.findTag("_chem_comp_bond.atom_id_2");
        const col_order = bond_loop.findTag("_chem_comp_bond.value_order");
        const col_arom = bond_loop.findTag("_chem_comp_bond.pdbx_aromatic_flag");

        const nrows = bond_loop.length();
        for (0..nrows) |row| {
            const comp_id = cif.asString(bond_loop.val(row, col_comp) orelse continue);

            const atom1 = if (col_a1) |c|
                cif.asString(bond_loop.val(row, c) orelse continue)
            else
                continue;

            const atom2 = if (col_a2) |c|
                cif.asString(bond_loop.val(row, c) orelse continue)
            else
                continue;

            const order: ccd_mod.BondOrder = if (col_order) |c|
                if (bond_loop.val(row, c)) |v| ccd_mod.BondOrder.fromString(cif.asString(v)) else .unknown
            else
                .unknown;

            const aromatic: bool = if (col_arom) |c|
                if (bond_loop.val(row, c)) |v| std.ascii.eqlIgnoreCase(cif.asString(v), "Y") else false
            else
                false;

            const gop = try bond_groups.getOrPut(comp_id);
            if (!gop.found_existing) gop.value_ptr.* = .empty;
            try gop.value_ptr.append(allocator, RawBond{
                .atom1 = atom1,
                .atom2 = atom2,
                .order = order,
                .aromatic = aromatic,
            });
        }
    }

    // ── Phase 3: build ComponentDict ─────────────────────────────────────────

    var dict = ccd_mod.ComponentDict{
        .components = std.StringHashMap(ccd_mod.Component).init(allocator),
        .allocator = allocator,
    };
    errdefer dict.deinit();

    // Collect all unique comp_ids from both groups.
    var all_ids = std.StringHashMap(void).init(allocator);
    defer all_ids.deinit();
    {
        var it = atom_groups.keyIterator();
        while (it.next()) |k| try all_ids.put(k.*, {});
    }
    {
        var it = bond_groups.keyIterator();
        while (it.next()) |k| try all_ids.put(k.*, {});
    }

    var id_it = all_ids.keyIterator();
    while (id_it.next()) |id_ptr| {
        const comp_id = id_ptr.*;

        // Atoms slice (may be empty if no atom loop).
        const raw_atoms: []const ccd_mod.CompAtom = if (atom_groups.get(comp_id)) |list|
            list.items
        else
            &[_]ccd_mod.CompAtom{};

        const owned_atoms = try allocator.dupe(ccd_mod.CompAtom, raw_atoms);
        errdefer allocator.free(owned_atoms);

        // Resolve bonds to index-based CompBond.
        var resolved_bonds = std.ArrayListUnmanaged(ccd_mod.CompBond).empty;
        defer resolved_bonds.deinit(allocator);

        if (bond_groups.get(comp_id)) |raw_bond_list| {
            for (raw_bond_list.items) |rb| {
                const idx1 = findAtomIdx(owned_atoms, rb.atom1) orelse continue;
                const idx2 = findAtomIdx(owned_atoms, rb.atom2) orelse continue;
                try resolved_bonds.append(allocator, ccd_mod.CompBond{
                    .atom_idx_1 = idx1,
                    .atom_idx_2 = idx2,
                    .order = rb.order,
                    .aromatic = rb.aromatic,
                });
            }
        }

        const owned_bonds = try allocator.dupe(ccd_mod.CompBond, resolved_bonds.items);
        errdefer allocator.free(owned_bonds);

        // Single allocation shared as both HashMap key and Component.comp_id.
        const key = try allocator.dupe(u8, comp_id);
        errdefer allocator.free(key);

        const owned_type = try allocator.dupe(u8, "");
        errdefer allocator.free(owned_type);

        try dict.components.put(key, ccd_mod.Component{
            .comp_id = key,
            .comp_type = owned_type,
            .atoms = owned_atoms,
            .bonds = owned_bonds,
        });
    }

    return dict;
}

/// For residues containing inter-residue bonded atoms, flag CCD leaving atoms
/// as bonded_inter_residue too. Leaving atoms (pdbx_leaving_atom_flag=Y) should
/// not exist in a properly bonded structure, so H placement on them is incorrect.
pub fn flagLeavingAtoms(
    mdl: *Model,
    inline_dict: ?*const ccd_mod.ComponentDict,
    ccd_dict: ?*const ccd_mod.ComponentDict,
) void {
    for (mdl.residues.items) |res| {
        const atoms = mdl.atoms.items[res.atom_start..res.atom_end];

        // Check if this residue has any inter-residue bonded atoms
        var has_bonded = false;
        for (atoms) |atom| {
            if (atom.flags.bonded_inter_residue) {
                has_bonded = true;
                break;
            }
        }
        if (!has_bonded) continue;

        // Look up CCD component (inline priority)
        const comp_id = res.compIdSlice();
        const component = if (inline_dict) |d| d.get(comp_id) else null;
        const effective = component orelse if (ccd_dict) |d| d.get(comp_id) else null;
        const comp = effective orelse continue;

        // Flag leaving atoms in the model
        for (comp.atoms) |comp_atom| {
            if (!comp_atom.leaving) continue;
            const leaving_name = comp_atom.nameSlice();
            // Find this atom in the model residue and mutate via direct index
            for (res.atom_start..res.atom_end) |i| {
                if (std.mem.eql(u8, mdl.atoms.items[i].nameSlice(), leaving_name)) {
                    mdl.atoms.items[i].flags.bonded_inter_residue = true;
                }
            }
        }
    }
}

// ── Tests ─────────────────────────────────────────────────────────────────────

const testing = std.testing;

test "parseInlineComponents returns ComponentDict" {
    const source = @embedFile("../test_data/inline_comp.cif");
    var doc = cif.readString(testing.allocator, source) catch unreachable;
    defer doc.deinit();
    const block = &doc.blocks.items[0];

    var dict = try parseInlineComponents(testing.allocator, block);
    defer if (dict) |*d| d.deinit();

    try testing.expect(dict != null);
    const ala = dict.?.get("ALA");
    try testing.expect(ala != null);
    try testing.expectEqual(@as(usize, 11), ala.?.atoms.len);
    try testing.expectEqual(@as(usize, 10), ala.?.bonds.len);
}

test "parseInlineComponents returns null when no inline data" {
    const source = @embedFile("../test_data/tiny.cif");
    var doc = cif.readString(testing.allocator, source) catch unreachable;
    defer doc.deinit();
    const block = &doc.blocks.items[0];

    var dict = try parseInlineComponents(testing.allocator, block);
    defer if (dict) |*d| d.deinit();

    try testing.expect(dict == null);
}

test "flagLeavingAtoms flags CCD leaving atoms on bonded residues" {
    const mmcif = @import("../mmcif.zig");
    const source = @embedFile("../test_data/leaving_atom.cif");
    var doc = cif.readString(testing.allocator, source) catch unreachable;
    defer doc.deinit();
    const block = &doc.blocks.items[0];

    var mdl = try mmcif.parseModel(testing.allocator, source);
    defer mdl.deinit();

    // Parse inline components (has leaving atom info)
    var inline_dict = try parseInlineComponents(testing.allocator, block);
    defer if (inline_dict) |*d| d.deinit();

    // Parse struct_conn -- flags C(0) and N(2) as bonded
    var lookup = try mmcif.buildAtomLookup(testing.allocator, block);
    defer lookup.deinit();
    try mmcif.parseStructConn(&mdl, block, &lookup);

    // Verify: C(0) bonded, OXT(1) NOT bonded yet, N(2) bonded
    try testing.expect(mdl.atoms.items[0].flags.bonded_inter_residue); // C
    try testing.expect(!mdl.atoms.items[1].flags.bonded_inter_residue); // OXT - not yet
    try testing.expect(mdl.atoms.items[2].flags.bonded_inter_residue); // N

    // Flag leaving atoms
    flagLeavingAtoms(&mdl, if (inline_dict) |*d| d else null, null);

    // Now OXT should also be flagged (it's a CCD leaving atom in a bonded residue)
    try testing.expect(mdl.atoms.items[1].flags.bonded_inter_residue); // OXT - now flagged
}
