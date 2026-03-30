//! CCD-derived hydrogen placement plan generation.
//!
//! Derives PlacementPlan from CCD component topology. For each hydrogen atom
//! in the component that is NOT already present in the model, determine the
//! placement type from hybridization (inferred from bond orders).
//! Used as fallback when no hardcoded plan exists (standard, nucleotide, modified).

const std = @import("std");
const ccd = @import("../ccd.zig");
const standard = @import("standard.zig");
const PlacementPlan = standard.PlacementPlan;
const PlacementType = standard.PlacementType;
const MoverHint = standard.MoverHint;
const element_mod = @import("../element.zig");

// ---------------------------------------------------------------------------
// Bond analysis types
// ---------------------------------------------------------------------------

const Hybridization = enum { sp3, sp2, sp, unknown };

const BondInfo = struct {
    total_bonds: u8,
    heavy_neighbor_count: u8,
    h_neighbor_count: u8,
    has_double: bool,
    has_triple: bool,
    has_aromatic: bool,
    hybridization: Hybridization,
};

// ---------------------------------------------------------------------------
// Public API
// ---------------------------------------------------------------------------

/// Derive placement plans for non-standard residue hydrogens from CCD component definition.
/// Only generates plans for H atoms not already present in existing_atom_names.
pub fn derivePlans(
    allocator: std.mem.Allocator,
    component: *const ccd.Component,
    existing_atom_names: []const [4]u8,
) ![]PlacementPlan {
    var plans = std.ArrayListUnmanaged(PlacementPlan){};
    errdefer plans.deinit(allocator);

    for (component.atoms, 0..) |atom, atom_idx| {
        // Skip non-hydrogen atoms
        if (atom.element_symbol[0] != 'H') continue;

        // Skip if already present in model
        if (nameExists(atom.name, existing_atom_names)) continue;

        // Find the heavy atom this H is bonded to
        const heavy_idx = findBondedHeavyAtom(component, @intCast(atom_idx)) orelse continue;
        const heavy_atom = component.atoms[heavy_idx];

        // Count bonds on heavy atom and determine hybridization
        const bond_info = analyzeBonds(component, @intCast(heavy_idx));

        // Determine placement type from hybridization
        const plan = deriveSinglePlan(component, atom, @intCast(heavy_idx), heavy_atom, bond_info) orelse continue;
        try plans.append(allocator, plan);
    }

    const result = try allocator.dupe(PlacementPlan, plans.items);
    plans.deinit(allocator);
    return result;
}

// ---------------------------------------------------------------------------
// Internal helpers
// ---------------------------------------------------------------------------

/// Check if atom name is in existing names list.
fn nameExists(name: [4]u8, list: []const [4]u8) bool {
    for (list) |existing| {
        if (std.mem.eql(u8, &name, &existing)) return true;
    }
    return false;
}

/// Find the non-H atom bonded to this H atom.
fn findBondedHeavyAtom(component: *const ccd.Component, h_idx: u16) ?u16 {
    for (component.bonds) |bond| {
        if (bond.atom_idx_1 == h_idx) {
            const other = bond.atom_idx_2;
            if (component.atoms[other].element_symbol[0] != 'H') return other;
        }
        if (bond.atom_idx_2 == h_idx) {
            const other = bond.atom_idx_1;
            if (component.atoms[other].element_symbol[0] != 'H') return other;
        }
    }
    return null;
}

/// Analyze bonds on the given atom index.
fn analyzeBonds(component: *const ccd.Component, atom_idx: u16) BondInfo {
    var info = BondInfo{
        .total_bonds = 0,
        .heavy_neighbor_count = 0,
        .h_neighbor_count = 0,
        .has_double = false,
        .has_triple = false,
        .has_aromatic = false,
        .hybridization = .unknown,
    };

    for (component.bonds) |bond| {
        const neighbor_idx: ?u16 = if (bond.atom_idx_1 == atom_idx)
            bond.atom_idx_2
        else if (bond.atom_idx_2 == atom_idx)
            bond.atom_idx_1
        else
            null;

        if (neighbor_idx) |ni| {
            info.total_bonds += 1;
            if (component.atoms[ni].element_symbol[0] == 'H') {
                info.h_neighbor_count += 1;
            } else {
                info.heavy_neighbor_count += 1;
            }
            switch (bond.order) {
                .double => info.has_double = true,
                .triple => info.has_triple = true,
                .aromatic => info.has_aromatic = true,
                else => {},
            }
        }
    }

    // Determine hybridization
    if (info.has_triple or (info.has_double and countDoubleBonds(component, atom_idx) >= 2)) {
        info.hybridization = .sp;
    } else if (info.has_double or info.has_aromatic) {
        info.hybridization = .sp2;
    } else {
        info.hybridization = .sp3;
    }

    return info;
}

/// Count the number of double bonds on an atom.
fn countDoubleBonds(component: *const ccd.Component, atom_idx: u16) u8 {
    var count: u8 = 0;
    for (component.bonds) |bond| {
        if ((bond.atom_idx_1 == atom_idx or bond.atom_idx_2 == atom_idx) and bond.order == .double) {
            count += 1;
        }
    }
    return count;
}

/// Select placement type based on hybridization and neighbor configuration.
fn deriveSinglePlan(
    component: *const ccd.Component,
    h_atom: ccd.CompAtom,
    heavy_idx: u16,
    heavy_atom: ccd.CompAtom,
    bond_info: BondInfo,
) ?PlacementPlan {
    // Determine bond length from element type
    const bond_len: f32 = switch (heavy_atom.element_symbol[0]) {
        'O' => 0.97,
        'N' => 1.02,
        'S' => 1.33,
        else => 1.10, // C-H default
    };

    // Determine atom type for the H
    const atom_type: element_mod.AtomType = switch (heavy_atom.element_symbol[0]) {
        'N', 'O', 'S' => .Hpol,
        else => if (bond_info.has_aromatic) .Har else .H,
    };

    // Determine mover hint
    const mover_hint: MoverHint = switch (heavy_atom.element_symbol[0]) {
        'O' => .rotate,
        'S' => .rotate,
        else => .none,
    };

    switch (bond_info.hybridization) {
        .sp3 => {
            if (bond_info.heavy_neighbor_count >= 3) {
                // HXR3: tetrahedral with 3 heavy neighbors
                const refs = findHeavyNeighborNames(component, heavy_idx, 3) orelse return null;
                return PlacementPlan{
                    .h_name = h_atom.name,
                    .placement_type = .hxr3,
                    .connected = .{ refs[0], refs[1], refs[2] },
                    .n_connected = 3,
                    .bond_len = bond_len,
                    .atom_type = atom_type,
                };
            } else if (bond_info.heavy_neighbor_count == 2) {
                // H2XR2: 2 heavy neighbors, dihedral-controlled
                const refs = findHeavyNeighborNames(component, heavy_idx, 2) orelse return null;
                const dihedral = estimateDihedral(component, h_atom, heavy_idx, refs[0]);
                return PlacementPlan{
                    .h_name = h_atom.name,
                    .placement_type = .h2xr2,
                    .connected = .{ heavy_atom.name, refs[0], blank },
                    .n_connected = 2,
                    .bond_len = bond_len,
                    .angle = 109.5,
                    .dihedral = dihedral,
                    .atom_type = atom_type,
                };
            } else if (bond_info.heavy_neighbor_count == 1) {
                // H3XR: dihedral-controlled (e.g., methyl)
                const refs = findHeavyNeighborNames(component, heavy_idx, 1) orelse return null;
                // Find a second reference for dihedral (neighbor of the heavy neighbor)
                const second_ref = findSecondReference(component, heavy_idx, refs[0]) orelse return null;
                const dihedral = estimateDihedralFromIdeal(component, h_atom, heavy_atom, refs[0]);
                return PlacementPlan{
                    .h_name = h_atom.name,
                    .placement_type = .h3xr,
                    .connected = .{ heavy_atom.name, refs[0], second_ref },
                    .n_connected = 3,
                    .bond_len = bond_len,
                    .angle = 109.5,
                    .dihedral = dihedral,
                    .atom_type = atom_type,
                    .mover_hint = if (bond_info.h_neighbor_count >= 3) .rotate_methyl else mover_hint,
                };
            }
            return null;
        },
        .sp2 => {
            // Planar placement
            if (bond_info.heavy_neighbor_count >= 2) {
                const refs = findHeavyNeighborNames(component, heavy_idx, 2) orelse return null;
                return PlacementPlan{
                    .h_name = h_atom.name,
                    .placement_type = .hxr2_planar,
                    .connected = .{ refs[0], refs[1], blank },
                    .n_connected = 2,
                    .bond_len = bond_len,
                    .atom_type = atom_type,
                };
            } else if (bond_info.heavy_neighbor_count == 1) {
                // Only one heavy neighbor on sp2 — use dihedral placement
                const refs = findHeavyNeighborNames(component, heavy_idx, 1) orelse return null;
                const second_ref = findSecondReference(component, heavy_idx, refs[0]) orelse return null;
                const dihedral = estimateDihedralFromIdeal(component, h_atom, heavy_atom, refs[0]);
                return PlacementPlan{
                    .h_name = h_atom.name,
                    .placement_type = .h3xr,
                    .connected = .{ heavy_atom.name, refs[0], second_ref },
                    .n_connected = 3,
                    .bond_len = bond_len,
                    .angle = 120.0,
                    .dihedral = dihedral,
                    .atom_type = atom_type,
                };
            }
            return null;
        },
        .sp => {
            // Linear placement (HXY)
            if (bond_info.heavy_neighbor_count >= 1) {
                const refs = findHeavyNeighborNames(component, heavy_idx, 1) orelse return null;
                return PlacementPlan{
                    .h_name = h_atom.name,
                    .placement_type = .hxy,
                    .connected = .{ heavy_atom.name, refs[0], blank },
                    .n_connected = 1,
                    .bond_len = bond_len,
                    .atom_type = atom_type,
                };
            }
            return null;
        },
        .unknown => return null,
    }
}

const blank: [4]u8 = .{ ' ', ' ', ' ', ' ' };

/// Find up to N heavy-atom neighbor names for the given atom.
fn findHeavyNeighborNames(component: *const ccd.Component, atom_idx: u16, comptime max: u8) ?[max][4]u8 {
    var result: [max][4]u8 = undefined;
    var count: u8 = 0;

    for (component.bonds) |bond| {
        const neighbor_idx: ?u16 = if (bond.atom_idx_1 == atom_idx)
            bond.atom_idx_2
        else if (bond.atom_idx_2 == atom_idx)
            bond.atom_idx_1
        else
            null;

        if (neighbor_idx) |ni| {
            if (component.atoms[ni].element_symbol[0] != 'H') {
                if (count < max) {
                    result[count] = component.atoms[ni].name;
                    count += 1;
                    if (count == max) return result;
                }
            }
        }
    }

    if (count >= max) return result;
    return null;
}

/// Find a reference atom bonded to the neighbor of the heavy atom (for dihedral reference).
fn findSecondReference(component: *const ccd.Component, heavy_idx: u16, first_ref_name: [4]u8) ?[4]u8 {
    // Find index of first_ref
    var first_ref_idx: ?u16 = null;
    for (component.atoms, 0..) |a, i| {
        if (std.mem.eql(u8, &a.name, &first_ref_name)) {
            first_ref_idx = @intCast(i);
            break;
        }
    }
    const ref_idx = first_ref_idx orelse return null;

    // Find a neighbor of ref_idx that is not heavy_idx and not H
    for (component.bonds) |bond| {
        const neighbor_idx: ?u16 = if (bond.atom_idx_1 == ref_idx)
            bond.atom_idx_2
        else if (bond.atom_idx_2 == ref_idx)
            bond.atom_idx_1
        else
            null;

        if (neighbor_idx) |ni| {
            if (ni != heavy_idx and component.atoms[ni].element_symbol[0] != 'H') {
                return component.atoms[ni].name;
            }
        }
    }
    return null;
}

/// Estimate dihedral angle for H placement from ideal coordinates.
/// TODO: Implement actual dihedral computation from CCD ideal coordinates.
/// Currently returns a fixed heuristic value (120.0 degrees).
fn estimateDihedral(component: *const ccd.Component, h_atom: ccd.CompAtom, heavy_idx: u16, ref_name: [4]u8) f32 {
    _ = component;
    _ = h_atom;
    _ = heavy_idx;
    _ = ref_name;
    return 120.0;
}

/// Estimate dihedral from ideal coordinates of the H and heavy atoms.
/// TODO: Implement actual dihedral computation from CCD ideal coordinates.
/// Currently returns a fixed heuristic value (180.0 degrees).
fn estimateDihedralFromIdeal(component: *const ccd.Component, h_atom: ccd.CompAtom, heavy_atom: ccd.CompAtom, ref_name: [4]u8) f32 {
    _ = component;
    _ = h_atom;
    _ = heavy_atom;
    _ = ref_name;
    return 180.0;
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

const testing = std.testing;

test "nameExists finds matching names" {
    const list = [_][4]u8{
        .{ 'N', ' ', ' ', ' ' },
        .{ 'C', 'A', ' ', ' ' },
    };
    try testing.expect(nameExists(.{ 'N', ' ', ' ', ' ' }, &list));
    try testing.expect(nameExists(.{ 'C', 'A', ' ', ' ' }, &list));
    try testing.expect(!nameExists(.{ 'H', ' ', ' ', ' ' }, &list));
}

test "derivePlans on empty component returns empty" {
    const comp = ccd.Component{
        .comp_id = "TST",
        .comp_type = "non-polymer",
        .atoms = &.{},
        .bonds = &.{},
    };
    const plans = try derivePlans(testing.allocator, &comp, &.{});
    defer testing.allocator.free(plans);
    try testing.expectEqual(@as(usize, 0), plans.len);
}

test "derivePlans skips existing H atoms" {
    // Build a simple component: C bonded to H
    var atoms = [_]ccd.CompAtom{
        .{ .name = .{ 'C', ' ', ' ', ' ' }, .name_len = 1, .element_symbol = .{ 'C', ' ' } },
        .{ .name = .{ 'H', '1', ' ', ' ' }, .name_len = 2, .element_symbol = .{ 'H', ' ' } },
    };
    var bonds = [_]ccd.CompBond{
        .{ .atom_idx_1 = 0, .atom_idx_2 = 1, .order = .single },
    };
    const comp = ccd.Component{
        .comp_id = "TST",
        .comp_type = "non-polymer",
        .atoms = &atoms,
        .bonds = &bonds,
    };

    // H1 already exists
    const existing = [_][4]u8{.{ 'H', '1', ' ', ' ' }};
    const plans = try derivePlans(testing.allocator, &comp, &existing);
    defer testing.allocator.free(plans);
    try testing.expectEqual(@as(usize, 0), plans.len);
}

test "analyzeBonds detects sp2" {
    var atoms = [_]ccd.CompAtom{
        .{ .name = .{ 'C', '1', ' ', ' ' }, .name_len = 2, .element_symbol = .{ 'C', ' ' } },
        .{ .name = .{ 'C', '2', ' ', ' ' }, .name_len = 2, .element_symbol = .{ 'C', ' ' } },
        .{ .name = .{ 'O', ' ', ' ', ' ' }, .name_len = 1, .element_symbol = .{ 'O', ' ' } },
        .{ .name = .{ 'H', ' ', ' ', ' ' }, .name_len = 1, .element_symbol = .{ 'H', ' ' } },
    };
    var bonds = [_]ccd.CompBond{
        .{ .atom_idx_1 = 0, .atom_idx_2 = 1, .order = .single },
        .{ .atom_idx_1 = 0, .atom_idx_2 = 2, .order = .double },
        .{ .atom_idx_1 = 0, .atom_idx_2 = 3, .order = .single },
    };
    const comp = ccd.Component{
        .comp_id = "TST",
        .comp_type = "non-polymer",
        .atoms = &atoms,
        .bonds = &bonds,
    };
    const info = analyzeBonds(&comp, 0);
    try testing.expectEqual(Hybridization.sp2, info.hybridization);
    try testing.expectEqual(@as(u8, 2), info.heavy_neighbor_count);
    try testing.expectEqual(@as(u8, 1), info.h_neighbor_count);
}

test "analyzeBonds detects sp3" {
    var atoms = [_]ccd.CompAtom{
        .{ .name = .{ 'C', ' ', ' ', ' ' }, .name_len = 1, .element_symbol = .{ 'C', ' ' } },
        .{ .name = .{ 'N', ' ', ' ', ' ' }, .name_len = 1, .element_symbol = .{ 'N', ' ' } },
        .{ .name = .{ 'H', '1', ' ', ' ' }, .name_len = 2, .element_symbol = .{ 'H', ' ' } },
        .{ .name = .{ 'H', '2', ' ', ' ' }, .name_len = 2, .element_symbol = .{ 'H', ' ' } },
        .{ .name = .{ 'H', '3', ' ', ' ' }, .name_len = 2, .element_symbol = .{ 'H', ' ' } },
    };
    var bonds = [_]ccd.CompBond{
        .{ .atom_idx_1 = 0, .atom_idx_2 = 1, .order = .single },
        .{ .atom_idx_1 = 0, .atom_idx_2 = 2, .order = .single },
        .{ .atom_idx_1 = 0, .atom_idx_2 = 3, .order = .single },
        .{ .atom_idx_1 = 0, .atom_idx_2 = 4, .order = .single },
    };
    const comp = ccd.Component{
        .comp_id = "TST",
        .comp_type = "non-polymer",
        .atoms = &atoms,
        .bonds = &bonds,
    };
    const info = analyzeBonds(&comp, 0);
    try testing.expectEqual(Hybridization.sp3, info.hybridization);
    try testing.expectEqual(@as(u8, 1), info.heavy_neighbor_count);
    try testing.expectEqual(@as(u8, 3), info.h_neighbor_count);
}

test "findBondedHeavyAtom returns correct index" {
    var atoms = [_]ccd.CompAtom{
        .{ .name = .{ 'C', ' ', ' ', ' ' }, .name_len = 1, .element_symbol = .{ 'C', ' ' } },
        .{ .name = .{ 'H', ' ', ' ', ' ' }, .name_len = 1, .element_symbol = .{ 'H', ' ' } },
    };
    var bonds = [_]ccd.CompBond{
        .{ .atom_idx_1 = 0, .atom_idx_2 = 1, .order = .single },
    };
    const comp = ccd.Component{
        .comp_id = "TST",
        .comp_type = "non-polymer",
        .atoms = &atoms,
        .bonds = &bonds,
    };
    const result = findBondedHeavyAtom(&comp, 1);
    try testing.expectEqual(@as(?u16, 0), result);
}
