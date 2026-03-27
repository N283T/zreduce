//! Unified hydrogen placement entry point.
//!
//! Adds hydrogens to a Model using standard plans for known residues (20 AA)
//! and CCD-derived plans for HET groups.

const std = @import("std");
const model_mod = @import("../model.zig");
const Model = model_mod.Model;
const Atom = model_mod.Atom;
const Residue = model_mod.Residue;
const ccd_mod = @import("../ccd.zig");
const ComponentDict = ccd_mod.ComponentDict;
const standard = @import("standard.zig");
const het = @import("het.zig");
const geometry = @import("geometry.zig");
const math_mod = @import("../math.zig");
const element = @import("../element.zig");

const Vec3f32 = math_mod.Vec3(f32);

pub const PlacementResult = struct {
    n_placed: u32 = 0,
    n_skipped: u32 = 0,
    n_residues: u32 = 0,
};

/// Add hydrogens to the model.
/// Uses standard plans for known residues (20 AA), CCD-derived plans for HET groups.
/// New atoms are appended to the end of model.atoms. Each new atom carries its residue_idx.
pub fn addHydrogens(
    mdl: *Model,
    ccd_dict: ?*const ComponentDict,
) !PlacementResult {
    var result = PlacementResult{};

    const n_residues = mdl.residues.items.len;
    for (0..n_residues) |res_idx| {
        const res = mdl.residues.items[res_idx];
        const comp_id = res.compIdSlice();

        if (standard.getPlans(comp_id)) |plans| {
            for (plans) |plan| {
                if (try executePlan(mdl, res, @intCast(res_idx), &plan)) {
                    result.n_placed += 1;
                } else {
                    result.n_skipped += 1;
                }
            }
            result.n_residues += 1;
        } else if (ccd_dict) |dict| {
            if (dict.get(comp_id)) |component| {
                const existing = try collectAtomNames(mdl.allocator, mdl, res);
                defer mdl.allocator.free(existing);
                const plans = try het.derivePlans(mdl.allocator, &component, existing);
                defer mdl.allocator.free(plans);

                for (plans) |plan| {
                    if (try executePlan(mdl, res, @intCast(res_idx), &plan)) {
                        result.n_placed += 1;
                    } else {
                        result.n_skipped += 1;
                    }
                }
                result.n_residues += 1;
            }
        }
    }

    return result;
}

// ---------------------------------------------------------------------------
// Plan execution
// ---------------------------------------------------------------------------

/// Execute a single placement plan: find reference atoms, compute H position, add to model.
/// Returns true if placed, false if skipped (missing reference atoms).
fn executePlan(mdl: *Model, res: Residue, res_idx: u32, plan: *const standard.PlacementPlan) !bool {
    switch (plan.placement_type) {
        .hxr3 => {
            // connected[0]=center, connected[1..2]=two known neighbors
            // Need to find 3rd heavy-atom neighbor of center
            const center_pos = findAtomPos(mdl, res, plan.connected[0]) orelse return false;
            const n1_pos = findAtomPos(mdl, res, plan.connected[1]) orelse return false;
            const n2_pos = findAtomPos(mdl, res, plan.connected[2]) orelse return false;
            const n3_pos = findThirdNeighbor(mdl, res, plan.connected[0], plan.connected[1], plan.connected[2]) orelse return false;

            const h_pos = geometry.placeHXR3(
                center_pos.cast(f64),
                n1_pos.cast(f64),
                n2_pos.cast(f64),
                n3_pos.cast(f64),
                plan.bond_len,
            );
            try appendHydrogen(mdl, h_pos.cast(f32), plan, res_idx);
            return true;
        },

        .h2xr2 => {
            // connected[0]=center, connected[1]=reference neighbor
            const center_pos = findAtomPos(mdl, res, plan.connected[0]) orelse return false;
            const n1_pos = findAtomPos(mdl, res, plan.connected[1]) orelse return false;
            const n2_pos = findOtherNeighbor(mdl, res, plan.connected[0], plan.connected[1]) orelse return false;

            const h_pos = geometry.placeH2XR2(
                center_pos.cast(f64),
                n1_pos.cast(f64),
                n2_pos.cast(f64),
                plan.bond_len,
                plan.angle,
                plan.dihedral,
            );
            try appendHydrogen(mdl, h_pos.cast(f32), plan, res_idx);
            return true;
        },

        .h3xr => {
            // connected[0]=a1 (center), connected[1]=a2, connected[2]=a3
            const a1_pos = findAtomPos(mdl, res, plan.connected[0]) orelse return false;
            const a2_pos = findAtomPos(mdl, res, plan.connected[1]) orelse return false;
            const a3_pos = if (!isBlank(plan.connected[2]))
                findAtomPos(mdl, res, plan.connected[2]) orelse return false
            else
                findOtherNeighbor(mdl, res, plan.connected[1], plan.connected[0]) orelse return false;

            const h_pos = geometry.placeH3XR(
                a1_pos.cast(f64),
                a2_pos.cast(f64),
                a3_pos.cast(f64),
                plan.bond_len,
                plan.angle,
                plan.dihedral,
            );
            try appendHydrogen(mdl, h_pos.cast(f32), plan, res_idx);
            return true;
        },

        .hxr2_planar => {
            // connected[0] and connected[1] are neighbors; center is the atom between them
            const n1_pos = findAtomPos(mdl, res, plan.connected[0]) orelse return false;
            const n2_pos = findAtomPos(mdl, res, plan.connected[1]) orelse return false;
            const center_pos = findAtomBetween(mdl, res, plan.connected[0], plan.connected[1]) orelse return false;

            const h_pos = geometry.placeHXR2Planar(
                center_pos.cast(f64),
                n1_pos.cast(f64),
                n2_pos.cast(f64),
                plan.bond_len,
                plan.fudge,
            );
            try appendHydrogen(mdl, h_pos.cast(f32), plan, res_idx);
            return true;
        },

        .hxr2_frac => {
            const a1_pos = findAtomPos(mdl, res, plan.connected[0]) orelse return false;
            const a2_pos = findAtomPos(mdl, res, plan.connected[1]) orelse return false;
            const a3_pos = findAtomPos(mdl, res, plan.connected[2]) orelse return false;

            const h_pos = geometry.placeHXR2Frac(
                a1_pos.cast(f64),
                a2_pos.cast(f64),
                a3_pos.cast(f64),
                plan.bond_len,
                plan.fudge,
            );
            try appendHydrogen(mdl, h_pos.cast(f32), plan, res_idx);
            return true;
        },

        .hxy => {
            const center_pos = findAtomPos(mdl, res, plan.connected[0]) orelse return false;
            const neighbor_pos = findAtomPos(mdl, res, plan.connected[1]) orelse return false;

            const h_pos = geometry.placeHXY(
                center_pos.cast(f64),
                neighbor_pos.cast(f64),
                plan.bond_len,
            );
            try appendHydrogen(mdl, h_pos.cast(f32), plan, res_idx);
            return true;
        },
    }
}

// ---------------------------------------------------------------------------
// Name comparison helpers
// ---------------------------------------------------------------------------

/// Compare a 4-char PDB-padded name with a model atom's nameSlice().
/// Trims leading/trailing spaces from the PDB name for comparison.
fn nameMatch(pdb_name: [4]u8, atom_name_slice: []const u8) bool {
    // Trim leading spaces
    var start: usize = 0;
    while (start < 4 and pdb_name[start] == ' ') start += 1;
    // Trim trailing spaces
    var end: usize = 4;
    while (end > start and pdb_name[end - 1] == ' ') end -= 1;
    const trimmed_len = end - start;
    if (trimmed_len != atom_name_slice.len) return false;
    for (start..end, 0..) |i, j| {
        if (pdb_name[i] != atom_name_slice[j]) return false;
    }
    return true;
}

/// Check if a 4-char name is blank (all spaces).
fn isBlank(name: [4]u8) bool {
    return name[0] == ' ' and name[1] == ' ' and name[2] == ' ' and name[3] == ' ';
}

// ---------------------------------------------------------------------------
// Atom lookup helpers
// ---------------------------------------------------------------------------

/// Find an atom by 4-char PDB name within a residue. Returns its position.
fn findAtomPos(mdl: *const Model, res: Residue, name: [4]u8) ?Vec3f32 {
    if (isBlank(name)) return null;
    const atoms = mdl.atoms.items[res.atom_start..res.atom_end];
    for (atoms) |a| {
        if (nameMatch(name, a.nameSlice())) {
            return a.pos;
        }
    }
    // Note: we only search the original residue atom range. Newly appended
    // H atoms are not referenced by plans (plans only reference heavy atoms).
    return null;
}

/// Find a heavy-atom neighbor of `center_name` that is NOT `exclude_name`.
/// Uses distance-based bonding (within 1.9 A of center).
fn findOtherNeighbor(mdl: *const Model, res: Residue, center_name: [4]u8, exclude_name: [4]u8) ?Vec3f32 {
    const center_pos = findAtomPos(mdl, res, center_name) orelse return null;
    const bond_cutoff: f32 = 1.9;

    const atoms = mdl.atoms.items[res.atom_start..res.atom_end];
    for (atoms) |a| {
        if (a.is_hydrogen) continue;
        const aname = a.nameSlice();
        if (nameMatch(center_name, aname)) continue;
        if (nameMatch(exclude_name, aname)) continue;
        if (a.pos.distance(center_pos) < bond_cutoff) {
            return a.pos;
        }
    }
    return null;
}

/// Find the 3rd heavy-atom neighbor of center (for HXR3 placement).
/// Excludes the two already-known neighbors.
fn findThirdNeighbor(mdl: *const Model, res: Residue, center_name: [4]u8, n1_name: [4]u8, n2_name: [4]u8) ?Vec3f32 {
    const center_pos = findAtomPos(mdl, res, center_name) orelse return null;
    const bond_cutoff: f32 = 1.9;

    const atoms = mdl.atoms.items[res.atom_start..res.atom_end];
    for (atoms) |a| {
        if (a.is_hydrogen) continue;
        const aname = a.nameSlice();
        if (nameMatch(center_name, aname)) continue;
        if (nameMatch(n1_name, aname)) continue;
        if (nameMatch(n2_name, aname)) continue;
        if (a.pos.distance(center_pos) < bond_cutoff) {
            return a.pos;
        }
    }
    return null;
}

/// Find an atom that is bonded to BOTH named atoms (the atom "between" them).
/// Used for planar placement where the center atom is implicit.
fn findAtomBetween(mdl: *const Model, res: Residue, name1: [4]u8, name2: [4]u8) ?Vec3f32 {
    const pos1 = findAtomPos(mdl, res, name1) orelse return null;
    const pos2 = findAtomPos(mdl, res, name2) orelse return null;
    const bond_cutoff: f32 = 1.9;

    const atoms = mdl.atoms.items[res.atom_start..res.atom_end];
    for (atoms) |a| {
        if (a.is_hydrogen) continue;
        const aname = a.nameSlice();
        if (nameMatch(name1, aname)) continue;
        if (nameMatch(name2, aname)) continue;
        if (a.pos.distance(pos1) < bond_cutoff and a.pos.distance(pos2) < bond_cutoff) {
            return a.pos;
        }
    }
    return null;
}

/// Collect existing atom names in a residue as [4]u8 arrays.
/// Caller must free the returned slice with the provided allocator.
fn collectAtomNames(allocator: std.mem.Allocator, mdl: *const Model, res: Residue) ![][4]u8 {
    const atoms = mdl.atoms.items[res.atom_start..res.atom_end];
    const result = try allocator.alloc([4]u8, atoms.len);
    for (atoms, 0..) |a, i| {
        result[i] = a.name;
    }
    return result;
}

/// Append a new hydrogen atom to the model.
fn appendHydrogen(mdl: *Model, pos: Vec3f32, plan: *const standard.PlacementPlan, res_idx: u32) !void {
    var atom = Atom{
        .pos = pos,
        .element_type = plan.atom_type,
        .residue_idx = res_idx,
        .is_hydrogen = true,
        .is_added = true,
        .vdw_radius = plan.atom_type.info().explicit_radius,
    };
    // Set name from plan h_name (trimmed of spaces)
    var start: usize = 0;
    while (start < 4 and plan.h_name[start] == ' ') start += 1;
    var end: usize = 4;
    while (end > start and plan.h_name[end - 1] == ' ') end -= 1;
    atom.setName(plan.h_name[start..end]);
    try mdl.atoms.append(mdl.allocator, atom);
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

const testing = std.testing;
const mmcif = @import("../mmcif.zig");

test "nameMatch trims PDB-style names correctly" {
    try testing.expect(nameMatch(.{ ' ', 'N', ' ', ' ' }, "N"));
    try testing.expect(nameMatch(.{ ' ', 'C', 'A', ' ' }, "CA"));
    try testing.expect(nameMatch(.{ 'H', 'G', '1', '1' }, "HG11"));
    try testing.expect(!nameMatch(.{ ' ', 'N', ' ', ' ' }, "CA"));
}

test "isBlank detects blank names" {
    try testing.expect(isBlank(.{ ' ', ' ', ' ', ' ' }));
    try testing.expect(!isBlank(.{ 'N', ' ', ' ', ' ' }));
}

test "place hydrogens on ALA" {
    const source = @embedFile("../test_data/tiny.cif");
    var mdl = try mmcif.parseModel(testing.allocator, source);
    defer mdl.deinit();

    const initial_count = mdl.atoms.items.len;
    try testing.expectEqual(@as(usize, 5), initial_count);

    const result = try addHydrogens(&mdl, null);

    // ALA should get H atoms added
    try testing.expect(result.n_placed > 0);
    try testing.expect(mdl.atoms.items.len > initial_count);
    try testing.expectEqual(@as(u32, 1), result.n_residues);

    // Find HA atom and check bond length to CA (~1.10 A)
    const ca_pos = mdl.atoms.items[1].pos; // CA is index 1
    var found_ha = false;
    for (mdl.atoms.items) |atom| {
        if (std.mem.eql(u8, atom.nameSlice(), "HA")) {
            found_ha = true;
            const dist = atom.pos.distance(ca_pos);
            try testing.expect(dist > 0.8);
            try testing.expect(dist < 1.4);
            break;
        }
    }
    try testing.expect(found_ha);
}

test "placed atoms have correct metadata" {
    const source = @embedFile("../test_data/tiny.cif");
    var mdl = try mmcif.parseModel(testing.allocator, source);
    defer mdl.deinit();

    _ = try addHydrogens(&mdl, null);

    // Check that newly added atoms have correct flags
    for (mdl.atoms.items[5..]) |atom| {
        try testing.expect(atom.is_hydrogen);
        try testing.expect(atom.is_added);
        try testing.expectEqual(@as(u32, 0), atom.residue_idx);
    }
}

test "PlacementResult tracks counts" {
    const source = @embedFile("../test_data/tiny.cif");
    var mdl = try mmcif.parseModel(testing.allocator, source);
    defer mdl.deinit();

    const result = try addHydrogens(&mdl, null);

    // Total placed + skipped should equal number of plans for ALA (5)
    try testing.expectEqual(@as(u32, 5), result.n_placed + result.n_skipped);
}
