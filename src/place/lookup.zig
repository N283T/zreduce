//! Atom lookup utilities for hydrogen placement.
//!
//! Provides name matching, atom position queries, bond-aware neighbor
//! lookups, and the ParentMeta helper used throughout the place/ modules.

const std = @import("std");
const model_mod = @import("../model.zig");
const Model = model_mod.Model;
const Atom = model_mod.Atom;
const Residue = model_mod.Residue;
const topology = @import("topology.zig");
const math_mod = @import("../math.zig");

const Vec3f32 = math_mod.Vec3(f32);

// ---------------------------------------------------------------------------
// Name comparison helpers
// ---------------------------------------------------------------------------

/// Compare a 4-char PDB-padded name with a model atom's nameSlice().
/// Trims leading/trailing spaces from the PDB name for comparison.
pub fn nameMatch(pdb_name: [4]u8, atom_name_slice: []const u8) bool {
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
pub fn isBlank(name: [4]u8) bool {
    return name[0] == ' ' and name[1] == ' ' and name[2] == ' ' and name[3] == ' ';
}

// ---------------------------------------------------------------------------
// Atom lookup helpers
// ---------------------------------------------------------------------------

/// Find an atom by name in the previous residue within the same chain.
/// Returns null if there is no previous residue (first in chain or chain break).
pub fn findPrevResAtomPos(mdl: *const Model, res_idx: u32, name: [4]u8, target_altloc: u8) ?Vec3f32 {
    if (res_idx == 0) return null;
    const cur_res = mdl.residues.items[res_idx];
    if (cur_res.is_chain_break_before) return null;
    const prev_res = mdl.residues.items[res_idx - 1];
    if (prev_res.chain_idx != cur_res.chain_idx) return null;
    return findAtomPos(mdl, prev_res, name, target_altloc);
}

/// Find an atom by 4-char PDB name within a residue. Returns its position.
/// When target_altloc is specified, prefers an atom with matching altloc,
/// but falls back to an atom with blank (' ') altloc if no exact match.
pub fn findAtomPos(mdl: *const Model, res: Residue, name: [4]u8, target_altloc: u8) ?Vec3f32 {
    if (isBlank(name)) return null;
    const atoms = mdl.atoms.items[res.atom_start..res.atom_end];
    var blank_match: ?Vec3f32 = null;
    for (atoms) |a| {
        if (nameMatch(name, a.nameSlice())) {
            if (a.altloc == target_altloc) return a.pos;
            if (a.altloc == ' ' and blank_match == null) blank_match = a.pos;
        }
    }
    // Note: we only search the original residue atom range. Newly appended
    // H atoms are not referenced by plans (plans only reference heavy atoms).
    return blank_match;
}

/// Find an atom by 4-char PDB name within a residue. Returns the full Atom.
/// When target_altloc is specified, prefers an atom with matching altloc,
/// but falls back to an atom with blank (' ') altloc if no exact match.
pub fn findAtom(mdl: *const Model, res: Residue, name: [4]u8, target_altloc: u8) ?Atom {
    if (isBlank(name)) return null;
    const atoms = mdl.atoms.items[res.atom_start..res.atom_end];
    var blank_match: ?Atom = null;
    for (atoms) |a| {
        if (nameMatch(name, a.nameSlice())) {
            if (a.altloc == target_altloc) return a;
            if (a.altloc == ' ' and blank_match == null) blank_match = a;
        }
    }
    return blank_match;
}

/// Check if an atom with the given name and altloc already exists in a residue.
/// Used to prevent duplicate hydrogen placement.
/// Only searches the original residue atom range, not newly appended H atoms.
/// altloc ' ' (blank) matches only blank; 'A' matches only 'A'.
pub fn existsInResidue(mdl: *const Model, res: Residue, name: [4]u8, altloc: u8) bool {
    const atoms = mdl.atoms.items[res.atom_start..res.atom_end];
    for (atoms) |a| {
        if (nameMatch(name, a.nameSlice()) and a.altloc == altloc) {
            return true;
        }
    }
    return false;
}

/// Find a heavy-atom neighbor of `center_name` that is NOT `exclude_name`.
/// Uses distance-based bonding (within 1.9 A of center).
pub fn findOtherNeighbor(mdl: *const Model, res: Residue, center_name: [4]u8, exclude_name: [4]u8, target_altloc: u8) ?Vec3f32 {
    const center_pos = findAtomPos(mdl, res, center_name, target_altloc) orelse return null;
    const bond_cutoff: f32 = 1.9;

    const atoms = mdl.atoms.items[res.atom_start..res.atom_end];
    for (atoms) |a| {
        if (a.is_hydrogen) continue;
        if (a.altloc != target_altloc and a.altloc != ' ') continue;
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
pub fn findThirdNeighbor(mdl: *const Model, res: Residue, center_name: [4]u8, n1_name: [4]u8, n2_name: [4]u8, target_altloc: u8) ?Vec3f32 {
    const center_pos = findAtomPos(mdl, res, center_name, target_altloc) orelse return null;
    const bond_cutoff: f32 = 1.9;

    const atoms = mdl.atoms.items[res.atom_start..res.atom_end];
    for (atoms) |a| {
        if (a.is_hydrogen) continue;
        if (a.altloc != target_altloc and a.altloc != ' ') continue;
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
pub fn findAtomBetween(mdl: *const Model, res: Residue, name1: [4]u8, name2: [4]u8, target_altloc: u8) ?Vec3f32 {
    const pos1 = findAtomPos(mdl, res, name1, target_altloc) orelse return null;
    const pos2 = findAtomPos(mdl, res, name2, target_altloc) orelse return null;
    const bond_cutoff: f32 = 1.9;

    const atoms = mdl.atoms.items[res.atom_start..res.atom_end];
    for (atoms) |a| {
        if (a.is_hydrogen) continue;
        if (a.altloc != target_altloc and a.altloc != ' ') continue;
        const aname = a.nameSlice();
        if (nameMatch(name1, aname)) continue;
        if (nameMatch(name2, aname)) continue;
        if (a.pos.distance(pos1) < bond_cutoff and a.pos.distance(pos2) < bond_cutoff) {
            return a.pos;
        }
    }
    return null;
}

// ---------------------------------------------------------------------------
// Bond-aware neighbor queries
// ---------------------------------------------------------------------------

/// Given a bond entry, return the partner of `atom_name`, or null if not involved.
pub fn bondPartner(bond: topology.BondEntry, atom_name: [4]u8) ?[4]u8 {
    if (std.mem.eql(u8, &bond.a1, &atom_name)) return bond.a2;
    if (std.mem.eql(u8, &bond.a2, &atom_name)) return bond.a1;
    return null;
}

/// Find a bonded neighbor of `center_name` that is NOT `exclude_name`, using topology.
/// Returns the first match; result depends on bond table ordering.
pub fn findBondedNeighbor(
    mdl: *const Model,
    res: Residue,
    bonds: []const topology.BondEntry,
    center_name: [4]u8,
    exclude_name: [4]u8,
    target_altloc: u8,
) ?Vec3f32 {
    for (bonds) |bond| {
        const partner = bondPartner(bond, center_name) orelse continue;
        if (std.mem.eql(u8, &partner, &exclude_name)) continue;
        if (findAtomPos(mdl, res, partner, target_altloc)) |pos| return pos;
    }
    return null;
}

/// Find the 3rd bonded neighbor of `center_name`, excluding `n1_name` and `n2_name`.
pub fn findThirdBondedNeighbor(
    mdl: *const Model,
    res: Residue,
    bonds: []const topology.BondEntry,
    center_name: [4]u8,
    n1_name: [4]u8,
    n2_name: [4]u8,
    target_altloc: u8,
) ?Vec3f32 {
    for (bonds) |bond| {
        const partner = bondPartner(bond, center_name) orelse continue;
        if (std.mem.eql(u8, &partner, &n1_name)) continue;
        if (std.mem.eql(u8, &partner, &n2_name)) continue;
        if (findAtomPos(mdl, res, partner, target_altloc)) |pos| return pos;
    }
    return null;
}

/// Find an atom bonded to BOTH `name1` and `name2` using topology.
pub fn findBondedAtomBetween(
    mdl: *const Model,
    res: Residue,
    bonds: []const topology.BondEntry,
    name1: [4]u8,
    name2: [4]u8,
    target_altloc: u8,
) ?Vec3f32 {
    // Collect atoms bonded to name1 (max 8 — no standard AA atom has more)
    var bonded_to_1: [8][4]u8 = undefined;
    var count_1: usize = 0;
    for (bonds) |bond| {
        const partner = bondPartner(bond, name1) orelse continue;
        if (count_1 < bonded_to_1.len) {
            bonded_to_1[count_1] = partner;
            count_1 += 1;
        }
    }

    // Check which are also bonded to name2
    for (bonds) |bond| {
        const partner = bondPartner(bond, name2) orelse continue;
        for (bonded_to_1[0..count_1]) |candidate| {
            if (std.mem.eql(u8, &partner, &candidate)) {
                if (findAtomPos(mdl, res, partner, target_altloc)) |pos| return pos;
            }
        }
    }
    return null;
}

/// Metadata inherited from the parent heavy atom.
pub const ParentMeta = struct {
    altloc: u8 = ' ',
    occupancy: f32 = 1.0,
    b_factor: f32 = 0.0,

    /// Extract metadata from an atom.
    pub fn fromAtom(a: Atom) ParentMeta {
        return .{
            .altloc = a.altloc,
            .occupancy = a.occupancy,
            .b_factor = a.b_factor,
        };
    }
};

/// Pad a short atom name (e.g. "H1") to a 4-char PDB-padded name.
/// name must be at most 4 bytes.
pub fn padName(name: []const u8) [4]u8 {
    std.debug.assert(name.len <= 4);
    var padded: [4]u8 = .{ ' ', ' ', ' ', ' ' };
    for (name, 0..) |c, i| padded[i] = c;
    return padded;
}

/// Trim leading/trailing spaces from a 4-char plan name.
pub fn trimPlanName(name: *const [4]u8) []const u8 {
    var start: usize = 0;
    while (start < 4 and name[start] == ' ') start += 1;
    var end: usize = 4;
    while (end > start and name[end - 1] == ' ') end -= 1;
    return name[start..end];
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

test "findAtom returns full atom with metadata" {
    const source = @embedFile("../test_data/tiny.cif");
    var mdl = try mmcif.parseModel(testing.allocator, source);
    defer mdl.deinit();

    const res = mdl.residues.items[0];

    // CA should be found with correct metadata
    const ca = findAtom(&mdl, res, .{ ' ', 'C', 'A', ' ' }, ' ');
    try testing.expect(ca != null);
    try testing.expectEqual(@as(f32, 1.0), ca.?.occupancy);
    try testing.expectEqual(@as(f32, 10.0), ca.?.b_factor);
    try testing.expectEqual(@as(u8, ' '), ca.?.altloc);

    // Non-existent atom returns null
    const xx = findAtom(&mdl, res, .{ ' ', 'X', 'X', ' ' }, ' ');
    try testing.expect(xx == null);
}

test "existsInResidue checks name and altloc" {
    const source = @embedFile("../test_data/tiny.cif");
    var mdl = try mmcif.parseModel(testing.allocator, source);
    defer mdl.deinit();

    const res = mdl.residues.items[0];

    // "N" with altloc ' ' exists in tiny.cif
    try testing.expect(existsInResidue(&mdl, res, .{ ' ', 'N', ' ', ' ' }, ' '));
    // "CA" with altloc ' ' exists
    try testing.expect(existsInResidue(&mdl, res, .{ ' ', 'C', 'A', ' ' }, ' '));
    // "HA" does not exist
    try testing.expect(!existsInResidue(&mdl, res, .{ ' ', 'H', 'A', ' ' }, ' '));
    // "N" with altloc 'A' does not exist (tiny.cif atoms have altloc=' ')
    try testing.expect(!existsInResidue(&mdl, res, .{ ' ', 'N', ' ', ' ' }, 'A'));
}

test "findAtomPos with altloc prefers matching conformer" {
    const source = @embedFile("../test_data/ala_altloc.cif");
    var mdl = try mmcif.parseModel(testing.allocator, source);
    defer mdl.deinit();

    const res = mdl.residues.items[0];

    // CA with altloc 'A' should return conformer A position (0, 0, 0)
    const ca_a = findAtomPos(&mdl, res, .{ ' ', 'C', 'A', ' ' }, 'A');
    try testing.expect(ca_a != null);
    try testing.expectApproxEqAbs(@as(f32, 0.0), ca_a.?.x, 1e-3);

    // CA with altloc 'B' should return conformer B position (0.1, 0.1, 0.1)
    const ca_b = findAtomPos(&mdl, res, .{ ' ', 'C', 'A', ' ' }, 'B');
    try testing.expect(ca_b != null);
    try testing.expectApproxEqAbs(@as(f32, 0.1), ca_b.?.x, 1e-3);

    // N has blank altloc — should be found by any target_altloc (fallback)
    const n_a = findAtomPos(&mdl, res, .{ ' ', 'N', ' ', ' ' }, 'A');
    try testing.expect(n_a != null);
    const n_b = findAtomPos(&mdl, res, .{ ' ', 'N', ' ', ' ' }, 'B');
    try testing.expect(n_b != null);
    try testing.expectApproxEqAbs(n_a.?.x, n_b.?.x, 1e-6);
}

test "findBondedNeighbor returns correct neighbor from topology" {
    const source = @embedFile("../test_data/tiny.cif");
    var mdl = try mmcif.parseModel(testing.allocator, source);
    defer mdl.deinit();

    const res = mdl.residues.items[0];
    const bonds = topology.getBonds("ALA").?;

    // CA is bonded to N, C, CB. Excluding N, should find C or CB.
    const result = findBondedNeighbor(&mdl, res, bonds, .{ ' ', 'C', 'A', ' ' }, .{ ' ', 'N', ' ', ' ' }, ' ');
    try testing.expect(result != null);

    // Verify it's C or CB (both bonded to CA, excluding N)
    const c_pos = findAtomPos(&mdl, res, .{ ' ', 'C', ' ', ' ' }, ' ');
    const cb_pos = findAtomPos(&mdl, res, .{ ' ', 'C', 'B', ' ' }, ' ');
    const pos = result.?;
    const is_c = c_pos != null and pos.x == c_pos.?.x and pos.y == c_pos.?.y and pos.z == c_pos.?.z;
    const is_cb = cb_pos != null and pos.x == cb_pos.?.x and pos.y == cb_pos.?.y and pos.z == cb_pos.?.z;
    try testing.expect(is_c or is_cb);
}

test "findThirdBondedNeighbor finds the third bonded atom" {
    const source = @embedFile("../test_data/tiny.cif");
    var mdl = try mmcif.parseModel(testing.allocator, source);
    defer mdl.deinit();

    const res = mdl.residues.items[0];
    const bonds = topology.getBonds("ALA").?;

    // CA bonded to N, C, CB. Excluding N and C → CB.
    const result = findThirdBondedNeighbor(&mdl, res, bonds, .{ ' ', 'C', 'A', ' ' }, .{ ' ', 'N', ' ', ' ' }, .{ ' ', 'C', ' ', ' ' }, ' ');
    try testing.expect(result != null);
    const cb_pos = findAtomPos(&mdl, res, .{ ' ', 'C', 'B', ' ' }, ' ').?;
    try testing.expectEqual(cb_pos.x, result.?.x);
    try testing.expectEqual(cb_pos.y, result.?.y);
    try testing.expectEqual(cb_pos.z, result.?.z);
}

test "findBondedAtomBetween finds atom bonded to both" {
    const source = @embedFile("../test_data/tiny.cif");
    var mdl = try mmcif.parseModel(testing.allocator, source);
    defer mdl.deinit();

    const res = mdl.residues.items[0];
    const bonds = topology.getBonds("ALA").?;

    // N and C are both bonded to CA
    const result = findBondedAtomBetween(&mdl, res, bonds, .{ ' ', 'N', ' ', ' ' }, .{ ' ', 'C', ' ', ' ' }, ' ');
    try testing.expect(result != null);
    const ca_pos = findAtomPos(&mdl, res, .{ ' ', 'C', 'A', ' ' }, ' ').?;
    try testing.expectEqual(ca_pos.x, result.?.x);
    try testing.expectEqual(ca_pos.y, result.?.y);
    try testing.expectEqual(ca_pos.z, result.?.z);
}

test "bond-based query finds stretched CB that distance-based misses" {
    const source = @embedFile("../test_data/ala_stretched.cif");
    var mdl = try mmcif.parseModel(testing.allocator, source);
    defer mdl.deinit();

    const res = mdl.residues.items[0];
    const bonds = topology.getBonds("ALA").?;

    // Distance-based: CB is >1.9A from CA, should NOT be found
    const dist_result = findThirdNeighbor(&mdl, res, .{ ' ', 'C', 'A', ' ' }, .{ ' ', 'N', ' ', ' ' }, .{ ' ', 'C', ' ', ' ' }, ' ');
    try testing.expect(dist_result == null);

    // Bond-based: CB is bonded to CA in topology, SHOULD be found
    const bond_result = findThirdBondedNeighbor(&mdl, res, bonds, .{ ' ', 'C', 'A', ' ' }, .{ ' ', 'N', ' ', ' ' }, .{ ' ', 'C', ' ', ' ' }, ' ');
    try testing.expect(bond_result != null);
}
