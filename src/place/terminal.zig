//! N-terminal and 3'-terminal hydrogen placement.
//!
//! Handles NH3+/NH2+ placement on peptide N-termini and HO3' on
//! nucleotide 3'-termini.

const std = @import("std");

/// Controls how the N-terminal backbone amine is protonated.
///
/// - `auto` (default): Place NH3+ (or NH2+ on PRO) only at real chain-first
///   residues. Residues after an observed chain break keep a single backbone
///   amide H. Matches ChimeraX addh and the behavior documented in #118.
/// - `aggressive`: Place NH3+/NH2+ on both real chain-first residues and
///   residues following a chain break. Matches reduce2's `first_in_chain`
///   mode and likely explains the additional NH3 placed by reduce2 on
///   AlphaFold models with chain breaks (#118 lists the reduce2 default as
///   an open question; see #114).
/// - `neutral`: Place a neutral NH2 (two H, no positive charge flag) at
///   non-PRO real N-termini. PRO is an exception: its backbone N is a
///   secondary amine (bonded to CA and CD), which under physiological
///   conditions is protonated to NH2+. PRO therefore keeps the default
///   NH2+ placement and the associated positive charge flag even in
///   `neutral` mode — making a full neutral secondary amine would require
///   a single-H placement path that does not exist today.
///
/// Residues after an observed chain break are **only** treated as N-termini
/// under `aggressive`; `auto` and `neutral` leave them with a single break
/// amide H regardless of mode.
pub const NtermMode = enum {
    auto,
    aggressive,
    neutral,

    /// Parse a mode string as used on the command line. Returns null on
    /// unrecognized input so the caller can emit a context-aware error.
    pub fn fromString(s: []const u8) ?NtermMode {
        if (std.ascii.eqlIgnoreCase(s, "auto")) return .auto;
        if (std.ascii.eqlIgnoreCase(s, "aggressive")) return .aggressive;
        if (std.ascii.eqlIgnoreCase(s, "neutral")) return .neutral;
        return null;
    }

    /// True iff residues after an observed chain break should be treated as
    /// N-termini (skip the single amide H and place NH3+/NH2+ instead).
    pub fn treatsBreakAsNterm(self: NtermMode) bool {
        return self == .aggressive;
    }

    /// True iff the backbone N of an N-terminal residue should carry the
    /// formal positive charge flag. The `is_pro` parameter captures the PRO
    /// exception: `neutral` mode still lets PRO be protonated because the
    /// placer keeps the NH2+ geometry for PRO (see dispatch in `placer.zig`).
    pub fn placesPositiveCharge(self: NtermMode, is_pro: bool) bool {
        return switch (self) {
            .auto, .aggressive => true,
            .neutral => is_pro,
        };
    }
};

const model_mod = @import("../model.zig");
const Model = model_mod.Model;
const Atom = model_mod.Atom;
const Residue = model_mod.Residue;
const math_mod = @import("../math.zig");
const element = @import("../element.zig");
const bond_policy = @import("bond_policy.zig");
const geometry = @import("geometry.zig");
const standard = @import("standard.zig");
const lookup = @import("lookup.zig");

const Vec3f32 = math_mod.Vec3(f32);
const ParentMeta = lookup.ParentMeta;
const findAtom = lookup.findAtom;
const findAtomPos = lookup.findAtomPos;
const existsInResidue = lookup.existsInResidue;
const padName = lookup.padName;
const nameMatch = lookup.nameMatch;

/// Check if a placement plan is for a backbone amide H (single NH or N-terminal H1/H2/H3).
/// Used by N-terminal skip logic — matches " H  ", " H1 ", " H2 ", " H3 " but only
/// when the parent atom (connected[0]) is backbone nitrogen " N  ".
/// This prevents false matches on nucleotide ring H (e.g. guanine H1 on N1 via C6).
pub fn isBackboneH(plan: *const standard.PlacementPlan) bool {
    const parent = plan.connected[0];
    if (!(parent[0] == ' ' and parent[1] == 'N' and parent[2] == ' ' and parent[3] == ' ')) return false;
    const h = plan.h_name;
    if (h[0] == ' ' and h[1] == 'H' and h[2] == ' ' and h[3] == ' ') return true;
    if (h[0] == ' ' and h[1] == 'H' and h[2] == '1' and h[3] == ' ') return true;
    if (h[0] == ' ' and h[1] == 'H' and h[2] == '2' and h[3] == ' ') return true;
    if (h[0] == ' ' and h[1] == 'H' and h[2] == '3' and h[3] == ' ') return true;
    return false;
}

/// Check if a plan is specifically the single backbone amide H (" H  ").
/// Used by peptide-plane placement — must NOT match N-terminal H1/H2/H3.
pub fn isBackboneAmideH(plan: *const standard.PlacementPlan) bool {
    const h = plan.h_name;
    return h[0] == ' ' and h[1] == 'H' and h[2] == ' ' and h[3] == ' ';
}

/// Append an N-terminal H atom to the model.
pub fn appendNtermH(mdl: *Model, h_pos: Vec3f32, name: []const u8, res_idx: u32, meta: ParentMeta, mover_hint: standard.MoverHint) !void {
    const hpol_info = element.AtomType.Hpol.info();
    var atom = Atom{
        .pos = h_pos,
        .element_type = .Hpol,
        .residue_idx = res_idx,
        .is_hydrogen = true,
        .is_added = true,
        .vdw_radius = hpol_info.explicit_radius,
        .altloc = meta.altloc,
        .occupancy = meta.occupancy,
        .b_factor = meta.b_factor,
        .flags = .{ .donor = true },
        .mover_hint = mover_hint,
    };
    atom.setName(name);
    try mdl.atoms.append(mdl.allocator, atom);
}

/// Result of a terminal placement attempt.
///
/// - `placed`: number of H atoms actually appended to the model
/// - `skipped`: number of plans skipped because an equivalent H already exists
/// - `missing_ref`: number of placements aborted because a required reference
///   atom (N, CA, C, CD, O3', C3', C4') could not be resolved at the target
///   altloc. Callers should tally this into `PlacementResult.n_skipped_missing_ref`
///   so the final run summary surfaces the shortfall instead of silently
///   dropping it (see #254).
pub const NtermResult = struct {
    placed: u32,
    skipped: u32,
    missing_ref: u32 = 0,
};

/// Place NH3+ hydrogens (H1, H2, H3) on the N-terminal residue.
/// Uses h3xr (dihedral-controlled) placement around the N-CA bond.
pub fn placeNtermNH3(mdl: *Model, res: Residue, res_idx: u32, target_altloc: u8, mode: bond_policy.BondLengthMode) !NtermResult {
    const n_atom = findAtom(mdl, res, .{ ' ', 'N', ' ', ' ' }, target_altloc) orelse return .{ .placed = 0, .skipped = 0, .missing_ref = 1 };
    const ca_pos = findAtomPos(mdl, res, .{ ' ', 'C', 'A', ' ' }, target_altloc) orelse return .{ .placed = 0, .skipped = 0, .missing_ref = 1 };
    const c_pos = findAtomPos(mdl, res, .{ ' ', 'C', ' ', ' ' }, target_altloc) orelse return .{ .placed = 0, .skipped = 0, .missing_ref = 1 };

    var meta = ParentMeta.fromAtom(n_atom);
    if (target_altloc != ' ') meta.altloc = target_altloc;

    const n64 = n_atom.pos.cast(f64);
    const ca64 = math_mod.Vec3(f64){ .x = ca_pos.x, .y = ca_pos.y, .z = ca_pos.z };
    const c64 = math_mod.Vec3(f64){ .x = c_pos.x, .y = c_pos.y, .z = c_pos.z };

    const bond_len: f64 = bond_policy.adjustedBondLength(mode, 1.00, n_atom.element_type, .Hpol);
    const angle_deg: f64 = 109.5;
    const dihedrals = [3]f64{ 180.0, 60.0, -60.0 };
    const names = [3][]const u8{ "H1", "H2", "H3" };

    var placed: u32 = 0;
    var skipped: u32 = 0;
    for (names, dihedrals) |name, dihedral| {
        if (existsInResidue(mdl, res, padName(name), meta.altloc)) {
            skipped += 1;
            continue;
        }

        const h64 = geometry.placeH3XR(n64, ca64, c64, bond_len, angle_deg, dihedral);
        const h_pos = Vec3f32{ .x = @floatCast(h64.x), .y = @floatCast(h64.y), .z = @floatCast(h64.z) };
        try appendNtermH(mdl, h_pos, name, res_idx, meta, .rotate_nh3);
        placed += 1;
    }
    return .{ .placed = placed, .skipped = skipped };
}

/// Place NH2+ hydrogens (H2, H3) on N-terminal PRO.
/// PRO has CD bonded to N (secondary amine), so only 2 H positions.
pub fn placeNtermNH2Pro(mdl: *Model, res: Residue, res_idx: u32, target_altloc: u8, mode: bond_policy.BondLengthMode) !NtermResult {
    const n_atom = findAtom(mdl, res, .{ ' ', 'N', ' ', ' ' }, target_altloc) orelse return .{ .placed = 0, .skipped = 0, .missing_ref = 1 };
    const ca_pos = findAtomPos(mdl, res, .{ ' ', 'C', 'A', ' ' }, target_altloc) orelse return .{ .placed = 0, .skipped = 0, .missing_ref = 1 };
    const cd_pos = findAtomPos(mdl, res, .{ ' ', 'C', 'D', ' ' }, target_altloc) orelse return .{ .placed = 0, .skipped = 0, .missing_ref = 1 };

    var meta = ParentMeta.fromAtom(n_atom);
    if (target_altloc != ' ') meta.altloc = target_altloc;

    // PRO N is sp3 with 2 heavy-atom neighbors (CA, CD) and 2 H.
    // Use h2xr2 (two H on atom with 2 heavy neighbors).
    const n64 = n_atom.pos.cast(f64);
    const ca64 = math_mod.Vec3(f64){ .x = ca_pos.x, .y = ca_pos.y, .z = ca_pos.z };
    const cd64 = math_mod.Vec3(f64){ .x = cd_pos.x, .y = cd_pos.y, .z = cd_pos.z };

    const bond_len: f64 = bond_policy.adjustedBondLength(mode, 1.00, n_atom.element_type, .Hpol);
    const angle_deg: f64 = 109.5;
    const names = [2][]const u8{ "H2", "H3" };
    const dihedrals = [2]f64{ 120.0, -120.0 };

    var placed: u32 = 0;
    var skipped: u32 = 0;
    for (names, dihedrals) |name, dihedral| {
        if (existsInResidue(mdl, res, padName(name), meta.altloc)) {
            skipped += 1;
            continue;
        }

        const h64 = geometry.placeH2XR2(n64, ca64, cd64, bond_len, angle_deg, dihedral);
        const h_pos = Vec3f32{ .x = @floatCast(h64.x), .y = @floatCast(h64.y), .z = @floatCast(h64.z) };
        try appendNtermH(mdl, h_pos, name, res_idx, meta, .rotate);
        placed += 1;
    }
    return .{ .placed = placed, .skipped = skipped };
}

/// Place neutral NH2 hydrogens (H2, H3) on a non-PRO N-terminal residue.
/// Same h3xr geometry as NH3+ (tetrahedral 109.5° angle, reused for code
/// simplicity even though a real primary amine is slightly different) but
/// skips the anti position. The resulting C-CA-N-H dihedrals are +60° and
/// -60° (both gauche to the backbone carbonyl C), leaving the 180° slot
/// empty for the nitrogen lone pair. The backbone N is treated as a
/// primary amine with one heavy neighbor, two H, and one lone pair.
pub fn placeNtermNH2Neutral(mdl: *Model, res: Residue, res_idx: u32, target_altloc: u8, mode: bond_policy.BondLengthMode) !NtermResult {
    const n_atom = findAtom(mdl, res, .{ ' ', 'N', ' ', ' ' }, target_altloc) orelse return .{ .placed = 0, .skipped = 0, .missing_ref = 1 };
    const ca_pos = findAtomPos(mdl, res, .{ ' ', 'C', 'A', ' ' }, target_altloc) orelse return .{ .placed = 0, .skipped = 0, .missing_ref = 1 };
    const c_pos = findAtomPos(mdl, res, .{ ' ', 'C', ' ', ' ' }, target_altloc) orelse return .{ .placed = 0, .skipped = 0, .missing_ref = 1 };

    var meta = ParentMeta.fromAtom(n_atom);
    if (target_altloc != ' ') meta.altloc = target_altloc;

    const n64 = n_atom.pos.cast(f64);
    const ca64 = math_mod.Vec3(f64){ .x = ca_pos.x, .y = ca_pos.y, .z = ca_pos.z };
    const c64 = math_mod.Vec3(f64){ .x = c_pos.x, .y = c_pos.y, .z = c_pos.z };

    const bond_len: f64 = bond_policy.adjustedBondLength(mode, 1.00, n_atom.element_type, .Hpol);
    const angle_deg: f64 = 109.5;
    // H2/H3 are placed at the two gauche C-CA-N-H dihedrals (+60°, -60°),
    // skipping the anti (180°) slot so the N lone pair sits trans to C'.
    // Same H naming convention as `placeNtermNH2Pro` for internal consistency.
    const names = [2][]const u8{ "H2", "H3" };
    const dihedrals = [2]f64{ 60.0, -60.0 };

    var placed: u32 = 0;
    var skipped: u32 = 0;
    for (names, dihedrals) |name, dihedral| {
        if (existsInResidue(mdl, res, padName(name), meta.altloc)) {
            skipped += 1;
            continue;
        }

        const h64 = geometry.placeH3XR(n64, ca64, c64, bond_len, angle_deg, dihedral);
        const h_pos = Vec3f32{ .x = @floatCast(h64.x), .y = @floatCast(h64.y), .z = @floatCast(h64.z) };
        // Use .rotate so the neutral NH2 rotates freely around the N-CA axis
        // (no dedicated neutral-amine mover exists; .rotate_nh3 is for NH3+).
        try appendNtermH(mdl, h_pos, name, res_idx, meta, .rotate);
        placed += 1;
    }
    return .{ .placed = placed, .skipped = skipped };
}

/// Place HO3' on 3' terminal nucleotide residues.
/// O3' is sp3 with 2 heavy-atom neighbors (C3', and normally the next P — absent at 3' terminus).
/// Uses h3xr geometry: O3' center, C3' and C4' as references, tetrahedral angle.
pub fn place3primeOH(mdl: *Model, res: Residue, res_idx: u32, target_altloc: u8, mode: bond_policy.BondLengthMode) !NtermResult {
    const o3_atom = findAtom(mdl, res, .{ ' ', 'O', '3', '\'' }, target_altloc) orelse return .{ .placed = 0, .skipped = 0, .missing_ref = 1 };
    const c3_pos = findAtomPos(mdl, res, .{ ' ', 'C', '3', '\'' }, target_altloc) orelse return .{ .placed = 0, .skipped = 0, .missing_ref = 1 };
    const c4_pos = findAtomPos(mdl, res, .{ ' ', 'C', '4', '\'' }, target_altloc) orelse return .{ .placed = 0, .skipped = 0, .missing_ref = 1 };

    var meta = ParentMeta.fromAtom(o3_atom);
    if (target_altloc != ' ') meta.altloc = target_altloc;

    const ho3_name = padName("HO3'");
    if (existsInResidue(mdl, res, ho3_name, meta.altloc)) {
        return .{ .placed = 0, .skipped = 1 };
    }

    const o3_64 = o3_atom.pos.cast(f64);
    const c3_64 = math_mod.Vec3(f64){ .x = c3_pos.x, .y = c3_pos.y, .z = c3_pos.z };
    const c4_64 = math_mod.Vec3(f64){ .x = c4_pos.x, .y = c4_pos.y, .z = c4_pos.z };

    const bond_len: f64 = bond_policy.adjustedBondLength(mode, 0.97, o3_atom.element_type, .Hpol);
    const h64 = geometry.placeH3XR(o3_64, c3_64, c4_64, bond_len, 109.5, 180.0);
    const h_pos = Vec3f32{ .x = @floatCast(h64.x), .y = @floatCast(h64.y), .z = @floatCast(h64.z) };

    const hpol_info = element.AtomType.Hpol.info();
    var atom = Atom{
        .pos = h_pos,
        .element_type = .Hpol,
        .residue_idx = res_idx,
        .is_hydrogen = true,
        .is_added = true,
        .vdw_radius = hpol_info.explicit_radius,
        .altloc = meta.altloc,
        .occupancy = meta.occupancy,
        .b_factor = meta.b_factor,
        .flags = .{ .donor = true },
        .mover_hint = .rotate,
    };
    atom.setName("HO3'");
    try mdl.atoms.append(mdl.allocator, atom);
    return .{ .placed = 1, .skipped = 0 };
}

/// Check if a residue has a nucleotide sugar-phosphate backbone (has O3' atom).
pub fn isNucleotideResidue(mdl: *const Model, res: Residue) bool {
    const atoms = mdl.atoms.items[res.atom_start..res.atom_end];
    for (atoms) |a| {
        if (nameMatch(.{ ' ', 'O', '3', '\'' }, a.nameSlice())) return true;
    }
    return false;
}

/// Check if a hydrogen name is a phosphate H (HOP2, HOP3).
/// In CCD, these are H atoms on phosphate oxygens OP2/OP3.
pub fn isPhosphateH(h_name: [4]u8) bool {
    return (h_name[0] == 'H' and h_name[1] == 'O' and h_name[2] == 'P');
}

// ── Tests ─────────────────────────────────────────────────────────────────────

const testing = @import("std").testing;
const mmcif = @import("../mmcif.zig");

test "placeNtermNH3 sets rotate_nh3 mover_hint on all 3 H atoms" {
    const source = @embedFile("../test_data/tiny.cif");
    var mdl = try mmcif.parseModel(testing.allocator, source);
    defer mdl.deinit();

    const res = mdl.residues.items[0];
    // tiny.cif is ALA N-terminal — needs N, CA, C atoms
    _ = try placeNtermNH3(&mdl, res, 0, ' ', .neutron);

    var nh3_count: u32 = 0;
    for (mdl.atoms.items) |atom| {
        if (!atom.is_added) continue;
        const name = atom.nameSlice();
        const is_nterm = std.mem.eql(u8, name, "H1") or
            std.mem.eql(u8, name, "H2") or
            std.mem.eql(u8, name, "H3");
        if (!is_nterm) continue;
        try testing.expectEqual(standard.MoverHint.rotate_nh3, atom.mover_hint);
        nh3_count += 1;
    }
    try testing.expectEqual(@as(u32, 3), nh3_count);
}

test "NtermMode.fromString parses canonical names case-insensitively" {
    try testing.expectEqual(@as(?NtermMode, .auto), NtermMode.fromString("auto"));
    try testing.expectEqual(@as(?NtermMode, .auto), NtermMode.fromString("AUTO"));
    try testing.expectEqual(@as(?NtermMode, .aggressive), NtermMode.fromString("aggressive"));
    try testing.expectEqual(@as(?NtermMode, .aggressive), NtermMode.fromString("Aggressive"));
    try testing.expectEqual(@as(?NtermMode, .neutral), NtermMode.fromString("neutral"));
    try testing.expectEqual(@as(?NtermMode, .neutral), NtermMode.fromString("NEUTRAL"));
    try testing.expectEqual(@as(?NtermMode, null), NtermMode.fromString("none"));
    try testing.expectEqual(@as(?NtermMode, null), NtermMode.fromString(""));
    try testing.expectEqual(@as(?NtermMode, null), NtermMode.fromString("garbage"));
}

test "placeNtermNH2Neutral places exactly 2 H atoms with gauche dihedrals" {
    // Use ala_stretched.cif (non-collinear N-CA-C backbone); tiny.cif has
    // a collinear backbone that degenerates h3xr dihedral placement.
    const source = @embedFile("../test_data/ala_stretched.cif");
    var mdl = try mmcif.parseModel(testing.allocator, source);
    defer mdl.deinit();

    const res = mdl.residues.items[0];
    const result = try placeNtermNH2Neutral(&mdl, res, 0, ' ', .neutron);
    try testing.expectEqual(@as(u32, 2), result.placed);
    try testing.expectEqual(@as(u32, 0), result.skipped);

    // Collect heavy-atom reference positions for geometric validation.
    var n_pos: ?Vec3f32 = null;
    var ca_pos: ?Vec3f32 = null;
    var c_pos: ?Vec3f32 = null;
    var h2_pos: ?Vec3f32 = null;
    var h3_pos: ?Vec3f32 = null;
    for (mdl.atoms.items) |atom| {
        const name = atom.nameSlice();
        if (std.mem.eql(u8, name, "N")) n_pos = atom.pos;
        if (std.mem.eql(u8, name, "CA")) ca_pos = atom.pos;
        if (std.mem.eql(u8, name, "C")) c_pos = atom.pos;
        if (atom.is_added and std.mem.eql(u8, name, "H2")) h2_pos = atom.pos;
        if (atom.is_added and std.mem.eql(u8, name, "H3")) h3_pos = atom.pos;
    }
    try testing.expect(n_pos != null and ca_pos != null and c_pos != null);
    try testing.expect(h2_pos != null and h3_pos != null);

    // Both H must be placed, non-NaN, and distinct from each other.
    try testing.expect(!std.math.isNan(h2_pos.?.x));
    try testing.expect(!std.math.isNan(h3_pos.?.x));
    try testing.expect(h2_pos.?.distance(h3_pos.?) > 0.1);

    // N-H bond length ≈ 1.00 Å (neutron mode, no isotope scaling).
    const nh2_len = n_pos.?.distance(h2_pos.?);
    const nh3_len = n_pos.?.distance(h3_pos.?);
    try testing.expectApproxEqAbs(@as(f32, 1.0), nh2_len, 0.02);
    try testing.expectApproxEqAbs(@as(f32, 1.0), nh3_len, 0.02);

    // H-N-CA angle ≈ 109.5° (tetrahedral).
    const n = n_pos.?;
    const ca = ca_pos.?;
    const nh2_vec = h2_pos.?.sub(n).normalize();
    const nh3_vec = h3_pos.?.sub(n).normalize();
    const nca_vec = ca.sub(n).normalize();
    const cos_h2 = nh2_vec.dot(nca_vec);
    const cos_h3 = nh3_vec.dot(nca_vec);
    const angle_h2 = std.math.radiansToDegrees(std.math.acos(std.math.clamp(cos_h2, -1.0, 1.0)));
    const angle_h3 = std.math.radiansToDegrees(std.math.acos(std.math.clamp(cos_h3, -1.0, 1.0)));
    try testing.expectApproxEqAbs(@as(f32, 109.5), angle_h2, 1.0);
    try testing.expectApproxEqAbs(@as(f32, 109.5), angle_h3, 1.0);

    var h_count: u32 = 0;
    for (mdl.atoms.items) |atom| {
        if (!atom.is_added) continue;
        const name = atom.nameSlice();
        const is_neutral_nh2 = std.mem.eql(u8, name, "H2") or std.mem.eql(u8, name, "H3");
        if (!is_neutral_nh2) continue;
        // Neutral NH2 rotates but is not an NH3+ rotator group.
        try testing.expectEqual(standard.MoverHint.rotate, atom.mover_hint);
        h_count += 1;
    }
    try testing.expectEqual(@as(u32, 2), h_count);

    // No H1 should be placed in neutral mode (skipping the anti dihedral slot).
    for (mdl.atoms.items) |atom| {
        if (!atom.is_added) continue;
        try testing.expect(!std.mem.eql(u8, atom.nameSlice(), "H1"));
    }
}

test "placeNtermNH2Pro sets rotate mover_hint on both H atoms" {
    // Build a minimal PRO-like model with N, CA, CD atoms
    var mdl = try mmcif.parseModel(testing.allocator,
        \\data_PRO_TEST
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
        \\_atom_site.auth_asym_id
        \\ATOM 1 N N   PRO A 1 1.000 2.000 3.000 1.00 10.0 . A
        \\ATOM 2 C CA  PRO A 1 2.000 3.000 4.000 1.00 10.0 . A
        \\ATOM 3 C C   PRO A 1 3.000 4.000 5.000 1.00 10.0 . A
        \\ATOM 4 O O   PRO A 1 4.000 5.000 6.000 1.00 10.0 . A
        \\ATOM 5 C CD  PRO A 1 0.500 3.000 3.500 1.00 10.0 . A
        \\#
    );
    defer mdl.deinit();

    const res = mdl.residues.items[0];
    _ = try placeNtermNH2Pro(&mdl, res, 0, ' ', .neutron);

    var h_count: u32 = 0;
    for (mdl.atoms.items) |atom| {
        if (!atom.is_added) continue;
        const name = atom.nameSlice();
        const is_pro_nterm = std.mem.eql(u8, name, "H2") or std.mem.eql(u8, name, "H3");
        if (!is_pro_nterm) continue;
        try testing.expectEqual(standard.MoverHint.rotate, atom.mover_hint);
        h_count += 1;
    }
    try testing.expectEqual(@as(u32, 2), h_count);
}

/// Check if a residue has atom named "P" (phosphorus — nucleotide backbone).
/// Matches both setName("P") format {'P',' ',' ',' '} and PDB-padded {' ','P',' ',' '}.
pub fn hasPhosphorusAtom(mdl: *const Model, res: Residue) bool {
    const atoms = mdl.atoms.items[res.atom_start..res.atom_end];
    for (atoms) |a| {
        if (nameMatch(.{ ' ', 'P', ' ', ' ' }, a.nameSlice())) return true;
    }
    return false;
}
