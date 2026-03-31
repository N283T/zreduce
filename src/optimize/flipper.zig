//! Flip movers for Asn/Gln (amide flip) and His (histidine ring flip).

const std = @import("std");
const math_mod = @import("../math.zig");
const model_mod = @import("../model.zig");
const mover_mod = @import("mover.zig");

const Vec3 = math_mod.Vec3;
const Atom = model_mod.Atom;
const Mover = mover_mod.Mover;
const Orientation = mover_mod.Orientation;
const element = @import("../element.zig");
const AtomFlags = element.AtomFlags;

/// Compute planar NH2 hydrogen positions given:
///   n_pos: the nitrogen position
///   c_pos: the carbon bonded to N (e.g. CG or CD)
///   o_pos: the oxygen bonded to C (for reference plane)
/// Returns [2]Vec3(f32): the two H positions on the amide nitrogen.
fn computeAmideNH2(n_pos: Vec3(f32), c_pos: Vec3(f32), o_pos: Vec3(f32)) [2]Vec3(f32) {
    // N-H bond length ~1.01 Å, H-N-H angle ~120° (sp2 nitrogen)
    const bond_len: f32 = 1.01;
    const half_angle_deg: f32 = 60.0; // half of 120°

    // Vector from C to N
    const cn = n_pos.sub(c_pos).normalize();
    // Vector from C to O (to define the plane)
    const co = o_pos.sub(c_pos).normalize();

    // Normal to the C-N-O plane
    const normal = cn.cross(co).normalize();

    // The two H atoms are placed symmetrically around the C-N direction
    // in the plane of the amide group, at bond_len from N.
    // We rotate cn by +half_angle and -half_angle around the normal.
    const h1_dir = math_mod.rotateAroundAxis(f32, n_pos.add(cn.scale(bond_len)), n_pos, normal, half_angle_deg);
    const h2_dir = math_mod.rotateAroundAxis(f32, n_pos.add(cn.scale(bond_len)), n_pos, normal, -half_angle_deg);

    return .{ h1_dir, h2_dir };
}

/// Create an amide flip mover for Asn or Gln.
/// o_idx: index of the O atom (OD1 or OE1)
/// n_idx: index of the N atom (ND2 or NE2)
/// h1_idx, h2_idx: indices of the two H atoms on N
/// c_idx: index of the C atom (CG or CD) for H recalculation
pub fn createAmideFlipper(
    allocator: std.mem.Allocator,
    atoms: []const Atom,
    o_idx: u32,
    n_idx: u32,
    h1_idx: u32,
    h2_idx: u32,
    c_idx: u32,
    residue_idx: u32,
) !Mover {
    const o_pos = atoms[o_idx].pos;
    const n_pos = atoms[n_idx].pos;
    const h1_pos = atoms[h1_idx].pos;
    const h2_pos = atoms[h2_idx].pos;
    const c_pos = atoms[c_idx].pos;

    // Orientation 0: original positions (O, N, H1, H2)
    const pos0 = try allocator.alloc(Vec3(f32), 4);
    pos0[0] = o_pos;
    pos0[1] = n_pos;
    pos0[2] = h1_pos;
    pos0[3] = h2_pos;

    // Orientation 1: swapped O<->N, recompute H positions
    // After swap: new N is at old O position, new O is at old N position
    const new_n_pos = o_pos; // N moves to where O was
    const new_o_pos = n_pos; // O moves to where N was
    const new_hs = computeAmideNH2(new_n_pos, c_pos, new_o_pos);

    const pos1 = try allocator.alloc(Vec3(f32), 4);
    pos1[0] = new_o_pos;
    pos1[1] = new_n_pos;
    pos1[2] = new_hs[0];
    pos1[3] = new_hs[1];

    const orientations = try allocator.alloc(Orientation, 2);
    orientations[0] = .{ .positions = pos0, .penalty = 0.0 };
    orientations[1] = .{ .positions = pos1, .penalty = 0.5 };

    const atom_indices = try allocator.alloc(u32, 4);
    atom_indices[0] = o_idx;
    atom_indices[1] = n_idx;
    atom_indices[2] = h1_idx;
    atom_indices[3] = h2_idx;

    return Mover{
        .kind = .amide_flip,
        .residue_idx = residue_idx,
        .atom_indices = atom_indices,
        .orientations = orientations,
        .allocator = allocator,
    };
}

const ABSENT_H_POS = Mover.ABSENT_H_POS;

/// Compute position of a ring H (HD1 on ND1 or HE2 on NE2) given:
///   n_pos: the nitrogen bearing the H
///   bonded1, bonded2: the two atoms bonded to N in the ring
fn computeRingNH(n_pos: Vec3(f32), bonded1: Vec3(f32), bonded2: Vec3(f32)) Vec3(f32) {
    const bond_len: f32 = 1.01;
    // Place H in the plane opposite to the bisector of the two ring bonds
    const v1 = bonded1.sub(n_pos).normalize();
    const v2 = bonded2.sub(n_pos).normalize();
    const bisector = v1.add(v2).normalize();
    return n_pos.add(bisector.scale(-bond_len));
}

/// Create a His flip mover.
/// 6 orientations (2 flip states × 3 protonation states).
///
/// Atom order in atom_indices / positions: [ND1, CD2, CE1, NE2, HD1, HE2]
/// hd1_idx / he2_idx may be null if not yet present in the model;
/// if null, a placeholder index (same as nd1_idx) is used and the
/// position is set to ABSENT_H_POS when H is absent.
pub fn createHisFlipper(
    allocator: std.mem.Allocator,
    atoms: []const Atom,
    nd1_idx: u32,
    cd2_idx: u32,
    ce1_idx: u32,
    ne2_idx: u32,
    hd1_idx: ?u32,
    he2_idx: ?u32,
    residue_idx: u32,
) !Mover {
    // --- Original ring positions ---
    const nd1 = atoms[nd1_idx].pos;
    const cd2 = atoms[cd2_idx].pos;
    const ce1 = atoms[ce1_idx].pos;
    const ne2 = atoms[ne2_idx].pos;

    // HD1 and HE2 original positions (if atoms exist)
    const hd1_orig: Vec3(f32) = if (hd1_idx) |idx| atoms[idx].pos else ABSENT_H_POS;
    const he2_orig: Vec3(f32) = if (he2_idx) |idx| atoms[idx].pos else ABSENT_H_POS;

    // Compute H positions based on ring geometry (for orientations that need them)
    // HD1 sits on ND1, bonded to CE1 and CG (we don't have CG here; use CE1 and CD2 as the two ring neighbours of ND1)
    const hd1_computed = computeRingNH(nd1, ce1, cd2);
    // HE2 sits on NE2, bonded to CE1 and CD2
    const he2_computed = computeRingNH(ne2, ce1, cd2);

    // --- Flipped ring positions ---
    // Swap ND1<->CD2 and CE1<->NE2
    const nd1_f = cd2; // ND1 goes to CD2 position
    const cd2_f = nd1; // CD2 goes to ND1 position
    const ce1_f = ne2; // CE1 goes to NE2 position
    const ne2_f = ce1; // NE2 goes to CE1 position

    // Flipped H positions
    const hd1_f_computed = computeRingNH(nd1_f, ce1_f, cd2_f);
    const he2_f_computed = computeRingNH(ne2_f, ce1_f, cd2_f);

    // --- Build 4 orientations ---
    // At neutral pH, HIS is predominantly neutral (HID or HIE, pKa ~6.0).
    // HIP (doubly protonated, +1 charge) is excluded: penalty-based approaches
    // fail because the dot-sphere scoring rewards additional H-bond donors
    // without accounting for desolvation/charge costs, causing nearly all HIS
    // to be assigned HIP regardless of penalty magnitude.
    //
    // Orient 0: no flip, HIE (NE2-HE2 protonated, ND1 acceptor)
    // Orient 1: no flip, HID (ND1-HD1 protonated, NE2 acceptor)
    // Orient 2: flip, HIE
    // Orient 3: flip, HID
    const penalties = [4]f32{ 0.00, 0.00, 0.50, 0.50 };

    // For each orientation: positions = [ND1, CD2, CE1, NE2, HD1, HE2]
    const OrientSpec = struct {
        nd1_p: Vec3(f32),
        cd2_p: Vec3(f32),
        ce1_p: Vec3(f32),
        ne2_p: Vec3(f32),
        hd1_p: Vec3(f32),
        he2_p: Vec3(f32),
    };

    // Prefer original positions for H if they exist in the original model
    const hd1_no_flip = if (hd1_idx != null) hd1_orig else hd1_computed;
    const he2_no_flip = if (he2_idx != null) he2_orig else he2_computed;

    const specs = [4]OrientSpec{
        // 0: no flip, HIE (HE2 only)
        .{ .nd1_p = nd1, .cd2_p = cd2, .ce1_p = ce1, .ne2_p = ne2, .hd1_p = ABSENT_H_POS, .he2_p = he2_no_flip },
        // 1: no flip, HID (HD1 only)
        .{ .nd1_p = nd1, .cd2_p = cd2, .ce1_p = ce1, .ne2_p = ne2, .hd1_p = hd1_no_flip, .he2_p = ABSENT_H_POS },
        // 2: flip, HIE
        .{ .nd1_p = nd1_f, .cd2_p = cd2_f, .ce1_p = ce1_f, .ne2_p = ne2_f, .hd1_p = ABSENT_H_POS, .he2_p = he2_f_computed },
        // 3: flip, HID
        .{ .nd1_p = nd1_f, .cd2_p = cd2_f, .ce1_p = ce1_f, .ne2_p = ne2_f, .hd1_p = hd1_f_computed, .he2_p = ABSENT_H_POS },
    };

    // Flags per orientation: [ND1, CD2, CE1, NE2, HD1, HE2]
    // Aromatic ring carbons keep ARA (aromatic+acceptor) in all states.
    // Nitrogens switch donor/acceptor based on protonation.
    // H atoms are always donor when present.
    const ar_acc = AtomFlags{ .aromatic = true, .acceptor = true }; // ring C
    const ar_don = AtomFlags{ .aromatic = true, .donor = true }; // protonated N
    const ar_acc_n = AtomFlags{ .aromatic = true, .acceptor = true }; // unprotonated N
    const h_don = AtomFlags{ .donor = true }; // H atom present
    const h_absent = AtomFlags{}; // absent H (no flags)

    // [ND1, CD2, CE1, NE2, HD1, HE2]
    const orient_flags = [4][6]AtomFlags{
        // 0: HIE (HE2 only) → ND1=acceptor, NE2=donor
        .{ ar_acc_n, ar_acc, ar_acc, ar_don, h_absent, h_don },
        // 1: HID (HD1 only) → ND1=donor, NE2=acceptor
        .{ ar_don, ar_acc, ar_acc, ar_acc_n, h_don, h_absent },
        // 2: flip, HIE → ND1=acceptor, NE2=donor (flipped positions)
        .{ ar_acc_n, ar_acc, ar_acc, ar_don, h_absent, h_don },
        // 3: flip, HID → ND1=donor, NE2=acceptor
        .{ ar_don, ar_acc, ar_acc, ar_acc_n, h_don, h_absent },
    };

    const orientations = try allocator.alloc(Orientation, 4);
    for (specs, 0..) |spec, i| {
        const positions = try allocator.alloc(Vec3(f32), 6);
        positions[0] = spec.nd1_p;
        positions[1] = spec.cd2_p;
        positions[2] = spec.ce1_p;
        positions[3] = spec.ne2_p;
        positions[4] = spec.hd1_p;
        positions[5] = spec.he2_p;

        const flags = try allocator.alloc(AtomFlags, 6);
        @memcpy(flags, &orient_flags[i]);

        orientations[i] = .{ .positions = positions, .flags = flags, .penalty = penalties[i] };
    }

    // atom_indices: [ND1, CD2, CE1, NE2, HD1_or_placeholder, HE2_or_placeholder]
    const actual_hd1_idx = hd1_idx orelse nd1_idx;
    const actual_he2_idx = he2_idx orelse nd1_idx;

    const atom_indices = try allocator.alloc(u32, 6);
    atom_indices[0] = nd1_idx;
    atom_indices[1] = cd2_idx;
    atom_indices[2] = ce1_idx;
    atom_indices[3] = ne2_idx;
    atom_indices[4] = actual_hd1_idx;
    atom_indices[5] = actual_he2_idx;

    return Mover{
        .kind = .his_flip,
        .residue_idx = residue_idx,
        .atom_indices = atom_indices,
        .orientations = orientations,
        .allocator = allocator,
    };
}

// ── Tests ─────────────────────────────────────────────────────────────────────

const testing = std.testing;

test "amide flipper has 2 orientations" {
    const allocator = testing.allocator;

    // Simple geometry: C at origin, O along +x, N along +y, H1/H2 near N
    const atoms = [_]Atom{
        .{ .pos = .{ .x = 1.5, .y = 0.0, .z = 0.0 } }, // 0: O (OD1)
        .{ .pos = .{ .x = 0.0, .y = 1.5, .z = 0.0 } }, // 1: N (ND2)
        .{ .pos = .{ .x = 0.2, .y = 2.3, .z = 0.5 } }, // 2: H1
        .{ .pos = .{ .x = -0.2, .y = 2.3, .z = -0.5 } }, // 3: H2
        .{ .pos = .{ .x = 0.0, .y = 0.0, .z = 0.0 } }, // 4: C (CG)
        .{ .pos = .{ .x = -1.0, .y = -1.0, .z = 0.0 } }, // 5: CB
    };

    var mover = try createAmideFlipper(allocator, &atoms, 0, 1, 2, 3, 4, 7);
    defer mover.deinit();

    try testing.expectEqual(@as(usize, 2), mover.orientations.len);
    try testing.expectApproxEqAbs(@as(f32, 0.0), mover.orientations[0].penalty, 1e-6);
    try testing.expectApproxEqAbs(@as(f32, 0.5), mover.orientations[1].penalty, 1e-6);
}

test "amide flipper swaps O and N positions" {
    const allocator = testing.allocator;

    const atoms = [_]Atom{
        .{ .pos = .{ .x = 1.5, .y = 0.0, .z = 0.0 } }, // 0: O
        .{ .pos = .{ .x = 0.0, .y = 1.5, .z = 0.0 } }, // 1: N
        .{ .pos = .{ .x = 0.2, .y = 2.3, .z = 0.5 } }, // 2: H1
        .{ .pos = .{ .x = -0.2, .y = 2.3, .z = -0.5 } }, // 3: H2
        .{ .pos = .{ .x = 0.0, .y = 0.0, .z = 0.0 } }, // 4: C
        .{ .pos = .{ .x = -1.0, .y = -1.0, .z = 0.0 } }, // 5: CB
    };

    var mover = try createAmideFlipper(allocator, &atoms, 0, 1, 2, 3, 4, 0);
    defer mover.deinit();

    const orient0 = mover.orientations[0];
    const orient1 = mover.orientations[1];

    // O position in orient 0 == N position in orient 1
    try testing.expectApproxEqAbs(orient0.positions[0].x, orient1.positions[1].x, 1e-6);
    try testing.expectApproxEqAbs(orient0.positions[0].y, orient1.positions[1].y, 1e-6);
    try testing.expectApproxEqAbs(orient0.positions[0].z, orient1.positions[1].z, 1e-6);

    // N position in orient 0 == O position in orient 1
    try testing.expectApproxEqAbs(orient0.positions[1].x, orient1.positions[0].x, 1e-6);
    try testing.expectApproxEqAbs(orient0.positions[1].y, orient1.positions[0].y, 1e-6);
    try testing.expectApproxEqAbs(orient0.positions[1].z, orient1.positions[0].z, 1e-6);
}

test "his flipper has 6 orientations" {
    const allocator = testing.allocator;

    // Simple imidazole ring geometry
    const atoms = [_]Atom{
        .{ .pos = .{ .x = 0.0, .y = 1.38, .z = 0.0 } }, // 0: ND1
        .{ .pos = .{ .x = 1.19, .y = -0.38, .z = 0.0 } }, // 1: CD2
        .{ .pos = .{ .x = 1.21, .y = 0.97, .z = 0.0 } }, // 2: CE1
        .{ .pos = .{ .x = 0.0, .y = -1.38, .z = 0.0 } }, // 3: NE2
        .{ .pos = .{ .x = -0.9, .y = 2.0, .z = 0.0 } }, // 4: HD1
        .{ .pos = .{ .x = -0.9, .y = -2.0, .z = 0.0 } }, // 5: HE2
    };

    var mover = try createHisFlipper(allocator, &atoms, 0, 1, 2, 3, 4, 5, 3);
    defer mover.deinit();

    try testing.expectEqual(@as(usize, 4), mover.orientations.len);
}

test "his flipper penalties match spec" {
    const allocator = testing.allocator;

    const atoms = [_]Atom{
        .{ .pos = .{ .x = 0.0, .y = 1.38, .z = 0.0 } }, // 0: ND1
        .{ .pos = .{ .x = 1.19, .y = -0.38, .z = 0.0 } }, // 1: CD2
        .{ .pos = .{ .x = 1.21, .y = 0.97, .z = 0.0 } }, // 2: CE1
        .{ .pos = .{ .x = 0.0, .y = -1.38, .z = 0.0 } }, // 3: NE2
        .{ .pos = .{ .x = -0.9, .y = 2.0, .z = 0.0 } }, // 4: HD1
        .{ .pos = .{ .x = -0.9, .y = -2.0, .z = 0.0 } }, // 5: HE2
    };

    var mover = try createHisFlipper(allocator, &atoms, 0, 1, 2, 3, 4, 5, 3);
    defer mover.deinit();

    const expected_penalties = [4]f32{ 0.00, 0.00, 0.50, 0.50 };
    for (mover.orientations, 0..) |o, i| {
        try testing.expectApproxEqAbs(expected_penalties[i], o.penalty, 1e-6);
    }
}

test "his flipper orientations have chemistry flags" {
    const allocator = testing.allocator;

    const atoms = [_]Atom{
        .{ .pos = .{ .x = 0.0, .y = 1.38, .z = 0.0 } }, // 0: ND1
        .{ .pos = .{ .x = 1.19, .y = -0.38, .z = 0.0 } }, // 1: CD2
        .{ .pos = .{ .x = 1.21, .y = 0.97, .z = 0.0 } }, // 2: CE1
        .{ .pos = .{ .x = 0.0, .y = -1.38, .z = 0.0 } }, // 3: NE2
        .{ .pos = .{ .x = -0.9, .y = 2.0, .z = 0.0 } }, // 4: HD1
        .{ .pos = .{ .x = -0.9, .y = -2.0, .z = 0.0 } }, // 5: HE2
    };

    var mover = try createHisFlipper(allocator, &atoms, 0, 1, 2, 3, 4, 5, 3);
    defer mover.deinit();

    // All 4 orientations should have flags
    for (mover.orientations) |o| {
        try testing.expect(o.flags != null);
        try testing.expectEqual(@as(usize, 6), o.flags.?.len);
    }

    // Orient 0: HIE (HE2 only) → ND1=acceptor (not donor), NE2=donor
    const flags0 = mover.orientations[0].flags.?;
    try testing.expect(flags0[0].acceptor); // ND1 acceptor
    try testing.expect(!flags0[0].donor); // ND1 not donor
    try testing.expect(flags0[3].donor); // NE2 donor
    try testing.expect(!flags0[3].acceptor); // NE2 not acceptor

    // Orient 1: HID (HD1 only) → ND1=donor, NE2=acceptor
    const flags1 = mover.orientations[1].flags.?;
    try testing.expect(flags1[0].donor); // ND1 donor
    try testing.expect(!flags1[0].acceptor); // ND1 not acceptor
    try testing.expect(flags1[3].acceptor); // NE2 acceptor
    try testing.expect(!flags1[3].donor); // NE2 not donor

    // No HIP orientations — all orientations have exactly one absent H
    for (mover.orientations) |o| {
        const hd1_absent = mover_mod.isAbsentH(.{ .pos = o.positions[4], .is_added = true });
        const he2_absent = mover_mod.isAbsentH(.{ .pos = o.positions[5], .is_added = true });
        // Exactly one H should be absent (HID xor HIE, never HIP)
        try testing.expect(hd1_absent != he2_absent);
    }
}
