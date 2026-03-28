//! Rotator movers: OH/SH (single H), NH3+, and methyl (CH3) rotators.
//! Each creates a Mover with discrete orientations by rotating hydrogens
//! around a bond axis.

const std = @import("std");
const math_mod = @import("../math.zig");
const model_mod = @import("../model.zig");
const mover_mod = @import("mover.zig");

const Vec3 = math_mod.Vec3;
const Atom = model_mod.Atom;
const Mover = mover_mod.Mover;
const Orientation = mover_mod.Orientation;

/// Create a SingleH rotator mover (OH, SH, SeH).
/// h_idx: index of the hydrogen atom to rotate
/// center_idx: index of the atom the H is bonded to (e.g., O in O-H)
/// axis_idx: index of the atom defining the rotation axis (e.g., C in C-O-H)
/// Produces 12 orientations at 30° intervals. No penalty for any orientation.
pub fn createSingleHRotator(
    allocator: std.mem.Allocator,
    atoms: []const Atom,
    h_idx: u32,
    center_idx: u32,
    axis_idx: u32,
    residue_idx: u32,
) !Mover {
    const center = atoms[center_idx].pos;
    const axis_atom = atoms[axis_idx].pos;
    const axis_dir = center.sub(axis_atom);
    const h_pos = atoms[h_idx].pos;

    const orientations = try allocator.alloc(Orientation, 12);
    for (0..12) |i| {
        const angle: f32 = @as(f32, @floatFromInt(i)) * 30.0;
        const rotated = math_mod.rotateAroundAxis(f32, h_pos, center, axis_dir, angle);
        const positions = try allocator.alloc(Vec3(f32), 1);
        positions[0] = rotated;
        orientations[i] = .{ .positions = positions, .penalty = 0.0 };
    }

    const atom_indices = try allocator.alloc(u32, 1);
    atom_indices[0] = h_idx;

    return Mover{
        .kind = .single_h_rotator,
        .residue_idx = residue_idx,
        .atom_indices = atom_indices,
        .orientations = orientations,
        .allocator = allocator,
        .center_idx = center_idx,
        .axis_idx = axis_idx,
    };
}

/// Create an NH3+ rotator mover.
/// h_indices: indices of the 3 hydrogen atoms
/// center_idx: the N atom index
/// axis_idx: the atom defining rotation axis (CA usually)
/// Produces 3 orientations at 0°, 60°, -60° offsets.
/// Non-default orientations have a small penalty of 0.05.
pub fn createNH3Rotator(
    allocator: std.mem.Allocator,
    atoms: []const Atom,
    h_indices: [3]u32,
    center_idx: u32,
    axis_idx: u32,
    residue_idx: u32,
) !Mover {
    const center = atoms[center_idx].pos;
    const axis_atom = atoms[axis_idx].pos;
    const axis_dir = center.sub(axis_atom);

    const angle_offsets = [3]f32{ 0.0, 60.0, -60.0 };
    const penalties = [3]f32{ 0.0, 0.05, 0.05 };

    const orientations = try allocator.alloc(Orientation, 3);
    for (0..3) |i| {
        const positions = try allocator.alloc(Vec3(f32), 3);
        for (0..3) |j| {
            const h_pos = atoms[h_indices[j]].pos;
            positions[j] = math_mod.rotateAroundAxis(f32, h_pos, center, axis_dir, angle_offsets[i]);
        }
        orientations[i] = .{ .positions = positions, .penalty = penalties[i] };
    }

    const atom_indices = try allocator.alloc(u32, 3);
    for (0..3) |i| {
        atom_indices[i] = h_indices[i];
    }

    return Mover{
        .kind = .nh3_rotator,
        .residue_idx = residue_idx,
        .atom_indices = atom_indices,
        .orientations = orientations,
        .allocator = allocator,
        .center_idx = center_idx,
        .axis_idx = axis_idx,
    };
}

/// Create a methyl rotator mover (CH3).
/// h_indices: indices of the 3 hydrogen atoms
/// center_idx: the C atom index
/// axis_idx: the atom defining rotation axis
/// Produces 3 orientations at 0°, 60°, -60° offsets. No penalty.
pub fn createMethylRotator(
    allocator: std.mem.Allocator,
    atoms: []const Atom,
    h_indices: [3]u32,
    center_idx: u32,
    axis_idx: u32,
    residue_idx: u32,
) !Mover {
    const center = atoms[center_idx].pos;
    const axis_atom = atoms[axis_idx].pos;
    const axis_dir = center.sub(axis_atom);

    const angle_offsets = [3]f32{ 0.0, 60.0, -60.0 };

    const orientations = try allocator.alloc(Orientation, 3);
    for (0..3) |i| {
        const positions = try allocator.alloc(Vec3(f32), 3);
        for (0..3) |j| {
            const h_pos = atoms[h_indices[j]].pos;
            positions[j] = math_mod.rotateAroundAxis(f32, h_pos, center, axis_dir, angle_offsets[i]);
        }
        orientations[i] = .{ .positions = positions, .penalty = 0.0 };
    }

    const atom_indices = try allocator.alloc(u32, 3);
    for (0..3) |i| {
        atom_indices[i] = h_indices[i];
    }

    return Mover{
        .kind = .methyl_rotator,
        .residue_idx = residue_idx,
        .atom_indices = atom_indices,
        .orientations = orientations,
        .allocator = allocator,
        .center_idx = center_idx,
        .axis_idx = axis_idx,
    };
}

// ── Tests ─────────────────────────────────────────────────────────────────────

const testing = std.testing;

test "single H rotator has 12 orientations" {
    const allocator = testing.allocator;

    // Geometry: axis at (0,0,0), center (O) at (1,0,0), H at (2,0,0)
    const center_pos = Vec3(f32){ .x = 1.0, .y = 0.0, .z = 0.0 };
    const atoms = [_]Atom{
        .{ .pos = .{ .x = 2.0, .y = 0.0, .z = 0.0 } }, // H at index 0
        .{ .pos = center_pos }, // center (O) at index 1
        .{ .pos = .{ .x = 0.0, .y = 0.0, .z = 0.0 } }, // axis atom at index 2
    };

    var mover = try createSingleHRotator(allocator, &atoms, 0, 1, 2, 0);
    defer mover.deinit();

    try testing.expectEqual(@as(usize, 12), mover.orientations.len);
    for (mover.orientations) |o| {
        const dist = o.positions[0].distance(center_pos);
        try testing.expectApproxEqAbs(@as(f32, 1.0), dist, 0.01);
    }
}

test "single H rotator has no penalties" {
    const allocator = testing.allocator;

    const atoms = [_]Atom{
        .{ .pos = .{ .x = 2.0, .y = 0.0, .z = 0.0 } },
        .{ .pos = .{ .x = 1.0, .y = 0.0, .z = 0.0 } },
        .{ .pos = .{ .x = 0.0, .y = 0.0, .z = 0.0 } },
    };

    var mover = try createSingleHRotator(allocator, &atoms, 0, 1, 2, 0);
    defer mover.deinit();

    for (mover.orientations) |o| {
        try testing.expectApproxEqAbs(@as(f32, 0.0), o.penalty, 0.001);
    }
}

test "NH3 rotator has 3 orientations with 3 atoms each" {
    const allocator = testing.allocator;

    // N at origin, CA at (-1,0,0), 3 H atoms placed around N
    const atoms = [_]Atom{
        .{ .pos = .{ .x = 0.0, .y = 1.0, .z = 0.0 } }, // H0 at index 0
        .{ .pos = .{ .x = 0.866, .y = -0.5, .z = 0.0 } }, // H1 at index 1
        .{ .pos = .{ .x = -0.866, .y = -0.5, .z = 0.0 } }, // H2 at index 2
        .{ .pos = .{ .x = 0.0, .y = 0.0, .z = 0.0 } }, // N at index 3
        .{ .pos = .{ .x = -1.0, .y = 0.0, .z = 0.0 } }, // CA at index 4
    };

    var mover = try createNH3Rotator(allocator, &atoms, .{ 0, 1, 2 }, 3, 4, 0);
    defer mover.deinit();

    try testing.expectEqual(@as(usize, 3), mover.orientations.len);
    try testing.expectEqual(@as(usize, 3), mover.orientations[0].positions.len);
    try testing.expectEqual(@as(usize, 3), mover.orientations[1].positions.len);
    try testing.expectEqual(@as(usize, 3), mover.orientations[2].positions.len);
}

test "NH3 rotator first orientation has zero penalty, others have 0.05" {
    const allocator = testing.allocator;

    const atoms = [_]Atom{
        .{ .pos = .{ .x = 0.0, .y = 1.0, .z = 0.0 } },
        .{ .pos = .{ .x = 0.866, .y = -0.5, .z = 0.0 } },
        .{ .pos = .{ .x = -0.866, .y = -0.5, .z = 0.0 } },
        .{ .pos = .{ .x = 0.0, .y = 0.0, .z = 0.0 } },
        .{ .pos = .{ .x = -1.0, .y = 0.0, .z = 0.0 } },
    };

    var mover = try createNH3Rotator(allocator, &atoms, .{ 0, 1, 2 }, 3, 4, 0);
    defer mover.deinit();

    try testing.expectApproxEqAbs(@as(f32, 0.0), mover.orientations[0].penalty, 0.001);
    try testing.expectApproxEqAbs(@as(f32, 0.05), mover.orientations[1].penalty, 0.001);
    try testing.expectApproxEqAbs(@as(f32, 0.05), mover.orientations[2].penalty, 0.001);
}

test "methyl rotator has 3 orientations" {
    const allocator = testing.allocator;

    const atoms = [_]Atom{
        .{ .pos = .{ .x = 0.0, .y = 1.0, .z = 0.0 } },
        .{ .pos = .{ .x = 0.866, .y = -0.5, .z = 0.0 } },
        .{ .pos = .{ .x = -0.866, .y = -0.5, .z = 0.0 } },
        .{ .pos = .{ .x = 0.0, .y = 0.0, .z = 0.0 } }, // C at index 3
        .{ .pos = .{ .x = -1.0, .y = 0.0, .z = 0.0 } }, // axis atom at index 4
    };

    var mover = try createMethylRotator(allocator, &atoms, .{ 0, 1, 2 }, 3, 4, 0);
    defer mover.deinit();

    try testing.expectEqual(@as(usize, 3), mover.orientations.len);
    try testing.expectEqual(@as(usize, 3), mover.atom_indices.len);
}

test "methyl rotator has no penalties" {
    const allocator = testing.allocator;

    const atoms = [_]Atom{
        .{ .pos = .{ .x = 0.0, .y = 1.0, .z = 0.0 } },
        .{ .pos = .{ .x = 0.866, .y = -0.5, .z = 0.0 } },
        .{ .pos = .{ .x = -0.866, .y = -0.5, .z = 0.0 } },
        .{ .pos = .{ .x = 0.0, .y = 0.0, .z = 0.0 } },
        .{ .pos = .{ .x = -1.0, .y = 0.0, .z = 0.0 } },
    };

    var mover = try createMethylRotator(allocator, &atoms, .{ 0, 1, 2 }, 3, 4, 0);
    defer mover.deinit();

    for (mover.orientations) |o| {
        try testing.expectApproxEqAbs(@as(f32, 0.0), o.penalty, 0.001);
    }
}

test "rotator atom distances preserved after rotation" {
    const allocator = testing.allocator;

    // H at (2,0,0), center at (1,0,0), axis at (0,0,0)
    // H is 1.0 away from center; this should be preserved for all orientations
    const center_pos = Vec3(f32){ .x = 1.0, .y = 0.0, .z = 0.0 };
    const atoms = [_]Atom{
        .{ .pos = .{ .x = 2.0, .y = 0.0, .z = 0.0 } },
        .{ .pos = center_pos },
        .{ .pos = .{ .x = 0.0, .y = 0.0, .z = 0.0 } },
    };

    var mover = try createSingleHRotator(allocator, &atoms, 0, 1, 2, 5);
    defer mover.deinit();

    try testing.expectEqual(@as(u32, 5), mover.residue_idx);
    for (mover.orientations) |o| {
        const d = o.positions[0].distance(center_pos);
        try testing.expectApproxEqAbs(@as(f32, 1.0), d, 0.01);
    }
}
