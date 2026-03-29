//! Dot-sphere bump/H-bond scoring for probe-based contact evaluation.
//! Direct port of the original reduce scoring logic (probe.cpp / ScoreContact).

const std = @import("std");
const math_mod = @import("../math.zig");
const Vec3 = math_mod.Vec3;
const element = @import("../element.zig");
const dot_sphere = @import("dot_sphere.zig");
const DotSphere = dot_sphere.DotSphere;

pub const ScoringParams = struct {
    gap_scale: f32 = 0.25, // Contact score gaussian width
    bump_weight: f32 = 10.0, // Clash penalty multiplier
    hb_weight: f32 = 4.0, // H-bond reward multiplier
    min_reg_hb_gap: f32 = 0.6, // Regular H-bond gap threshold (Å)
    min_charged_hb_gap: f32 = 0.8, // Charged H-bond gap threshold (Å)
    bad_bump_gap_cut: f32 = 0.4, // Bad bump classification threshold (Å)
    dot_density: f32 = 16.0, // Dots per Å²
    probe_radius: f32 = 0.0, // Probe sphere radius
};

pub const ScoreResult = struct {
    total: f32 = 0.0,
    contact_sub: f32 = 0.0,
    bump_sub: f32 = 0.0,
    hb_sub: f32 = 0.0,
    n_bad_bumps: u32 = 0,
};

/// Score a single atom against its environment using dot sphere probes.
/// atom_pos: position of the atom being scored
/// atom_radius: VDW radius of the atom
/// atom_flags: donor/acceptor flags of the atom
/// neighbors: positions of nearby atoms
/// neighbor_radii: VDW radii of nearby atoms
/// neighbor_flags: flags of nearby atoms
/// sphere: pre-generated dot sphere for this atom's radius
/// params: scoring parameters
pub fn scoreAtom(
    atom_pos: Vec3(f32),
    atom_radius: f32,
    atom_flags: element.AtomFlags,
    neighbors: []const Vec3(f32),
    neighbor_radii: []const f32,
    neighbor_flags: []const element.AtomFlags,
    sphere: *const DotSphere,
    params: ScoringParams,
) ScoreResult {
    _ = atom_radius; // radius is encoded in the sphere points already
    var result = ScoreResult{};
    const dot_scale: f32 = 1.0 / @as(f32, @floatFromInt(sphere.points.len));

    for (sphere.points) |dot| {
        // Scale dot to atom surface (sphere is at atom_radius already)
        const probe_pos = Vec3(f32){
            .x = atom_pos.x + dot.x,
            .y = atom_pos.y + dot.y,
            .z = atom_pos.z + dot.z,
        };

        // Find minimum gap to any neighbor
        var min_gap: f32 = std.math.floatMax(f32);
        var min_idx: usize = 0;

        for (neighbors, 0..) |npos, i| {
            const dist = probe_pos.distance(npos);
            const gap = dist - neighbor_radii[i] - params.probe_radius;
            if (gap < min_gap) {
                min_gap = gap;
                min_idx = i;
            }
        }

        if (min_gap >= std.math.floatMax(f32) * 0.5) continue; // no neighbors

        // Classify and score this dot
        if (min_gap >= 0.0) {
            // Contact (van der Waals contact)
            const score = math_mod.fastExp(-(min_gap / params.gap_scale) * (min_gap / params.gap_scale));
            result.contact_sub += score * dot_scale;
            result.total += score * dot_scale;
        } else {
            // Overlap — check for H-bond
            const is_hb = isHBond(atom_flags, neighbor_flags[min_idx], min_gap, params);
            if (is_hb) {
                const score = params.hb_weight * (-0.5 * min_gap);
                result.hb_sub += score * dot_scale;
                result.total += score * dot_scale;
            } else {
                // Clash (bump)
                const score = -params.bump_weight * (-0.5 * min_gap);
                result.bump_sub += score * dot_scale;
                result.total += score * dot_scale;
                if (-min_gap >= params.bad_bump_gap_cut) {
                    result.n_bad_bumps += 1;
                }
            }
        }
    }

    return result;
}

/// Check if an overlap qualifies as an H-bond.
/// Requires one atom to be a donor (H) and the other an acceptor,
/// and the overlap must be within the H-bond gap threshold.
pub fn isHBond(flags_a: element.AtomFlags, flags_b: element.AtomFlags, gap: f32, params: ScoringParams) bool {
    const a_donor = flags_a.donor;
    const b_donor = flags_b.donor;
    const a_acc = flags_a.acceptor;
    const b_acc = flags_b.acceptor;

    // One must be donor, other must be acceptor
    if (!((a_donor and b_acc) or (b_donor and a_acc))) return false;

    // Check gap threshold
    const charged = (flags_a.positive or flags_a.negative) and (flags_b.positive or flags_b.negative);
    const threshold = if (charged) params.min_charged_hb_gap else params.min_reg_hb_gap;

    return (-gap) <= threshold;
}

// ── Tests ─────────────────────────────────────────────────────────────────────

const testing = std.testing;

test "score no overlap gives positive contact" {
    // Two atoms at VDW contact distance should give positive score
    // Atom at origin with r=1.70, neighbor at distance 3.10 (slight contact)
    var sphere = try DotSphere.generate(testing.allocator, 1.70, 16.0);
    defer sphere.deinit();

    const atom_pos = Vec3(f32){ .x = 0.0, .y = 0.0, .z = 0.0 };
    const neighbor_pos = [_]Vec3(f32){.{ .x = 3.20, .y = 0.0, .z = 0.0 }}; // slight contact (gap ~0.1Å)
    const neighbor_radii = [_]f32{1.40};
    const neighbor_flags = [_]element.AtomFlags{.{}};

    const result = scoreAtom(atom_pos, 1.70, .{}, &neighbor_pos, &neighbor_radii, &neighbor_flags, &sphere, .{});
    try testing.expect(result.total > 0);
    try testing.expect(result.bump_sub == 0);
}

test "score clash gives negative bump" {
    // Overlapping atoms should produce negative bump score
    var sphere = try DotSphere.generate(testing.allocator, 1.70, 16.0);
    defer sphere.deinit();

    const atom_pos = Vec3(f32){ .x = 0.0, .y = 0.0, .z = 0.0 };
    const neighbor_pos = [_]Vec3(f32){.{ .x = 2.0, .y = 0.0, .z = 0.0 }}; // heavy overlap
    const neighbor_radii = [_]f32{1.70};
    const neighbor_flags = [_]element.AtomFlags{.{}};

    const result = scoreAtom(atom_pos, 1.70, .{}, &neighbor_pos, &neighbor_radii, &neighbor_flags, &sphere, .{});
    try testing.expect(result.total < 0);
    try testing.expect(result.bump_sub < 0);
}

test "score H-bond overlap gives positive" {
    // Donor-acceptor overlap within H-bond threshold should be positive
    var sphere = try DotSphere.generate(testing.allocator, 1.05, 16.0);
    defer sphere.deinit();

    const atom_pos = Vec3(f32){ .x = 0.0, .y = 0.0, .z = 0.0 };
    const neighbor_pos = [_]Vec3(f32){.{ .x = 2.15, .y = 0.0, .z = 0.0 }}; // slight overlap, H-bond range
    const neighbor_radii = [_]f32{1.40};
    const donor_flags = element.AtomFlags{ .donor = true };
    const acc_flags = [_]element.AtomFlags{element.AtomFlags{ .acceptor = true }};

    const result = scoreAtom(atom_pos, 1.05, donor_flags, &neighbor_pos, &neighbor_radii, &acc_flags, &sphere, .{});
    try testing.expect(result.hb_sub > 0);
}
