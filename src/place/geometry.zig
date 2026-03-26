//! Hydrogen placement geometry functions (Types 1-6).
//! Computes hydrogen atom positions based on reference atom positions
//! using different geometric rules.

const std = @import("std");
const testing = std.testing;
const math_mod = @import("../math.zig");
const Vec3 = math_mod.Vec3;

/// Type 1 (HXR3): Tetrahedral — H opposite to 3 neighbors (sp3, 3 neighbors present).
pub fn placeHXR3(center: Vec3(f64), n1: Vec3(f64), n2: Vec3(f64), n3: Vec3(f64), bond_len: f64) Vec3(f64) {
    const v1 = n1.sub(center).normalize();
    const v2 = n2.sub(center).normalize();
    const v3 = n3.sub(center).normalize();
    const dir = v1.add(v2).add(v3).scaleTo(-bond_len);
    return center.add(dir);
}

/// Type 2 (H2XR2): Two H on sp2 atom. Computes one H position;
/// call with different dihedral for the second H.
pub fn placeH2XR2(center: Vec3(f64), n1: Vec3(f64), n2: Vec3(f64), bond_len: f64, angle_deg: f64, dihedral_deg: f64) Vec3(f64) {
    const v1 = n1.sub(center).normalize();
    const v2 = n2.sub(center).normalize();
    const between = center.add(v1.add(v2).scale(0.5));
    return placeH3XR(center, between, n1, bond_len, angle_deg, dihedral_deg);
}

/// Type 3 (H3XR): Dihedral-controlled placement.
/// a1=center atom, a2=bonded atom, a3=reference for dihedral.
/// theta_deg=bond angle, phi_deg=dihedral angle.
pub fn placeH3XR(a1: Vec3(f64), a2: Vec3(f64), a3: Vec3(f64), bond_len: f64, theta_deg: f64, phi_deg: f64) Vec3(f64) {
    const v21 = a1.sub(a2).normalize();
    const v23 = a3.sub(a2).normalize();
    const norm = v21.cross(v23).scaleTo(bond_len);
    const pos4 = math_mod.rotateAroundAxis(f64, a1.add(norm), a2, a1.sub(a2), phi_deg - 90.0);
    const v14 = pos4.sub(a1).normalize();
    const v12 = a2.sub(a1).normalize();
    const pos5 = v14.cross(v12).add(a1);
    return math_mod.rotateAroundAxis(f64, pos4, a1, pos5.sub(a1), 90.0 - theta_deg);
}

/// Type 4 (HXR2): Planar bisector with fudge factor.
/// Places H in the plane of center-n1-n2, opposite to neighbors.
/// fudge biases the direction (0.0 = exact bisector).
pub fn placeHXR2Planar(center: Vec3(f64), n1: Vec3(f64), n2: Vec3(f64), bond_len: f64, fudge: f64) Vec3(f64) {
    const v1 = n1.sub(center).normalize();
    const v2 = n2.sub(center).normalize();
    const t = 0.5 + fudge;
    const interp = v1.scale(1.0 - t).add(v2.scale(t));
    return center.add(interp.scaleTo(-bond_len));
}

/// Type 5 (HXR2): Fractional angle distribution.
/// Used for backbone NH placement.
pub fn placeHXR2Frac(a1: Vec3(f64), a2: Vec3(f64), a3: Vec3(f64), bond_len: f64, fract: f64) Vec3(f64) {
    const v12 = a2.sub(a1).scaleTo(bond_len);
    const v13 = a3.sub(a1).scaleTo(bond_len);
    const pos4 = v12.cross(v13).add(a1);
    const cnca_angle = math_mod.angle(f64, a2, a1, a3);
    const hnca_angle = fract * (360.0 - cnca_angle);
    return math_mod.rotateAroundAxis(f64, a1.add(v12), pos4, a1.sub(pos4), hnca_angle);
}

/// Type 6 (HXY): Linear extension along center→away_from direction.
pub fn placeHXY(center: Vec3(f64), neighbor_atom: Vec3(f64), bond_len: f64) Vec3(f64) {
    const dir = center.sub(neighbor_atom).scaleTo(bond_len);
    return center.add(dir);
}

// ── Tests ─────────────────────────────────────────────────────────────────────

test "type1 HXR3 tetrahedral" {
    const center = Vec3(f64){ .x = 0.0, .y = 0.0, .z = 0.0 };
    const n1 = Vec3(f64){ .x = 1.0, .y = 0.0, .z = 0.0 };
    const n2 = Vec3(f64){ .x = 0.0, .y = 1.0, .z = 0.0 };
    const n3 = Vec3(f64){ .x = 0.0, .y = 0.0, .z = 1.0 };
    const h = placeHXR3(center, n1, n2, n3, 1.10);
    // H should be opposite to the 3 neighbors
    try testing.expect(h.x < 0);
    try testing.expect(h.y < 0);
    try testing.expect(h.z < 0);
    // Bond length should be 1.10
    try testing.expectApproxEqAbs(h.distance(center), 1.10, 1e-3);
}

test "type1 HXR3 symmetric tetrahedral bond length" {
    // All neighbors equidistant along axes: result should point equally negative
    const center = Vec3(f64){ .x = 1.0, .y = 1.0, .z = 1.0 };
    const n1 = Vec3(f64){ .x = 2.0, .y = 1.0, .z = 1.0 };
    const n2 = Vec3(f64){ .x = 1.0, .y = 2.0, .z = 1.0 };
    const n3 = Vec3(f64){ .x = 1.0, .y = 1.0, .z = 2.0 };
    const h = placeHXR3(center, n1, n2, n3, 1.08);
    try testing.expectApproxEqAbs(h.distance(center), 1.08, 1e-6);
    // H should be in negative direction from center
    try testing.expect(h.x < center.x);
    try testing.expect(h.y < center.y);
    try testing.expect(h.z < center.z);
}

test "type3 H3XR dihedral placement bond length" {
    // Test that placeH3XR produces a bond of the requested length
    const a1 = Vec3(f64){ .x = 0.0, .y = 0.0, .z = 0.0 };
    const a2 = Vec3(f64){ .x = 1.5, .y = 0.0, .z = 0.0 };
    const a3 = Vec3(f64){ .x = 1.5, .y = 1.5, .z = 0.0 };
    const h = placeH3XR(a1, a2, a3, 1.09, 109.5, 180.0);
    try testing.expectApproxEqAbs(h.distance(a1), 1.09, 1e-3);
}

test "type3 H3XR different dihedral produces different position" {
    const a1 = Vec3(f64){ .x = 0.0, .y = 0.0, .z = 0.0 };
    const a2 = Vec3(f64){ .x = 1.5, .y = 0.0, .z = 0.0 };
    const a3 = Vec3(f64){ .x = 1.5, .y = 1.5, .z = 0.0 };
    const h1 = placeH3XR(a1, a2, a3, 1.09, 109.5, 60.0);
    const h2 = placeH3XR(a1, a2, a3, 1.09, 109.5, 180.0);
    // Different dihedrals should produce different positions
    try testing.expect(h1.distance(h2) > 0.01);
}

test "type4 HXR2Planar bisector bond length" {
    // Symmetric case: two neighbors at 90 degrees, H should point along bisector
    const center = Vec3(f64){ .x = 0.0, .y = 0.0, .z = 0.0 };
    const n1 = Vec3(f64){ .x = 1.0, .y = 0.0, .z = 0.0 };
    const n2 = Vec3(f64){ .x = 0.0, .y = 1.0, .z = 0.0 };
    const h = placeHXR2Planar(center, n1, n2, 1.0, 0.0);
    try testing.expectApproxEqAbs(h.distance(center), 1.0, 1e-6);
    // H should point in the negative bisector direction
    try testing.expect(h.x < 0);
    try testing.expect(h.y < 0);
    // z should be approximately zero (planar)
    try testing.expectApproxEqAbs(h.z, 0.0, 1e-6);
}

test "type4 HXR2Planar fudge biases direction" {
    const center = Vec3(f64){ .x = 0.0, .y = 0.0, .z = 0.0 };
    const n1 = Vec3(f64){ .x = 1.0, .y = 0.0, .z = 0.0 };
    const n2 = Vec3(f64){ .x = 0.0, .y = 1.0, .z = 0.0 };
    const h0 = placeHXR2Planar(center, n1, n2, 1.0, 0.0);
    const h_fudge = placeHXR2Planar(center, n1, n2, 1.0, 0.2);
    // Fudge should produce a different direction
    try testing.expect(h0.distance(h_fudge) > 1e-6);
}

test "type6 HXY linear" {
    const center = Vec3(f64){ .x = 0.0, .y = 0.0, .z = 0.0 };
    const neighbor_atom = Vec3(f64){ .x = -1.5, .y = 0.0, .z = 0.0 };
    const h = placeHXY(center, neighbor_atom, 1.0);
    try testing.expectApproxEqAbs(h.x, 1.0, 1e-6);
    try testing.expectApproxEqAbs(h.y, 0.0, 1e-6);
}

test "type6 HXY bond length preserved" {
    const center = Vec3(f64){ .x = 1.0, .y = 2.0, .z = 3.0 };
    const neighbor = Vec3(f64){ .x = 0.0, .y = 0.0, .z = 0.0 };
    const h = placeHXY(center, neighbor, 1.08);
    try testing.expectApproxEqAbs(h.distance(center), 1.08, 1e-6);
}

test "type2 H2XR2 bond length and distinct positions" {
    const center = Vec3(f64){ .x = 0.0, .y = 0.0, .z = 0.0 };
    const n1 = Vec3(f64){ .x = 1.5, .y = 0.0, .z = 0.0 };
    const n2 = Vec3(f64){ .x = 0.0, .y = 1.5, .z = 0.0 };
    const h1 = placeH2XR2(center, n1, n2, 1.10, 109.5, 120.0);
    const h2 = placeH2XR2(center, n1, n2, 1.10, 109.5, -120.0);
    // Both should have correct bond length
    try testing.expectApproxEqAbs(h1.distance(center), 1.10, 0.05);
    try testing.expectApproxEqAbs(h2.distance(center), 1.10, 0.05);
    // Different dihedrals should produce distinct positions
    try testing.expect(h1.distance(h2) > 0.1);
}

test "type5 HXR2Frac bond length" {
    const a1 = Vec3(f64){ .x = 0.0, .y = 0.0, .z = 0.0 };
    const a2 = Vec3(f64){ .x = 1.5, .y = 0.0, .z = 0.0 };
    const a3 = Vec3(f64){ .x = 0.0, .y = 1.5, .z = 0.0 };
    const h = placeHXR2Frac(a1, a2, a3, 1.02, 0.5);
    try testing.expectApproxEqAbs(h.distance(a1), 1.02, 0.05);
}
