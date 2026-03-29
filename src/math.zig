//! Math utilities for molecular geometry calculations.
//! Provides Vec3(T), rotation, dihedral, and angle functions.

const std = @import("std");
const math = std.math;

/// Generic 3D vector type.
pub fn Vec3(comptime T: type) type {
    return struct {
        const Self = @This();
        x: T,
        y: T,
        z: T,

        pub const zero = Self{ .x = 0, .y = 0, .z = 0 };

        pub fn add(self: Self, other: Self) Self {
            const a: @Vector(4, T) = .{ self.x, self.y, self.z, 0 };
            const b: @Vector(4, T) = .{ other.x, other.y, other.z, 0 };
            const r = a + b;
            return .{ .x = r[0], .y = r[1], .z = r[2] };
        }

        pub fn sub(self: Self, other: Self) Self {
            const a: @Vector(4, T) = .{ self.x, self.y, self.z, 0 };
            const b: @Vector(4, T) = .{ other.x, other.y, other.z, 0 };
            const r = a - b;
            return .{ .x = r[0], .y = r[1], .z = r[2] };
        }

        pub fn scale(self: Self, s: T) Self {
            return .{ .x = self.x * s, .y = self.y * s, .z = self.z * s };
        }

        /// Scale to given length. Returns zero if current length is zero (< 1e-10).
        pub fn scaleTo(self: Self, len: T) Self {
            const l = self.length();
            if (l < 1e-10) return Self.zero;
            return self.scale(len / l);
        }

        pub fn dot(self: Self, other: Self) T {
            const a: @Vector(4, T) = .{ self.x, self.y, self.z, 0 };
            const b: @Vector(4, T) = .{ other.x, other.y, other.z, 0 };
            return @reduce(.Add, a * b);
        }

        pub fn cross(self: Self, other: Self) Self {
            return .{
                .x = self.y * other.z - self.z * other.y,
                .y = self.z * other.x - self.x * other.z,
                .z = self.x * other.y - self.y * other.x,
            };
        }

        pub fn length(self: Self) T {
            return math.sqrt(self.dot(self));
        }

        pub fn distance(self: Self, other: Self) T {
            return self.sub(other).length();
        }

        pub fn normalize(self: Self) Self {
            const l = self.length();
            if (l < 1e-10) return Self.zero;
            return self.scale(1.0 / l);
        }

        pub fn negate(self: Self) Self {
            return .{ .x = -self.x, .y = -self.y, .z = -self.z };
        }

        /// Cast to Vec3(U) — supports f64 <-> f32 conversions.
        pub fn cast(self: Self, comptime U: type) Vec3(U) {
            return .{
                .x = @floatCast(self.x),
                .y = @floatCast(self.y),
                .z = @floatCast(self.z),
            };
        }
    };
}

/// Rotate point around axis (origin + direction) by angle in degrees.
/// Uses Rodrigues' rotation formula:
///   p' = p*cos + (u×p)*sin + u*(u·p)*(1-cos)
pub fn rotateAroundAxis(
    comptime T: type,
    point: Vec3(T),
    origin: Vec3(T),
    axis_dir: Vec3(T),
    degrees: T,
) Vec3(T) {
    const rad = degrees * math.pi / 180.0;
    const cos_a = math.cos(rad);
    const sin_a = math.sin(rad);

    // Translate point to origin-relative coordinates
    const p = point.sub(origin);
    const u = axis_dir.normalize();

    // Rodrigues: p' = p*cos + (u×p)*sin + u*(u·p)*(1-cos)
    const term1 = p.scale(cos_a);
    const term2 = u.cross(p).scale(sin_a);
    const term3 = u.scale(u.dot(p) * (1.0 - cos_a));

    return term1.add(term2).add(term3).add(origin);
}

/// Compute dihedral angle in degrees for 4 points a-b-c-d.
/// Uses the formula with normal vectors to planes a-b-c and b-c-d.
pub fn dihedralAngle(comptime T: type, a: Vec3(T), b: Vec3(T), c: Vec3(T), d: Vec3(T)) T {
    const b1 = b.sub(a);
    const b2 = c.sub(b);
    const b3 = d.sub(c);

    const n1 = b1.cross(b2);
    const n2 = b2.cross(b3);
    const m1 = n1.cross(b2.normalize());

    const y = m1.dot(n2);
    const x = n1.dot(n2);

    return math.atan2(y, x) * 180.0 / math.pi;
}

/// Compute angle in degrees at vertex b for points a-b-c.
/// Clamps dot product to [-1, 1] to avoid NaN from acos.
pub fn angle(comptime T: type, a: Vec3(T), b: Vec3(T), c: Vec3(T)) T {
    const ba = a.sub(b).normalize();
    const bc = c.sub(b).normalize();
    const cos_a = math.clamp(ba.dot(bc), -1.0, 1.0);
    return math.acos(cos_a) * 180.0 / math.pi;
}

// ── Tests ─────────────────────────────────────────────────────────────────────

test "Vec3 add" {
    const V = Vec3(f64);
    const a = V{ .x = 1, .y = 2, .z = 3 };
    const b = V{ .x = 4, .y = 5, .z = 6 };
    const result = a.add(b);
    try std.testing.expectApproxEqAbs(5.0, result.x, 1e-9);
    try std.testing.expectApproxEqAbs(7.0, result.y, 1e-9);
    try std.testing.expectApproxEqAbs(9.0, result.z, 1e-9);
}

test "Vec3 cross product" {
    const V = Vec3(f64);
    const x_axis = V{ .x = 1, .y = 0, .z = 0 };
    const y_axis = V{ .x = 0, .y = 1, .z = 0 };
    const result = x_axis.cross(y_axis);
    try std.testing.expectApproxEqAbs(0.0, result.x, 1e-9);
    try std.testing.expectApproxEqAbs(0.0, result.y, 1e-9);
    try std.testing.expectApproxEqAbs(1.0, result.z, 1e-9);
}

test "Vec3 normalize" {
    const V = Vec3(f64);
    const v = V{ .x = 3, .y = 4, .z = 0 };
    const n = v.normalize();
    try std.testing.expectApproxEqAbs(1.0, n.length(), 1e-9);
    try std.testing.expectApproxEqAbs(0.6, n.x, 1e-9);
    try std.testing.expectApproxEqAbs(0.8, n.y, 1e-9);
}

test "rotateAroundAxis 90 degrees around Z" {
    const V = Vec3(f64);
    const point = V{ .x = 1, .y = 0, .z = 0 };
    const origin = V.zero;
    const z_axis = V{ .x = 0, .y = 0, .z = 1 };
    const result = rotateAroundAxis(f64, point, origin, z_axis, 90.0);
    try std.testing.expectApproxEqAbs(0.0, result.x, 1e-9);
    try std.testing.expectApproxEqAbs(1.0, result.y, 1e-9);
    try std.testing.expectApproxEqAbs(0.0, result.z, 1e-9);
}

test "rotateAroundAxis 180 degrees around Z" {
    const V = Vec3(f64);
    const point = V{ .x = 1, .y = 0, .z = 0 };
    const origin = V.zero;
    const z_axis = V{ .x = 0, .y = 0, .z = 1 };
    const result = rotateAroundAxis(f64, point, origin, z_axis, 180.0);
    try std.testing.expectApproxEqAbs(-1.0, result.x, 1e-9);
    try std.testing.expectApproxEqAbs(0.0, result.y, 1e-9);
    try std.testing.expectApproxEqAbs(0.0, result.z, 1e-9);
}

test "dihedral angle trans (~180 degrees)" {
    // Classic trans dihedral using a zigzag geometry:
    // a=(0,1,0), b=(0,0,0), c=(1,0,0), d=(1,-1,0) -> 180 degrees
    const V = Vec3(f64);
    const a = V{ .x = 0, .y = 1, .z = 0 };
    const b = V{ .x = 0, .y = 0, .z = 0 };
    const c = V{ .x = 1, .y = 0, .z = 0 };
    const d = V{ .x = 1, .y = -1, .z = 0 };
    const result = dihedralAngle(f64, a, b, c, d);
    // Trans dihedral should be ~180 degrees
    try std.testing.expectApproxEqAbs(180.0, @abs(result), 1e-6);
}

test "angle 90 degrees" {
    const V = Vec3(f64);
    // Right angle: a=(1,0,0), b=(0,0,0), c=(0,1,0) -> 90 degrees at b
    const a = V{ .x = 1, .y = 0, .z = 0 };
    const b = V.zero;
    const c = V{ .x = 0, .y = 1, .z = 0 };
    const result = angle(f64, a, b, c);
    try std.testing.expectApproxEqAbs(90.0, result, 1e-9);
}

/// Fast approximation of exp(x) for x in [-87, 0] range.
/// Uses Schraudolph's integer-cast method. Max relative error ~6%.
/// Suitable for contact scoring where relative ranking matters, not absolute values.
pub fn fastExp(x: f32) f32 {
    const clamped = @max(x, -87.0);
    const v: i32 = @intFromFloat(12102203.0 * clamped + 1065353216.0);
    return @bitCast(@as(u32, @intCast(@max(v, 0))));
}

test "fastExp approximation within 7% for scoring range" {
    const test_values = [_]f32{ 0.0, -0.5, -1.0, -2.0, -4.0 };
    for (test_values) |x| {
        const exact = @exp(x);
        const approx = fastExp(x);
        const rel_err = @abs(approx - exact) / exact;
        try std.testing.expect(rel_err < 0.07);
    }
}

test "fastExp clamp boundary: x = -87 is on the boundary" {
    // x = -87 is at the clamp boundary; result should match exp(-87) within 7%.
    const exact = @exp(@as(f32, -87.0));
    const approx = fastExp(-87.0);
    const rel_err = @abs(approx - exact) / exact;
    try std.testing.expect(rel_err < 0.07);
}

test "fastExp clamp: x below -87 is clamped to exp(-87)" {
    // x = -100 is below the clamp boundary; fastExp should return the same as fastExp(-87).
    const clamped = fastExp(-87.0);
    const below = fastExp(-100.0);
    // Both should produce the same bit-identical result since -100 is clamped to -87.
    try std.testing.expectEqual(clamped, below);
}

test "fastExp positive input: x = 0.5 documents behavior" {
    // fastExp is designed for x in [-87, 0]; positive x is outside the intended range.
    // Document that it returns a value greater than 1.0 (exp(0.5) > 1).
    // The approximation accuracy for positive x is not guaranteed, but it must not crash.
    const result = fastExp(0.5);
    try std.testing.expect(result > 1.0);
}

test "fastExp sweep: max relative error stays within 7% across [-80, 0]" {
    // Sweep 160 evenly-spaced points in [-80, 0] and verify the worst-case error.
    var x: f32 = -80.0;
    while (x <= 0.0) : (x += 0.5) {
        const exact = @exp(x);
        const approx = fastExp(x);
        const rel_err = @abs(approx - exact) / exact;
        try std.testing.expect(rel_err < 0.07);
    }
}
