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
            return .{ .x = self.x + other.x, .y = self.y + other.y, .z = self.z + other.z };
        }

        pub fn sub(self: Self, other: Self) Self {
            return .{ .x = self.x - other.x, .y = self.y - other.y, .z = self.z - other.z };
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
            return self.x * other.x + self.y * other.y + self.z * other.z;
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
