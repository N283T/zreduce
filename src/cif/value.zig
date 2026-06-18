//! CIF value helpers: null-checking and type conversions.

const std = @import("std");

/// Returns true if the value is the CIF null sentinel ("." or "?").
pub fn isNull(v: []const u8) bool {
    return std.mem.eql(u8, v, ".") or std.mem.eql(u8, v, "?");
}

/// Returns "" for null values, otherwise returns v unchanged.
pub fn asString(v: []const u8) []const u8 {
    if (isNull(v)) return "";
    return v;
}

/// Parse v as f32. Handles uncertainty suffixes like "1.234(5)".
/// Returns null for CIF null values or parse failures.
pub fn asFloat(v: []const u8) ?f32 {
    if (isNull(v)) return null;
    // Strip uncertainty suffix "(N...)" if present
    const src = if (std.mem.findScalar(u8, v, '(')) |paren_idx|
        v[0..paren_idx]
    else
        v;
    return std.fmt.parseFloat(f32, src) catch null;
}

/// Like asFloat but returns `default` instead of null.
pub fn asFloatOr(v: []const u8, default: f32) f32 {
    return asFloat(v) orelse default;
}

/// Parse v as integer type T. Returns null for CIF null values or parse failures.
pub fn asInt(comptime T: type, v: []const u8) ?T {
    if (isNull(v)) return null;
    return std.fmt.parseInt(T, v, 10) catch null;
}

/// Like asInt but returns `default` instead of null.
pub fn asIntOr(comptime T: type, v: []const u8, default: T) T {
    return asInt(T, v) orelse default;
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

const testing = std.testing;

test "isNull" {
    try testing.expect(isNull("."));
    try testing.expect(isNull("?"));
    try testing.expect(!isNull("1.0"));
    try testing.expect(!isNull("hello"));
}

test "asFloat with uncertainty suffix" {
    const result = asFloat("1.234(5)");
    try testing.expect(result != null);
    try testing.expectApproxEqAbs(@as(f32, 1.234), result.?, 1e-4);
}

test "asFloat null values" {
    try testing.expectEqual(@as(?f32, null), asFloat("."));
    try testing.expectEqual(@as(?f32, null), asFloat("?"));
}

test "asFloat plain" {
    const result = asFloat("3.14");
    try testing.expect(result != null);
    try testing.expectApproxEqAbs(@as(f32, 3.14), result.?, 1e-4);
}

test "asInt" {
    try testing.expectEqual(@as(?i32, 42), asInt(i32, "42"));
    try testing.expectEqual(@as(?i32, null), asInt(i32, "."));
    try testing.expectEqual(@as(?i32, null), asInt(i32, "abc"));
}

test "asString" {
    try testing.expectEqualStrings("", asString("."));
    try testing.expectEqualStrings("", asString("?"));
    try testing.expectEqualStrings("hello", asString("hello"));
}
