//! CIF value formatting, element symbol conversion, and fixed-point float output helpers.
//!
//! Extracted from mmcif_writer.zig to keep that module focused on
//! atom-site and document-level output.

const std = @import("std");
const element = @import("../element.zig");
const AtomType = element.AtomType;

/// Convert AtomType back to a 1-2 char element symbol string.
pub fn elementSymbol(atom_type: AtomType) []const u8 {
    return switch (atom_type) {
        .H, .Har, .Hpol, .Ha_p, .HOd => "H",
        .C, .Car, .C_eq_O => "C",
        .N, .Nacc => "N",
        .O => "O",
        .P => "P",
        .S => "S",
        .Se => "Se",
        .F => "F",
        .Cl => "Cl",
        .Br => "Br",
        .I => "I",
        .Li => "Li",
        .Na => "Na",
        .Mg => "Mg",
        .K => "K",
        .Ca => "Ca",
        .Mn => "Mn",
        .Fe => "Fe",
        .Co => "Co",
        .Ni => "Ni",
        .Cu => "Cu",
        .Zn => "Zn",
        .As => "As",
        .Rb => "Rb",
        .Sr => "Sr",
        .Mo => "Mo",
        .Ag => "Ag",
        .Cd => "Cd",
        .Sn => "Sn",
        .Cs => "Cs",
        .Ba => "Ba",
        .W => "W",
        .Pt => "Pt",
        .Au => "Au",
        .Hg => "Hg",
        .Pb => "Pb",
        .U => "U",
        .unknown => "X",
    };
}

/// Write an atom name, trimming leading and trailing spaces.
/// Atom names in CIF are typically unquoted even when they contain leading spaces.
pub fn writeAtomName(writer: anytype, name: []const u8) !void {
    // Trim trailing spaces
    var end = name.len;
    while (end > 0 and name[end - 1] == ' ') end -= 1;
    // Trim leading spaces
    var start: usize = 0;
    while (start < end and name[start] == ' ') start += 1;
    if (start >= end) {
        try writer.writeByte('.');
    } else {
        try writer.writeAll(name[start..end]);
    }
}

/// Write a CIF value, quoting if it contains spaces, quotes, or special characters.
/// Note: '.' and '?' are written unquoted as CIF null/unknown markers.
/// This is correct for round-tripping parsed CIF values where the parser
/// already stripped quotes from actual data values.
///
/// When `in_loop` is true, the semicolon text-field form is avoided because
/// it requires the ';' to appear at the start of a line, which is not
/// guaranteed mid-row. Instead, the value is wrapped in double-quotes with
/// any internal '"' replaced by '\'' (single-quote). This is slightly lossy
/// for the exotic case where a value contains both quote types, but produces
/// valid CIF that parses correctly.
pub fn writeCifValue(writer: anytype, val: []const u8) !void {
    return writeCifValueImpl(writer, val, false);
}

/// Loop-safe variant: never emits semicolon text fields.
pub fn writeCifValueInLoop(writer: anytype, val: []const u8) !void {
    return writeCifValueImpl(writer, val, true);
}

fn writeCifValueImpl(writer: anytype, val: []const u8, in_loop: bool) !void {
    if (val.len == 0) {
        try writer.writeByte('.');
        return;
    }
    // Check if quoting is needed
    var needs_quote = false;
    var has_single = false;
    var has_double = false;
    var has_newline = false;
    for (val) |c| {
        if (c == ' ' or c == '\t') needs_quote = true;
        if (c == '\'') has_single = true;
        if (c == '"') has_double = true;
        if (c == '\n' or c == '\r') has_newline = true;
    }
    // Starts with special char?
    if (val[0] == '_' or val[0] == '#' or val[0] == '$' or val[0] == ';' or
        val[0] == '[' or val[0] == ']' or val[0] == '{' or val[0] == '}') needs_quote = true;
    // Could be confused with CIF keyword?
    if (std.ascii.startsWithIgnoreCase(val, "data_") or
        std.ascii.startsWithIgnoreCase(val, "save_") or
        std.ascii.eqlIgnoreCase(val, "loop_") or
        std.ascii.eqlIgnoreCase(val, "stop_") or
        std.ascii.eqlIgnoreCase(val, "global_"))
    {
        needs_quote = true;
    }
    if (has_newline) needs_quote = true;

    if (!needs_quote) {
        try writer.writeAll(val);
    } else if (!has_single) {
        try writer.writeByte('\'');
        for (val) |c| {
            if (c == '\n' or c == '\r') {
                try writer.writeByte(' ');
            } else {
                try writer.writeByte(c);
            }
        }
        try writer.writeByte('\'');
    } else if (!has_double) {
        try writer.writeByte('"');
        for (val) |c| {
            if (c == '\n' or c == '\r') {
                try writer.writeByte(' ');
            } else {
                try writer.writeByte(c);
            }
        }
        try writer.writeByte('"');
    } else if (in_loop) {
        // Both quote types in a loop context: semicolon text fields are not safe
        // here because CIF requires ';' to appear at column 1 (start of line),
        // but we may be mid-row. Wrap in double-quotes and replace internal '"'
        // with '\'' so the output is valid CIF. Newlines are replaced with spaces.
        try writer.writeByte('"');
        for (val) |c| {
            if (c == '\n' or c == '\r') {
                try writer.writeByte(' ');
            } else if (c == '"') {
                try writer.writeByte('\'');
            } else {
                try writer.writeByte(c);
            }
        }
        try writer.writeByte('"');
    } else {
        // Both quote types outside a loop — use semicolon text field (lossless).
        // Semicolon text fields must start on a new line with ';' as the first char.
        try writer.writeAll("\n;");
        try writer.writeAll(val);
        if (val.len > 0 and val[val.len - 1] != '\n') try writer.writeByte('\n');
        try writer.writeAll(";\n");
    }
}

/// Write a CIF pair tag-value. For pairs (not in loops), semicolon text
/// fields are allowed since they appear on their own lines.
pub fn writePairCifValue(writer: anytype, val: []const u8) !void {
    if (val.len == 0) {
        try writer.writeByte('.');
        return;
    }
    var has_newline = false;
    for (val) |c| {
        if (c == '\n' or c == '\r') {
            has_newline = true;
            break;
        }
    }
    if (has_newline) {
        try writer.writeAll("\n;");
        try writer.writeAll(val);
        if (val[val.len - 1] != '\n') try writer.writeByte('\n');
        try writer.writeAll(";\n");
    } else {
        try writeCifValue(writer, val);
    }
}

/// Write a CIF pair tag-value, quoting the value if needed.
pub fn writePairValue(writer: anytype, tag: []const u8, val: []const u8) !void {
    try writer.writeAll(tag);
    var has_newline = false;
    for (val) |c| {
        if (c == '\n' or c == '\r') {
            has_newline = true;
            break;
        }
    }

    if (has_newline) {
        try writePairCifValue(writer, val);
        return;
    }

    try writer.writeByte(' ');
    try writePairCifValue(writer, val);
    try writer.writeByte('\n');
}

/// Write an altloc identifier ('.' for blank).
pub fn writeAltId(writer: anytype, altloc: u8) !void {
    if (altloc == ' ') {
        try writer.writeByte('.');
    } else {
        try writer.writeByte(altloc);
    }
}

// ── Fixed-point float helpers ────────────────────────────────────────────────

/// Write a float with exactly 3 decimal places using integer arithmetic.
/// Avoids std.fmt.formatFloat overhead for the common coordinate case.
pub fn writeFixedFloat3(writer: anytype, val: f32) !void {
    // Guard against NaN, Inf, and values too large for integer conversion.
    if (std.math.isNan(val) or std.math.isInf(val) or @abs(val) > 1.0e15) {
        return writer.print("{d:.3}", .{val});
    }
    if (val < 0) {
        const abs_val = -val;
        const scaled: u64 = @intFromFloat(@round(abs_val * 1000.0));
        if (scaled == 0) return writeFixedFloat3(writer, 0.0);
        try writer.writeByte('-');
        return writeFixedFloat3(writer, abs_val);
    }
    const scaled: u64 = @intFromFloat(@round(val * 1000.0));
    const int_part = scaled / 1000;
    const frac_part: u32 = @intCast(scaled % 1000);
    try writer.print("{d}", .{int_part});
    try writer.writeByte('.');
    if (frac_part < 10) {
        try writer.writeAll("00");
    } else if (frac_part < 100) {
        try writer.writeByte('0');
    }
    try writer.print("{d}", .{frac_part});
}

/// Write a float with exactly 2 decimal places using integer arithmetic.
/// Avoids std.fmt.formatFloat overhead for occupancy/b-factor output.
pub fn writeFixedFloat2(writer: anytype, val: f32) !void {
    // Guard against NaN, Inf, and values too large for integer conversion.
    if (std.math.isNan(val) or std.math.isInf(val) or @abs(val) > 1.0e15) {
        return writer.print("{d:.2}", .{val});
    }
    if (val < 0) {
        const abs_val = -val;
        const scaled: u64 = @intFromFloat(@round(abs_val * 100.0));
        if (scaled == 0) return writeFixedFloat2(writer, 0.0);
        try writer.writeByte('-');
        return writeFixedFloat2(writer, abs_val);
    }
    const scaled: u64 = @intFromFloat(@round(val * 100.0));
    const int_part = scaled / 100;
    const frac_part: u32 = @intCast(scaled % 100);
    try writer.print("{d}", .{int_part});
    try writer.writeByte('.');
    if (frac_part < 10) {
        try writer.writeByte('0');
    }
    try writer.print("{d}", .{frac_part});
}

// ── Tests ────────────────────────────────────────────────────────────────────

const testing = std.testing;

test "elementSymbol returns correct symbols" {
    try testing.expectEqualStrings("H", elementSymbol(.H));
    try testing.expectEqualStrings("H", elementSymbol(.Hpol));
    try testing.expectEqualStrings("C", elementSymbol(.Car));
    try testing.expectEqualStrings("N", elementSymbol(.Nacc));
    try testing.expectEqualStrings("Fe", elementSymbol(.Fe));
    try testing.expectEqualStrings("Cl", elementSymbol(.Cl));
}

test "writeFixedFloat3 formats coordinates correctly" {
    var buf: [32]u8 = undefined;
    var w: std.Io.Writer = .fixed(&buf);

    // Zero
    w = .fixed(&buf);
    try writeFixedFloat3(&w, 0.0);
    try testing.expectEqualStrings("0.000", w.buffered());

    // Positive integer-valued
    w = .fixed(&buf);
    try writeFixedFloat3(&w, 12.0);
    try testing.expectEqualStrings("12.000", w.buffered());

    // Positive with decimals
    w = .fixed(&buf);
    try writeFixedFloat3(&w, 1.5);
    try testing.expectEqualStrings("1.500", w.buffered());

    // Rounding: 123.4567 -> 123.457
    w = .fixed(&buf);
    try writeFixedFloat3(&w, 123.4567);
    try testing.expectEqualStrings("123.457", w.buffered());

    // Negative value
    w = .fixed(&buf);
    try writeFixedFloat3(&w, -4.321);
    try testing.expectEqualStrings("-4.321", w.buffered());

    // Small negative near zero: -0.001
    w = .fixed(&buf);
    try writeFixedFloat3(&w, -0.001);
    try testing.expectEqualStrings("-0.001", w.buffered());

    // Leading zeros in fractional part: 1.005
    w = .fixed(&buf);
    try writeFixedFloat3(&w, 1.005);
    try testing.expectEqualStrings("1.005", w.buffered());

    // Two leading zeros in fractional: 1.001
    w = .fixed(&buf);
    try writeFixedFloat3(&w, 1.001);
    try testing.expectEqualStrings("1.001", w.buffered());
}

test "writeFixedFloat2 formats occupancy/b-factor correctly" {
    var buf: [32]u8 = undefined;
    var w: std.Io.Writer = .fixed(&buf);

    // Zero
    w = .fixed(&buf);
    try writeFixedFloat2(&w, 0.0);
    try testing.expectEqualStrings("0.00", w.buffered());

    // 1.0 (common occupancy)
    w = .fixed(&buf);
    try writeFixedFloat2(&w, 1.0);
    try testing.expectEqualStrings("1.00", w.buffered());

    // 0.5 (half occupancy)
    w = .fixed(&buf);
    try writeFixedFloat2(&w, 0.5);
    try testing.expectEqualStrings("0.50", w.buffered());

    // Rounding: 10.456 -> 10.46
    w = .fixed(&buf);
    try writeFixedFloat2(&w, 10.456);
    try testing.expectEqualStrings("10.46", w.buffered());

    // Negative value
    w = .fixed(&buf);
    try writeFixedFloat2(&w, -3.14);
    try testing.expectEqualStrings("-3.14", w.buffered());

    // Leading zero in fractional part: 20.05
    w = .fixed(&buf);
    try writeFixedFloat2(&w, 20.05);
    try testing.expectEqualStrings("20.05", w.buffered());

    // Large b-factor
    w = .fixed(&buf);
    try writeFixedFloat2(&w, 99.99);
    try testing.expectEqualStrings("99.99", w.buffered());
}

test "writeFixedFloat3 does not panic on NaN, Inf, or overflow" {
    var buf: [64]u8 = undefined;
    var w: std.Io.Writer = .fixed(&buf);

    // NaN falls back to std.fmt (must not panic)
    w = .fixed(&buf);
    try writeFixedFloat3(&w, std.math.nan(f32));
    try testing.expect(w.buffered().len > 0);

    // Positive infinity falls back
    w = .fixed(&buf);
    try writeFixedFloat3(&w, std.math.inf(f32));
    try testing.expect(w.buffered().len > 0);

    // Negative infinity falls back
    w = .fixed(&buf);
    try writeFixedFloat3(&w, -std.math.inf(f32));
    try testing.expect(w.buffered().len > 0);

    // Very large value (1e20) falls back
    w = .fixed(&buf);
    try writeFixedFloat3(&w, 1.0e20);
    try testing.expect(w.buffered().len > 0);
}

test "writeFixedFloat2 does not panic on NaN, Inf, or overflow" {
    var buf: [64]u8 = undefined;
    var w: std.Io.Writer = .fixed(&buf);

    // NaN falls back to std.fmt (must not panic)
    w = .fixed(&buf);
    try writeFixedFloat2(&w, std.math.nan(f32));
    try testing.expect(w.buffered().len > 0);

    // Positive infinity falls back
    w = .fixed(&buf);
    try writeFixedFloat2(&w, std.math.inf(f32));
    try testing.expect(w.buffered().len > 0);

    // Negative infinity falls back
    w = .fixed(&buf);
    try writeFixedFloat2(&w, -std.math.inf(f32));
    try testing.expect(w.buffered().len > 0);

    // Very large value (1e20) falls back
    w = .fixed(&buf);
    try writeFixedFloat2(&w, 1.0e20);
    try testing.expect(w.buffered().len > 0);
}

test "writeCifValue quotes special-char-prefixed values" {
    // Test that values starting with [, ], {, } get quoted
    const cases = [_][]const u8{ "[bracket", "]close", "{brace", "}close" };
    for (cases) |input| {
        var buf: std.Io.Writer.Allocating = .init(testing.allocator);
        defer buf.deinit();
        try writeCifValue(&buf.writer, input);
        const output = buf.writer.buffered();
        // All should be quoted (start with ' or ")
        try testing.expect(output[0] == '\'' or output[0] == '"');
    }
    // Plain value should NOT be quoted
    {
        var buf: std.Io.Writer.Allocating = .init(testing.allocator);
        defer buf.deinit();
        try writeCifValue(&buf.writer, "hello");
        try testing.expectEqualStrings("hello", buf.writer.buffered());
    }
    // CIF null markers must remain bare for round-tripping preserved data.
    for ([_][]const u8{ ".", "?" }) |marker| {
        var buf: std.Io.Writer.Allocating = .init(testing.allocator);
        defer buf.deinit();
        try writeCifValue(&buf.writer, marker);
        try testing.expectEqualStrings(marker, buf.writer.buffered());
    }
}

test "writeCifValue uses semicolon field when value has both quote types" {
    const val = "it's a \"test\"";
    var buf: std.Io.Writer.Allocating = .init(testing.allocator);
    defer buf.deinit();
    try writeCifValue(&buf.writer, val);
    // Should start with newline+semicolon and end with newline+semicolon+newline
    try testing.expect(std.mem.startsWith(u8, buf.writer.buffered(), "\n;"));
    try testing.expect(std.mem.endsWith(u8, buf.writer.buffered(), ";\n"));
    // The original value must be preserved verbatim
    try testing.expect(std.mem.find(u8, buf.writer.buffered(), val) != null);
}

test "writeCifValueInLoop avoids semicolon field when value has both quote types" {
    const val = "it's a \"test\"";
    var buf: std.Io.Writer.Allocating = .init(testing.allocator);
    defer buf.deinit();
    try writeCifValueInLoop(&buf.writer, val);
    // Must NOT produce a semicolon text field (no leading newline+semicolon)
    try testing.expect(!std.mem.startsWith(u8, buf.writer.buffered(), "\n;"));
    // Must be wrapped in double-quotes
    try testing.expect(std.mem.startsWith(u8, buf.writer.buffered(), "\""));
    try testing.expect(std.mem.endsWith(u8, buf.writer.buffered(), "\""));
    // Internal double-quote should have been replaced (no raw " inside)
    const inner = buf.writer.buffered()[1 .. buf.writer.buffered().len - 1];
    try testing.expect(std.mem.find(u8, inner, "\"") == null);
}

test "writeFixedFloat3 no negative zero" {
    var buf: [32]u8 = undefined;
    var w: std.Io.Writer = .fixed(&buf);

    // -0.0004 rounds to 0 → must output "0.000", not "-0.000"
    w = .fixed(&buf);
    try writeFixedFloat3(&w, -0.0004);
    try testing.expectEqualStrings("0.000", w.buffered());

    // -0.0 → "0.000"
    w = .fixed(&buf);
    try writeFixedFloat3(&w, -0.0);
    try testing.expectEqualStrings("0.000", w.buffered());

    // A genuine small negative should still carry the minus sign
    w = .fixed(&buf);
    try writeFixedFloat3(&w, -0.001);
    try testing.expectEqualStrings("-0.001", w.buffered());
}

test "writeFixedFloat2 no negative zero" {
    var buf: [32]u8 = undefined;
    var w: std.Io.Writer = .fixed(&buf);

    // -0.004 rounds to 0 → must output "0.00", not "-0.00"
    w = .fixed(&buf);
    try writeFixedFloat2(&w, -0.004);
    try testing.expectEqualStrings("0.00", w.buffered());

    // -0.0 → "0.00"
    w = .fixed(&buf);
    try writeFixedFloat2(&w, -0.0);
    try testing.expectEqualStrings("0.00", w.buffered());

    // A genuine small negative should still carry the minus sign
    w = .fixed(&buf);
    try writeFixedFloat2(&w, -0.01);
    try testing.expectEqualStrings("-0.01", w.buffered());
}
