//! CIF value formatting and fixed-point float output helpers.
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

/// Write an atom name, trimming trailing spaces.
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
    var fbs = std.io.fixedBufferStream(&buf);
    const w = fbs.writer();

    // Zero
    fbs.reset();
    try writeFixedFloat3(w, 0.0);
    try testing.expectEqualStrings("0.000", fbs.getWritten());

    // Positive integer-valued
    fbs.reset();
    try writeFixedFloat3(w, 12.0);
    try testing.expectEqualStrings("12.000", fbs.getWritten());

    // Positive with decimals
    fbs.reset();
    try writeFixedFloat3(w, 1.5);
    try testing.expectEqualStrings("1.500", fbs.getWritten());

    // Rounding: 123.4567 -> 123.457
    fbs.reset();
    try writeFixedFloat3(w, 123.4567);
    try testing.expectEqualStrings("123.457", fbs.getWritten());

    // Negative value
    fbs.reset();
    try writeFixedFloat3(w, -4.321);
    try testing.expectEqualStrings("-4.321", fbs.getWritten());

    // Small negative near zero: -0.001
    fbs.reset();
    try writeFixedFloat3(w, -0.001);
    try testing.expectEqualStrings("-0.001", fbs.getWritten());

    // Leading zeros in fractional part: 1.005
    fbs.reset();
    try writeFixedFloat3(w, 1.005);
    try testing.expectEqualStrings("1.005", fbs.getWritten());

    // Two leading zeros in fractional: 1.001
    fbs.reset();
    try writeFixedFloat3(w, 1.001);
    try testing.expectEqualStrings("1.001", fbs.getWritten());
}

test "writeFixedFloat2 formats occupancy/b-factor correctly" {
    var buf: [32]u8 = undefined;
    var fbs = std.io.fixedBufferStream(&buf);
    const w = fbs.writer();

    // Zero
    fbs.reset();
    try writeFixedFloat2(w, 0.0);
    try testing.expectEqualStrings("0.00", fbs.getWritten());

    // 1.0 (common occupancy)
    fbs.reset();
    try writeFixedFloat2(w, 1.0);
    try testing.expectEqualStrings("1.00", fbs.getWritten());

    // 0.5 (half occupancy)
    fbs.reset();
    try writeFixedFloat2(w, 0.5);
    try testing.expectEqualStrings("0.50", fbs.getWritten());

    // Rounding: 10.456 -> 10.46
    fbs.reset();
    try writeFixedFloat2(w, 10.456);
    try testing.expectEqualStrings("10.46", fbs.getWritten());

    // Negative value
    fbs.reset();
    try writeFixedFloat2(w, -3.14);
    try testing.expectEqualStrings("-3.14", fbs.getWritten());

    // Leading zero in fractional part: 20.05
    fbs.reset();
    try writeFixedFloat2(w, 20.05);
    try testing.expectEqualStrings("20.05", fbs.getWritten());

    // Large b-factor
    fbs.reset();
    try writeFixedFloat2(w, 99.99);
    try testing.expectEqualStrings("99.99", fbs.getWritten());
}

test "writeFixedFloat3 does not panic on NaN, Inf, or overflow" {
    var buf: [64]u8 = undefined;
    var fbs = std.io.fixedBufferStream(&buf);
    const w = fbs.writer();

    // NaN falls back to std.fmt (must not panic)
    fbs.reset();
    try writeFixedFloat3(w, std.math.nan(f32));
    try testing.expect(fbs.getWritten().len > 0);

    // Positive infinity falls back
    fbs.reset();
    try writeFixedFloat3(w, std.math.inf(f32));
    try testing.expect(fbs.getWritten().len > 0);

    // Negative infinity falls back
    fbs.reset();
    try writeFixedFloat3(w, -std.math.inf(f32));
    try testing.expect(fbs.getWritten().len > 0);

    // Very large value (1e20) falls back
    fbs.reset();
    try writeFixedFloat3(w, 1.0e20);
    try testing.expect(fbs.getWritten().len > 0);
}

test "writeFixedFloat2 does not panic on NaN, Inf, or overflow" {
    var buf: [64]u8 = undefined;
    var fbs = std.io.fixedBufferStream(&buf);
    const w = fbs.writer();

    // NaN falls back to std.fmt (must not panic)
    fbs.reset();
    try writeFixedFloat2(w, std.math.nan(f32));
    try testing.expect(fbs.getWritten().len > 0);

    // Positive infinity falls back
    fbs.reset();
    try writeFixedFloat2(w, std.math.inf(f32));
    try testing.expect(fbs.getWritten().len > 0);

    // Negative infinity falls back
    fbs.reset();
    try writeFixedFloat2(w, -std.math.inf(f32));
    try testing.expect(fbs.getWritten().len > 0);

    // Very large value (1e20) falls back
    fbs.reset();
    try writeFixedFloat2(w, 1.0e20);
    try testing.expect(fbs.getWritten().len > 0);
}

test "writeCifValue quotes special-char-prefixed values" {
    // Test that values starting with [, ], {, } get quoted
    const cases = [_][]const u8{ "[bracket", "]close", "{brace", "}close" };
    for (cases) |input| {
        var buf = std.ArrayList(u8).empty;
        defer buf.deinit(testing.allocator);
        try writeCifValue(buf.writer(testing.allocator), input);
        const output = buf.items;
        // All should be quoted (start with ' or ")
        try testing.expect(output[0] == '\'' or output[0] == '"');
    }
    // Plain value should NOT be quoted
    {
        var buf = std.ArrayList(u8).empty;
        defer buf.deinit(testing.allocator);
        try writeCifValue(buf.writer(testing.allocator), "hello");
        try testing.expectEqualStrings("hello", buf.items);
    }
    // CIF null markers must remain bare for round-tripping preserved data.
    for ([_][]const u8{ ".", "?" }) |marker| {
        var buf = std.ArrayList(u8).empty;
        defer buf.deinit(testing.allocator);
        try writeCifValue(buf.writer(testing.allocator), marker);
        try testing.expectEqualStrings(marker, buf.items);
    }
}

test "writeCifValue uses semicolon field when value has both quote types" {
    const val = "it's a \"test\"";
    var buf = std.ArrayList(u8).empty;
    defer buf.deinit(testing.allocator);
    try writeCifValue(buf.writer(testing.allocator), val);
    // Should start with newline+semicolon and end with newline+semicolon+newline
    try testing.expect(std.mem.startsWith(u8, buf.items, "\n;"));
    try testing.expect(std.mem.endsWith(u8, buf.items, ";\n"));
    // The original value must be preserved verbatim
    try testing.expect(std.mem.indexOf(u8, buf.items, val) != null);
}

test "writeCifValueInLoop avoids semicolon field when value has both quote types" {
    const val = "it's a \"test\"";
    var buf = std.ArrayList(u8).empty;
    defer buf.deinit(testing.allocator);
    try writeCifValueInLoop(buf.writer(testing.allocator), val);
    // Must NOT produce a semicolon text field (no leading newline+semicolon)
    try testing.expect(!std.mem.startsWith(u8, buf.items, "\n;"));
    // Must be wrapped in double-quotes
    try testing.expect(std.mem.startsWith(u8, buf.items, "\""));
    try testing.expect(std.mem.endsWith(u8, buf.items, "\""));
    // Internal double-quote should have been replaced (no raw " inside)
    const inner = buf.items[1 .. buf.items.len - 1];
    try testing.expect(std.mem.indexOf(u8, inner, "\"") == null);
}

test "writeFixedFloat3 no negative zero" {
    var buf: [32]u8 = undefined;
    var fbs = std.io.fixedBufferStream(&buf);
    const w = fbs.writer();

    // -0.0004 rounds to 0 → must output "0.000", not "-0.000"
    fbs.reset();
    try writeFixedFloat3(w, -0.0004);
    try testing.expectEqualStrings("0.000", fbs.getWritten());

    // -0.0 → "0.000"
    fbs.reset();
    try writeFixedFloat3(w, -0.0);
    try testing.expectEqualStrings("0.000", fbs.getWritten());

    // A genuine small negative should still carry the minus sign
    fbs.reset();
    try writeFixedFloat3(w, -0.001);
    try testing.expectEqualStrings("-0.001", fbs.getWritten());
}

test "writeFixedFloat2 no negative zero" {
    var buf: [32]u8 = undefined;
    var fbs = std.io.fixedBufferStream(&buf);
    const w = fbs.writer();

    // -0.004 rounds to 0 → must output "0.00", not "-0.00"
    fbs.reset();
    try writeFixedFloat2(w, -0.004);
    try testing.expectEqualStrings("0.00", fbs.getWritten());

    // -0.0 → "0.00"
    fbs.reset();
    try writeFixedFloat2(w, -0.0);
    try testing.expectEqualStrings("0.00", fbs.getWritten());

    // A genuine small negative should still carry the minus sign
    fbs.reset();
    try writeFixedFloat2(w, -0.01);
    try testing.expectEqualStrings("-0.01", fbs.getWritten());
}
