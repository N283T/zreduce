const std = @import("std");
const model_mod = @import("../model.zig");
const mover_mod = @import("mover.zig");

const Model = model_mod.Model;
const Mover = mover_mod.Mover;
const MoverKind = mover_mod.MoverKind;

pub const Selector = struct {
    chain_id: []const u8,
    auth_seq_id: i32,
    ins_code: u8 = ' ',

    pub fn matches(self: Selector, mdl: *const Model, mover: *const Mover) bool {
        const res = mdl.residues.items[mover.residue_idx];
        const chain = mdl.chains.items[res.chain_idx];
        if (res.auth_seq_id != self.auth_seq_id) return false;
        if (res.ins_code != self.ins_code) return false;
        return std.mem.eql(u8, chain.authSlice(), self.chain_id) or
            std.mem.eql(u8, chain.labelSlice(), self.chain_id);
    }
};

pub const Entry = struct {
    selector: Selector,
    comp_id: [3]u8,
    comp_id_len: u2,
    target: []const u8,
    value: Value,

    pub fn compIdSlice(self: *const Entry) []const u8 {
        return self.comp_id[0..self.comp_id_len];
    }
};

pub const Value = union(enum) {
    orientation: u16,
    amide: enum { original, flip },
    his: enum { hie, hid, hie_flip, hid_flip },
};

pub const FixOverrides = struct {
    allocator: std.mem.Allocator,
    entries: []Entry,

    pub fn deinit(self: *FixOverrides) void {
        for (self.entries) |entry| {
            self.allocator.free(entry.selector.chain_id);
            self.allocator.free(entry.target);
        }
        self.allocator.free(self.entries);
        self.entries = &.{};
    }

    pub fn warnUnmatched(self: *const FixOverrides, mdl: *const Model, movers: []const Mover) void {
        for (self.entries) |entry| {
            var matched = false;
            for (movers) |*m| {
                if (!entry.selector.matches(mdl, m)) continue;
                const res = mdl.residues.items[m.residue_idx];
                if (!std.mem.eql(u8, entry.compIdSlice(), res.compIdSlice())) continue;
                const target = moverTargetName(mdl, m);
                if (std.mem.eql(u8, target, entry.target)) {
                    matched = true;
                    break;
                }
            }
            if (!matched) {
                std.debug.print("  warning: fix override for {s}:{d} {s} {s} did not match any mover\n", .{
                    entry.selector.chain_id,
                    entry.selector.auth_seq_id,
                    entry.compIdSlice(),
                    entry.target,
                });
            }
        }
    }
};

pub fn parseFile(allocator: std.mem.Allocator, path: []const u8) !FixOverrides {
    const file = try std.fs.cwd().openFile(path, .{});
    defer file.close();
    const source = try file.readToEndAlloc(allocator, 1024 * 1024);
    defer allocator.free(source);
    return parseString(allocator, source);
}

pub fn parseString(allocator: std.mem.Allocator, source: []const u8) !FixOverrides {
    var entries = std.ArrayListUnmanaged(Entry){};
    errdefer {
        for (entries.items) |entry| {
            allocator.free(entry.selector.chain_id);
            allocator.free(entry.target);
        }
        entries.deinit(allocator);
    }

    var lines = std.mem.splitScalar(u8, source, '\n');
    var line_num: u32 = 0;
    while (lines.next()) |line_raw| {
        line_num += 1;
        const line = std.mem.trim(u8, line_raw, " \t\r");
        if (line.len == 0 or line[0] == '#') continue;

        var toks = std.mem.tokenizeAny(u8, line, " \t");
        const selector_tok = toks.next() orelse {
            std.debug.print("Warning: fix override line {d}: expected 'chain:seq comp target value'\n", .{line_num});
            continue;
        };
        const comp_tok = toks.next() orelse {
            std.debug.print("Warning: fix override line {d}: missing comp_id\n", .{line_num});
            continue;
        };
        const target_tok = toks.next() orelse {
            std.debug.print("Warning: fix override line {d}: missing target\n", .{line_num});
            continue;
        };
        const value_tok = toks.next() orelse {
            std.debug.print("Warning: fix override line {d}: missing value\n", .{line_num});
            continue;
        };
        if (toks.next() != null) {
            std.debug.print("Warning: fix override line {d}: unexpected extra tokens\n", .{line_num});
            continue;
        }
        if (comp_tok.len == 0 or comp_tok.len > 3) {
            std.debug.print("Warning: fix override line {d}: invalid comp_id '{s}'\n", .{ line_num, comp_tok });
            continue;
        }

        const selector = parseSelector(allocator, selector_tok) catch {
            std.debug.print("Warning: fix override line {d}: invalid selector '{s}'\n", .{ line_num, selector_tok });
            continue;
        };
        errdefer allocator.free(selector.chain_id);
        const target = try allocator.dupe(u8, target_tok);
        errdefer allocator.free(target);
        const value = parseValue(target_tok, value_tok) catch {
            std.debug.print("Warning: fix override line {d}: invalid value '{s}' for target '{s}'\n", .{ line_num, value_tok, target_tok });
            continue;
        };

        var comp_id: [3]u8 = .{ ' ', ' ', ' ' };
        const comp_id_len: u2 = @intCast(comp_tok.len);
        for (0..comp_id_len) |i| comp_id[i] = std.ascii.toUpper(comp_tok[i]);

        try entries.append(allocator, .{
            .selector = selector,
            .comp_id = comp_id,
            .comp_id_len = comp_id_len,
            .target = target,
            .value = value,
        });
    }

    return .{ .allocator = allocator, .entries = try entries.toOwnedSlice(allocator) };
}

pub fn applyFixes(overrides: *const FixOverrides, mdl: *const Model, movers: []Mover) void {
    for (movers) |*m| {
        const res = mdl.residues.items[m.residue_idx];
        var i = overrides.entries.len;
        while (i > 0) {
            i -= 1;
            const entry = &overrides.entries[i];
            if (!entry.selector.matches(mdl, m)) continue;
            if (!std.mem.eql(u8, entry.compIdSlice(), res.compIdSlice())) continue;
            if (!std.mem.eql(u8, moverTargetName(mdl, m), entry.target)) continue;
            if (entryOrientation(m, entry.value)) |idx| {
                m.lockToOrientation(idx);
            } else {
                std.debug.print("  warning: fix override for {s}:{d} {s} {s} is not valid for this mover\n", .{
                    entry.selector.chain_id,
                    entry.selector.auth_seq_id,
                    entry.compIdSlice(),
                    entry.target,
                });
            }
            break;
        }
    }
}

pub fn dumpMovers(writer: anytype, mdl: *const Model, movers: []const Mover) !void {
    try writer.writeAll(
        "# selector comp_id target current states\n",
    );
    for (movers) |*m| {
        const res = mdl.residues.items[m.residue_idx];
        const chain = mdl.chains.items[res.chain_idx];
        const target = moverTargetName(mdl, m);
        try writer.print("{s}:{d}", .{ chain.authSlice(), res.auth_seq_id });
        if (res.ins_code != ' ') try writer.print(":{c}", .{res.ins_code});
        try writer.print(" {s} {s} {d} ", .{ res.compIdSlice(), target, m.best_orientation });
        try dumpStates(writer, m);
        try writer.writeByte('\n');
    }
}

fn dumpStates(writer: anytype, m: *const Mover) !void {
    switch (m.kind) {
        .amide_flip => try writer.writeAll("ORIGINAL|FLIP"),
        .his_flip => try writer.writeAll("HIE|HID|HIE_FLIP|HID_FLIP"),
        else => try writer.print("0..{d}", .{m.nOrientations() - 1}),
    }
}

pub fn moverTargetName(mdl: *const Model, m: *const Mover) []const u8 {
    return switch (m.kind) {
        .amide_flip => "amide",
        .his_flip => "his",
        .single_h_rotator, .nh3_rotator, .methyl_rotator, .aromatic_methyl => blk: {
            const center_idx = m.center_idx orelse break :blk "unknown";
            break :blk mdl.atoms.items[center_idx].nameSlice();
        },
    };
}

fn entryOrientation(m: *const Mover, value: Value) ?u16 {
    return switch (value) {
        .orientation => |idx| if (idx < m.nOrientations()) idx else null,
        .amide => |amide| switch (m.kind) {
            .amide_flip => switch (amide) {
                .original => 0,
                .flip => 1,
            },
            else => null,
        },
        .his => |his| switch (m.kind) {
            .his_flip => switch (his) {
                .hie => 0,
                .hid => 1,
                .hie_flip => 2,
                .hid_flip => 3,
            },
            else => null,
        },
    };
}

fn parseSelector(allocator: std.mem.Allocator, token: []const u8) !Selector {
    var parts = std.mem.splitScalar(u8, token, ':');
    const chain_tok = parts.next() orelse return error.InvalidFixOverride;
    const seq_tok = parts.next() orelse return error.InvalidFixOverride;
    const ins_tok = parts.next();
    if (parts.next() != null) return error.InvalidFixOverride;
    const auth_seq_id = std.fmt.parseInt(i32, seq_tok, 10) catch return error.InvalidFixOverride;
    var ins_code: u8 = ' ';
    if (ins_tok) |tok| {
        if (tok.len == 1) ins_code = tok[0] else if (!std.mem.eql(u8, tok, ".")) return error.InvalidFixOverride;
    }
    return .{
        .chain_id = try allocator.dupe(u8, chain_tok),
        .auth_seq_id = auth_seq_id,
        .ins_code = ins_code,
    };
}

fn parseValue(target_tok: []const u8, value_tok: []const u8) !Value {
    if (std.ascii.eqlIgnoreCase(target_tok, "amide")) {
        if (std.ascii.eqlIgnoreCase(value_tok, "ORIGINAL")) return .{ .amide = .original };
        if (std.ascii.eqlIgnoreCase(value_tok, "FLIP")) return .{ .amide = .flip };
        return error.InvalidFixOverride;
    }
    if (std.ascii.eqlIgnoreCase(target_tok, "his")) {
        if (std.ascii.eqlIgnoreCase(value_tok, "HIE")) return .{ .his = .hie };
        if (std.ascii.eqlIgnoreCase(value_tok, "HID")) return .{ .his = .hid };
        if (std.ascii.eqlIgnoreCase(value_tok, "HIE_FLIP")) return .{ .his = .hie_flip };
        if (std.ascii.eqlIgnoreCase(value_tok, "HID_FLIP")) return .{ .his = .hid_flip };
        return error.InvalidFixOverride;
    }
    const idx = std.fmt.parseInt(u16, value_tok, 10) catch return error.InvalidFixOverride;
    return .{ .orientation = idx };
}

const testing = std.testing;
const mmcif = @import("../mmcif.zig");
const place = @import("../place.zig");
const optimize = @import("optimize.zig");

test "parse fix override file" {
    var overrides = try parseString(testing.allocator,
        \\A:1 ASN amide FLIP
        \\A:2 HIS his HID_FLIP
        \\B:9 SER OG 6
    );
    defer overrides.deinit();

    try testing.expectEqual(@as(usize, 3), overrides.entries.len);
    try testing.expectEqualStrings("amide", overrides.entries[0].target);
    try testing.expect(overrides.entries[0].value == .amide);
    try testing.expectEqual(@as(u16, 6), overrides.entries[2].value.orientation);
}

test "invalid fix lines are ignored" {
    var overrides = try parseString(testing.allocator,
        \\A:1 ASN amide FLIP
        \\bad line
        \\A:2 HIST his HIE
    );
    defer overrides.deinit();

    try testing.expectEqual(@as(usize, 1), overrides.entries.len);
}

test "apply fix locks his mover state" {
    const source = @embedFile("../test_data/his.cif");
    var mdl = try mmcif.parseModel(testing.allocator, source);
    defer mdl.deinit();
    mdl.residues.items[0].auth_seq_id = 1;

    place.applyChemistry(&mdl);
    _ = try place.addHydrogens(&mdl, null, null);
    const gen = try optimize.generateMovers(testing.allocator, &mdl, false, null, null, null);
    const movers = gen.movers;
    defer {
        for (movers) |*m| m.deinit();
        testing.allocator.free(movers);
    }

    var overrides = try parseString(testing.allocator,
        \\A:1 HIS his HID_FLIP
    );
    defer overrides.deinit();
    applyFixes(&overrides, &mdl, movers);

    try testing.expectEqual(@as(usize, 1), movers.len);
    try testing.expect(movers[0].is_fixed);
    try testing.expectEqual(@as(u16, 3), movers[0].best_orientation);
}

test "apply fix locks amide and rotator movers" {
    {
        const source = @embedFile("../test_data/asn.cif");
        var mdl = try mmcif.parseModel(testing.allocator, source);
        defer mdl.deinit();
        mdl.residues.items[0].auth_seq_id = 1;

        place.applyChemistry(&mdl);
        _ = try place.addHydrogens(&mdl, null, null);
        const gen = try optimize.generateMovers(testing.allocator, &mdl, false, null, null, null);
        const movers = gen.movers;
        defer {
            for (movers) |*m| m.deinit();
            testing.allocator.free(movers);
        }

        var overrides = try parseString(testing.allocator,
            \\A:1 ASN amide FLIP
        );
        defer overrides.deinit();
        applyFixes(&overrides, &mdl, movers);

        try testing.expect(movers[0].is_fixed);
        try testing.expectEqual(@as(u16, 1), movers[0].best_orientation);
    }

    {
        const source = @embedFile("../test_data/tiny.cif");
        var mdl = try mmcif.parseModel(testing.allocator, source);
        defer mdl.deinit();
        mdl.residues.items[0].auth_seq_id = 1;

        place.applyChemistry(&mdl);
        _ = try place.addHydrogens(&mdl, null, null);
        const gen = try optimize.generateMovers(testing.allocator, &mdl, false, null, null, null);
        const movers = gen.movers;
        defer {
            for (movers) |*m| m.deinit();
            testing.allocator.free(movers);
        }

        var overrides = try parseString(testing.allocator,
            \\A:1 ALA CB 2
        );
        defer overrides.deinit();
        applyFixes(&overrides, &mdl, movers);

        try testing.expect(movers[0].is_fixed);
        try testing.expectEqual(@as(u16, 2), movers[0].best_orientation);
    }
}

test "dump movers lists symbolic states" {
    const source = @embedFile("../test_data/his.cif");
    var mdl = try mmcif.parseModel(testing.allocator, source);
    defer mdl.deinit();
    mdl.residues.items[0].auth_seq_id = 1;

    place.applyChemistry(&mdl);
    _ = try place.addHydrogens(&mdl, null, null);
    const gen = try optimize.generateMovers(testing.allocator, &mdl, false, null, null, null);
    const movers = gen.movers;
    defer {
        for (movers) |*m| m.deinit();
        testing.allocator.free(movers);
    }

    var buf: std.ArrayListUnmanaged(u8) = .empty;
    defer buf.deinit(testing.allocator);
    try dumpMovers(buf.writer(testing.allocator), &mdl, movers);
    try testing.expect(std.mem.indexOf(u8, buf.items, "A:1 HIS his 0 HIE|HID|HIE_FLIP|HID_FLIP") != null);
}
