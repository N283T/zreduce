//! Atom struct representing a single atom in the molecular model.

const math = @import("../math.zig");
const element = @import("../element.zig");

pub const AtomFlags = element.AtomFlags;

pub const Atom = struct {
    pos: math.Vec3(f32),
    name: [4]u8 = .{ ' ', ' ', ' ', ' ' },
    name_len: u4 = 0,
    element_type: element.AtomType = .unknown,
    residue_idx: u32 = 0,
    altloc: u8 = ' ',
    occupancy: f32 = 1.0,
    b_factor: f32 = 0.0,
    is_hydrogen: bool = false,
    is_added: bool = false,
    vdw_radius: f32 = 1.70,
    flags: AtomFlags = .{},
    serial: u32 = 0,

    pub fn nameSlice(self: *const Atom) []const u8 {
        return self.name[0..self.name_len];
    }

    pub fn setName(self: *Atom, name: []const u8) void {
        const len: u4 = @intCast(@min(name.len, 4));
        self.name = .{ ' ', ' ', ' ', ' ' };
        for (0..len) |i| self.name[i] = name[i];
        self.name_len = len;
    }
};

// ── Tests ─────────────────────────────────────────────────────────────────────

const std = @import("std");

test "Atom setName and nameSlice" {
    var a = Atom{
        .pos = .{ .x = 0, .y = 0, .z = 0 },
    };
    a.setName("CA");
    try std.testing.expectEqualStrings("CA", a.nameSlice());
    try std.testing.expectEqual(@as(u4, 2), a.name_len);
}

test "Atom setName truncates at 4" {
    var a = Atom{
        .pos = .{ .x = 0, .y = 0, .z = 0 },
    };
    a.setName("ABCDE");
    try std.testing.expectEqualStrings("ABCD", a.nameSlice());
    try std.testing.expectEqual(@as(u4, 4), a.name_len);
}
