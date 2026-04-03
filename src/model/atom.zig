//! Atom struct representing a single atom in the molecular model.

const math = @import("../math.zig");
const element = @import("../element.zig");
const standard = @import("../place/standard.zig");
const fixed_string = @import("fixed_string.zig");

pub const AtomFlags = element.AtomFlags;
pub const FixedString = fixed_string.FixedString;

pub const Atom = struct {
    pos: math.Vec3(f32),
    name: FixedString(4) = .{},
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
    mover_hint: standard.MoverHint = .none,

    pub fn nameSlice(self: *const Atom) []const u8 {
        return self.name.slice();
    }

    pub fn setName(self: *Atom, n: []const u8) void {
        self.name.set(n);
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
    try std.testing.expectEqual(@as(u3, 2), a.name.len);
}

test "Atom setName truncates at 4" {
    var a = Atom{
        .pos = .{ .x = 0, .y = 0, .z = 0 },
    };
    a.setName("ABCDE");
    try std.testing.expectEqualStrings("ABCD", a.nameSlice());
    try std.testing.expectEqual(@as(u3, 4), a.name.len);
}
