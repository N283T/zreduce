//! Chain struct representing a single chain in the molecular model.

const std = @import("std");
const fixed_string = @import("fixed_string.zig");

const FixedString = fixed_string.FixedString;

pub const Chain = struct {
    label_asym_id: FixedString(4) = .{},
    auth_asym_id: FixedString(4) = .{},
    entity_id: FixedString(4) = .{},
    residue_start: u32 = 0,
    residue_end: u32 = 0,

    pub fn labelSlice(self: *const Chain) []const u8 {
        return self.label_asym_id.slice();
    }

    pub fn setLabelAsymId(self: *Chain, id: []const u8) void {
        self.label_asym_id.set(id);
    }

    pub fn authSlice(self: *const Chain) []const u8 {
        return self.auth_asym_id.slice();
    }

    pub fn setAuthAsymId(self: *Chain, id: []const u8) void {
        self.auth_asym_id.set(id);
    }

    pub fn entityIdSlice(self: *const Chain) []const u8 {
        return self.entity_id.slice();
    }

    pub fn setEntityId(self: *Chain, id: []const u8) void {
        self.entity_id.set(id);
    }

    pub fn residueCount(self: *const Chain) u32 {
        return self.residue_end - self.residue_start;
    }
};

// ── Tests ─────────────────────────────────────────────────────────────────────

test "Chain labelSlice" {
    var c = Chain{};
    c.setLabelAsymId("A");
    try std.testing.expectEqualStrings("A", c.labelSlice());
}

test "Chain setAuthAsymId and authSlice" {
    var c = Chain{};
    c.setAuthAsymId("AA");
    try std.testing.expectEqualStrings("AA", c.authSlice());
}

test "Chain entity_id string" {
    var c = Chain{};
    c.setEntityId("2");
    try std.testing.expectEqualStrings("2", c.entityIdSlice());
}
