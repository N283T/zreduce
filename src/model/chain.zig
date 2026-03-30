//! Chain struct representing a single chain in the molecular model.

const std = @import("std");

pub const Chain = struct {
    label_asym_id: [4]u8 = .{ ' ', ' ', ' ', ' ' },
    label_asym_id_len: u4 = 0,
    auth_asym_id: [4]u8 = .{ ' ', ' ', ' ', ' ' },
    auth_asym_id_len: u4 = 0,
    entity_id: [4]u8 = .{ ' ', ' ', ' ', ' ' },
    entity_id_len: u4 = 0,
    residue_start: u32 = 0,
    residue_end: u32 = 0,

    pub fn labelSlice(self: *const Chain) []const u8 {
        return self.label_asym_id[0..@min(@as(usize, self.label_asym_id_len), 4)];
    }

    pub fn setLabelAsymId(self: *Chain, id: []const u8) void {
        const len: u4 = @intCast(@min(id.len, 4));
        self.label_asym_id = .{ ' ', ' ', ' ', ' ' };
        for (0..len) |i| self.label_asym_id[i] = id[i];
        self.label_asym_id_len = len;
    }

    pub fn authSlice(self: *const Chain) []const u8 {
        return self.auth_asym_id[0..@min(@as(usize, self.auth_asym_id_len), 4)];
    }

    pub fn setAuthAsymId(self: *Chain, id: []const u8) void {
        const len: u4 = @intCast(@min(id.len, 4));
        self.auth_asym_id = .{ ' ', ' ', ' ', ' ' };
        for (0..len) |i| self.auth_asym_id[i] = id[i];
        self.auth_asym_id_len = len;
    }

    pub fn entityIdSlice(self: *const Chain) []const u8 {
        return self.entity_id[0..@min(@as(usize, self.entity_id_len), 4)];
    }

    pub fn setEntityId(self: *Chain, id: []const u8) void {
        const len: u4 = @intCast(@min(id.len, 4));
        self.entity_id = .{ ' ', ' ', ' ', ' ' };
        for (0..len) |i| self.entity_id[i] = id[i];
        self.entity_id_len = len;
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
