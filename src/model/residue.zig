//! Residue struct representing a single residue in the molecular model.

const std = @import("std");

pub const EntityType = enum { polymer, non_polymer, water, unknown };

pub const Residue = struct {
    comp_id: [3]u8 = .{ ' ', ' ', ' ' },
    comp_id_len: u3 = 0,
    chain_idx: u16 = 0,
    seq_id: i32 = 0,
    ins_code: u8 = ' ',
    atom_start: u32 = 0,
    atom_end: u32 = 0,
    entity_type: EntityType = .unknown,

    pub fn compIdSlice(self: *const Residue) []const u8 {
        return self.comp_id[0..self.comp_id_len];
    }

    pub fn setCompId(self: *Residue, id: []const u8) void {
        const len: u3 = @intCast(@min(id.len, 3));
        self.comp_id = .{ ' ', ' ', ' ' };
        for (0..len) |i| self.comp_id[i] = id[i];
        self.comp_id_len = len;
    }

    pub fn atomCount(self: *const Residue) u32 {
        return self.atom_end - self.atom_start;
    }
};

// ── Tests ─────────────────────────────────────────────────────────────────────

test "Residue setCompId and compIdSlice" {
    var r = Residue{};
    r.setCompId("ALA");
    try std.testing.expectEqualStrings("ALA", r.compIdSlice());
    try std.testing.expectEqual(@as(u3, 3), r.comp_id_len);
}

test "Residue atomCount" {
    const r = Residue{ .atom_start = 10, .atom_end = 18 };
    try std.testing.expectEqual(@as(u32, 8), r.atomCount());
}
