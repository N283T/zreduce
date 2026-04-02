//! Residue struct representing a single residue in the molecular model.

const std = @import("std");
const fixed_string = @import("fixed_string.zig");

const FixedString = fixed_string.FixedString;

pub const EntityType = enum { polymer, non_polymer, branched, water, unknown };

pub const Residue = struct {
    comp_id: FixedString(5) = .{},
    chain_idx: u32 = 0,
    seq_id: i32 = 0,
    auth_seq_id: i32 = 0,
    ins_code: u8 = ' ',
    atom_start: u32 = 0,
    atom_end: u32 = 0,
    entity_type: EntityType = .unknown,
    is_chain_break_before: bool = false,

    pub fn compIdSlice(self: *const Residue) []const u8 {
        return self.comp_id.slice();
    }

    pub fn setCompId(self: *Residue, id: []const u8) void {
        self.comp_id.set(id);
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
    try std.testing.expectEqual(@as(u3, 3), r.comp_id.len);
    r.setCompId("ATP");
    try std.testing.expectEqualStrings("ATP", r.compIdSlice());
    r.setCompId("BGLA");
    try std.testing.expectEqualStrings("BGLA", r.compIdSlice());
    try std.testing.expectEqual(@as(u3, 4), r.comp_id.len);
    r.setCompId("BGLAN");
    try std.testing.expectEqualStrings("BGLAN", r.compIdSlice());
    try std.testing.expectEqual(@as(u3, 5), r.comp_id.len);
}

test "Residue atomCount" {
    const r = Residue{ .atom_start = 10, .atom_end = 18 };
    try std.testing.expectEqual(@as(u32, 8), r.atomCount());
}

test "EntityType has branched variant" {
    const e: EntityType = .branched;
    try std.testing.expect(e == .branched);
}

test "Residue auth_seq_id default" {
    const r = Residue{};
    try std.testing.expectEqual(@as(i32, 0), r.auth_seq_id);
}
