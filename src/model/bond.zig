//! Bond struct representing a covalent bond between two atoms.

const std = @import("std");

pub const BondOrder = enum(u3) {
    single,
    double,
    triple,
    aromatic,
    delocalized,
    unknown,

    pub fn fromString(s: []const u8) BondOrder {
        if (std.mem.eql(u8, s, "SING")) return .single;
        if (std.mem.eql(u8, s, "DOUB")) return .double;
        if (std.mem.eql(u8, s, "TRIP")) return .triple;
        if (std.mem.eql(u8, s, "AROM")) return .aromatic;
        if (std.mem.eql(u8, s, "DELO")) return .delocalized;
        return .unknown;
    }
};

pub const BondSource = enum(u3) {
    component_template,
    struct_conn,
    polymer_backbone,
    branch_link,
    inferred,
};

pub const Bond = struct {
    atom_1: u32,
    atom_2: u32,
    order: BondOrder = .single,
    source: BondSource = .inferred,
};

// ── Tests ─────────────────────────────────────────────────────────────────────

test "BondOrder fromString" {
    try std.testing.expectEqual(BondOrder.single, BondOrder.fromString("SING"));
    try std.testing.expectEqual(BondOrder.double, BondOrder.fromString("DOUB"));
    try std.testing.expectEqual(BondOrder.triple, BondOrder.fromString("TRIP"));
    try std.testing.expectEqual(BondOrder.aromatic, BondOrder.fromString("AROM"));
    try std.testing.expectEqual(BondOrder.delocalized, BondOrder.fromString("DELO"));
    try std.testing.expectEqual(BondOrder.unknown, BondOrder.fromString("XXXX"));
}
