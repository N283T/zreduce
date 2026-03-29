//! Aggregate model containing all molecular data.

const std = @import("std");
const Allocator = std.mem.Allocator;

const atom_mod = @import("atom.zig");
const residue_mod = @import("residue.zig");
const chain_mod = @import("chain.zig");
const bond_mod = @import("bond.zig");

pub const Atom = atom_mod.Atom;
pub const Residue = residue_mod.Residue;
pub const Chain = chain_mod.Chain;
pub const Bond = bond_mod.Bond;

pub const Model = struct {
    atoms: std.ArrayListUnmanaged(Atom),
    residues: std.ArrayListUnmanaged(Residue),
    chains: std.ArrayListUnmanaged(Chain),
    bonds: std.ArrayListUnmanaged(Bond),
    allocator: Allocator,
    n_unobs_atoms: u32 = 0,

    pub fn init(allocator: Allocator) Model {
        return .{
            .atoms = .empty,
            .residues = .empty,
            .chains = .empty,
            .bonds = .empty,
            .allocator = allocator,
        };
    }

    pub fn deinit(self: *Model) void {
        self.atoms.deinit(self.allocator);
        self.residues.deinit(self.allocator);
        self.chains.deinit(self.allocator);
        self.bonds.deinit(self.allocator);
    }

    /// Return a slice of atoms belonging to the given residue.
    pub fn residueAtoms(self: *const Model, res: Residue) []const Atom {
        return self.atoms.items[res.atom_start..res.atom_end];
    }

    /// Return a mutable slice of atoms belonging to the given residue.
    pub fn residueAtomsMut(self: *Model, res: Residue) []Atom {
        return self.atoms.items[res.atom_start..res.atom_end];
    }

    /// Find the index of an atom by name within a residue. Returns null if not found.
    pub fn findAtomInResidue(self: *const Model, res: Residue, name: []const u8) ?u32 {
        const atoms = self.residueAtoms(res);
        for (atoms, 0..) |*a, i| {
            if (std.mem.eql(u8, a.nameSlice(), name)) {
                return @intCast(res.atom_start + i);
            }
        }
        return null;
    }
};

// ── Tests ─────────────────────────────────────────────────────────────────────

test "Model init and deinit" {
    const allocator = std.testing.allocator;
    var m = Model.init(allocator);
    defer m.deinit();
    try std.testing.expectEqual(@as(usize, 0), m.atoms.items.len);
}

test "Model residueAtoms" {
    const allocator = std.testing.allocator;
    var m = Model.init(allocator);
    defer m.deinit();

    // Add two atoms
    var a1 = Atom{ .pos = .{ .x = 1, .y = 0, .z = 0 } };
    a1.setName("N");
    var a2 = Atom{ .pos = .{ .x = 2, .y = 0, .z = 0 } };
    a2.setName("CA");
    try m.atoms.append(m.allocator, a1);
    try m.atoms.append(m.allocator, a2);

    const res = Residue{ .atom_start = 0, .atom_end = 2 };
    const atoms = m.residueAtoms(res);
    try std.testing.expectEqual(@as(usize, 2), atoms.len);
}

test "Model findAtomInResidue" {
    const allocator = std.testing.allocator;
    var m = Model.init(allocator);
    defer m.deinit();

    var a1 = Atom{ .pos = .{ .x = 1, .y = 0, .z = 0 } };
    a1.setName("N");
    var a2 = Atom{ .pos = .{ .x = 2, .y = 0, .z = 0 } };
    a2.setName("CA");
    try m.atoms.append(m.allocator, a1);
    try m.atoms.append(m.allocator, a2);

    const res = Residue{ .atom_start = 0, .atom_end = 2 };
    const idx = m.findAtomInResidue(res, "CA");
    try std.testing.expectEqual(@as(?u32, 1), idx);
    try std.testing.expectEqual(@as(?u32, null), m.findAtomInResidue(res, "CB"));
}
