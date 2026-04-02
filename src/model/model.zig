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
    /// Number of atoms present before hydrogen placement. Atoms at indices
    /// [0, original_atom_count) are original; atoms at [original_atom_count, len)
    /// are added hydrogens that carry a `residue_idx` back-pointer.
    original_atom_count: u32 = 0,

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

    /// Remove all hydrogen atoms from the model and rebuild indices.
    /// Returns the number of hydrogen atoms removed.
    pub fn stripHydrogens(self: *Model) u32 {
        var write: u32 = 0;
        var removed: u32 = 0;

        // Compact atoms, skipping hydrogens, and update residue atom_start/atom_end.
        for (self.residues.items) |*res| {
            const new_start = write;
            var read = res.atom_start;
            while (read < res.atom_end) : (read += 1) {
                if (self.atoms.items[read].is_hydrogen) {
                    removed += 1;
                } else {
                    if (write != read) {
                        self.atoms.items[write] = self.atoms.items[read];
                    }
                    write += 1;
                }
            }
            res.atom_start = new_start;
            res.atom_end = write;
        }

        self.atoms.items.len = write;

        // Bonds reference atom indices — remove bonds involving removed atoms.
        // Since we compacted, we need to rebuild bond indices. For simplicity,
        // clear all bonds (they will be re-parsed from struct_conn after strip).
        self.bonds.items.len = 0;

        return removed;
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

test "stripHydrogens removes H and rebuilds indices" {
    const allocator = std.testing.allocator;
    var m = Model.init(allocator);
    defer m.deinit();

    // Residue 0: N, H, CA (H should be removed)
    var a_n = Atom{ .pos = .{ .x = 1, .y = 0, .z = 0 } };
    a_n.setName("N");
    var a_h = Atom{ .pos = .{ .x = 1.5, .y = 0, .z = 0 }, .is_hydrogen = true };
    a_h.setName("H");
    var a_ca = Atom{ .pos = .{ .x = 2, .y = 0, .z = 0 } };
    a_ca.setName("CA");
    try m.atoms.append(m.allocator, a_n);
    try m.atoms.append(m.allocator, a_h);
    try m.atoms.append(m.allocator, a_ca);

    // Residue 1: C, O (no H)
    var a_c = Atom{ .pos = .{ .x = 3, .y = 0, .z = 0 } };
    a_c.setName("C");
    var a_o = Atom{ .pos = .{ .x = 4, .y = 0, .z = 0 } };
    a_o.setName("O");
    try m.atoms.append(m.allocator, a_c);
    try m.atoms.append(m.allocator, a_o);

    var r0 = Residue{ .atom_start = 0, .atom_end = 3 };
    r0.setCompId("ALA");
    var r1 = Residue{ .atom_start = 3, .atom_end = 5 };
    r1.setCompId("GLY");
    try m.residues.append(m.allocator, r0);
    try m.residues.append(m.allocator, r1);

    const removed = m.stripHydrogens();
    try std.testing.expectEqual(@as(u32, 1), removed);
    try std.testing.expectEqual(@as(usize, 4), m.atoms.items.len);

    // Residue 0: N, CA (H removed)
    try std.testing.expectEqual(@as(u32, 0), m.residues.items[0].atom_start);
    try std.testing.expectEqual(@as(u32, 2), m.residues.items[0].atom_end);
    try std.testing.expectEqualStrings("N", m.atoms.items[0].nameSlice());
    try std.testing.expectEqualStrings("CA", m.atoms.items[1].nameSlice());

    // Residue 1: C, O (unchanged content, updated indices)
    try std.testing.expectEqual(@as(u32, 2), m.residues.items[1].atom_start);
    try std.testing.expectEqual(@as(u32, 4), m.residues.items[1].atom_end);
    try std.testing.expectEqualStrings("C", m.atoms.items[2].nameSlice());
    try std.testing.expectEqualStrings("O", m.atoms.items[3].nameSlice());
}

test "stripHydrogens handles all-hydrogen residue" {
    const allocator = std.testing.allocator;
    var m = Model.init(allocator);
    defer m.deinit();

    var a_h1 = Atom{ .pos = .{ .x = 1, .y = 0, .z = 0 }, .is_hydrogen = true };
    a_h1.setName("H1");
    var a_h2 = Atom{ .pos = .{ .x = 2, .y = 0, .z = 0 }, .is_hydrogen = true };
    a_h2.setName("H2");
    try m.atoms.append(m.allocator, a_h1);
    try m.atoms.append(m.allocator, a_h2);

    var r0 = Residue{ .atom_start = 0, .atom_end = 2 };
    r0.setCompId("HOH");
    try m.residues.append(m.allocator, r0);

    const removed = m.stripHydrogens();
    try std.testing.expectEqual(@as(u32, 2), removed);
    try std.testing.expectEqual(@as(usize, 0), m.atoms.items.len);
    // Empty residue: atom_start == atom_end
    try std.testing.expectEqual(m.residues.items[0].atom_start, m.residues.items[0].atom_end);
}

test "stripHydrogens on empty model" {
    const allocator = std.testing.allocator;
    var m = Model.init(allocator);
    defer m.deinit();

    const removed = m.stripHydrogens();
    try std.testing.expectEqual(@as(u32, 0), removed);
    try std.testing.expectEqual(@as(usize, 0), m.atoms.items.len);
}
