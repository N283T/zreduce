//! Hardcoded heavy-atom bond lists for the 20 standard amino acids.
//!
//! Atom names follow PDB convention (4-character, space-padded).
//! Only heavy (non-hydrogen) atoms are included.

const std = @import("std");

// ---------------------------------------------------------------------------
// Public types
// ---------------------------------------------------------------------------

pub const BondEntry = struct {
    a1: [4]u8,
    a2: [4]u8,
};

// ---------------------------------------------------------------------------
// Helpers
// ---------------------------------------------------------------------------

const n = @import("lookup.zig").padName;

fn b(comptime a: []const u8, comptime bname: []const u8) BondEntry {
    return .{ .a1 = n(a), .a2 = n(bname) };
}

// ---------------------------------------------------------------------------
// Backbone bonds (shared by all standard residues)
// ---------------------------------------------------------------------------

const backbone_bonds = [_]BondEntry{
    b(" N  ", " CA "),
    b(" CA ", " C  "),
    b(" C  ", " O  "),
};

// ---------------------------------------------------------------------------
// ALA — 4 bonds (3 backbone + CA-CB)
// ---------------------------------------------------------------------------

const ala_bonds = backbone_bonds ++ [_]BondEntry{
    b(" CA ", " CB "),
};

// ---------------------------------------------------------------------------
// GLY — 3 bonds (backbone only)
// ---------------------------------------------------------------------------

const gly_bonds = backbone_bonds ++ [_]BondEntry{};

// ---------------------------------------------------------------------------
// VAL — 6 bonds
// ---------------------------------------------------------------------------

const val_bonds = backbone_bonds ++ [_]BondEntry{
    b(" CA ", " CB "),
    b(" CB ", "CG1 "),
    b(" CB ", "CG2 "),
};

// ---------------------------------------------------------------------------
// LEU — 7 bonds
// ---------------------------------------------------------------------------

const leu_bonds = backbone_bonds ++ [_]BondEntry{
    b(" CA ", " CB "),
    b(" CB ", " CG "),
    b(" CG ", "CD1 "),
    b(" CG ", "CD2 "),
};

// ---------------------------------------------------------------------------
// ILE — 7 bonds
// ---------------------------------------------------------------------------

const ile_bonds = backbone_bonds ++ [_]BondEntry{
    b(" CA ", " CB "),
    b(" CB ", "CG1 "),
    b(" CB ", "CG2 "),
    b("CG1 ", "CD1 "),
};

// ---------------------------------------------------------------------------
// PRO — 7 bonds (ring: CD-N closes the ring)
// ---------------------------------------------------------------------------

const pro_bonds = backbone_bonds ++ [_]BondEntry{
    b(" CA ", " CB "),
    b(" CB ", " CG "),
    b(" CG ", " CD "),
    b(" CD ", " N  "),
};

// ---------------------------------------------------------------------------
// PHE — 11 bonds
// ---------------------------------------------------------------------------

const phe_bonds = backbone_bonds ++ [_]BondEntry{
    b(" CA ", " CB "),
    b(" CB ", " CG "),
    b(" CG ", "CD1 "),
    b(" CG ", "CD2 "),
    b("CD1 ", "CE1 "),
    b("CD2 ", "CE2 "),
    b("CE1 ", " CZ "),
    b("CE2 ", " CZ "),
};

// ---------------------------------------------------------------------------
// TYR — 12 bonds (PHE + CZ-OH)
// ---------------------------------------------------------------------------

const tyr_bonds = backbone_bonds ++ [_]BondEntry{
    b(" CA ", " CB "),
    b(" CB ", " CG "),
    b(" CG ", "CD1 "),
    b(" CG ", "CD2 "),
    b("CD1 ", "CE1 "),
    b("CD2 ", "CE2 "),
    b("CE1 ", " CZ "),
    b("CE2 ", " CZ "),
    b(" CZ ", " OH "),
};

// ---------------------------------------------------------------------------
// TRP — 15 bonds (indole bicyclic system)
// 5-ring: CG-CD1-NE1-CE2-CD2 (with CG-CD2 closing)
// 6-ring: CD2-CE3-CZ3-CH2-CZ2-CE2 (CE2-CD2 already in 5-ring)
// ---------------------------------------------------------------------------

const trp_bonds = backbone_bonds ++ [_]BondEntry{
    b(" CA ", " CB "),
    b(" CB ", " CG "),
    b(" CG ", "CD1 "),
    b(" CG ", "CD2 "),
    b("CD1 ", "NE1 "),
    b("NE1 ", "CE2 "),
    b("CE2 ", "CD2 "),
    b("CE2 ", "CZ2 "),
    b("CD2 ", "CE3 "),
    b("CE3 ", "CZ3 "),
    b("CZ3 ", "CH2 "),
    b("CH2 ", "CZ2 "),
};

// ---------------------------------------------------------------------------
// SER — 5 bonds
// ---------------------------------------------------------------------------

const ser_bonds = backbone_bonds ++ [_]BondEntry{
    b(" CA ", " CB "),
    b(" CB ", " OG "),
};

// ---------------------------------------------------------------------------
// THR — 6 bonds
// ---------------------------------------------------------------------------

const thr_bonds = backbone_bonds ++ [_]BondEntry{
    b(" CA ", " CB "),
    b(" CB ", "OG1 "),
    b(" CB ", "CG2 "),
};

// ---------------------------------------------------------------------------
// CYS — 5 bonds
// ---------------------------------------------------------------------------

const cys_bonds = backbone_bonds ++ [_]BondEntry{
    b(" CA ", " CB "),
    b(" CB ", " SG "),
};

// ---------------------------------------------------------------------------
// MET — 7 bonds
// ---------------------------------------------------------------------------

const met_bonds = backbone_bonds ++ [_]BondEntry{
    b(" CA ", " CB "),
    b(" CB ", " CG "),
    b(" CG ", " SD "),
    b(" SD ", " CE "),
};

// ---------------------------------------------------------------------------
// ASP — 7 bonds
// ---------------------------------------------------------------------------

const asp_bonds = backbone_bonds ++ [_]BondEntry{
    b(" CA ", " CB "),
    b(" CB ", " CG "),
    b(" CG ", "OD1 "),
    b(" CG ", "OD2 "),
};

// ---------------------------------------------------------------------------
// GLU — 8 bonds
// ---------------------------------------------------------------------------

const glu_bonds = backbone_bonds ++ [_]BondEntry{
    b(" CA ", " CB "),
    b(" CB ", " CG "),
    b(" CG ", " CD "),
    b(" CD ", "OE1 "),
    b(" CD ", "OE2 "),
};

// ---------------------------------------------------------------------------
// ASN — 7 bonds
// ---------------------------------------------------------------------------

const asn_bonds = backbone_bonds ++ [_]BondEntry{
    b(" CA ", " CB "),
    b(" CB ", " CG "),
    b(" CG ", "OD1 "),
    b(" CG ", "ND2 "),
};

// ---------------------------------------------------------------------------
// GLN — 8 bonds
// ---------------------------------------------------------------------------

const gln_bonds = backbone_bonds ++ [_]BondEntry{
    b(" CA ", " CB "),
    b(" CB ", " CG "),
    b(" CG ", " CD "),
    b(" CD ", "OE1 "),
    b(" CD ", "NE2 "),
};

// ---------------------------------------------------------------------------
// LYS — 8 bonds
// ---------------------------------------------------------------------------

const lys_bonds = backbone_bonds ++ [_]BondEntry{
    b(" CA ", " CB "),
    b(" CB ", " CG "),
    b(" CG ", " CD "),
    b(" CD ", " CE "),
    b(" CE ", " NZ "),
};

// ---------------------------------------------------------------------------
// ARG — 10 bonds
// ---------------------------------------------------------------------------

const arg_bonds = backbone_bonds ++ [_]BondEntry{
    b(" CA ", " CB "),
    b(" CB ", " CG "),
    b(" CG ", " CD "),
    b(" CD ", " NE "),
    b(" NE ", " CZ "),
    b(" CZ ", "NH1 "),
    b(" CZ ", "NH2 "),
};

// ---------------------------------------------------------------------------
// HIS — 10 bonds (imidazole ring)
// ---------------------------------------------------------------------------

const his_bonds = backbone_bonds ++ [_]BondEntry{
    b(" CA ", " CB "),
    b(" CB ", " CG "),
    b(" CG ", "ND1 "),
    b(" CG ", "CD2 "),
    b("ND1 ", "CE1 "),
    b("CD2 ", "NE2 "),
    b("CE1 ", "NE2 "),
};

// ---------------------------------------------------------------------------
// Lookup
// ---------------------------------------------------------------------------

pub fn getBonds(comp_id: []const u8) ?[]const BondEntry {
    const map = std.StaticStringMap([]const BondEntry).initComptime(.{
        .{ "ALA", &ala_bonds },
        .{ "GLY", &gly_bonds },
        .{ "VAL", &val_bonds },
        .{ "LEU", &leu_bonds },
        .{ "ILE", &ile_bonds },
        .{ "PRO", &pro_bonds },
        .{ "PHE", &phe_bonds },
        .{ "TYR", &tyr_bonds },
        .{ "TRP", &trp_bonds },
        .{ "SER", &ser_bonds },
        .{ "THR", &thr_bonds },
        .{ "CYS", &cys_bonds },
        .{ "MET", &met_bonds },
        .{ "ASP", &asp_bonds },
        .{ "GLU", &glu_bonds },
        .{ "ASN", &asn_bonds },
        .{ "GLN", &gln_bonds },
        .{ "LYS", &lys_bonds },
        .{ "ARG", &arg_bonds },
        .{ "HIS", &his_bonds },
    });
    return map.get(comp_id);
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

test "ALA has 4 heavy-atom bonds" {
    const bonds = getBonds("ALA");
    try std.testing.expect(bonds != null);
    try std.testing.expectEqual(@as(usize, 4), bonds.?.len);
}

test "GLY has 3 heavy-atom bonds" {
    const bonds = getBonds("GLY");
    try std.testing.expect(bonds != null);
    try std.testing.expectEqual(@as(usize, 3), bonds.?.len);
}

test "VAL has 6 heavy-atom bonds" {
    const bonds = getBonds("VAL");
    try std.testing.expect(bonds != null);
    try std.testing.expectEqual(@as(usize, 6), bonds.?.len);
}

test "LEU has 7 heavy-atom bonds" {
    const bonds = getBonds("LEU");
    try std.testing.expect(bonds != null);
    try std.testing.expectEqual(@as(usize, 7), bonds.?.len);
}

test "ILE has 7 heavy-atom bonds" {
    const bonds = getBonds("ILE");
    try std.testing.expect(bonds != null);
    try std.testing.expectEqual(@as(usize, 7), bonds.?.len);
}

test "PRO has 7 heavy-atom bonds" {
    const bonds = getBonds("PRO");
    try std.testing.expect(bonds != null);
    try std.testing.expectEqual(@as(usize, 7), bonds.?.len);
}

test "PHE has 11 heavy-atom bonds" {
    const bonds = getBonds("PHE");
    try std.testing.expect(bonds != null);
    try std.testing.expectEqual(@as(usize, 11), bonds.?.len);
}

test "TYR has 12 heavy-atom bonds" {
    const bonds = getBonds("TYR");
    try std.testing.expect(bonds != null);
    try std.testing.expectEqual(@as(usize, 12), bonds.?.len);
}

test "TRP has 15 heavy-atom bonds" {
    const bonds = getBonds("TRP");
    try std.testing.expect(bonds != null);
    try std.testing.expectEqual(@as(usize, 15), bonds.?.len);
}

test "SER has 5 heavy-atom bonds" {
    const bonds = getBonds("SER");
    try std.testing.expect(bonds != null);
    try std.testing.expectEqual(@as(usize, 5), bonds.?.len);
}

test "THR has 6 heavy-atom bonds" {
    const bonds = getBonds("THR");
    try std.testing.expect(bonds != null);
    try std.testing.expectEqual(@as(usize, 6), bonds.?.len);
}

test "CYS has 5 heavy-atom bonds" {
    const bonds = getBonds("CYS");
    try std.testing.expect(bonds != null);
    try std.testing.expectEqual(@as(usize, 5), bonds.?.len);
}

test "MET has 7 heavy-atom bonds" {
    const bonds = getBonds("MET");
    try std.testing.expect(bonds != null);
    try std.testing.expectEqual(@as(usize, 7), bonds.?.len);
}

test "ASP has 7 heavy-atom bonds" {
    const bonds = getBonds("ASP");
    try std.testing.expect(bonds != null);
    try std.testing.expectEqual(@as(usize, 7), bonds.?.len);
}

test "GLU has 8 heavy-atom bonds" {
    const bonds = getBonds("GLU");
    try std.testing.expect(bonds != null);
    try std.testing.expectEqual(@as(usize, 8), bonds.?.len);
}

test "ASN has 7 heavy-atom bonds" {
    const bonds = getBonds("ASN");
    try std.testing.expect(bonds != null);
    try std.testing.expectEqual(@as(usize, 7), bonds.?.len);
}

test "GLN has 8 heavy-atom bonds" {
    const bonds = getBonds("GLN");
    try std.testing.expect(bonds != null);
    try std.testing.expectEqual(@as(usize, 8), bonds.?.len);
}

test "LYS has 8 heavy-atom bonds" {
    const bonds = getBonds("LYS");
    try std.testing.expect(bonds != null);
    try std.testing.expectEqual(@as(usize, 8), bonds.?.len);
}

test "ARG has 10 heavy-atom bonds" {
    const bonds = getBonds("ARG");
    try std.testing.expect(bonds != null);
    try std.testing.expectEqual(@as(usize, 10), bonds.?.len);
}

test "HIS has 10 heavy-atom bonds" {
    const bonds = getBonds("HIS");
    try std.testing.expect(bonds != null);
    try std.testing.expectEqual(@as(usize, 10), bonds.?.len);
}

test "unknown residue returns null" {
    try std.testing.expect(getBonds("XYZ") == null);
}

test "all 20 amino acids have bonds" {
    const aas = [_][]const u8{
        "ALA", "GLY", "VAL", "LEU", "ILE", "PRO", "PHE", "TYR",
        "TRP", "SER", "THR", "CYS", "MET", "ASP", "GLU", "ASN",
        "GLN", "LYS", "ARG", "HIS",
    };
    for (aas) |aa| {
        try std.testing.expect(getBonds(aa) != null);
    }
}
