//! Hardcoded hydrogen placement plans for the 20 standard amino acids.
//!
//! Each amino acid has a list of PlacementPlan entries describing how to
//! place hydrogen atoms based on the heavy-atom geometry. Atom names follow
//! PDB convention (4-character, space-padded).

const std = @import("std");
const element = @import("../element.zig");

// ---------------------------------------------------------------------------
// Public types
// ---------------------------------------------------------------------------

pub const PlacementType = enum(u3) {
    hxr3, // tetrahedral (sp3, 3 neighbors)
    h2xr2, // two H on sp2
    h3xr, // dihedral-controlled
    hxr2_planar, // planar bisector
    hxr2_frac, // fractional angle
    hxy, // linear
};

pub const MoverHint = enum(u3) {
    none,
    rotate, // OH, SH
    rotate_nh3, // NH3+
    rotate_methyl, // CH3
    flip_amide, // Asn/Gln
    flip_his, // His
};

pub const PlacementPlan = struct {
    h_name: [4]u8,
    placement_type: PlacementType,
    connected: [3][4]u8, // reference atom names (up to 3)
    n_connected: u2, // number of reference atoms used
    bond_len: f32,
    angle: f32 = 0.0,
    dihedral: f32 = 0.0,
    fudge: f32 = 0.0,
    atom_type: element.AtomType,
    mover_hint: MoverHint = .none,
};

// ---------------------------------------------------------------------------
// Helpers
// ---------------------------------------------------------------------------

fn n(comptime s: []const u8) [4]u8 {
    var buf: [4]u8 = .{ ' ', ' ', ' ', ' ' };
    for (s, 0..) |c, i| {
        if (i >= 4) break;
        buf[i] = c;
    }
    return buf;
}

const blank = n("    ");

// Shorthand constructors for common patterns

fn hxr3(comptime h: []const u8, comptime a: []const u8, comptime b: []const u8, comptime c: []const u8) PlacementPlan {
    return .{
        .h_name = n(h),
        .placement_type = .hxr3,
        .connected = .{ n(a), n(b), n(c) },
        .n_connected = 3,
        .bond_len = 1.10,
        .atom_type = .H,
    };
}

fn h2xr2(comptime h: []const u8, comptime a: []const u8, comptime b: []const u8, dihedral: f32) PlacementPlan {
    return .{
        .h_name = n(h),
        .placement_type = .h2xr2,
        .connected = .{ n(a), n(b), blank },
        .n_connected = 2,
        .bond_len = 1.10,
        .angle = 109.5,
        .dihedral = dihedral,
        .atom_type = .H,
    };
}

fn methyl(comptime h: []const u8, comptime center: []const u8, comptime bonded: []const u8, dihedral: f32) PlacementPlan {
    return .{
        .h_name = n(h),
        .placement_type = .h3xr,
        .connected = .{ n(center), n(bonded), blank },
        .n_connected = 3,
        .bond_len = 1.10,
        .angle = 109.5,
        .dihedral = dihedral,
        .atom_type = .H,
        .mover_hint = .rotate_methyl,
    };
}

fn planarH(comptime h: []const u8, comptime a: []const u8, comptime b: []const u8) PlacementPlan {
    return .{
        .h_name = n(h),
        .placement_type = .hxr2_planar,
        .connected = .{ n(a), n(b), blank },
        .n_connected = 2,
        .bond_len = 1.10,
        .atom_type = .Har,
    };
}

fn planarPol(comptime h: []const u8, comptime a: []const u8, comptime b: []const u8, bond_len: f32) PlacementPlan {
    return .{
        .h_name = n(h),
        .placement_type = .hxr2_planar,
        .connected = .{ n(a), n(b), blank },
        .n_connected = 2,
        .bond_len = bond_len,
        .atom_type = .Hpol,
    };
}

const backbone_h = PlacementPlan{
    .h_name = n(" H  "),
    .placement_type = .h3xr,
    .connected = .{ n(" N  "), n(" CA "), n(" C  ") },
    .n_connected = 3,
    .bond_len = 1.02,
    .angle = 119.0,
    .dihedral = 180.0,
    .atom_type = .Hpol,
};

// ---------------------------------------------------------------------------
// ALA (5 H)
// ---------------------------------------------------------------------------

const ala_plans = [_]PlacementPlan{
    hxr3(" HA ", " CA ", " N  ", " C  "),
    methyl(" HB1", " CB ", " CA ", 180.0),
    methyl(" HB2", " CB ", " CA ", 60.0),
    methyl(" HB3", " CB ", " CA ", -60.0),
    backbone_h,
};

// ---------------------------------------------------------------------------
// GLY (3 H)
// ---------------------------------------------------------------------------

const gly_plans = [_]PlacementPlan{
    h2xr2(" HA2", " CA ", " N  ", 120.0),
    h2xr2(" HA3", " CA ", " N  ", -120.0),
    backbone_h,
};

// ---------------------------------------------------------------------------
// VAL (8 H)
// ---------------------------------------------------------------------------

const val_plans = [_]PlacementPlan{
    hxr3(" HA ", " CA ", " N  ", " C  "),
    hxr3(" HB ", " CB ", " CA ", "CG1 "),
    methyl("HG11", "CG1 ", " CB ", 180.0),
    methyl("HG12", "CG1 ", " CB ", 60.0),
    methyl("HG13", "CG1 ", " CB ", -60.0),
    methyl("HG21", "CG2 ", " CB ", 180.0),
    methyl("HG22", "CG2 ", " CB ", 60.0),
    methyl("HG23", "CG2 ", " CB ", -60.0),
    backbone_h,
};

// ---------------------------------------------------------------------------
// LEU (9 H)
// ---------------------------------------------------------------------------

const leu_plans = [_]PlacementPlan{
    hxr3(" HA ", " CA ", " N  ", " C  "),
    h2xr2(" HB2", " CB ", " CA ", 120.0),
    h2xr2(" HB3", " CB ", " CA ", -120.0),
    hxr3(" HG ", " CG ", " CB ", "CD1 "),
    methyl("HD11", "CD1 ", " CG ", 180.0),
    methyl("HD12", "CD1 ", " CG ", 60.0),
    methyl("HD13", "CD1 ", " CG ", -60.0),
    methyl("HD21", "CD2 ", " CG ", 180.0),
    methyl("HD22", "CD2 ", " CG ", 60.0),
    methyl("HD23", "CD2 ", " CG ", -60.0),
    backbone_h,
};

// ---------------------------------------------------------------------------
// ILE (9 H)
// ---------------------------------------------------------------------------

const ile_plans = [_]PlacementPlan{
    hxr3(" HA ", " CA ", " N  ", " C  "),
    hxr3(" HB ", " CB ", " CA ", "CG1 "),
    h2xr2("HG12", "CG1 ", " CB ", 120.0),
    h2xr2("HG13", "CG1 ", " CB ", -120.0),
    methyl("HG21", "CG2 ", " CB ", 180.0),
    methyl("HG22", "CG2 ", " CB ", 60.0),
    methyl("HG23", "CG2 ", " CB ", -60.0),
    methyl("HD11", "CD1 ", "CG1 ", 180.0),
    methyl("HD12", "CD1 ", "CG1 ", 60.0),
    methyl("HD13", "CD1 ", "CG1 ", -60.0),
    backbone_h,
};

// ---------------------------------------------------------------------------
// PRO (5 H — no backbone H)
// ---------------------------------------------------------------------------

const pro_plans = [_]PlacementPlan{
    hxr3(" HA ", " CA ", " N  ", " C  "),
    h2xr2(" HB2", " CB ", " CA ", 120.0),
    h2xr2(" HB3", " CB ", " CA ", -120.0),
    h2xr2(" HG2", " CG ", " CB ", 120.0),
    h2xr2(" HG3", " CG ", " CB ", -120.0),
    h2xr2(" HD2", " CD ", " CG ", 120.0),
    h2xr2(" HD3", " CD ", " CG ", -120.0),
};

// ---------------------------------------------------------------------------
// PHE (7 H)
// ---------------------------------------------------------------------------

const phe_plans = [_]PlacementPlan{
    hxr3(" HA ", " CA ", " N  ", " C  "),
    h2xr2(" HB2", " CB ", " CA ", 120.0),
    h2xr2(" HB3", " CB ", " CA ", -120.0),
    planarH(" HD1", " CG ", "CE1 "),
    planarH(" HD2", " CG ", "CE2 "),
    planarH(" HE1", "CD1 ", " CZ "),
    planarH(" HE2", "CD2 ", " CZ "),
    planarH(" HZ ", "CE1 ", "CE2 "),
    backbone_h,
};

// ---------------------------------------------------------------------------
// TYR (7 H + OH)
// ---------------------------------------------------------------------------

const tyr_plans = [_]PlacementPlan{
    hxr3(" HA ", " CA ", " N  ", " C  "),
    h2xr2(" HB2", " CB ", " CA ", 120.0),
    h2xr2(" HB3", " CB ", " CA ", -120.0),
    planarH(" HD1", " CG ", "CE1 "),
    planarH(" HD2", " CG ", "CE2 "),
    planarH(" HE1", "CD1 ", " CZ "),
    planarH(" HE2", "CD2 ", " CZ "),
    PlacementPlan{
        .h_name = n(" HH "),
        .placement_type = .hxy,
        .connected = .{ n(" OH "), n(" CZ "), blank },
        .n_connected = 1,
        .bond_len = 0.97,
        .atom_type = .Hpol,
        .mover_hint = .rotate,
    },
    backbone_h,
};

// ---------------------------------------------------------------------------
// TRP (8 H)
// ---------------------------------------------------------------------------

const trp_plans = [_]PlacementPlan{
    hxr3(" HA ", " CA ", " N  ", " C  "),
    h2xr2(" HB2", " CB ", " CA ", 120.0),
    h2xr2(" HB3", " CB ", " CA ", -120.0),
    planarH(" HD1", " CG ", "NE1 "),
    planarPol(" HE1", "CD1 ", "CE2 ", 1.02),
    planarH(" HE3", "CD2 ", "CZ3 "),
    planarH(" HZ2", "CE2 ", "CH2 "),
    planarH(" HZ3", "CE3 ", "CH2 "),
    planarH(" HH2", "CZ2 ", "CZ3 "),
    backbone_h,
};

// ---------------------------------------------------------------------------
// SER (3 H)
// ---------------------------------------------------------------------------

const ser_plans = [_]PlacementPlan{
    hxr3(" HA ", " CA ", " N  ", " C  "),
    h2xr2(" HB2", " CB ", " CA ", 120.0),
    h2xr2(" HB3", " CB ", " CA ", -120.0),
    PlacementPlan{
        .h_name = n(" HG "),
        .placement_type = .hxy,
        .connected = .{ n(" OG "), n(" CB "), blank },
        .n_connected = 1,
        .bond_len = 0.97,
        .atom_type = .Hpol,
        .mover_hint = .rotate,
    },
    backbone_h,
};

// ---------------------------------------------------------------------------
// THR (4 H)
// ---------------------------------------------------------------------------

const thr_plans = [_]PlacementPlan{
    hxr3(" HA ", " CA ", " N  ", " C  "),
    hxr3(" HB ", " CB ", " CA ", "OG1 "),
    PlacementPlan{
        .h_name = n(" HG1"),
        .placement_type = .hxy,
        .connected = .{ n("OG1 "), n(" CB "), blank },
        .n_connected = 1,
        .bond_len = 0.97,
        .atom_type = .Hpol,
        .mover_hint = .rotate,
    },
    methyl("HG21", "CG2 ", " CB ", 180.0),
    methyl("HG22", "CG2 ", " CB ", 60.0),
    methyl("HG23", "CG2 ", " CB ", -60.0),
    backbone_h,
};

// ---------------------------------------------------------------------------
// CYS (3 H)
// ---------------------------------------------------------------------------

const cys_plans = [_]PlacementPlan{
    hxr3(" HA ", " CA ", " N  ", " C  "),
    h2xr2(" HB2", " CB ", " CA ", 120.0),
    h2xr2(" HB3", " CB ", " CA ", -120.0),
    PlacementPlan{
        .h_name = n(" HG "),
        .placement_type = .hxy,
        .connected = .{ n(" SG "), n(" CB "), blank },
        .n_connected = 1,
        .bond_len = 1.33,
        .atom_type = .Hpol,
        .mover_hint = .rotate,
    },
    backbone_h,
};

// ---------------------------------------------------------------------------
// MET (6 H)
// ---------------------------------------------------------------------------

const met_plans = [_]PlacementPlan{
    hxr3(" HA ", " CA ", " N  ", " C  "),
    h2xr2(" HB2", " CB ", " CA ", 120.0),
    h2xr2(" HB3", " CB ", " CA ", -120.0),
    h2xr2(" HG2", " CG ", " CB ", 120.0),
    h2xr2(" HG3", " CG ", " CB ", -120.0),
    methyl(" HE1", " CE ", " SD ", 180.0),
    methyl(" HE2", " CE ", " SD ", 60.0),
    methyl(" HE3", " CE ", " SD ", -60.0),
    backbone_h,
};

// ---------------------------------------------------------------------------
// ASP (2 H)
// ---------------------------------------------------------------------------

const asp_plans = [_]PlacementPlan{
    hxr3(" HA ", " CA ", " N  ", " C  "),
    h2xr2(" HB2", " CB ", " CA ", 120.0),
    h2xr2(" HB3", " CB ", " CA ", -120.0),
    backbone_h,
};

// ---------------------------------------------------------------------------
// GLU (4 H)
// ---------------------------------------------------------------------------

const glu_plans = [_]PlacementPlan{
    hxr3(" HA ", " CA ", " N  ", " C  "),
    h2xr2(" HB2", " CB ", " CA ", 120.0),
    h2xr2(" HB3", " CB ", " CA ", -120.0),
    h2xr2(" HG2", " CG ", " CB ", 120.0),
    h2xr2(" HG3", " CG ", " CB ", -120.0),
    backbone_h,
};

// ---------------------------------------------------------------------------
// ASN (4 H)
// ---------------------------------------------------------------------------

const asn_plans = [_]PlacementPlan{
    hxr3(" HA ", " CA ", " N  ", " C  "),
    h2xr2(" HB2", " CB ", " CA ", 120.0),
    h2xr2(" HB3", " CB ", " CA ", -120.0),
    PlacementPlan{
        .h_name = n("HD21"),
        .placement_type = .h3xr,
        .connected = .{ n("ND2 "), n(" CG "), n(" CB ") },
        .n_connected = 3,
        .bond_len = 1.02,
        .angle = 119.0,
        .dihedral = 180.0,
        .atom_type = .Hpol,
        .mover_hint = .flip_amide,
    },
    PlacementPlan{
        .h_name = n("HD22"),
        .placement_type = .h3xr,
        .connected = .{ n("ND2 "), n(" CG "), n(" CB ") },
        .n_connected = 3,
        .bond_len = 1.02,
        .angle = 119.0,
        .dihedral = 0.0,
        .atom_type = .Hpol,
        .mover_hint = .flip_amide,
    },
    backbone_h,
};

// ---------------------------------------------------------------------------
// GLN (6 H)
// ---------------------------------------------------------------------------

const gln_plans = [_]PlacementPlan{
    hxr3(" HA ", " CA ", " N  ", " C  "),
    h2xr2(" HB2", " CB ", " CA ", 120.0),
    h2xr2(" HB3", " CB ", " CA ", -120.0),
    h2xr2(" HG2", " CG ", " CB ", 120.0),
    h2xr2(" HG3", " CG ", " CB ", -120.0),
    PlacementPlan{
        .h_name = n("HE21"),
        .placement_type = .h3xr,
        .connected = .{ n("NE2 "), n(" CD "), n(" CG ") },
        .n_connected = 3,
        .bond_len = 1.02,
        .angle = 119.0,
        .dihedral = 180.0,
        .atom_type = .Hpol,
        .mover_hint = .flip_amide,
    },
    PlacementPlan{
        .h_name = n("HE22"),
        .placement_type = .h3xr,
        .connected = .{ n("NE2 "), n(" CD "), n(" CG ") },
        .n_connected = 3,
        .bond_len = 1.02,
        .angle = 119.0,
        .dihedral = 0.0,
        .atom_type = .Hpol,
        .mover_hint = .flip_amide,
    },
    backbone_h,
};

// ---------------------------------------------------------------------------
// LYS (9 H)
// ---------------------------------------------------------------------------

const lys_plans = [_]PlacementPlan{
    hxr3(" HA ", " CA ", " N  ", " C  "),
    h2xr2(" HB2", " CB ", " CA ", 120.0),
    h2xr2(" HB3", " CB ", " CA ", -120.0),
    h2xr2(" HG2", " CG ", " CB ", 120.0),
    h2xr2(" HG3", " CG ", " CB ", -120.0),
    h2xr2(" HD2", " CD ", " CG ", 120.0),
    h2xr2(" HD3", " CD ", " CG ", -120.0),
    h2xr2(" HE2", " CE ", " CD ", 120.0),
    h2xr2(" HE3", " CE ", " CD ", -120.0),
    PlacementPlan{
        .h_name = n(" HZ1"),
        .placement_type = .h3xr,
        .connected = .{ n(" NZ "), n(" CE "), n(" CD ") },
        .n_connected = 3,
        .bond_len = 1.05,
        .angle = 109.5,
        .dihedral = 180.0,
        .atom_type = .Hpol,
        .mover_hint = .rotate_nh3,
    },
    PlacementPlan{
        .h_name = n(" HZ2"),
        .placement_type = .h3xr,
        .connected = .{ n(" NZ "), n(" CE "), n(" CD ") },
        .n_connected = 3,
        .bond_len = 1.05,
        .angle = 109.5,
        .dihedral = 60.0,
        .atom_type = .Hpol,
        .mover_hint = .rotate_nh3,
    },
    PlacementPlan{
        .h_name = n(" HZ3"),
        .placement_type = .h3xr,
        .connected = .{ n(" NZ "), n(" CE "), n(" CD ") },
        .n_connected = 3,
        .bond_len = 1.05,
        .angle = 109.5,
        .dihedral = -60.0,
        .atom_type = .Hpol,
        .mover_hint = .rotate_nh3,
    },
    backbone_h,
};

// ---------------------------------------------------------------------------
// ARG (7 H)
// ---------------------------------------------------------------------------

const arg_plans = [_]PlacementPlan{
    hxr3(" HA ", " CA ", " N  ", " C  "),
    h2xr2(" HB2", " CB ", " CA ", 120.0),
    h2xr2(" HB3", " CB ", " CA ", -120.0),
    h2xr2(" HG2", " CG ", " CB ", 120.0),
    h2xr2(" HG3", " CG ", " CB ", -120.0),
    h2xr2(" HD2", " CD ", " CG ", 120.0),
    h2xr2(" HD3", " CD ", " CG ", -120.0),
    planarPol(" HE ", " CD ", " CZ ", 1.02),
    PlacementPlan{
        .h_name = n("HH11"),
        .placement_type = .h3xr,
        .connected = .{ n("NH1 "), n(" CZ "), n(" NE ") },
        .n_connected = 3,
        .bond_len = 1.02,
        .angle = 119.0,
        .dihedral = 180.0,
        .atom_type = .Hpol,
    },
    PlacementPlan{
        .h_name = n("HH12"),
        .placement_type = .h3xr,
        .connected = .{ n("NH1 "), n(" CZ "), n(" NE ") },
        .n_connected = 3,
        .bond_len = 1.02,
        .angle = 119.0,
        .dihedral = 0.0,
        .atom_type = .Hpol,
    },
    PlacementPlan{
        .h_name = n("HH21"),
        .placement_type = .h3xr,
        .connected = .{ n("NH2 "), n(" CZ "), n(" NE ") },
        .n_connected = 3,
        .bond_len = 1.02,
        .angle = 119.0,
        .dihedral = 180.0,
        .atom_type = .Hpol,
    },
    PlacementPlan{
        .h_name = n("HH22"),
        .placement_type = .h3xr,
        .connected = .{ n("NH2 "), n(" CZ "), n(" NE ") },
        .n_connected = 3,
        .bond_len = 1.02,
        .angle = 119.0,
        .dihedral = 0.0,
        .atom_type = .Hpol,
    },
    backbone_h,
};

// ---------------------------------------------------------------------------
// HIS (4 H + ring H)
// ---------------------------------------------------------------------------

const his_plans = [_]PlacementPlan{
    hxr3(" HA ", " CA ", " N  ", " C  "),
    h2xr2(" HB2", " CB ", " CA ", 120.0),
    h2xr2(" HB3", " CB ", " CA ", -120.0),
    planarH(" HD2", " CG ", "NE2 "),
    planarH(" HE1", "ND1 ", "NE2 "),
    PlacementPlan{
        .h_name = n(" HD1"),
        .placement_type = .hxr2_planar,
        .connected = .{ n("CE1 "), n(" CG "), blank },
        .n_connected = 2,
        .bond_len = 1.02,
        .atom_type = .Hpol,
        .mover_hint = .flip_his,
    },
    PlacementPlan{
        .h_name = n(" HE2"),
        .placement_type = .hxr2_planar,
        .connected = .{ n("CE1 "), n("CD2 "), blank },
        .n_connected = 2,
        .bond_len = 1.02,
        .atom_type = .Hpol,
        .mover_hint = .flip_his,
    },
    backbone_h,
};

// ---------------------------------------------------------------------------
// Lookup
// ---------------------------------------------------------------------------

pub fn getPlans(comp_id: []const u8) ?[]const PlacementPlan {
    const map = std.StaticStringMap([]const PlacementPlan).initComptime(.{
        .{ "ALA", &ala_plans },
        .{ "GLY", &gly_plans },
        .{ "VAL", &val_plans },
        .{ "LEU", &leu_plans },
        .{ "ILE", &ile_plans },
        .{ "PRO", &pro_plans },
        .{ "PHE", &phe_plans },
        .{ "TYR", &tyr_plans },
        .{ "TRP", &trp_plans },
        .{ "SER", &ser_plans },
        .{ "THR", &thr_plans },
        .{ "CYS", &cys_plans },
        .{ "MET", &met_plans },
        .{ "ASP", &asp_plans },
        .{ "GLU", &glu_plans },
        .{ "ASN", &asn_plans },
        .{ "GLN", &gln_plans },
        .{ "LYS", &lys_plans },
        .{ "ARG", &arg_plans },
        .{ "HIS", &his_plans },
    });
    return map.get(comp_id);
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

test "ALA plans" {
    const plans = getPlans("ALA");
    try std.testing.expect(plans != null);
    try std.testing.expectEqual(@as(usize, 5), plans.?.len); // HA + 3xHB + H
}

test "GLY plans" {
    const plans = getPlans("GLY");
    try std.testing.expect(plans != null);
    try std.testing.expectEqual(@as(usize, 3), plans.?.len); // HA2 + HA3 + H
}

test "PRO has no backbone H" {
    const plans = getPlans("PRO");
    try std.testing.expect(plans != null);
    // PRO should NOT have backbone H
    for (plans.?) |p| {
        const h_name = p.h_name;
        try std.testing.expect(!(h_name[0] == ' ' and h_name[1] == 'H' and h_name[2] == ' ' and h_name[3] == ' '));
    }
}

test "all 20 amino acids have plans" {
    const aas = [_][]const u8{ "ALA", "GLY", "VAL", "LEU", "ILE", "PRO", "PHE", "TYR", "TRP", "SER", "THR", "CYS", "MET", "ASP", "GLU", "ASN", "GLN", "LYS", "ARG", "HIS" };
    for (aas) |aa| {
        try std.testing.expect(getPlans(aa) != null);
    }
}

test "VAL plan count" {
    const plans = getPlans("VAL");
    try std.testing.expect(plans != null);
    // HA + HB + 3xHG1 + 3xHG2 + H = 9
    try std.testing.expectEqual(@as(usize, 9), plans.?.len);
}

test "LEU plan count" {
    const plans = getPlans("LEU");
    try std.testing.expect(plans != null);
    // HA + 2xHB + HG + 3xHD1 + 3xHD2 + H = 11
    try std.testing.expectEqual(@as(usize, 11), plans.?.len);
}

test "ILE plan count" {
    const plans = getPlans("ILE");
    try std.testing.expect(plans != null);
    // HA + HB + 2xHG1 + 3xHG2 + 3xHD1 + H = 11
    try std.testing.expectEqual(@as(usize, 11), plans.?.len);
}

test "PRO plan count" {
    const plans = getPlans("PRO");
    try std.testing.expect(plans != null);
    // HA + 2xHB + 2xHG + 2xHD = 7
    try std.testing.expectEqual(@as(usize, 7), plans.?.len);
}

test "PHE plan count" {
    const plans = getPlans("PHE");
    try std.testing.expect(plans != null);
    // HA + 2xHB + HD1 + HD2 + HE1 + HE2 + HZ + H = 9
    try std.testing.expectEqual(@as(usize, 9), plans.?.len);
}

test "TYR plan count" {
    const plans = getPlans("TYR");
    try std.testing.expect(plans != null);
    // HA + 2xHB + HD1 + HD2 + HE1 + HE2 + HH + H = 9
    try std.testing.expectEqual(@as(usize, 9), plans.?.len);
}

test "TRP plan count" {
    const plans = getPlans("TRP");
    try std.testing.expect(plans != null);
    // HA + 2xHB + HD1 + HE1 + HE3 + HZ2 + HZ3 + HH2 + H = 10
    try std.testing.expectEqual(@as(usize, 10), plans.?.len);
}

test "SER plan count" {
    const plans = getPlans("SER");
    try std.testing.expect(plans != null);
    // HA + 2xHB + HG + H = 5
    try std.testing.expectEqual(@as(usize, 5), plans.?.len);
}

test "THR plan count" {
    const plans = getPlans("THR");
    try std.testing.expect(plans != null);
    // HA + HB + HG1 + 3xHG2 + H = 7
    try std.testing.expectEqual(@as(usize, 7), plans.?.len);
}

test "CYS plan count" {
    const plans = getPlans("CYS");
    try std.testing.expect(plans != null);
    // HA + 2xHB + HG + H = 5
    try std.testing.expectEqual(@as(usize, 5), plans.?.len);
}

test "MET plan count" {
    const plans = getPlans("MET");
    try std.testing.expect(plans != null);
    // HA + 2xHB + 2xHG + 3xHE + H = 9
    try std.testing.expectEqual(@as(usize, 9), plans.?.len);
}

test "ASP plan count" {
    const plans = getPlans("ASP");
    try std.testing.expect(plans != null);
    // HA + 2xHB + H = 4
    try std.testing.expectEqual(@as(usize, 4), plans.?.len);
}

test "GLU plan count" {
    const plans = getPlans("GLU");
    try std.testing.expect(plans != null);
    // HA + 2xHB + 2xHG + H = 6
    try std.testing.expectEqual(@as(usize, 6), plans.?.len);
}

test "ASN plan count" {
    const plans = getPlans("ASN");
    try std.testing.expect(plans != null);
    // HA + 2xHB + 2xHD2 + H = 6
    try std.testing.expectEqual(@as(usize, 6), plans.?.len);
}

test "GLN plan count" {
    const plans = getPlans("GLN");
    try std.testing.expect(plans != null);
    // HA + 2xHB + 2xHG + 2xHE2 + H = 8
    try std.testing.expectEqual(@as(usize, 8), plans.?.len);
}

test "LYS plan count" {
    const plans = getPlans("LYS");
    try std.testing.expect(plans != null);
    // HA + 2xHB + 2xHG + 2xHD + 2xHE + 3xHZ + H = 13
    try std.testing.expectEqual(@as(usize, 13), plans.?.len);
}

test "ARG plan count" {
    const plans = getPlans("ARG");
    try std.testing.expect(plans != null);
    // HA + 2xHB + 2xHG + 2xHD + HE + 2xHH1 + 2xHH2 + H = 13
    try std.testing.expectEqual(@as(usize, 13), plans.?.len);
}

test "HIS plan count" {
    const plans = getPlans("HIS");
    try std.testing.expect(plans != null);
    // HA + 2xHB + HD2 + HE1 + HD1 + HE2 + H = 8
    try std.testing.expectEqual(@as(usize, 8), plans.?.len);
}

test "unknown residue returns null" {
    try std.testing.expect(getPlans("XYZ") == null);
}
