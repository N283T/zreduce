//! Hardcoded hydrogen placement plans for common modified amino acids.
//! Derived from parent residue plans in standard.zig with specific atom changes.
//! Covers: MSE, SEP, TPO, CSO, PCA, PTR.

const std = @import("std");
const standard = @import("standard.zig");
const element = @import("../element.zig");
const PlacementPlan = standard.PlacementPlan;
const PlacementType = standard.PlacementType;
const MoverHint = standard.MoverHint;

// Atom name helper (same as standard.zig)
fn n(comptime s: []const u8) [4]u8 {
    var buf: [4]u8 = .{ ' ', ' ', ' ', ' ' };
    for (s, 0..) |c, i| {
        if (i >= 4) break;
        buf[i] = c;
    }
    return buf;
}

const blank = n("    ");

// Bond length constants from standard.zig (canonical CCD-derived values)
const c_h_sp3 = standard.c_h_sp3;
const c_h_arom = standard.c_h_arom;
const n_h_backbone = standard.n_h_backbone;

// Reuse standard.zig helper patterns inline

fn hxr3(comptime h: []const u8, comptime a: []const u8, comptime b: []const u8, comptime c: []const u8) PlacementPlan {
    return .{
        .h_name = n(h),
        .placement_type = .hxr3,
        .connected = .{ n(a), n(b), n(c) },
        .n_connected = 3,
        .bond_len = c_h_sp3,
        .atom_type = .H,
    };
}

fn h2xr2(comptime h: []const u8, comptime a: []const u8, comptime b: []const u8, dihedral: f32) PlacementPlan {
    return .{
        .h_name = n(h),
        .placement_type = .h2xr2,
        .connected = .{ n(a), n(b), blank },
        .n_connected = 2,
        .bond_len = c_h_sp3,
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
        .bond_len = c_h_sp3,
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
        .bond_len = c_h_arom,
        .atom_type = .Har,
    };
}

const backbone_h = PlacementPlan{
    .h_name = n(" H  "),
    .placement_type = .h3xr,
    .connected = .{ n(" N  "), n(" CA "), n(" C  ") },
    .n_connected = 3,
    .bond_len = n_h_backbone,
    .angle = 119.0,
    .dihedral = 180.0,
    .atom_type = .Hpol,
};

// ===========================================================================
// MSE: Selenomethionine (parent: MET, SD → SE)
// ===========================================================================

const mse_plans = [_]PlacementPlan{
    hxr3(" HA ", " CA ", " N  ", " C  "),
    h2xr2(" HB2", " CB ", " CA ", 120.0),
    h2xr2(" HB3", " CB ", " CA ", -120.0),
    h2xr2(" HG2", " CG ", " CB ", 120.0),
    h2xr2(" HG3", " CG ", " CB ", -120.0),
    methyl(" HE1", " CE ", " SE ", 180.0),
    methyl(" HE2", " CE ", " SE ", 60.0),
    methyl(" HE3", " CE ", " SE ", -60.0),
    backbone_h,
};

// ===========================================================================
// SEP: Phosphoserine (parent: SER, OG-H removed — phosphorylated)
// ===========================================================================

const sep_plans = [_]PlacementPlan{
    hxr3(" HA ", " CA ", " N  ", " C  "),
    h2xr2(" HB2", " CB ", " CA ", 120.0),
    h2xr2(" HB3", " CB ", " CA ", -120.0),
    backbone_h,
};

// ===========================================================================
// TPO: Phosphothreonine (parent: THR, OG1-H removed — phosphorylated)
// ===========================================================================

const tpo_plans = [_]PlacementPlan{
    hxr3(" HA ", " CA ", " N  ", " C  "),
    hxr3(" HB ", " CB ", " CA ", "OG1 "),
    methyl("HG21", "CG2 ", " CB ", 180.0),
    methyl("HG22", "CG2 ", " CB ", 60.0),
    methyl("HG23", "CG2 ", " CB ", -60.0),
    backbone_h,
};

// ===========================================================================
// CSO: S-hydroxycysteine (parent: CYS, SG-H removed, OD-H added)
// ===========================================================================

const cso_plans = [_]PlacementPlan{
    hxr3(" HA ", " CA ", " N  ", " C  "),
    h2xr2(" HB2", " CB ", " CA ", 120.0),
    h2xr2(" HB3", " CB ", " CA ", -120.0),
    PlacementPlan{
        .h_name = n(" HD "),
        .placement_type = .h3xr,
        .connected = .{ n(" OD "), n(" SG "), blank },
        .n_connected = 2,
        .bond_len = 0.97,
        .angle = 109.5,
        .dihedral = 180.0,
        .atom_type = .Hpol,
        .mover_hint = .rotate,
    },
    backbone_h,
};

// ===========================================================================
// PCA: Pyroglutamic acid (parent: GLU, cyclized N-CD, no backbone H)
// ===========================================================================

const pca_plans = [_]PlacementPlan{
    hxr3(" HA ", " CA ", " N  ", " C  "),
    h2xr2(" HB2", " CB ", " CA ", 120.0),
    h2xr2(" HB3", " CB ", " CA ", -120.0),
    h2xr2(" HG2", " CG ", " CB ", 120.0),
    h2xr2(" HG3", " CG ", " CB ", -120.0),
    // No backbone H — N is in the pyrrolidone ring (bonded to CD)
};

// ===========================================================================
// PTR: O-phosphotyrosine (parent: TYR, OH-H removed — phosphorylated)
// ===========================================================================

const ptr_plans = [_]PlacementPlan{
    hxr3(" HA ", " CA ", " N  ", " C  "),
    h2xr2(" HB2", " CB ", " CA ", 120.0),
    h2xr2(" HB3", " CB ", " CA ", -120.0),
    planarH(" HD1", " CG ", "CE1 "),
    planarH(" HD2", " CG ", "CE2 "),
    planarH(" HE1", "CD1 ", " CZ "),
    planarH(" HE2", "CD2 ", " CZ "),
    backbone_h,
};

// ===========================================================================
// Lookup
// ===========================================================================

pub fn getPlans(comp_id: []const u8) ?[]const PlacementPlan {
    const map = std.StaticStringMap([]const PlacementPlan).initComptime(.{
        .{ "MSE", &mse_plans },
        .{ "SEP", &sep_plans },
        .{ "TPO", &tpo_plans },
        .{ "CSO", &cso_plans },
        .{ "PCA", &pca_plans },
        .{ "PTR", &ptr_plans },
    });
    return map.get(comp_id);
}

// ===========================================================================
// Tests
// ===========================================================================

test "all 6 modified residues have plans" {
    const mods = [_][]const u8{ "MSE", "SEP", "TPO", "CSO", "PCA", "PTR" };
    for (mods) |m| {
        try std.testing.expect(getPlans(m) != null);
    }
}

test "plan counts" {
    try std.testing.expectEqual(@as(usize, 9), getPlans("MSE").?.len);
    try std.testing.expectEqual(@as(usize, 4), getPlans("SEP").?.len);
    try std.testing.expectEqual(@as(usize, 6), getPlans("TPO").?.len);
    try std.testing.expectEqual(@as(usize, 5), getPlans("CSO").?.len);
    try std.testing.expectEqual(@as(usize, 5), getPlans("PCA").?.len);
    try std.testing.expectEqual(@as(usize, 8), getPlans("PTR").?.len);
}

test "MSE has methyl rotator on CE" {
    const plans = getPlans("MSE").?;
    var methyl_count: usize = 0;
    for (plans) |p| {
        if (p.mover_hint == .rotate_methyl) methyl_count += 1;
    }
    try std.testing.expectEqual(@as(usize, 3), methyl_count);
    // Verify SE is referenced (not SD)
    for (plans) |p| {
        const h = std.mem.trim(u8, &p.h_name, " ");
        if (std.mem.eql(u8, h, "HE1")) {
            // connected[1] should be SE, not SD
            try std.testing.expectEqualStrings("SE", std.mem.trim(u8, &p.connected[1], " "));
        }
    }
}

test "CSO has OH rotator on OD" {
    const plans = getPlans("CSO").?;
    var found_hd = false;
    for (plans) |p| {
        const h = std.mem.trim(u8, &p.h_name, " ");
        if (std.mem.eql(u8, h, "HD")) {
            found_hd = true;
            try std.testing.expectEqual(MoverHint.rotate, p.mover_hint);
            try std.testing.expectEqual(element.AtomType.Hpol, p.atom_type);
        }
    }
    try std.testing.expect(found_hd);
}

test "PCA has no backbone H" {
    const plans = getPlans("PCA").?;
    for (plans) |p| {
        const h = std.mem.trim(u8, &p.h_name, " ");
        try std.testing.expect(!std.mem.eql(u8, h, "H"));
    }
}

test "SEP, TPO, PTR have no phosphate or hydroxyl OH H" {
    const phospho = [_][]const u8{ "SEP", "TPO", "PTR" };
    for (phospho) |comp| {
        for (getPlans(comp).?) |p| {
            const h = std.mem.trim(u8, &p.h_name, " ");
            // No phosphate H
            try std.testing.expect(!std.mem.eql(u8, h, "HOP2"));
            try std.testing.expect(!std.mem.eql(u8, h, "HOP3"));
            try std.testing.expect(!std.mem.eql(u8, h, "HO2P"));
            try std.testing.expect(!std.mem.eql(u8, h, "HO3P"));
            // No parent OH H (HG for SER, HG1 for THR, HH for TYR)
            try std.testing.expect(!std.mem.eql(u8, h, "HG"));
            try std.testing.expect(!std.mem.eql(u8, h, "HG1"));
            try std.testing.expect(!std.mem.eql(u8, h, "HH"));
        }
    }
}

test "TPO has methyl rotator on CG2" {
    const plans = getPlans("TPO").?;
    var methyl_count: usize = 0;
    for (plans) |p| {
        if (p.mover_hint == .rotate_methyl) methyl_count += 1;
    }
    try std.testing.expectEqual(@as(usize, 3), methyl_count);
}

test "PTR has aromatic ring H" {
    const plans = getPlans("PTR").?;
    var ar_count: usize = 0;
    for (plans) |p| {
        if (p.atom_type == .Har) ar_count += 1;
    }
    try std.testing.expectEqual(@as(usize, 4), ar_count); // HD1, HD2, HE1, HE2
}

test "unknown modified residue returns null" {
    try std.testing.expect(getPlans("XYZ") == null);
    try std.testing.expect(getPlans("ALA") == null);
    try std.testing.expect(getPlans("DA") == null);
}
