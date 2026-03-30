//! Hardcoded hydrogen placement plans for standard DNA/RNA nucleotides.
//! Covers: DA, DC, DG, DT (DNA) and A, C, G, U (RNA).

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

// ---------------------------------------------------------------------------
// Sugar H helpers
// ---------------------------------------------------------------------------

/// Tetrahedral H on sp3 carbon (3 neighbors known).
/// connected[0] = parent carbon, connected[1],[2] = two known neighbors.
/// The 3rd neighbor is found at runtime by the placer.
fn sugar_hxr3(comptime h: []const u8, comptime center: []const u8, comptime a: []const u8, comptime b: []const u8, comptime c_: []const u8) PlacementPlan {
    _ = c_; // 3rd neighbor found at runtime
    return .{
        .h_name = n(h),
        .placement_type = .hxr3,
        .connected = .{ n(center), n(a), n(b) },
        .n_connected = 3,
        .bond_len = 1.09,
        .atom_type = .H,
    };
}

/// Two H atoms on sp3 carbon (2 neighbors known), specified by dihedral angle sign.
fn sugar_h2xr2(comptime h: []const u8, comptime center: []const u8, comptime a: []const u8, comptime b: []const u8, dihedral: f32) PlacementPlan {
    return .{
        .h_name = n(h),
        .placement_type = .h2xr2,
        .connected = .{ n(center), n(a), n(b) },
        .n_connected = 2,
        .bond_len = 1.09,
        .angle = 109.5,
        .dihedral = dihedral,
        .atom_type = .H,
    };
}

/// OH rotator (dihedral-controlled, e.g., HO2' on RNA ribose).
fn oh_rotator(comptime h: []const u8, comptime o: []const u8, comptime c: []const u8, comptime c2: []const u8) PlacementPlan {
    return .{
        .h_name = n(h),
        .placement_type = .h3xr,
        .connected = .{ n(o), n(c), n(c2) },
        .n_connected = 2,
        .bond_len = 0.98,
        .angle = 109.5,
        .dihedral = 180.0,
        .atom_type = .Hpol,
        .mover_hint = .rotate,
    };
}

// ---------------------------------------------------------------------------
// Base H helpers
// ---------------------------------------------------------------------------

/// Aromatic C-H on planar ring (e.g., H2 on C2, H8 on C8).
fn aromatic_ch(comptime h: []const u8, comptime c: []const u8, comptime nb1: []const u8, comptime nb2: []const u8) PlacementPlan {
    return .{
        .h_name = n(h),
        .placement_type = .hxr2_planar,
        .connected = .{ n(c), n(nb1), n(nb2) },
        .n_connected = 2,
        .bond_len = 1.08,
        .atom_type = .Har,
    };
}

/// N-H on planar ring (e.g., H1 on N1 of guanine, H3 on N3 of uracil).
fn ring_nh(comptime h: []const u8, comptime nn: []const u8, comptime nb1: []const u8, comptime nb2: []const u8) PlacementPlan {
    return .{
        .h_name = n(h),
        .placement_type = .hxr2_planar,
        .connected = .{ n(nn), n(nb1), n(nb2) },
        .n_connected = 2,
        .bond_len = 1.02,
        .atom_type = .Hpol,
    };
}

/// sp2 NH2 hydrogen (e.g., H61/H62 on N6 of adenine).
fn amino_h(comptime h: []const u8, comptime nn: []const u8, comptime c: []const u8, comptime c2: []const u8, angle: f32) PlacementPlan {
    return .{
        .h_name = n(h),
        .placement_type = .h3xr,
        .connected = .{ n(nn), n(c), n(c2) },
        .n_connected = 2,
        .bond_len = 1.02,
        .angle = angle,
        .dihedral = 180.0,
        .atom_type = .Hpol,
    };
}

/// Thymine methyl H (dihedral-controlled, 3-fold rotator).
fn thymine_methyl(comptime h: []const u8, dihedral: f32) PlacementPlan {
    return .{
        .h_name = n(h),
        .placement_type = .h3xr,
        .connected = .{ n(" C7 "), n(" C5 "), n(" C4 ") },
        .n_connected = 2,
        .bond_len = 1.09,
        .angle = 109.5,
        .dihedral = dihedral,
        .atom_type = .H,
        .mover_hint = .rotate_methyl,
    };
}

// ===========================================================================
// DNA nucleotides
// ===========================================================================

// DA: 2'-deoxyadenosine (purine, glycosidic bond at N9)
const da_plans = [_]PlacementPlan{
    // Sugar (deoxyribose)
    sugar_hxr3(" H1'", " C1'", " O4'", " C2'", " N9 "),
    sugar_h2xr2(" H2'", " C2'", " C3'", " C1'", -126.5),
    sugar_h2xr2("H2''", " C2'", " C3'", " C1'", 126.5),
    sugar_hxr3(" H3'", " C3'", " C4'", " C2'", " O3'"),
    sugar_hxr3(" H4'", " C4'", " C5'", " C3'", " O4'"),
    sugar_h2xr2(" H5'", " C5'", " C4'", " O5'", -126.5),
    sugar_h2xr2("H5''", " C5'", " C4'", " O5'", 126.5),
    // Base (adenine)
    aromatic_ch(" H2 ", " C2 ", " N1 ", " N3 "),
    aromatic_ch(" H8 ", " C8 ", " N7 ", " N9 "),
    amino_h("H61 ", " N6 ", " C6 ", " C5 ", 120.0),
    amino_h("H62 ", " N6 ", " C6 ", " C5 ", -120.0),
};

// DC: 2'-deoxycytidine (pyrimidine, glycosidic bond at N1)
const dc_plans = [_]PlacementPlan{
    // Sugar (deoxyribose)
    sugar_hxr3(" H1'", " C1'", " O4'", " C2'", " N1 "),
    sugar_h2xr2(" H2'", " C2'", " C3'", " C1'", -126.5),
    sugar_h2xr2("H2''", " C2'", " C3'", " C1'", 126.5),
    sugar_hxr3(" H3'", " C3'", " C4'", " C2'", " O3'"),
    sugar_hxr3(" H4'", " C4'", " C5'", " C3'", " O4'"),
    sugar_h2xr2(" H5'", " C5'", " C4'", " O5'", -126.5),
    sugar_h2xr2("H5''", " C5'", " C4'", " O5'", 126.5),
    // Base (cytosine)
    aromatic_ch(" H5 ", " C5 ", " C4 ", " C6 "),
    aromatic_ch(" H6 ", " C6 ", " C5 ", " N1 "),
    amino_h("H41 ", " N4 ", " C4 ", " N3 ", 120.0),
    amino_h("H42 ", " N4 ", " C4 ", " N3 ", -120.0),
};

// DG: 2'-deoxyguanosine (purine, glycosidic bond at N9)
const dg_plans = [_]PlacementPlan{
    // Sugar (deoxyribose)
    sugar_hxr3(" H1'", " C1'", " O4'", " C2'", " N9 "),
    sugar_h2xr2(" H2'", " C2'", " C3'", " C1'", -126.5),
    sugar_h2xr2("H2''", " C2'", " C3'", " C1'", 126.5),
    sugar_hxr3(" H3'", " C3'", " C4'", " C2'", " O3'"),
    sugar_hxr3(" H4'", " C4'", " C5'", " C3'", " O4'"),
    sugar_h2xr2(" H5'", " C5'", " C4'", " O5'", -126.5),
    sugar_h2xr2("H5''", " C5'", " C4'", " O5'", 126.5),
    // Base (guanine)
    aromatic_ch(" H8 ", " C8 ", " N7 ", " N9 "),
    ring_nh(" H1 ", " N1 ", " C6 ", " C2 "),
    amino_h(" H21", " N2 ", " C2 ", " N1 ", 120.0),
    amino_h(" H22", " N2 ", " C2 ", " N1 ", -120.0),
};

// DT: 2'-deoxythymidine (pyrimidine, glycosidic bond at N1)
const dt_plans = [_]PlacementPlan{
    // Sugar (deoxyribose)
    sugar_hxr3(" H1'", " C1'", " O4'", " C2'", " N1 "),
    sugar_h2xr2(" H2'", " C2'", " C3'", " C1'", -126.5),
    sugar_h2xr2("H2''", " C2'", " C3'", " C1'", 126.5),
    sugar_hxr3(" H3'", " C3'", " C4'", " C2'", " O3'"),
    sugar_hxr3(" H4'", " C4'", " C5'", " C3'", " O4'"),
    sugar_h2xr2(" H5'", " C5'", " C4'", " O5'", -126.5),
    sugar_h2xr2("H5''", " C5'", " C4'", " O5'", 126.5),
    // Base (thymine)
    aromatic_ch(" H6 ", " C6 ", " C5 ", " N1 "),
    ring_nh(" H3 ", " N3 ", " C4 ", " C2 "),
    thymine_methyl("H71 ", 180.0),
    thymine_methyl("H72 ", -60.0),
    thymine_methyl("H73 ", 60.0),
};

// ===========================================================================
// RNA nucleotides
// ===========================================================================

// A: adenosine (purine, glycosidic bond at N9, ribose with 2'-OH)
const a_plans = [_]PlacementPlan{
    // Sugar (ribose -- H2' is single, HO2' is OH rotator)
    sugar_hxr3(" H1'", " C1'", " O4'", " C2'", " N9 "),
    sugar_hxr3(" H2'", " C2'", " C3'", " C1'", " O2'"),
    oh_rotator("HO2'", " O2'", " C2'", " C3'"),
    sugar_hxr3(" H3'", " C3'", " C4'", " C2'", " O3'"),
    sugar_hxr3(" H4'", " C4'", " C5'", " C3'", " O4'"),
    sugar_h2xr2(" H5'", " C5'", " C4'", " O5'", -126.5),
    sugar_h2xr2("H5''", " C5'", " C4'", " O5'", 126.5),
    // Base (adenine)
    aromatic_ch(" H2 ", " C2 ", " N1 ", " N3 "),
    aromatic_ch(" H8 ", " C8 ", " N7 ", " N9 "),
    amino_h("H61 ", " N6 ", " C6 ", " C5 ", 120.0),
    amino_h("H62 ", " N6 ", " C6 ", " C5 ", -120.0),
};

// C: cytidine (pyrimidine, glycosidic bond at N1, ribose with 2'-OH)
const c_plans = [_]PlacementPlan{
    sugar_hxr3(" H1'", " C1'", " O4'", " C2'", " N1 "),
    sugar_hxr3(" H2'", " C2'", " C3'", " C1'", " O2'"),
    oh_rotator("HO2'", " O2'", " C2'", " C3'"),
    sugar_hxr3(" H3'", " C3'", " C4'", " C2'", " O3'"),
    sugar_hxr3(" H4'", " C4'", " C5'", " C3'", " O4'"),
    sugar_h2xr2(" H5'", " C5'", " C4'", " O5'", -126.5),
    sugar_h2xr2("H5''", " C5'", " C4'", " O5'", 126.5),
    // Base (cytosine)
    aromatic_ch(" H5 ", " C5 ", " C4 ", " C6 "),
    aromatic_ch(" H6 ", " C6 ", " C5 ", " N1 "),
    amino_h("H41 ", " N4 ", " C4 ", " N3 ", 120.0),
    amino_h("H42 ", " N4 ", " C4 ", " N3 ", -120.0),
};

// G: guanosine (purine, glycosidic bond at N9, ribose with 2'-OH)
const g_plans = [_]PlacementPlan{
    sugar_hxr3(" H1'", " C1'", " O4'", " C2'", " N9 "),
    sugar_hxr3(" H2'", " C2'", " C3'", " C1'", " O2'"),
    oh_rotator("HO2'", " O2'", " C2'", " C3'"),
    sugar_hxr3(" H3'", " C3'", " C4'", " C2'", " O3'"),
    sugar_hxr3(" H4'", " C4'", " C5'", " C3'", " O4'"),
    sugar_h2xr2(" H5'", " C5'", " C4'", " O5'", -126.5),
    sugar_h2xr2("H5''", " C5'", " C4'", " O5'", 126.5),
    // Base (guanine)
    aromatic_ch(" H8 ", " C8 ", " N7 ", " N9 "),
    ring_nh(" H1 ", " N1 ", " C6 ", " C2 "),
    amino_h(" H21", " N2 ", " C2 ", " N1 ", 120.0),
    amino_h(" H22", " N2 ", " C2 ", " N1 ", -120.0),
};

// U: uridine (pyrimidine, glycosidic bond at N1, ribose with 2'-OH)
const u_plans = [_]PlacementPlan{
    sugar_hxr3(" H1'", " C1'", " O4'", " C2'", " N1 "),
    sugar_hxr3(" H2'", " C2'", " C3'", " C1'", " O2'"),
    oh_rotator("HO2'", " O2'", " C2'", " C3'"),
    sugar_hxr3(" H3'", " C3'", " C4'", " C2'", " O3'"),
    sugar_hxr3(" H4'", " C4'", " C5'", " C3'", " O4'"),
    sugar_h2xr2(" H5'", " C5'", " C4'", " O5'", -126.5),
    sugar_h2xr2("H5''", " C5'", " C4'", " O5'", 126.5),
    // Base (uracil)
    aromatic_ch(" H5 ", " C5 ", " C4 ", " C6 "),
    aromatic_ch(" H6 ", " C6 ", " C5 ", " N1 "),
    ring_nh(" H3 ", " N3 ", " C4 ", " C2 "),
};

// ===========================================================================
// Lookup
// ===========================================================================

pub fn getPlans(comp_id: []const u8) ?[]const PlacementPlan {
    const map = std.StaticStringMap([]const PlacementPlan).initComptime(.{
        .{ "DA", &da_plans },
        .{ "DC", &dc_plans },
        .{ "DG", &dg_plans },
        .{ "DT", &dt_plans },
        .{ "A", &a_plans },
        .{ "C", &c_plans },
        .{ "G", &g_plans },
        .{ "U", &u_plans },
    });
    return map.get(comp_id);
}

// ===========================================================================
// Tests
// ===========================================================================

test "all 8 nucleotides have plans" {
    const nucs = [_][]const u8{ "DA", "DC", "DG", "DT", "A", "C", "G", "U" };
    for (nucs) |nuc| {
        try std.testing.expect(getPlans(nuc) != null);
    }
}

test "DNA plan counts" {
    try std.testing.expectEqual(@as(usize, 11), getPlans("DA").?.len); // 7 sugar + 4 base
    try std.testing.expectEqual(@as(usize, 11), getPlans("DC").?.len); // 7 sugar + 4 base
    try std.testing.expectEqual(@as(usize, 11), getPlans("DG").?.len); // 7 sugar + 4 base
    try std.testing.expectEqual(@as(usize, 12), getPlans("DT").?.len); // 7 sugar + 5 base (methyl)
}

test "RNA plan counts" {
    try std.testing.expectEqual(@as(usize, 11), getPlans("A").?.len); // 7 sugar + 4 base
    try std.testing.expectEqual(@as(usize, 11), getPlans("C").?.len); // 7 sugar + 4 base
    try std.testing.expectEqual(@as(usize, 11), getPlans("G").?.len); // 7 sugar + 4 base
    try std.testing.expectEqual(@as(usize, 10), getPlans("U").?.len); // 7 sugar + 3 base
}

test "DA has adenine base H atoms" {
    const plans = getPlans("DA").?;
    var has_h2 = false;
    var has_h8 = false;
    var has_h61 = false;
    var has_h62 = false;
    for (plans) |p| {
        const h = std.mem.trim(u8, &p.h_name, " ");
        if (std.mem.eql(u8, h, "H2")) has_h2 = true;
        if (std.mem.eql(u8, h, "H8")) has_h8 = true;
        if (std.mem.eql(u8, h, "H61")) has_h61 = true;
        if (std.mem.eql(u8, h, "H62")) has_h62 = true;
    }
    try std.testing.expect(has_h2);
    try std.testing.expect(has_h8);
    try std.testing.expect(has_h61);
    try std.testing.expect(has_h62);
}

test "DT has methyl rotator" {
    const plans = getPlans("DT").?;
    var methyl_count: usize = 0;
    for (plans) |p| {
        if (p.mover_hint == .rotate_methyl) methyl_count += 1;
    }
    try std.testing.expectEqual(@as(usize, 3), methyl_count); // H71, H72, H73
}

test "RNA nucleotides have HO2' rotator" {
    const rna = [_][]const u8{ "A", "C", "G", "U" };
    for (rna) |nuc| {
        const plans = getPlans(nuc).?;
        var has_ho2 = false;
        for (plans) |p| {
            const h = std.mem.trimRight(u8, &p.h_name, " ");
            if (std.mem.eql(u8, h, "HO2'") and p.mover_hint == .rotate) has_ho2 = true;
        }
        try std.testing.expect(has_ho2);
    }
}

test "DNA nucleotides have H2'' (deoxyribose)" {
    const dna = [_][]const u8{ "DA", "DC", "DG", "DT" };
    for (dna) |nuc| {
        const plans = getPlans(nuc).?;
        var has_h2pp = false;
        for (plans) |p| {
            const h = std.mem.trimRight(u8, &p.h_name, " ");
            if (std.mem.eql(u8, h, "H2''")) has_h2pp = true;
        }
        try std.testing.expect(has_h2pp);
    }
}

test "no polymer linking H in any nucleotide" {
    const all = [_][]const u8{ "DA", "DC", "DG", "DT", "A", "C", "G", "U" };
    for (all) |nuc| {
        for (getPlans(nuc).?) |p| {
            const h = std.mem.trimRight(u8, &p.h_name, " ");
            try std.testing.expect(!std.mem.eql(u8, h, "HOP2"));
            try std.testing.expect(!std.mem.eql(u8, h, "HOP3"));
            try std.testing.expect(!std.mem.eql(u8, h, "HO3'"));
            try std.testing.expect(!std.mem.eql(u8, h, "HO5'"));
        }
    }
}

test "unknown nucleotide returns null" {
    try std.testing.expect(getPlans("XYZ") == null);
    try std.testing.expect(getPlans("ALA") == null);
}
