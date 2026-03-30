# Nucleotide Hydrogen Placement Templates Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add hardcoded hydrogen placement plans for 8 standard nucleotides (DA, DC, DG, DT, A, C, G, U) so they work without a CCD dictionary.

**Architecture:** New `src/place/nucleotide.zig` defines per-residue placement plans using the same `PlacementPlan` type from `standard.zig`. Placer falls back from standard → nucleotide → CCD. Sugar backbone plans are defined per-residue (not shared) because H1' glycosidic neighbor differs by base type.

**Tech Stack:** Zig 0.15, existing PlacementPlan/PlacementType/MoverHint types

---

### Task 1: Create nucleotide.zig with DA (deoxyadenosine) plans

**Files:**
- Create: `src/place/nucleotide.zig`

- [ ] **Step 1: Create nucleotide.zig with imports and helpers**

```zig
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

// -- Sugar H helpers --

/// Tetrahedral H on sp3 carbon (3 neighbors known)
fn sugar_hxr3(comptime h: []const u8, comptime center: []const u8, comptime a: []const u8, comptime b: []const u8, comptime c_: []const u8) PlacementPlan {
    return .{
        .h_name = n(h),
        .placement_type = .hxr3,
        .connected = .{ n(center), n(a), n(b) },
        .n_connected = 3,
        .bond_len = 1.09,
        .atom_type = .H,
    };
}

// Note: hxr3 in standard.zig takes (h, a, b, c) where a is the parent carbon.
// For consistency we use the same convention: connected[0] = parent, [1],[2] = other neighbors.
// But we need 4 args for center atom, so we use an explicit helper.

/// Two H atoms on sp3 carbon (2 neighbors known), specified by dihedral angle sign
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

/// OH rotator (dihedral-controlled, e.g., HO2' on RNA ribose)
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

// -- Base H helpers --

/// Aromatic C-H on planar ring (e.g., H2 on C2, H8 on C8)
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

/// N-H on planar ring (e.g., H1 on N1 of guanine, H3 on N3 of uracil)
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

/// sp2 NH2 hydrogen (e.g., H61/H62 on N6 of adenine)
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

/// Thymine methyl H (dihedral-controlled, 3-fold rotator)
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
```

- [ ] **Step 2: Add DA plans (deoxyadenosine)**

Append after the helpers:

```zig
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
```

- [ ] **Step 3: Add getPlans function and DA test**

```zig
pub fn getPlans(comp_id: []const u8) ?[]const PlacementPlan {
    const map = std.StaticStringMap([]const PlacementPlan).initComptime(.{
        .{ "DA", &da_plans },
    });
    return map.get(comp_id);
}

// ===========================================================================
// Tests
// ===========================================================================

test "DA plans" {
    const plans = getPlans("DA");
    try std.testing.expect(plans != null);
    try std.testing.expectEqual(@as(usize, 11), plans.?.len); // 7 sugar + 4 base
}

test "DA has adenine base H atoms" {
    const plans = getPlans("DA").?;
    var has_h2 = false;
    var has_h8 = false;
    var has_h61 = false;
    var has_h62 = false;
    for (plans) |p| {
        const h = std.mem.trimRight(u8, &p.h_name, " ");
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

test "DA has no polymer linking H (HOP2, HOP3, HO3')" {
    const plans = getPlans("DA").?;
    for (plans) |p| {
        const h = std.mem.trimRight(u8, &p.h_name, " ");
        try std.testing.expect(!std.mem.eql(u8, h, "HOP2"));
        try std.testing.expect(!std.mem.eql(u8, h, "HOP3"));
        try std.testing.expect(!std.mem.eql(u8, h, "HO3'"));
    }
}
```

- [ ] **Step 4: Run tests**

Run: `zig build test --summary all`
Expected: All existing tests pass + 3 new DA tests pass.

Note: the file won't be discovered by test runner yet — it needs to be imported. We'll add it to place.zig in Task 3. For now, verify the file compiles by adding a temporary import in the test command or waiting until Task 3.

Actually, since Zig discovers tests through imports and `nucleotide.zig` isn't imported yet, the tests won't run. Skip this step and verify in Task 3.

- [ ] **Step 5: Commit**

```bash
git add src/place/nucleotide.zig
git commit -m "feat: add nucleotide.zig with DA (deoxyadenosine) placement plans"
```

---

### Task 2: Add remaining 7 nucleotide plans (DC, DG, DT, A, C, G, U)

**Files:**
- Modify: `src/place/nucleotide.zig`

- [ ] **Step 1: Add DC plans (deoxycytidine)**

```zig
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
```

- [ ] **Step 2: Add DG plans (deoxyguanosine)**

```zig
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
```

- [ ] **Step 3: Add DT plans (deoxythymidine)**

```zig
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
```

- [ ] **Step 4: Add RNA nucleotide plans (A, C, G, U)**

```zig
// ===========================================================================
// RNA nucleotides
// ===========================================================================

// A: adenosine (purine, glycosidic bond at N9, ribose with 2'-OH)
const a_plans = [_]PlacementPlan{
    // Sugar (ribose — H2' is single, HO2' is OH rotator)
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
```

- [ ] **Step 5: Update getPlans to include all 8 residues**

```zig
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
```

- [ ] **Step 6: Add tests for all residues**

```zig
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
    try std.testing.expectEqual(@as(usize, 11), getPlans("A").?.len);  // 8 sugar + 4 base - 1 (HO2' counted in 8)
    // Actually: 7 sugar-C-H + HO2' = 8 sugar, + 4 base = 12? Let me recount...
    // H1', H2'(1), HO2', H3', H4', H5', H5'' = 7 atoms in sugar
    // + 4 base = 11
    try std.testing.expectEqual(@as(usize, 11), getPlans("C").?.len);
    try std.testing.expectEqual(@as(usize, 11), getPlans("G").?.len);
    try std.testing.expectEqual(@as(usize, 10), getPlans("U").?.len);  // 7 sugar + 3 base
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
```

- [ ] **Step 7: Commit**

```bash
git add src/place/nucleotide.zig
git commit -m "feat: add placement plans for all 8 standard nucleotides"
```

---

### Task 3: Integrate nucleotide plans into placer.zig

**Files:**
- Modify: `src/place/placer.zig:13,70`
- Modify: `src/place.zig` (add nucleotide export)

- [ ] **Step 1: Add nucleotide import to placer.zig**

At line 13 of `src/place/placer.zig`, after `const standard = @import("standard.zig");`, add:

```zig
const nucleotide = @import("nucleotide.zig");
```

- [ ] **Step 2: Add nucleotide fallback in placement logic**

At line 70, change:

```zig
if (standard.getPlans(comp_id)) |plans| {
```

The block after this handles standard residues. After its closing brace (the `else` that tries CCD), insert the nucleotide fallback. The structure becomes:

```zig
if (standard.getPlans(comp_id)) |plans| {
    // ... existing standard residue handling ...
} else if (nucleotide.getPlans(comp_id)) |plans| {
    // Nucleotide with hardcoded plans (same execution path as standard)
    const altlocs = collectAltlocs(mdl, res);
    for (plans) |*plan| {
        for (altlocs) |alt| {
            if (try executePlan(mdl, res, @intCast(res_idx), plan, null, alt)) {
                n_placed += 1;
            }
        }
    }
} else if (dict) |d| {
    // CCD fallback...
```

Note: nucleotide plans pass `null` for bonds (no topology.zig bond table for nucleotides — the placement plans encode all needed geometry directly). The `executePlan` function already handles `bonds == null` by falling back to distance-based neighbor finding.

- [ ] **Step 3: Add nucleotide export to place.zig**

In `src/place.zig`, add after the `standard` line:

```zig
pub const nucleotide = @import("place/nucleotide.zig");
```

- [ ] **Step 4: Run tests**

Run: `zig build test --summary all`
Expected: All tests pass including the new nucleotide tests (they are now discovered through the place.zig → nucleotide.zig import chain).

- [ ] **Step 5: Commit**

```bash
git add src/place/placer.zig src/place.zig
git commit -m "feat: integrate nucleotide templates into placer fallback chain"
```

---

### Task 4: End-to-end verification with real structures

**Files:**
- No code changes — verification only

- [ ] **Step 1: Build release and test with a DNA-containing structure**

```bash
zig build -Doptimize=ReleaseFast

# Test with fold_test2 if it has DNA/RNA, or use E. coli structures
# First, find a structure with nucleotides
ls /Users/nagaet/pdb/afdb/UP000000625_83333_ECOLI_v6/cif/ | head -5
```

Note: AlphaFold models are protein-only. We need a PDB structure with DNA/RNA. Check examples/data or test_data for any.

If no suitable test file exists, create a minimal one or use one from the PDB. The batch test with E. coli structures already passed (protein-only), so this test focuses on verifying nucleotide placement without CCD.

- [ ] **Step 2: Verify nucleotide H counts match between hardcoded and CCD**

If a structure with nucleotides is available:
```bash
# Without CCD (uses new hardcoded templates)
./zig-out/bin/zreduce run structure_with_dna.cif -o /tmp/no_ccd.cif

# With CCD (uses CCD-derived plans)
./zig-out/bin/zreduce run structure_with_dna.cif -d /Users/nagaet/pdb/data/monomers/components.cif.gz -o /tmp/with_ccd.cif

# Compare hydrogen counts
grep "^ATOM.*H " /tmp/no_ccd.cif | wc -l
grep "^ATOM.*H " /tmp/with_ccd.cif | wc -l
```

- [ ] **Step 3: Run full test suite**

Run: `zig build test --summary all`
Expected: All tests pass.

- [ ] **Step 4: Commit design doc and plan**

```bash
git add docs/specs/2026-03-30-nucleotide-templates-design.md docs/specs/2026-03-30-nucleotide-templates-plan.md
git commit -m "docs: add nucleotide templates design and implementation plan"
```
