# Nucleotide Hydrogen Placement Templates Design

## Context

Currently only 20 standard amino acids have hardcoded placement plans in `standard.zig`. DNA/RNA nucleotides (DA, DC, DG, DT, A, C, G, U) require the CCD dictionary for hydrogen placement. These 8 residues appear in 9,000-12,000+ PDB entries each.

Hardcoding removes the CCD dependency for the most common non-protein residues.

## Reference Sources

- Original reduce `StdResH.cpp`: has templates for all 8 nucleotides + shared "ribose phosphate backbone"
- CCD `chem_comp_atom` / `chem_comp_bond`: authoritative atom names and bond topology
- CCD individual files at `~/pdb/refdata/chem_comp/`

## Design

### New file: `src/place/nucleotide.zig`

Separate from `standard.zig` (which is already 777 lines). Same `PlacementPlan` type, same `getPlans()` pattern.

### Architecture

```
standard.zig:getPlans(comp_id)    -> amino acid plans (20 residues)
nucleotide.zig:getPlans(comp_id)  -> nucleotide plans (8 residues)
placer.zig                        -> tries standard, then nucleotide, then CCD fallback
```

### Placement Plan Groups

#### Shared Sugar Backbone (DNA)

All DNA residues (DA, DC, DG, DT) share:

| H | Parent | Type | Bond | Neighbors | Notes |
|---|--------|------|------|-----------|-------|
| H1' | C1' | hxr3 (tetrahedral) | 1.09 | O4', C2', N9/N1 | glycosidic carbon |
| H2' | C2' | h2xr2 (2H on sp3) | 1.09 | C3', C1' | deoxyribose, angle -126.5 |
| H2'' | C2' | h2xr2 | 1.09 | C3', C1' | deoxyribose, angle +126.5 |
| H3' | C3' | hxr3 | 1.09 | C4', C2', O3' | |
| H4' | C4' | hxr3 | 1.09 | C5', C3', O4' | |
| H5' | C5' | h2xr2 | 1.09 | C4', O5' | angle -126.5 |
| H5'' | C5' | h2xr2 | 1.09 | C4', O5' | angle +126.5 |

#### Shared Sugar Backbone (RNA)

RNA residues (A, C, G, U) share the same as DNA except at C2':

| H | Parent | Type | Bond | Neighbors | Notes |
|---|--------|------|------|-----------|-------|
| H1' | C1' | hxr3 | 1.09 | O4', C2', N9/N1 | |
| H2' | C2' | hxr3 | 1.09 | C3', C1', O2' | ribose, single H (O2' present) |
| HO2' | O2' | h3xr (dihedral) | 0.98 | C2', C3' | OH rotator, 109.5°, dihedral 180° |
| H3' | C3' | hxr3 | 1.09 | C4', C2', O3' | |
| H4' | C4' | hxr3 | 1.09 | C5', C3', O4' | |
| H5' | C5' | h2xr2 | 1.09 | C4', O5' | |
| H5'' | C5' | h2xr2 | 1.09 | C4', O5' | |

#### Base-Specific H Atoms

**Adenine (DA, A)** — purine, glycosidic bond at N9:

| H | Parent | Type | Bond | Neighbors |
|---|--------|------|------|-----------|
| H2 | C2 | hxr2_planar (aromatic) | 1.08 | N1, N3 |
| H8 | C8 | hxr2_planar (aromatic) | 1.08 | N7, N9 |
| H61 | N6 | h3xr (sp2 NH2) | 1.02 | C6, C5 | angle 120°, dihedral 180° |
| H62 | N6 | h3xr (sp2 NH2) | 1.02 | C6, C5 | angle -120°, dihedral 180° |

**Guanine (DG, G)** — purine, glycosidic bond at N9:

| H | Parent | Type | Bond | Neighbors |
|---|--------|------|------|-----------|
| H8 | C8 | hxr2_planar | 1.08 | N7, N9 |
| H1 | N1 | hxr2_planar | 1.02 | C6, C2 |
| H21 | N2 | h3xr (sp2 NH2) | 1.02 | C2, N1 | angle 120°, dihedral 180° |
| H22 | N2 | h3xr (sp2 NH2) | 1.02 | C2, N1 | angle -120°, dihedral 180° |

**Cytosine (DC, C)** — pyrimidine, glycosidic bond at N1:

| H | Parent | Type | Bond | Neighbors |
|---|--------|------|------|-----------|
| H5 | C5 | hxr2_planar | 1.08 | C4, C6 |
| H6 | C6 | hxr2_planar | 1.08 | C5, N1 |
| H41 | N4 | h3xr (sp2 NH2) | 1.02 | C4, N3 | angle 120°, dihedral 180° |
| H42 | N4 | h3xr (sp2 NH2) | 1.02 | C4, N3 | angle -120°, dihedral 180° |

**Thymine (DT)** — pyrimidine, glycosidic bond at N1:

| H | Parent | Type | Bond | Neighbors |
|---|--------|------|------|-----------|
| H6 | C6 | hxr2_planar | 1.08 | C5, N1 |
| H3 | N3 | hxr2_planar | 1.02 | C4, C2 |
| H71 | C7 | h3xr (methyl) | 1.09 | C5, C4 | dihedral 180°, methyl rotator |
| H72 | C7 | h3xr | 1.09 | C5, C4 | dihedral -60° |
| H73 | C7 | h3xr | 1.09 | C5, C4 | dihedral 60° |

**Uracil (U)** — pyrimidine, glycosidic bond at N1:

| H | Parent | Type | Bond | Neighbors |
|---|--------|------|------|-----------|
| H5 | C5 | hxr2_planar | 1.08 | C4, C6 |
| H6 | C6 | hxr2_planar | 1.08 | C5, N1 |
| H3 | N3 | hxr2_planar | 1.02 | C4, C2 |

### Excluded H Atoms (polymer linking)

These are present in CCD but NOT placed for polymer residues:
- HOP2, HOP3 — phosphate H (consumed in phosphodiester bond)
- HO3' — 3' terminal only (next residue's phosphate binds here)
- HO5' — 5' terminal only (phosphate binds here)

The `IFNOPO4` flag in reduce controls this. In zreduce, these are simply omitted from the plans (same approach as backbone amide H for PRO).

### MoverHints

| H | MoverHint | Reason |
|---|-----------|--------|
| HO2' (RNA) | `.rotate` | OH rotatable, 12 orientations |
| H71/H72/H73 (DT methyl) | `.rotate_methyl` | CH3 rotatable, 3 orientations |
| All others | `.none` | Fixed geometry |

### Integration: `placer.zig` Changes

Currently `placer.zig` calls `standard.getPlans(comp_id)`. Add nucleotide lookup:

```zig
const plans = standard.getPlans(comp_id) orelse
    nucleotide.getPlans(comp_id) orelse
    // CCD fallback...
```

### Bond Lengths (from reduce/CCD)

| Type | reduce X-ray | reduce nuclear | Used |
|------|-------------|----------------|------|
| C-H aromatic | 0.93 | 1.08 | 1.08 |
| C-H sp3 | 0.97 | 1.09 | 1.09 |
| N-H sp2 | 0.86 | 1.02 | 1.02 |
| O-H | 0.84 | 0.98 | 0.98 |

We use nuclear (riding) distances, consistent with existing amino acid plans.

### H1' Glycosidic N Atom

The third neighbor for H1' placement differs by base type:
- Purines (A, G, DA, DG): N9
- Pyrimidines (C, U, DC, DT): N1

This means sugar plans cannot be 100% shared — H1' needs to know the glycosidic N. Solution: define sugar plans per residue (small duplication, but simple and explicit).

### File Structure

```zig
// src/place/nucleotide.zig

// DNA sugar helper (H1' uses glycosidic_n parameter)
fn dna_sugar(glycosidic_n) -> [7]PlacementPlan  // H1', H2', H2'', H3', H4', H5', H5''

// RNA sugar helper
fn rna_sugar(glycosidic_n) -> [8]PlacementPlan  // H1', H2', HO2', H3', H4', H5', H5''

// Per-residue: sugar + base concatenated
const da_plans: []const PlacementPlan  // 7 sugar + 4 base = 11
const dc_plans: []const PlacementPlan  // 7 sugar + 4 base = 11
const dg_plans: []const PlacementPlan  // 7 sugar + 4 base = 11
const dt_plans: []const PlacementPlan  // 7 sugar + 5 base = 12 (methyl)
const a_plans: []const PlacementPlan   // 8 sugar + 4 base = 12
const c_plans: []const PlacementPlan   // 8 sugar + 4 base = 12
const g_plans: []const PlacementPlan   // 8 sugar + 4 base = 12
const u_plans: []const PlacementPlan   // 8 sugar + 3 base = 11

pub fn getPlans(comp_id) -> ?[]const PlacementPlan
```

### Tests

1. All 8 nucleotides return non-null plans
2. Plan count per residue matches expected
3. DA/A have H61, H62 (adenine NH2)
4. DG/G have H1 (imino), H21, H22 (amino)
5. DC/C have H41, H42 (amino)
6. DT has H71/H72/H73 (methyl) with rotate_methyl hint
7. RNA residues have HO2' with rotate hint
8. DNA residues have H2'' (deoxyribose 2nd H)
9. No residue has HOP2, HOP3, HO3' (excluded polymer H)

### Verification

After implementation, compare output of:
```bash
# With CCD (current)
zreduce run structure_with_dna.cif -d components.cif -o with_ccd.cif

# Without CCD (new hardcoded templates)
zreduce run structure_with_dna.cif -o without_ccd.cif

# Diff hydrogen counts and positions
```

Hydrogen names and counts should match. Positions may differ slightly due to bond length/angle differences between hardcoded and CCD-derived values.
