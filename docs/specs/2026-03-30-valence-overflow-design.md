# Valence Overflow Correction Design

**Date:** 2026-03-30
**Issue:** #64
**Status:** Approved

## Overview

When inter-residue bonds (from `_struct_conn` / `_pdbx_entity_branch_link`) connect to atoms in non-standard residues, the CCD template's bond orders may conflict, causing valence overflow. For example, a carbon with a C=O double bond in the CCD template that also forms an inter-residue covalent bond would have valence 5 (>4).

This affects `ccd_derive.analyzeBonds` which determines hybridization from template bond orders. The extra inter-residue bond is not counted, leading to incorrect hybridization (sp2 instead of sp3).

## Solution

Modify `analyzeBonds` in `ccd_derive.zig` to accept an `extra_bonds` count representing inter-residue bonds on the atom. This count is added to `total_bonds` and `heavy_neighbor_count`, and triggers double→single bond demotion when valence exceeds the element's maximum.

### Changes

1. **`analyzeBonds`**: Add `extra_bonds: u8` parameter. After counting template bonds, add `extra_bonds` to totals. If total valence exceeds max for the element, demote double→single bonds.

2. **`derivePlans`**: Add `inter_residue_atoms: []const [4]u8` parameter — names of atoms with `bonded_inter_residue` flag. For each heavy atom, count how many inter-residue bonds it has (from the list) and pass to `analyzeBonds`.

3. **`placer.zig`**: In the CCD-derived placement branch, collect atom names with `bonded_inter_residue` flag from the residue and pass to `derivePlans`.

### Valence Table

| Element | Max Valence |
|---------|-------------|
| H | 1 |
| C | 4 |
| N | 3 |
| O | 2 |
| S | 2 |
| P | 3 |

### Demotion Logic

When valence > max:
1. Find double bonds on the atom, demote to single (valence -= 1)
2. Find triple bonds, demote to double (valence -= 1)
3. Repeat until valence <= max or no more bonds to demote
4. If still over, adjust hybridization to sp3 (failsafe)
