# Issue #14: Bond-Graph-Based Reference Selection for Standard Residues

## Problem

Standard-residue hydrogen placement infers bonded neighbors using a fixed 1.9 A
distance cutoff (`findOtherNeighbor`, `findThirdNeighbor`, `findAtomBetween`).
This is fragile for distorted coordinates, crowded geometry, and alternate conformers.

HET placement already uses explicit bond topology from CCD components. Standard
residues need the same approach.

## Design

### Part A: Standard Residue Bond Topology Table

New file `src/place/topology.zig` provides hardcoded heavy-atom bond lists for
the 20 standard amino acids.

```
BondEntry = struct { a1: [4]u8, a2: [4]u8 }

fn getBonds(comp_id: []const u8) ?[]const BondEntry
```

- Covers all 20 standard amino acids
- Heavy atoms only (H bonds are implicit via PlacementPlan.connected[0])
- Backbone bonds (N-CA, CA-C, C-O) shared across all residues
- Side-chain bonds specific to each amino acid
- Returns null for unknown comp_id (triggers distance-based fallback)

### Part B: Bond-Aware Neighbor Query Functions

Replace the three distance-based functions in `placer.zig`:

1. `findBondedNeighbor(mdl, res, bonds, center_name, exclude_name) -> ?Vec3f32`
   - Replaces `findOtherNeighbor`
   - Searches `bonds` for entries matching `center_name`, returns position of
     the bonded partner that is not `exclude_name`

2. `findThirdBondedNeighbor(mdl, res, bonds, center_name, n1_name, n2_name) -> ?Vec3f32`
   - Replaces `findThirdNeighbor`
   - Searches `bonds` for entries matching `center_name`, returns position of
     the bonded partner that is not `n1_name` or `n2_name`

3. `findBondedAtomBetween(mdl, res, bonds, name1, name2) -> ?Vec3f32`
   - Replaces `findAtomBetween`
   - Searches `bonds` for an atom that is bonded to both `name1` and `name2`

Each function takes `[]const BondEntry` as parameter. When bonds are null
(non-standard residue without topology), the existing distance-based functions
serve as fallback.

### Part C: Integration

In `addHydrogens`, the standard-residue loop:
1. Calls `topology.getBonds(comp_id)` to get bond list
2. Passes bond list through to `executePlan`
3. `executePlan` dispatches to bond-aware or distance-based helpers depending
   on whether bonds are available

The `executePlan` signature gains an optional bonds parameter:
```
fn executePlan(mdl, res, res_idx, plan, bonds: ?[]const topology.BondEntry) !bool
```

### Part D: Distance Fallback

The original `findOtherNeighbor`, `findThirdNeighbor`, `findAtomBetween` remain
as private fallback functions. They are used only when `topology.getBonds` returns
null (unknown residue type reaching standard placement path, which should not
happen in practice but provides safety).

## Changed Files

| File | Changes |
|------|---------|
| `src/place/topology.zig` | New: amino acid bond topology table |
| `src/place/placer.zig` | Add bond-aware queries, update executePlan signature, wire topology |

## Test Plan

1. Topology table correctness: verify bond counts for ALA, GLY, TRP (simplest,
   smallest, largest side chains).
2. Bond-aware queries: unit tests with known residue geometry.
3. Distorted geometry: CIF fixture where 1.9 A cutoff would pick wrong neighbor
   but bond-based selection picks correctly.
4. Regression: existing placement tests continue to pass.

## Out of Scope

- Populating `Model.bonds` globally (reserved for inter-residue bonds, Issue #18)
- Parsing `_struct_conn` from mmCIF
- CCD-based topology for standard residues (CCD is optional)
