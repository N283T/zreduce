# Issue #12: Avoid Duplicate Standard-Residue Hydrogens and Preserve Metadata

## Problem

Standard-residue hydrogen placement has two correctness issues:

1. **Duplicate hydrogens**: Plans are always executed without checking if the target
   hydrogen already exists in the model. Partially protonated inputs get duplicate atoms.
2. **Lost metadata**: `appendHydrogen()` and `appendNtermH()` use hardcoded defaults
   (`occupancy=1.0`, `b_factor=0.0`, `altloc=' '`) instead of inheriting from the
   parent heavy atom.

HET placement already solves the duplicate problem via `collectAtomNames()` +
`nameExists()`, but standard-residue placement does not use this pattern.

## Design

### Part A: Duplicate Hydrogen Prevention

Apply the proven HET pattern to standard residues with altloc awareness.

**Duplicate detection key:** `(atom_name, altloc)` pair.

- `altloc=' '` existing H blocks only `altloc=' '` new H
- `altloc='A'` existing H blocks only `altloc='A'` new H
- Different altlocs are independent (not duplicates)

**Implementation:**

1. Before executing standard-residue plans, collect existing atom `(name, altloc)` pairs
   for the residue.
2. In `executePlan()`, after resolving the base atom and determining the target altloc,
   check if `(h_name, altloc)` already exists. If so, skip placement.
3. Apply the same check in N-terminal hydrogen placement (`placeNtermNH3`, `placeNtermNH2Pro`).

### Part B: Metadata Inheritance

Newly placed hydrogens inherit `altloc`, `occupancy`, and `b_factor` from the plan's
base atom (the heavy atom the hydrogen is bonded to).

**Implementation:**

1. Add a metadata struct or parameters to `appendHydrogen()` and `appendNtermH()` for
   `altloc`, `occupancy`, and `b_factor`.
2. In `executePlan()`, resolve the base atom via `findAtomPos()`, then extract its
   metadata. Pass it through to `appendHydrogen()`.
3. In `placeNtermNH3()` / `placeNtermNH2Pro()`, use the N atom's metadata for the
   placed H atoms.
4. The base atom is the first reference atom in each `PlacementPlan` — this is the
   atom the hydrogen is directly bonded to.

### Changed Files

| File | Changes |
|------|---------|
| `src/place/placer.zig` | Add duplicate check, extend `appendHydrogen`/`appendNtermH` signatures, pass metadata through `executePlan` and N-term placement |
| `src/place/placer.zig` (tests) | Add tests for duplicate prevention and metadata inheritance |
| `src/test_data/` | Add test CIF with pre-existing H and/or altlocs if needed |

### Test Plan

1. **Duplicate prevention**: Input with existing H atoms → no duplicate H added.
2. **Altloc-aware duplicates**: Input with altloc A/B H atoms → correct per-altloc handling.
3. **Metadata inheritance**: New H atoms inherit `occupancy`, `b_factor`, `altloc` from
   parent heavy atom.
4. **N-terminal metadata**: N-terminal H atoms inherit metadata from the N atom.
5. **Regression**: Existing ALA placement test continues to pass unchanged.

### Out of Scope

- Conformer-aware placement (Issue #13) — this issue focuses on metadata preservation,
  not on splitting residue contexts by conformer.
- Bond-graph-based neighbor inference (Issue #14).
