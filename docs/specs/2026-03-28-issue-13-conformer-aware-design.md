# Issue #13: Conformer-Aware Standard-Residue Placement

## Problem

Atom lookup functions (`findAtomPos`, `findAtom`) ignore altloc, returning the
first name match regardless of conformer. In multi-conformer residues this
causes H placement to use mixed-conformer geometry (e.g., N from altloc A with
CB from altloc B), producing geometrically invalid hydrogen positions.

## Design

### Part A: altloc-Aware Atom Lookup

Add `target_altloc: u8` parameter to `findAtomPos` and `findAtom`. Fallback rule:

1. Prefer atom matching `(name, target_altloc)`
2. Fall back to atom matching `(name, altloc=' ')` (shared/blank)
3. Return null if neither found

This matches PDB convention: backbone atoms with blank altloc are shared across
all conformers, while side-chain atoms with specific altloc belong to one conformer.

### Part B: Conformer-Per-Residue Placement Loop

In `addHydrogens`, for each residue:

1. Collect distinct non-blank altloc values from the residue's atoms
2. If none found (all blank): execute plans once with `target_altloc=' '`
3. If conformers exist: execute plans once per conformer with the corresponding
   `target_altloc` ('A', 'B', etc.)

Each conformer's placed H atoms inherit the conformer's altloc via the existing
`ParentMeta` mechanism (Issue #12).

### Part C: Function Signature Changes

- `findAtomPos(mdl, res, name)` → `findAtomPos(mdl, res, name, target_altloc)`
- `findAtom(mdl, res, name)` → `findAtom(mdl, res, name, target_altloc)`
- `executePlan(mdl, res, res_idx, plan, bonds)` →
  `executePlan(mdl, res, res_idx, plan, bonds, target_altloc)`
- Bond-aware functions (`findBondedNeighbor`, etc.) gain `target_altloc` and
  pass it through to `findAtomPos`
- Distance-based fallback functions also gain `target_altloc`

### Part D: N-Terminal Conformer Handling

`placeNtermNH3` and `placeNtermNH2Pro` also need `target_altloc` to look up
the N atom in the correct conformer. They already use `findAtom` for N, which
will gain the altloc parameter.

## Changed Files

| File | Changes |
|------|---------|
| `src/place/placer.zig` | All lookup functions gain altloc param, conformer loop in addHydrogens |
| `src/test_data/` | CIF fixture with altloc A/B conformers |

## Test Plan

1. Multi-conformer fixture: ALA with shared backbone N, alternate CA/CB for A and B
2. Conformer A H uses conformer A geometry
3. Conformer B H uses conformer B geometry
4. Shared backbone atoms accessible from both conformers
5. Single-conformer (blank altloc) residues work unchanged
6. No regressions on existing tests

## Out of Scope

- Insertion code (icode) handling
- Residue splitting by conformer in the parser
- HET group conformer handling
