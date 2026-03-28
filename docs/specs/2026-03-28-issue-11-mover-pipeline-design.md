# Issue #11: Mover Generation & Optimizer Pipeline Wiring (Rotators Only)

## Problem

After hydrogen placement, all standard-residue movable groups remain in their
initial template positions. The optimizer infrastructure exists (rotator/flipper
constructors, scorer, clique detection) but is not connected to the main pipeline.
`--no-opt` and `--no-flip` flags have no effect.

## Scope

Phase 1: Wire rotators only (rotate, rotate_methyl, rotate_nh3). Flippers
(amide, his) are deferred to Issue #17.

## Design

### Part A: Add `mover_hint` Field to Atom

Add `mover_hint: MoverHint = .none` to the Atom struct. Set it in
`appendHydrogen` from `plan.mover_hint`. This preserves the placement plan's
optimization intent through to the mover generation phase.

### Part B: Mover Generation

New file `src/optimize/mover_gen.zig` with:

```
pub fn generateMovers(allocator, model, no_flip: bool) ![]Mover
```

Algorithm:
1. Iterate all residues
2. For each residue, scan added H atoms for mover_hints
3. Group H atoms by (center_atom_name, mover_hint) — e.g., all 3 HB methyl H
   atoms on ALA share center CB
4. For each group, resolve atom indices and call the appropriate constructor:
   - `.rotate` → `createSingleHRotator(h_idx, center_idx, axis_idx)`
     - center = parent atom (e.g., OG for SER)
     - axis = atom bonded to center (e.g., CB for SER OG)
   - `.rotate_methyl` → `createMethylRotator(h_indices, center_idx, axis_idx)`
     - center = carbon with 3 H (e.g., CB for ALA)
     - axis = atom bonded to center (e.g., CA for ALA CB)
   - `.rotate_nh3` → `createNH3Rotator(h_indices, center_idx, axis_idx)`
     - center = NZ (for LYS)
     - axis = CE (for LYS)
   - `.flip_amide`, `.flip_his` → skip (Issue #17)
5. Return array of Movers

Center and axis atoms are resolved from the PlacementPlan's `connected` array,
which is stored in `standard.zig`. For each H atom, `connected[0]` is the center
(parent heavy atom) and `connected[1]` is the axis reference.

However, since we only have `mover_hint` on the atom (not the full plan), we
need to re-derive center/axis from the atom's name and residue context. This
can be done by looking up the standard plan for the residue and matching by
H atom name to find the corresponding connected atoms.

### Part C: Pipeline Wiring (main.zig)

After `addHydrogens` and before output writing:

```
if (!config.no_opt) {
    const movers = try generateMovers(allocator, &mdl, config.no_flip);
    defer { for (movers) |*m| m.deinit(); allocator.free(movers); }
    const opt_result = try optimizer.optimize(allocator, movers, &mdl, opt_config);
    // Print optimization summary
}
```

### Part D: Atom Index Resolution

Use `Model.findAtomInResidue(residue_idx, atom_name)` to resolve atom names
to global indices. This function already exists in the model.

For H atoms added by placement, they are appended after `res.atom_end`, so
we need to search the full atoms array for added H with matching residue_idx
and name.

## Changed Files

| File | Changes |
|------|---------|
| `src/model/atom.zig` | Add `mover_hint` field |
| `src/place/placer.zig` | Set `mover_hint` in `appendHydrogen` |
| `src/optimize/mover_gen.zig` | New: generateMovers function |
| `src/optimize.zig` | Export mover_gen |
| `src/main.zig` | Wire mover generation + optimizer call |

## Test Plan

1. `mover_hint` correctly set on placed H atoms
2. ALA produces 1 methyl_rotator for CB
3. SER produces 1 single_h_rotator for OG-HG
4. `generateMovers` returns correct count for multi-residue model
5. Optimizer execution changes at least one H position
6. Existing placement tests unaffected

## Out of Scope

- Amide flippers (Issue #17)
- His flippers (Issue #17)
- CellList-based scoring optimization
- JSON log population (separate enhancement)
