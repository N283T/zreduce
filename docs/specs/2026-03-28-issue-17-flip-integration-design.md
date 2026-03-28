# Issue #17: Amide/His Flip Mover Integration

## Problem

Amide (ASN/GLN) and histidine flip movers exist as constructors in `flipper.zig`
but are not wired into the mover generation pipeline. `mover_gen.zig` skips
`flip_amide` and `flip_his` hints with `continue`. HIS placement only creates
HE2 but not HD1, which the His flipper needs.

## Scope

Wire existing flipper constructors into `mover_gen.zig`. Add HD1 to HIS
placement plans. Respect `--no-flip` flag. State-dependent chemistry updates
(donor/acceptor toggling) deferred to Issue #26.

## Design

### Part A: Add HD1 to HIS Placement Plans

Add a `PlacementPlan` for HD1 in `standard.zig` HIS plans:
- `hxr2_planar` type (same as HE2), bonded to ND1
- connected: CE1 and CG (ring neighbors of ND1)
- `mover_hint = .flip_his`
- Atom type: Hpol

### Part B: Amide Flipper Generation in mover_gen.zig

For `.flip_amide` hint atoms (ASN HD21/HD22, GLN HE21/HE22):
- Deduplicate by (residue_idx, mover_hint) since both H atoms trigger
- Find the 2 H atoms, the N atom (ND2 or NE2), O atom (OD1 or OE1),
  and C atom (CG or CD) by looking up related plans
- Call `flipper.createAmideFlipper`

Atom resolution for amide:
- H atoms: the two atoms with `flip_amide` hint in the residue
- N atom: `connected[0]` from the H's plan (ND2 for ASN, NE2 for GLN)
- C atom: `connected[1]` from the H's plan (CG for ASN, CD for GLN)
- O atom: for ASN → OD1, for GLN → OE1. Derive from comp_id.

### Part C: His Flipper Generation in mover_gen.zig

For `.flip_his` hint atoms (HIS HE2, HD1):
- Deduplicate by (residue_idx, mover_hint)
- Find heavy atoms: ND1, CD2, CE1, NE2 by name
- Find H atoms: HD1 and HE2 (may be absent → null)
- Call `flipper.createHisFlipper`

### Part D: --no-flip Flag

In `generateMovers`, when `no_flip` is true, skip `flip_amide` and `flip_his`
hints (already structured as a switch case, just add the check).

## Changed Files

| File | Changes |
|------|---------|
| `src/place/standard.zig` | Add HD1 plan to HIS |
| `src/optimize/mover_gen.zig` | Implement amide and His flipper generation |

## Test Plan

1. HIS placement now produces both HD1 and HE2
2. ASN generates 1 amide_flip mover
3. HIS generates 1 his_flip mover with 6 orientations
4. `no_flip=true` suppresses flip movers
5. Existing tests pass unchanged (HIS plan count updated)

## Out of Scope

- State-dependent donor/acceptor toggling (Issue #26)
- Ion-aware lock-down behavior
- Amide fixup geometry improvements beyond current flipper
