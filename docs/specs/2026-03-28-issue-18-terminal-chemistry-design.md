# Issue #18: C-Terminal/OXT Handling & Terminal Charge Annotation

## Problem

N-terminal protonation exists but has no charge annotation. C-terminal handling
is completely absent: no OXT recognition, no C-terminal detection, no terminal
carboxylate charge annotation. This affects scoring quality for terminal residues.

Chain-break detection is out of scope (separate issue).

## Design

### Part A: C-Terminal Detection

In `addHydrogens`, detect C-terminal alongside existing N-terminal:

```
const is_cterm = (res_idx == chain.residue_end - 1);
```

### Part B: Terminal Chemistry Annotations

Add `getTerminalAnnotation(atom_name, is_nterm, is_cterm)` to `chemistry.zig`:

- **OXT** (C-terminal): O type + negative + acceptor
- **C-terminal backbone O**: add negative flag (carboxylate)
- **N-terminal backbone N**: add positive flag (NH3+)

Terminal annotations are applied as flag OR-merge on top of standard annotations,
not as replacement. This preserves existing donor/acceptor flags while adding
charge information.

### Part C: applyChemistry Extension

Extend `applyChemistry` to accept terminal state per residue. After applying
standard annotations, apply terminal annotations via flag merge:

```
atom.flags |= terminal_ann.flags  // OR merge, not replace
atom.element_type = terminal_ann.atom_type  // only if overridden
```

The function needs chain/residue context to determine terminal status. Either:
- Pass the full model (already available) and compute is_nterm/is_cterm internally
- The function already iterates residues, so it can check chain boundaries

### Part D: No H on OXT

OXT is an oxygen with no H in standard protein chemistry (deprotonated COO-).
Standard placement plans do not include OXT-bonded H, so no explicit suppression
is needed. The duplicate-H check (Issue #12) handles any pre-existing case.

## Changed Files

| File | Changes |
|------|---------|
| `src/place/chemistry.zig` | Add `getTerminalAnnotation`, flag merge constants |
| `src/place/placer.zig` | Extend `applyChemistry` with terminal detection |
| `src/test_data/` | CIF fixture with OXT atom |

## Test Plan

1. C-terminal detected correctly (last residue in chain)
2. OXT gets negative + acceptor flags
3. C-terminal backbone O gets negative flag added
4. N-terminal backbone N gets positive flag added
5. Standard annotations preserved after terminal annotation merge
6. No regressions

## Out of Scope

- Chain-break detection (separate issue)
- C-terminal protonation variants (e.g., COOH for low pH)
- HET group terminal handling
