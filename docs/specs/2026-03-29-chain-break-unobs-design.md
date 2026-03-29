# Chain-Break Detection + Unobserved Atom Handling

## Problem

N-terminal detection uses only "first residue in chain array". Sequence gaps
(unobserved residues in the middle of a chain) are invisible — the residue
after a gap should be treated as N-terminal (NH3+) but currently gets a
peptide NH instead.

Additionally, `_pdbx_unobs_or_zero_occ_atoms` is not parsed. While unobserved
atoms are already absent from `_atom_site` (so placement naturally skips them),
explicit tracking enables better diagnostic logging.

## Design

### Part A: Parse `_pdbx_poly_seq_scheme` for Chain-Break Detection

After `_atom_site` parsing in `parseModel`, find the `_pdbx_poly_seq_scheme`
loop and scan it per chain:

- Extract `asym_id`, `seq_id`, `auth_seq_num`
- When `auth_seq_num` is `?` or `.`, the residue is unobserved
- Track which observed `seq_id` values have a gap before them
- Set `is_chain_break_before = true` on the corresponding model Residue

If `_pdbx_poly_seq_scheme` is absent, no chain breaks are detected (backward
compatible).

### Part B: Parse `_pdbx_unobs_or_zero_occ_atoms` for Diagnostics

Parse the loop and store a set of `(asym_id, seq_id, atom_name)` tuples as
unobserved atoms on the Model. This is used for:

- Diagnostic logging in `generateMovers` when a mover is skipped because a
  required heavy atom is missing
- Future use: distinguishing "atom not in structure" from "atom unobserved
  but expected"

Stored as a simple `UnobservedAtomSet` on the Model (hash set or sorted array).

### Part C: Extend N-terminal / C-terminal Detection

In `addHydrogens` and `applyChemistry`:

```
is_nterm = (res_idx == chain.residue_start) or res.is_chain_break_before
is_cterm = (res_idx == chain.residue_end - 1) or
           (res_idx + 1 < n_residues and mdl.residues.items[res_idx + 1].is_chain_break_before)
```

A residue followed by a chain break is C-terminal. A residue after a chain
break is N-terminal.

## Changed Files

| File | Changes |
|------|---------|
| `src/model/residue.zig` | Add `is_chain_break_before: bool` field |
| `src/mmcif.zig` | Parse `_pdbx_poly_seq_scheme` and `_pdbx_unobs_or_zero_occ_atoms` |
| `src/model/model.zig` | Add unobserved atom set (optional) |
| `src/place/placer.zig` | Extend is_nterm/is_cterm with chain-break awareness |

## Test Plan

1. CIF with sequence gap → residue after gap gets NH3+
2. CIF with sequence gap → residue before gap gets C-terminal charge
3. CIF without `_pdbx_poly_seq_scheme` → backward compatible
4. Existing single-residue and multi-chain tests pass unchanged

## Out of Scope

- Insertion code handling (Issue #28)
- Using unobserved residue list to add placeholder residues
- Nucleic acid chain-break handling
