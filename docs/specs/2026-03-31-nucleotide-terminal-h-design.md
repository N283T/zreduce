# Issue #96: Nucleotide 3'/5' Terminal H Placement

## Problem

1. **Missing HO3'**: 3' terminal nucleotides lack HO3' (4 cases in fold_test2)
2. **Extra HOP2**: Modified nucleotides (CCD derive) place HOP2 on mid-chain P (4 cases)

## Root Cause

- Hardcoded nucleotide plans intentionally exclude terminal H (HO3', HO5', HOP2, HOP3) — mid-chain assumption
- No nucleotide terminal detection exists (unlike AA `is_nterm`/`is_cterm`)
- Phosphodiester bonds are implicit polymer bonds (not in struct_conn), so `bonded_inter_residue` flag is never set on P/O3'/O5'

## Design

### 1. 3' Terminal HO3' Placement

Add `place3primeOH` in `placer.zig`, analogous to `placeNtermNH3`:

- **Detection**: `is_cterm` logic — last residue in chain or next has `is_chain_break_before`
- **Applies to**: standard nucleotides (DA, DC, DG, DT, A, C, G, U) + modified nucleotides with O3'
- **Geometry**: h3xr on O3' (tetrahedral, bond_len=0.97, using C3'-O3' axis)
- **MoverHint**: `.rotator` for OH optimization

### 2. Mid-chain HOP2 Suppression (CCD derive)

For non-terminal nucleotide residues processed via CCD derive, skip H atoms bonded to P:

- **Detection**: residue is nucleotide-like (has P atom) AND not 5' terminal
- **Implementation**: add P-bonded H names (HOP2, HOP3) to skip list in `addHydrogens` when calling `derivePlans`, OR filter out P-parent plans post-derive for non-terminal
- **Simpler approach**: in `addHydrogens`, after CCD derive, check if plan's parent is P and residue is non-terminal → skip. This avoids changing `derivePlans` interface.

### 3. Nucleotide Detection

Need a way to identify nucleotide residues (not just standard ones):
- Standard: matched by `nucleotide.getPlans()`
- Modified: falls through to CCD derive — detect by presence of atoms like P, O3', O5', C3', C4', C5' (sugar-phosphate backbone)

Simple heuristic: residue has atom named "P" → nucleotide-like.

### Changes

1. **`placer.zig` — `addHydrogens`**:
   - After standard nucleotide plan execution: if `is_cterm`, call `place3primeOH()`
   - In CCD derive path: if residue has P atom and is not 5'-terminal, filter out plans where parent is P

2. **`placer.zig` — `place3primeOH`** (new function):
   - Find O3', C3', C4' atoms
   - Place H via h3xr geometry (tetrahedral)
   - Set MoverHint.rotator for optimization

3. **No changes to**: `nucleotide.zig`, `ccd_derive.zig`, `mmcif.zig`

### Non-goals

- 5' terminal HO5' placement (zreduce and ChimeraX both skip on standard nucleotides)
- 5' terminal HOP2/HOP3 placement (separate consideration)
- Modifying struct_conn parsing or `bonded_inter_residue` for phosphodiester bonds
