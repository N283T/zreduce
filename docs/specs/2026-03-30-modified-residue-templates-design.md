# Modified Residue Placement Templates Design

## Context

6 most common modified amino acids in the PDB require CCD dictionary for hydrogen placement. Hardcoding them removes the CCD dependency.

| comp_id | Name | Parent | PDB entries | Modification |
|---------|------|--------|-------------|--------------|
| MSE | Selenomethionine | MET | 10,254 | SD → SE |
| SEP | Phosphoserine | SER | 2,215 | OG-H → OG-PO3 |
| TPO | Phosphothreonine | THR | 1,753 | OG1-H → OG1-PO3 |
| CSO | S-hydroxycysteine | CYS | 1,209 | SG-H → SG-OH |
| PCA | Pyroglutamic acid | GLU | 975 | Cyclized (N-CD bond), no backbone H |
| PTR | O-phosphotyrosine | TYR | 905 | OH-H → OH-PO3 |

## Design

### New file: `src/place/modified.zig`

Same pattern as `nucleotide.zig` — uses `PlacementPlan` from `standard.zig`.

### Derivation from Parent Plans

Each modified residue's H plans are derived from the parent with specific changes:

**MSE (parent: MET)**
- HA (hxr3 CA, N, C, CB)
- HB2, HB3 (h2xr2 CB, CA, CG)
- HG2, HG3 (h2xr2 CG, CB, SE) — SD→SE name change
- HE1, HE2, HE3 (methyl CE, SE, CG) — SD→SE, rotate_methyl
- backbone H (h3xr N, CA, C)
- **Total: 9 plans** (same as MET minus the backbone_h being handled separately — actually 9 + backbone = 10. Wait, MET has 9 plans in standard.zig. Let me recheck.)

Actually, standard.zig MET plans include backbone_h. So MSE should too. Let me count MET:
- HA + 2xHB + 2xHG + 3xHE + H(backbone) = 9

MSE: identical layout, SE replaces SD:
- HA (hxr3 CA, N, C, CB)
- HB2, HB3 (h2xr2 CB, CA, CG)
- HG2, HG3 (h2xr2 CG, CB, SE)
- HE1, HE2, HE3 (methyl CE, SE, CG, rotate_methyl)
- H (backbone, h3xr N, CA, C)
- **Total: 9**

**SEP (parent: SER)**
SER has: HA + 2xHB + HG(OG rotator) + H(backbone) = 5
SEP: phosphorylation removes OG-H (HG), replaces with phosphate group.
- HA (hxr3 CA, N, C, CB)
- HB2, HB3 (h2xr2 CB, CA, OG)
- H (backbone)
- **Total: 4** (SER minus HG)

**TPO (parent: THR)**
THR has: HA + HB + HG1(OG1 rotator) + 3xHG2(methyl) + H(backbone) = 7
TPO: phosphorylation removes OG1-H (HG1).
- HA (hxr3 CA, N, C, CB)
- HB (hxr3 CB, CA, OG1, CG2)
- HG21, HG22, HG23 (methyl CG2, CB, CA, rotate_methyl)
- H (backbone)
- **Total: 6** (THR minus HG1)

**CSO (parent: CYS)**
CYS has: HA + 2xHB + HG(SG rotator) + H(backbone) = 5
CSO: S-hydroxylation. SG now bonded to OD. HG removed, HD added on OD.
- HA (hxr3 CA, N, C, CB)
- HB2, HB3 (h2xr2 CB, CA, SG)
- HD (h3xr OD, SG, CB, rotate — OH rotator)
- H (backbone)
- **Total: 5** (CYS: swap HG→HD on OD)

**PCA (parent: GLU)**
GLU has: HA + 2xHB + 2xHG + H(backbone) = 6
PCA: pyroglutamic acid — cyclized, N-CD bond. N has only 1 bond to backbone (no NH).
- HA (hxr3 CA, N, C, CB)
- HB2, HB3 (h2xr2 CB, CA, CG)
- HG2, HG3 (h2xr2 CG, CB, CD)
- **No backbone H** (N is in the ring, no NH)
- **Total: 5**

**PTR (parent: TYR)**
TYR has: HA + 2xHB + HD1 + HD2 + HE1 + HE2 + HH(OH rotator) + H(backbone) = 9
PTR: phosphorylation removes OH-H (HH).
- HA (hxr3 CA, N, C, CB)
- HB2, HB3 (h2xr2 CB, CA, CG)
- HD1, HD2 (aromatic planar)
- HE1, HE2 (aromatic planar)
- H (backbone)
- **Total: 8** (TYR minus HH)

### Integration

Same pattern as nucleotide.zig:
- `placer.zig`: fallback chain becomes standard → nucleotide → modified → CCD
- `mover_gen.zig`: `findPlanForH` searches modified plans too
- `place.zig`: exports modified module

### MoverHints

| Residue | H | Hint |
|---------|---|------|
| MSE | HE1/HE2/HE3 | .rotate_methyl |
| TPO | HG21/HG22/HG23 | .rotate_methyl |
| CSO | HD | .rotate |
| All | backbone H | .none |
| All others | | .none |

### Tests

1. All 6 residues return non-null plans
2. Plan counts per residue (MSE=9, SEP=4, TPO=6, CSO=5, PCA=5, PTR=8)
3. MSE has methyl rotator on CE (HE1/2/3)
4. CSO has OH rotator on OD (HD)
5. PCA has no backbone H
6. SEP/TPO/PTR have no phosphate OH (HOP2, HOP3 excluded)
7. Unknown returns null
