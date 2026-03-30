# Bond Graph Parsing and Integration Design

**Date:** 2026-03-30
**Issue:** #60
**Status:** Approved

## Overview

Parse `_struct_conn`, `_chem_comp_bond`/`_chem_comp_atom` (inline), and `_pdbx_entity_branch_link` from mmCIF structure files to build a bond graph. Use the bond graph to skip hydrogen placement on leaving atoms (e.g., disulfide SG-HG, glycosidic OH-H) and to replace distance-based fallback with explicit bond topology for non-standard residues.

## Architecture

```
mmCIF file
  ├─ _atom_site                → Model (existing)
  ├─ _chem_comp_atom           → InlineComponentDict (new)
  ├─ _chem_comp_bond           → InlineComponentDict (new)
  ├─ _struct_conn              → Model.bonds + AtomFlags (new)
  └─ _pdbx_entity_branch_link  → Model.bonds + AtomFlags (new)
                                      │
                            placer.zig: skip H if bonded_inter_residue flag set
```

## Components

### 1. Inline Component Dictionary (mmcif.zig)

Parse `_chem_comp_atom` and `_chem_comp_bond` loops embedded in structure files. Reuse existing `ccd.Component`, `ccd.CompAtom`, and `ccd.CompBond` structs to produce a `ComponentDict` compatible with the existing CCD interface.

**Priority:** Inline components take precedence over external CCD dictionary. CCD dictionary (`-d` flag) serves as fallback for components not present inline.

**Key columns — `_chem_comp_atom`:**
- `comp_id`, `atom_id`, `type_symbol`, `charge`

**Key columns — `_chem_comp_bond`:**
- `comp_id`, `atom_id_1`, `atom_id_2`, `value_order`, `pdbx_aromatic_flag`

### 2. _struct_conn Parser (mmcif.zig)

Parse inter-residue connections (disulfide, covalent, metallic coordination).

**Key columns:**
- `conn_type_id` — bond type (disulf, covale, metalc, etc.)
- `ptnr1_label_asym_id`, `ptnr1_label_comp_id`, `ptnr1_label_atom_id`, `ptnr1_label_seq_id`
- `ptnr2_label_asym_id`, `ptnr2_label_comp_id`, `ptnr2_label_atom_id`, `ptnr2_label_seq_id`
- `pdbx_value_order` — bond order (SING, DOUB, etc.)

**Processing:**
1. Parse loop rows into intermediate structs
2. Resolve atom identities to Model atom indices via lookup map
3. Append to `Model.bonds` with `source = .struct_conn`
4. Set `bonded_inter_residue` flag on both partner atoms

### 3. _pdbx_entity_branch_link Parser (mmcif.zig)

Parse glycan branch linkages with explicit leaving atom information.

**Key columns:**
- `entity_id`, `num_1`, `num_2`
- `comp_id_1`, `atom_id_1`, `leaving_atom_id_1`
- `comp_id_2`, `atom_id_2`, `leaving_atom_id_2`

**Processing:**
1. Parse loop rows
2. Resolve atom identities via entity/num/comp/atom lookup
3. Append to `Model.bonds` with `source = .branch_link`
4. Set `bonded_inter_residue` flag on leaving atoms

### 4. AtomFlags Extension (element.zig)

Add `bonded_inter_residue` bit to the existing `AtomFlags` packed struct. This flag indicates the atom participates in an inter-residue covalent bond and should not receive hydrogen placement.

### 5. Placer Integration (placer.zig)

In the hydrogen placement loop, check the parent heavy atom's `bonded_inter_residue` flag. If set, skip H placement for that atom. This handles:
- Disulfide bonds: CYS SG (skip HG)
- Glycosidic bonds: leaving oxygen (skip OH hydrogen)
- Other covalent modifications

### 6. Pipeline Changes (run.zig)

Updated execution order:

```
1. readFile()                    → mmCIF source text
2. cif.readString()              → CIF Document
3. mmcif.parseModel()            → Model
4. mmcif.parseInlineComponents() → InlineComponentDict (NEW)
5. mmcif.parseStructConn()       → Model.bonds + flags (NEW)
6. mmcif.parseBranchLinks()      → Model.bonds + flags (NEW)
7. place.applyChemistry()        → enrich atom annotations
8. place.addHydrogens()          → uses inline_dict ?? ccd_dict
9. optimize / validate / write
```

### 7. Dictionary Priority

```
For each residue comp_id:
  1. Check InlineComponentDict (from structure file _chem_comp_bond/_chem_comp_atom)
  2. If not found, check external CCD dictionary (if provided via -d flag)
  3. If neither available, fall back to distance-based neighbor resolution (existing)
```

## Atom Lookup Strategy

Build a HashMap during bond resolution: `(label_asym_id, label_seq_id, atom_name) → atom_index`. Construct from Model's chain/residue/atom hierarchy after `parseModel()`. Handle:
- Altloc: match primary conformer (blank or 'A')
- Branch entities: use entity_id + num for `_pdbx_entity_branch_link`

## topology.zig Deprecation

With inline `_chem_comp_bond` available from structure files, the hardcoded 20-AA bond tables in `topology.zig` become redundant. Mark as deprecated but retain for backward compatibility with files lacking `_chem_comp_bond`. Remove in a future cleanup pass.

## Testing Strategy

**Unit tests:**
- `_struct_conn` parsing with disulfide bond fixture
- `_pdbx_entity_branch_link` parsing with glycan fixture
- `_chem_comp_atom`/`_chem_comp_bond` inline parsing
- Atom lookup resolution
- `bonded_inter_residue` flag propagation

**Integration tests:**
- H placement skips on flagged atoms (disulfide SG, glycosidic O)
- Inline dict overrides CCD dict for same comp_id
- Fallback to CCD when inline missing

**Regression:**
- All existing 253+ tests pass unchanged
- Real-world validation with 5fyl (glycoprotein), disulfide-containing structures

## Reference Implementation

- `~/zig-cif-graph-parser/src/graph.zig` — Full Zig bond graph builder with atom lookup, struct_conn, branch_link parsing
- `mmcif-dict` CLI for dictionary definitions
- `pdb-mine` for real data queries
