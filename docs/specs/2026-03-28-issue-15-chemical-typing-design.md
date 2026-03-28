# Issue #15: Residue/Atom-Specific Chemical Typing for Standard Residues

## Problem

After mmCIF parsing, heavy atoms receive only element-based typing (all C → `.C`,
all N → `.N`, all O → `.O`). The `flags` field on atoms defaults to all-false.
This means:
- The scorer's `isHBond()` cannot identify H-bond donor/acceptor pairs on heavy atoms
- No distinction between carbonyl C (smaller radius) and sp3 C
- No aromatic, charge, or acceptor flags on standard residue heavy atoms

Placed hydrogens receive `atom_type` from plans but `flags` are not set (except
N-terminal H which explicitly sets `donor`).

## Design

### Part A: Chemical Annotation Table

New file `src/place/chemistry.zig` provides a lookup table keyed by
`(comp_id, atom_name)` returning chemical properties.

```
ChemAnnotation = struct {
    atom_type: element.AtomType,
    flags: element.AtomFlags,
}

fn getAnnotation(comp_id: []const u8, atom_name: [4]u8) ?ChemAnnotation
```

Annotations cover:
- **Backbone** (all residues): C → C_eq_O/acceptor, O → acceptor, N → donor
- **Charged atoms**: ASP OD1/OD2 → negative, GLU OE1/OE2 → negative,
  LYS NZ → positive, ARG NH1/NH2/NE/CZ → positive
- **Aromatic rings**: PHE/TYR/TRP/HIS ring carbons → Car/aromatic,
  HIS ND1/NE2 → Nacc/aromatic/acceptor
- **Acceptors**: SER OG, THR OG1, TYR OH, ASN OD1, GLN OE1, MET SD, CYS SG
- **Carbonyl-adjacent**: ASN CG, GLN CD → C_eq_O

### Part B: Apply Chemistry to Model

New public function in `placer.zig`:

```
pub fn applyChemistry(mdl: *Model) void
```

- Iterates all residues in the model
- For each heavy atom in a standard residue, looks up annotation
- Updates `element_type`, `flags`, and `vdw_radius` on the atom
- Called from `main.zig` before `addHydrogens`

### Part C: Fix Hydrogen Flags

Update `appendHydrogen` to set `flags` from `plan.atom_type.info().flags`.
Currently only N-terminal H gets flags set explicitly.

## Changed Files

| File | Changes |
|------|---------|
| `src/place/chemistry.zig` | New: annotation table and lookup |
| `src/place/placer.zig` | `applyChemistry` function, `appendHydrogen` flags fix |
| `src/place.zig` | Export chemistry module |
| `src/main.zig` | Call `applyChemistry` before `addHydrogens` |

## Test Plan

1. Backbone C typed as C_eq_O after applyChemistry
2. Backbone O has acceptor flag
3. ASP OD1/OD2 have negative + acceptor flags
4. LYS NZ has positive + donor flags
5. PHE ring carbons have aromatic flag and Car type
6. Placed H atoms have correct flags from element table
7. All existing tests pass unchanged

## Out of Scope

- Flip/fixup chemistry state changes (Issue #17)
- N-terminal / C-terminal charge annotation (Issue #18)
- HET group chemistry annotation
