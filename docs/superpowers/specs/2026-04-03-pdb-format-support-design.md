# PDB Format Input/Output Support

## Summary

Add PDB format reading and writing to zreduce, enabling `zreduce run input.pdb -o output.pdb`.
Currently zreduce supports mmCIF only. The internal model is format-agnostic, so the
hydrogen placement and optimization engines require no changes.

## Decisions

| Decision | Choice | Rationale |
|----------|--------|-----------|
| Output format | PDB in -> PDB out only | Keep it simple; mmCIF users use mmCIF path |
| gzip | Supported (.pdb.gz) | Consistent with existing mmCIF gzip support |
| USER MOD records | Not emitted | `--validate` already provides diagnostics |
| CONECT records | Dropped on output | Unreliable in PDB; downstream tools don't depend on them |
| Multi-model | MODEL 1 only | mmCIF side is also single-model; revisit later |
| Bond information | topology.zig + external CCD | Original reduce also doesn't parse SSBOND for placement |
| Batch support | Not in initial scope | Will work naturally once run.zig supports PDB |

## Architecture

### Format Detection

Determined by file extension:
- `.pdb`, `.pdb.gz`, `.ent`, `.ent.gz` -> PDB path
- `.cif`, `.cif.gz`, `.mmcif`, `.mmcif.gz` -> mmCIF path (existing)

### New Files

#### `src/pdb.zig` — PDB Parser (~500-700 lines)

Responsibilities:
- Parse ATOM/HETATM records (fixed-width columns) into `Model`
- Handle Hybrid-36 serial and sequence numbers
- TER records mark chain boundaries
- MODEL/ENDMDL: extract MODEL 1 only, skip others
- Collect all non-ATOM lines as raw text for passthrough output

PDB ATOM fixed-width format (columns 1-based):
```
 1- 6  Record type ("ATOM  " or "HETATM")
 7-11  Serial number (Hybrid-36)
13-16  Atom name (4 chars, space-padded)
17     Altloc
18-20  Residue name (3 chars)
22     Chain ID (single char for classic PDB)
23-26  Residue sequence number (Hybrid-36)
27     Insertion code
31-38  X coordinate (8.3f)
39-46  Y coordinate (8.3f)
47-54  Z coordinate (8.3f)
55-60  Occupancy (6.2f)
61-66  B-factor (6.2f)
77-78  Element symbol (right-justified)
79-80  Charge
```

Key mappings to internal Model:
- `label_asym_id` / `auth_asym_id` = chain ID (single char, same for both)
- `label_seq_id` / `auth_seq_id` = residue sequence number
- `entity_id` = derived from chain grouping (sequential assignment)
- `entity_type` = heuristic: HOH -> water, standard AA/NA -> polymer, else non_polymer
- `is_chain_break_before` = detected from sequence number gaps within a chain

#### `src/writer/pdb_writer.zig` — PDB Writer (~300-400 lines)

Responsibilities:
- Output preserved non-ATOM records (passthrough)
- Generate ATOM/HETATM lines from Model atoms (original + added H)
- Serial number renumbering (sequential, Hybrid-36 if >99999)
- Drop CONECT and MASTER records
- gzip output support (reuse existing gzip infrastructure)

Output ordering per residue:
1. Original heavy atoms (preserved order)
2. Added hydrogen atoms (appended after heavy atoms of same residue)

#### Changes to `src/run.zig`

- Add format detection based on file extension
- PDB path: `pdb.parseModel()` instead of `cif.readString()` + `mmcif.parseModel()`
- PDB path skips: `_struct_conn`, `_pdbx_entity_branch_link`, inline components
- PDB path uses: `pdb_writer` instead of `mmcif_writer`
- Shared path: `applyChemistry()`, `addHydrogens()`, `generateMovers()`, `optimize()` unchanged

#### Changes to `src/main.zig`

- Update help text to mention PDB support
- No CLI flag changes needed (format auto-detected from extension)

### What Does NOT Change

- `src/place/` — All hydrogen placement logic (format-agnostic)
- `src/optimize/` — All optimization logic (format-agnostic)
- `src/model/` — Internal model structures
- `src/cif/` — CIF tokenizer/parser
- `src/mmcif.zig` — mmCIF-specific parsing
- `src/writer/mmcif_writer.zig` — mmCIF output

### Passthrough Record Storage

The PDB parser stores non-ATOM records as a list of tagged entries:

```
PdbRecord = union(enum) {
    atom_site: void,       // placeholder for where atoms go
    raw_line: []const u8,  // any non-ATOM line (HEADER, REMARK, HELIX, etc.)
    ter: void,             // TER record (regenerated on output)
}
```

The writer iterates this list, emitting raw lines as-is and inserting
Model atoms where `atom_site` markers appear.

### Entity Type Heuristic

Since PDB lacks `_entity.type`, classify residues by comp_id:
- HOH -> `.water`
- Standard 20 amino acids + standard nucleotides -> `.polymer`
- Everything else -> `.non_polymer`

This is sufficient for hydrogen placement decisions.

## Verification

### mmCIF Regression Check

After implementation, verify mmCIF path is unaffected:
1. Run existing test suite (`zig build test --summary all`)
2. Benchmark on AF model (309-residue) — confirm ~0.03s unchanged
3. Benchmark on large structure (2339-residue) — confirm ~1.0s unchanged
4. Batch benchmark on E. coli set — confirm throughput unchanged
5. Diff mmCIF output against pre-change baseline for identical results

### PDB-Specific Tests

- Round-trip: known PDB -> zreduce -> compare H positions with original reduce output
- Passthrough: non-ATOM records preserved exactly
- gzip: `.pdb.gz` input and output
- Hybrid-36: serial numbers > 99999
- Multi-model: only MODEL 1 processed
- Edge cases: missing element column, short lines, HETATM ligands
