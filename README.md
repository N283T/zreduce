# zreduce

A from-scratch Zig implementation of the [reduce](https://github.com/rlabduke/reduce) hydrogen placement tool for macromolecular structures.

zreduce reads mmCIF structures, adds hydrogen atoms using geometric rules and CCD bond topology, then optimizes rotatable/flippable groups via clique-based search with probe dot-sphere scoring.

## Features

- **mmCIF-native** input and output (no PDB format dependency)
- **CCD-driven** hydrogen placement for non-standard residues (nucleotides, glycans, ligands, modified residues)
- **20 standard amino acids** with hardcoded placement plans and bond topology
- **6 geometry types**: tetrahedral, sp2, dihedral-controlled, planar bisector, fractional angle, linear
- **Conformer-aware**: altloc handling with per-conformer placement and optimization
- **Residue-specific chemistry**: donor/acceptor/aromatic/charge annotations for scoring
- **Chain-break detection**: from `_pdbx_poly_seq_scheme` for correct terminal chemistry
- **Dot-sphere scoring** with bump/H-bond classification (matching original reduce constants)
- **Clique-based optimization**: singleton, brute-force, iterative greedy, with fine angular search
- **Rotation movers**: OH/SH (12 orientations), NH3+ (3), methyl (3) — standard + CCD-derived
- **Flip movers**: Asn/Gln amide flip (2 orientations), His ring flip (6 orientations)
- **CellList spatial index**: O(N) scoring instead of O(N^2)
- **Multithreaded optimization**: parallel singleton and fine-search via thread pool
- **SIMD acceleration**: Vec3 operations via `@Vector(4, T)`, fast exp approximation
- **Batch processing**: parallel file-level processing with atomic work-stealing
- **Model validation**: sentinel detection, NaN/Inf coordinate checks
- **Zero-copy CIF parser** with streaming CCD support (~40K components)

## Quick Start

### Requirements

- [Zig](https://ziglang.org/) 0.14+ (tested with 0.15.2)

### Build

```bash
zig build                           # Debug
zig build -Doptimize=ReleaseFast    # Optimized (recommended for real structures)
```

### Run

```bash
# Single file processing
zreduce run input.cif -o output.cif

# With CCD dictionary (enables non-standard residue H placement + optimization)
zreduce run input.cif -d components.cif -o output.cif

# Placement only (skip optimization)
zreduce run input.cif -o output.cif --no-opt

# Disable flips (keep rotators only)
zreduce run input.cif -o output.cif --no-flip

# With validation diagnostics
zreduce run input.cif -o output.cif --validate

# Add water hydrogens with occupancy/B-factor filtering
zreduce run input.cif -o output.cif --water
zreduce run input.cif -o output.cif --water-phantom
zreduce run input.cif -o output.cif --water --water-occ-cutoff 0.5 --water-b-cutoff 30

# Write JSON optimization log
zreduce run input.cif -o output.cif --json log.json

# Force residue protonation states from a control file
zreduce run input.cif -o output.cif --protonation protonation.txt

# Dump mover IDs / allowed states, then force selected movers
zreduce run input.cif --dump-movers movers.txt --no-opt
zreduce run input.cif -o output.cif --fix fix.txt

# Batch processing (parallel)
zreduce batch input_dir/ -o output_dir/
zreduce batch input_dir/ -d components.cif --jsonl log.jsonl --protonation protonation.txt
zreduce batch input_dir/ -d components.cif --fix fix.txt
zreduce batch input_dir/ -j 4    # limit to 4 threads
```

### Test

```bash
zig build test --summary all    # 325 tests
```

## CLI

zreduce uses subcommands: `run` for single files, `batch` for directories.

### `zreduce run` — single file

| Flag | Description |
|------|-------------|
| `-h, --help` | Show help message |
| `-o, --output PATH` | Output mmCIF file (default: stdout) |
| `-d, --dict PATH` | Path to components.cif for non-standard residues |
| `--json PATH` | Write JSON optimization log |
| `--protonation PATH` | Residue protonation override file |
| `--fix PATH` | Force mover states from control file |
| `--dump-movers PATH` | Write available mover IDs/states to file |
| `--no-opt` | Skip optimization (placement only) |
| `--no-flip` | Disable Asn/Gln/His flips |
| `--validate` | Print detailed validation diagnostics |
| `--water` | Add water hydrogens |
| `--water-phantom` | Allow zero-occupancy phantom water hydrogens |
| `--water-occ-cutoff N` | Skip waters with occupancy below `N` |
| `--water-b-cutoff N` | Skip waters with B-factor above `N` |

### `zreduce batch` — directory

| Flag | Description |
|------|-------------|
| `-h, --help` | Show help message |
| `-o, --output PATH` | Output directory (default: `<input>_reduced/`) |
| `-d, --dict PATH` | CCD dictionary (loaded once, shared across files) |
| `-j, --threads N` | Thread count (default: auto-detect CPU count) |
| `--jsonl PATH` | Aggregated JSONL log file |
| `--protonation PATH` | Residue protonation override file |
| `--fix PATH` | Force mover states from control file |
| `--no-opt` | Skip optimization |
| `--no-flip` | Disable flips |
| `--quiet` | Suppress progress output |
| `--water` | Add water hydrogens |
| `--water-phantom` | Allow zero-occupancy phantom water hydrogens |
| `--water-occ-cutoff N` | Skip waters with occupancy below `N` |
| `--water-b-cutoff N` | Skip waters with B-factor above `N` |

### Global flags

| Flag | Description |
|------|-------------|
| `-h, --help` | Show help and subcommand list |
| `-V, --version` | Show version |

### Protonation Override File

One override per line:

```text
# chain:auth_seq[:ins_code] comp_id state
A:57 HIS HIE
A:102 ASP OD2
B:14 GLU DEPROTONATED
C:88 LYS NEUTRAL
D:5 CYS THIOLATE
```

Supported states:

- `HIS`: `AUTO`, `HID`, `HIE`, `HIP`
- `ASP`: `DEPROTONATED`, `OD1`, `OD2`
- `GLU`: `DEPROTONATED`, `OE1`, `OE2`
- `LYS`: `CHARGED`, `NEUTRAL`
- `CYS`: `THIOL`, `THIOLATE`

### Fix Override File

One override per line:

```text
# chain:auth_seq[:ins_code] comp_id target value
A:57 ASN amide FLIP
A:88 HIS his HID_FLIP
B:14 SER OG 6
```

Targets and values:

- `amide`: `ORIGINAL`, `FLIP`
- `his`: `HIE`, `HID`, `HIE_FLIP`, `HID_FLIP`
- rotators: center atom name plus coarse orientation index
  Example: `SER OG 6`, `ALA CB 2`, `LYS NZ 1`

## Architecture

```
mmCIF input
    |
    v
+------------------+
| CIF Parser       |  Zero-copy tokenizer + parser
+--------+---------+
         |
         v
+------------------+
| CCD Loader       |  Streaming parser for components.cif
+--------+---------+
         |
         v
+------------------+
| Model Build      |  Atom[], Residue[], Chain[]
| + Chemistry      |  Donor/acceptor/aromatic/charge annotations
| + Chain Breaks   |  From _pdbx_poly_seq_scheme
+--------+---------+
         |
         v
+------------------+
| H Placement      |  Conformer-aware, altloc-consistent
| Standard + CCD   |  Bond topology for neighbor resolution
+--------+---------+
         |
         v
+------------------+
| Mover Generation |  Standard plans + CCD topology fallback
+--------+---------+
         |
         v
+------------------+
| Optimizer        |  CellList scoring + clique search
| + Fine Search    |  Angular refinement around coarse best
+--------+---------+
         |
         v
+------------------+
| Validation       |  Sentinel, NaN/Inf checks
+--------+---------+
         |
         v
+------------------+
| Output Writer    |  mmCIF with H atoms + JSON log
+------------------+

Batch mode wraps the above pipeline with file-level parallelism:
  Directory scan → Thread pool → processFile per file → JSONL log
```

## Performance

Benchmarked on Apple Silicon (ReleaseFast):

| Structure | Residues | Movers | Time |
|-----------|----------|--------|------|
| AF small protein | 16 | 12 | 0.008s |
| AF medium protein | 309 | 299 | 0.03s |
| AF large protein | 1486 | 1320 | 0.39s |
| AF extra-large protein | 2339 | 2434 | 1.0s |

Batch processing: 4370 E. coli proteome structures in 72s (10.5x CPU utilization).

## Project Structure

```
src/
  main.zig              CLI entry point (subcommand dispatch: run/batch)
  run.zig               Single-file processing pipeline
  batch.zig             Batch processing (parallel, JSONL log)
  root.zig              Library re-exports
  validate.zig          Post-placement model validation
  math.zig              Vec3(T), rotation, dihedral
  element.zig           AtomType, VDW radii, AtomFlags
  cif/                  CIF parser subsystem
  mmcif.zig             atom_site + poly_seq_scheme + unobs atoms
  ccd.zig               CCD component dictionary
  model/                Molecular model structs
  place/                Hydrogen placement
    placer.zig          Unified placer (conformer-aware)
    standard.zig        20 AA placement plans
    het.zig             CCD-derived placement
    topology.zig        Bond topology tables
    chemistry.zig       Chemical annotations
  optimize/             Optimization engine
    optimizer.zig       Clique search + fine search + multithreaded optimization
    scoring.zig         CellList-based scoring with SoA layout
    mover_gen.zig       Mover generation (standard + CCD)
    scorer.zig          Dot-sphere scoring
    rotator.zig         Rotation movers
    flipper.zig         Flip movers
    mover.zig           Mover types
    clique.zig          Interaction graph
    dot_sphere.zig      Dot generation
  writer/               Output writers
examples/
  data/                 Input structures
  result/               zreduce output
```

## License

TBD

## References

- [Original reduce (C++)](https://github.com/rlabduke/reduce)
- [Word et al. (1999) J Mol Biol 285:1735-1747](https://doi.org/10.1006/jmbi.1998.2401)
