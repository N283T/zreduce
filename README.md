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
# Basic: place hydrogens on standard amino acids
zreduce input.cif -o output.cif

# With CCD dictionary (enables non-standard residue H placement + optimization)
zreduce input.cif -d components.cif -o output.cif

# Placement only (skip optimization)
zreduce input.cif -o output.cif --no-opt

# Disable flips (keep rotators only)
zreduce input.cif -o output.cif --no-flip

# With validation diagnostics
zreduce input.cif -o output.cif --validate

# Write JSON optimization log
zreduce input.cif -o output.cif --json log.json
```

### Test

```bash
zig build test --summary all    # 225+ tests
```

## CLI Options

| Flag | Description |
|------|-------------|
| `-h, --help` | Show help message |
| `-V, --version` | Show version |
| `-o, --output PATH` | Output mmCIF file (default: stdout) |
| `-d, --dict PATH` | Path to components.cif for non-standard residues |
| `--json PATH` | Write JSON optimization log |
| `--no-opt` | Skip optimization (placement only) |
| `--no-flip` | Disable Asn/Gln/His flips |
| `--validate` | Print detailed validation diagnostics |

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
```

## Performance

Benchmarked on Apple Silicon (ReleaseFast):

| Structure | Residues | Atoms | Movers | Time |
|-----------|----------|-------|--------|------|
| AF small protein | 16 | ~200 | 12 | <1s |
| AF medium protein | 309 | ~3K | 299 | 1s |
| AF large protein | 2339 | ~18K | 2434 | ~3s |
| Protein+DNA+RNA+glycan (CCD) | 507 | ~10K | 696 | 6s |

## Project Structure

```
src/
  main.zig              CLI entry point
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
    optimizer.zig       Clique search + fine search + CellList scoring
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
