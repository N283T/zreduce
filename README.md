# zreduce

A from-scratch Zig implementation of the [reduce](https://github.com/rlabduke/reduce) hydrogen placement tool for macromolecular structures.

zreduce reads mmCIF structures, adds hydrogen atoms using geometric rules and CCD bond topology, then optimizes rotatable/flippable groups via clique-based search with probe dot-sphere scoring.

## Features

- **mmCIF-native** input and output (no PDB format dependency)
- **CCD-driven** hydrogen placement for HET groups (ligands, modified residues)
- **20 standard amino acids** with hardcoded placement plans and mover hints
- **6 geometry types**: tetrahedral, sp2, dihedral-controlled, planar bisector, fractional angle, linear
- **Dot-sphere scoring** with bump/H-bond classification (matching original reduce constants)
- **Clique-based optimization**: singleton, brute-force, and greedy strategies
- **Rotation movers**: OH/SH (12 orientations), NH3+ (3), methyl (3)
- **Flip movers**: Asn/Gln amide flip (2 orientations), His ring flip (6 orientations)
- **Zero-copy CIF parser** with streaming CCD support (~40K components)
- **Spatial neighbor list** using counting-sort cell grid

## Quick Start

### Requirements

- [Zig](https://ziglang.org/) 0.14+ (tested with 0.15.2)
- zlib (for gzip CCD support, linked via system library)

### Build

```bash
zig build
```

### Run

```bash
# Basic: place hydrogens, write to stdout
zig build run -- input.cif

# Write to file
zig build run -- input.cif -o output.cif

# With CCD dictionary for HET groups
zig build run -- input.cif -d components.cif -o output.cif

# Write JSON log
zig build run -- input.cif -o output.cif --json log.json
```

### Test

```bash
zig build test
```

## CLI Options

| Flag | Description |
|------|-------------|
| `-h, --help` | Show help message |
| `-V, --version` | Show version |
| `-o, --output PATH` | Output mmCIF file (default: stdout) |
| `-d, --dict PATH` | Path to components.cif for CCD HET groups |
| `--json PATH` | Write JSON optimization log |
| `--no-opt` | Skip optimization (placement only) * |
| `--no-flip` | Disable Asn/Gln/His flips * |

\* Not yet implemented in v0.1.0 — optimization pipeline integration is planned for v0.2.0.

## Architecture

```
mmCIF input
    |
    v
+-------------+
| CIF Parser  |  Zero-copy tokenizer + parser
+------+------+
       |
       v
+-------------+
| CCD Loader  |  Streaming parser for components.cif
+------+------+
       |
       v
+-------------+
| Model Build |  Atom[], Residue[], Chain[], Bond[]
+------+------+
       |
       v
+-----------------+
| H Placement     |  Geometric placement (type 1-6)
| (Standard+HET)  |  Standard: hardcoded, HET: CCD-derived
+------+----------+
       |
       v
+-----------------+
| Optimizer       |  Dot-sphere scoring + clique search
+------+----------+
       |
       v
+-----------------+
| Output Writer   |  mmCIF with H atoms + JSON log
+-----------------+
```

## Project Structure

```
src/
  main.zig              CLI entry point
  root.zig              Library re-exports
  math.zig              Vec3(T), rotation, dihedral
  element.zig           AtomType, VDW radii, flags
  cif/                  CIF parser subsystem
  cif.zig               CIF module re-exports
  mmcif.zig             atom_site extraction
  ccd.zig               CCD component dictionary
  model/                Molecular model structs
  model.zig             Model module re-exports
  place/                Hydrogen placement
  place.zig             Place module re-exports
  optimize/             Optimization engine
  writer/               Output writers
  writer.zig            Writer module re-exports
  integration_test.zig  End-to-end tests
```

## Scoring Constants (from original reduce)

| Parameter | Value | Description |
|-----------|-------|-------------|
| dot_density | 16.0 dots/A^2 | Dot sphere density |
| probe_radius | 0.0 A | Probe sphere radius |
| gap_scale | 0.25 | Contact score gaussian width |
| bump_weight | 10.0 | Clash penalty multiplier |
| hb_weight | 4.0 | H-bond reward multiplier |
| bad_bump_gap_cut | 0.4 A | Bad bump threshold |

## Known Limitations (v0.1.0)

- Optimizer is not yet integrated into the main pipeline (placement only)
- CCD dihedral estimation uses fixed heuristic values (not computed from ideal coordinates)
- No PDB format support (mmCIF only)
- No Python bindings

## License

TBD

## References

- [Original reduce (C++)](https://github.com/rlabduke/reduce)
- [Word et al. (1999) J Mol Biol 285:1735-1747](https://doi.org/10.1006/jmbi.1998.2401)
