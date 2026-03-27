# Changelog

All notable changes to zreduce will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/),
and this project adheres to [Semantic Versioning](https://semver.org/).

## [0.1.0] - 2026-03-26

### Added

#### Phase 1: Foundation
- CIF tokenizer and parser with zero-copy design (Document/Block/Loop/Pair)
- CIF value helpers (null detection, float/int parsing with uncertainty suffix)
- Vec3 math library with rotation, dihedral angle, and angle computation
- Element table with comptime VDW radii and atom flags (45 atom types)
- Molecular model structs (Atom, Residue, Chain, Bond, Model)
- Spatial neighbor list using counting-sort cell grid (CellList)
- mmCIF atom_site extraction with chain/residue boundary tracking

#### Phase 2: Hydrogen Placement
- Type 1-6 hydrogen placement geometry functions
- Hardcoded placement plans for all 20 standard amino acids
- Streaming CCD component dictionary parser (~40K blocks without full Document tree)
- CCD-derived HET group placement via bond topology hybridization analysis
- Unified placer: standard plans for 20 AA, CCD fallback for HET groups

#### Phase 3: Optimization Engine
- Dot sphere generation with concentric ring algorithm (matching original reduce)
- Dot-sphere bump/H-bond scorer with contact/clash/H-bond classification
- OH/SH rotation movers (12 orientations at 30 degree intervals)
- NH3+ rotation movers (3 orientations)
- Methyl rotation movers (3 orientations)
- Asn/Gln amide flip movers (2 orientations, penalty 0.0/0.5)
- His ring flip movers (6 orientations: 2 flip states x 3 protonation states)
- Interaction graph construction and connected component detection
- Clique optimizer: singleton, brute-force (up to 100K states), greedy fallback

#### Phase 4: Output + CLI
- mmCIF writer with atom_site loop and custom _zreduce_log category
- JSON log writer for optimization results
- Full CLI pipeline: parse args, read mmCIF, load CCD, place H, write output
- CLI flags: -o (output), -d (CCD dict), --json (log), --no-opt, --no-flip
- End-to-end integration tests (placement, multi-chain, bond lengths, round-trip)

### Known Limitations
- Optimizer not yet integrated into main pipeline (placement only in v0.1.0)
- CCD dihedral estimation uses fixed heuristic values
- `--no-opt` and `--no-flip` flags accepted but not yet functional
- No PDB format support
