# Changelog

All notable changes to zreduce will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/),
and this project adheres to [Semantic Versioning](https://semver.org/).

## [0.2.0] - 2026-03-29

### Added

#### Placement Quality
- Duplicate hydrogen prevention with (name, altloc) duplicate key (#12)
- Parent atom metadata inheritance (altloc, occupancy, b_factor) for placed H (#12)
- Bond-graph-based neighbor inference replacing distance-only 1.9A cutoff (#14)
- Bond topology tables for all 20 standard amino acids (#14)
- Residue/atom-specific chemical typing (donor, acceptor, aromatic, charge) (#15)
- Terminal chemistry: C-terminal/OXT handling, N-terminal positive charge (#18)
- Conformer-aware placement with per-altloc atom lookup and blank fallback (#13)
- Insertion code (`pdbx_PDB_ins_code`) in residue boundary detection (#28)

#### Optimization Pipeline
- Mover generation from placed atoms using MoverHint metadata (#11)
- Standard residue rotator movers: OH/SH, methyl, NH3+ (#11)
- Amide flipper integration for ASN/GLN (#17)
- Histidine ring flipper integration with 6-state protonation model (#17)
- HD1 placement plan for HIS (previously only HE2) (#17)
- State-dependent donor/acceptor flags per His flip orientation (#26)
- Fine angular search after coarse optimization (#16)
- Iterative greedy for large cliques (replaces simple independent greedy) (#16)
- CellList spatial index for O(N) scoring, ~13x speedup on large structures (#35)
- CCD-based rotator generation for non-standard residues (nucleotides, glycans, ligands) (#36)
- `--no-opt` and `--no-flip` flags now functional

#### mmCIF Parsing
- Chain-break detection from `_pdbx_poly_seq_scheme` with gap-aware N/C-terminal chemistry
- `_pdbx_unobs_or_zero_occ_atoms` count for diagnostics
- Correct altloc output for added hydrogen atoms in mmCIF writer
- CIF value quoting for `[`, `]`, `{`, `}`, bare `?` and `.`

#### Quality & Diagnostics
- Model validation: sentinel position, NaN, Inf coordinate detection
- `--validate` CLI flag for detailed diagnostics
- Mover generation logging with skip count and per-skip warnings
- `isAbsentH` unified sentinel detection across codebase

### Fixed
- Flipper absent H atoms (sentinel position) excluded from output
- CIF quoting for SMILES containing `[` (e.g., `[Na+]`)
- Terminal annotation flag merge (OR) instead of type replacement
- Conformer iteration override of meta.altloc for shared backbone atoms

### Changed
- `applyChemistry` now detects chain-break terminals in addition to chain boundaries
- `executePlan` uses bond topology with distance-based fallback
- `generateMovers` accepts CCD dictionary for non-standard residue mover creation
- `scoreMover` uses CellList neighbor queries instead of all-atom iteration

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
- CCD-derived placement via bond topology hybridization analysis
- Unified placer: standard plans for 20 AA, CCD fallback for non-standard residues

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
