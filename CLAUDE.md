# CLAUDE.md — zreduce project instructions

## What is zreduce?

A Zig implementation of the [reduce](https://github.com/rlabduke/reduce) hydrogen placement tool for macromolecular structures (mmCIF format). Reads mmCIF, adds hydrogen atoms, optimizes orientations via clique-based scoring.

## Build & Test

```bash
zig build                           # Debug build
zig build -Doptimize=ReleaseFast    # Optimized build
zig build test --summary all        # Run all tests (253+)
```

## Run

```bash
# Single file
./zig-out/bin/zreduce run input.cif -o output.cif
./zig-out/bin/zreduce run input.cif -d components.cif -o output.cif  # with CCD
./zig-out/bin/zreduce run input.cif -o output.cif --validate         # with diagnostics

# Batch processing (parallel)
./zig-out/bin/zreduce batch input_dir/ -o output_dir/
./zig-out/bin/zreduce batch input_dir/ -d components.cif --jsonl log.jsonl
./zig-out/bin/zreduce batch input_dir/ -j 4  # limit to 4 threads
```

## Project structure

```
src/
  main.zig              CLI entry point (subcommand dispatch: run/batch)
  run.zig               Single-file processing pipeline
  batch.zig             Batch processing (parallel, JSONL log)
  root.zig              Library re-exports
  validate.zig          Post-placement model validation
  math.zig              Vec3, rotation, dihedral
  element.zig           AtomType, VDW radii, AtomFlags, mergeFlags
  cif/                  CIF parser (zero-copy tokenizer)
  mmcif.zig             _atom_site + _pdbx_poly_seq_scheme + _pdbx_unobs_or_zero_occ_atoms
  ccd.zig               CCD component dictionary (streaming parser)
  model/                Atom, Residue, Chain, Bond, CellList
  place/                Hydrogen placement
    placer.zig          Main placement logic (conformer-aware, altloc)
    standard.zig        20 AA placement plans + MoverHint
    het.zig             CCD-derived placement with hybridization analysis
    topology.zig        Bond topology tables for 20 AAs
    chemistry.zig       Residue/atom-specific chemical annotations
  optimize/             Optimization engine
    optimizer.zig       Clique-based search + fine angular search + multithreaded optimization
    scoring.zig         CellList-based scoring with SoA layout + centroid early-exit
    scorer.zig          Dot-sphere bump/H-bond scoring
    mover.zig           Mover struct, Orientation, isAbsentH
    rotator.zig         OH/SH, NH3+, methyl rotators + fine orientations
    flipper.zig         Asn/Gln amide flip, His ring flip
    mover_gen.zig       Mover generation from placed atoms (standard + CCD)
    clique.zig          Interaction graph + connected components
    dot_sphere.zig      Concentric ring dot generation
  writer/               mmCIF + JSON output
  test_data/            Test CIF fixtures
examples/
  data/                 Input structures (AF models, fold_test2)
  result/               zreduce output
```

## Key conventions

- **Atom names**: PDB 4-char space-padded (`[4]u8`). Use `nameSlice()` for trimmed view.
- **Altloc**: `' '` = blank/shared. Conformer-aware: `findAtomPos(mdl, res, name, target_altloc)`.
- **MoverHint**: Stored on Atom from PlacementPlan. Drives mover generation.
- **ParentMeta**: Inherits altloc/occupancy/b_factor from parent heavy atom to placed H.
- **Bond topology**: `topology.zig` for standard AAs, CCD for non-standard.
- **Chemistry**: `chemistry.zig` for donor/acceptor/aromatic/charge annotations.
- **Terminal detection**: chain boundary + `is_chain_break_before` from `_pdbx_poly_seq_scheme`.
- **Sentinel**: `ABSENT_H_POS = (1000,1000,1000)` for flipper absent H. Use `isAbsentH()`.

## Testing

- `zig build test` runs all tests (~225)
- Test fixtures in `src/test_data/`: tiny.cif, multi_chain.cif, ala_with_h.cif, ala_altloc.cif, ala_stretched.cif, ala_cterm.cif, gap_chain.cif, ins_code.cif, asn.cif, his.cif
- Real-world test data in `examples/data/` (AF models, fold_test2 with DNA/RNA/glycan)

## Performance

- ReleaseFast build for benchmarks
- CellList spatial index + centroid early-exit + SoA layout for scoring
- SIMD Vec3 + fast exp approximation
- Multithreaded singleton/fine-search optimization
- ~0.03s for 309-residue, ~1.0s for 2339-residue (Apple Silicon)
- Batch: 4370 E. coli structures in 72s (file-level parallelism)

## Open areas

- CCD dihedral estimation uses fixed heuristics (not computed from ideal coords)
- Multi-model JSON log omits per-model mover details (movers freed per model)
- State-dependent chemistry updates for flip movers are per-orientation flags only
