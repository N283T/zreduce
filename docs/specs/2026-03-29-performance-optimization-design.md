# Performance Optimization Design

## Context

zreduce processes 2339-residue structures in 3.1 seconds (ReleaseFast). Profiling shows 88% of time in `scoreMover` (optimizer), 11% in mmCIF writer. Target: 0.2-0.3 seconds (10-15x speedup).

### Profile Breakdown (AF-P76347, 2339 residues, 2434 movers)

| Component | Time % | Root Cause |
|-----------|--------|------------|
| scoreMover mover-vs-mover loop | ~60% | O(M) scan of all movers per scoring call |
| scorePair math (sqrt, exp) | ~25% | Scalar distance computation, libc exp() |
| CellList neighbor query | ~5% | Already spatial-indexed, minor |
| mmCIF writer | ~11% | std.fmt float formatting |

### Current Bottleneck Detail

`scoreMover` iterates all 2434 movers for every atom of the current mover being scored. With ~12 orientations per singleton and 2434 singletons, this is ~2434 x 12 x 2434 x 2 = ~142M distance checks in the mover-vs-mover loop alone. Most of these are between movers that are far apart and will score 0.0.

## Design

### Phase 1: Mover-vs-Mover Early-Exit via Residue Centroid

**Goal:** Eliminate 95%+ of mover-vs-mover distance computations.

**Approach:**
- Precompute a representative position (centroid) for each mover from its parent heavy atom positions at optimize start
- In `scoreMover`, before iterating another mover's atom_indices, check distance between mover centroids
- Skip if centroid distance > `search_radius + max_displacement_margin` (search_radius=5.0, margin=1.5 for H bond length)
- Store centroid array in `ScoreContext` as `mover_centroids: []Vec3(f32)`

**Why centroid, not CellList for movers:**
- Mover atoms move with each orientation trial. CellList would need per-orientation updates (the problem we just fixed).
- Centroids are based on parent heavy atoms which do not move during optimization, so they are computed once.
- H atoms are at most ~1.5A from parent, so centroid + 1.5A margin is conservative and sufficient.

**Changes:**
- `buildScoreContext`: compute and store `mover_centroids[]` from mover parent atoms
- `scoreMover`: add centroid distance check before mover-vs-mover inner loop
- No API changes, no test behavior changes (same results, just faster)

**Expected:** 3.1s -> ~1.3s

### Phase 2: Multithreaded Singleton and Fine Search

**Goal:** Utilize all CPU cores for independent mover optimization.

**Approach:**
- Singleton movers and fine-search movers are independent: each writes only to its own `atom_indices` (disjoint, verified by debug assert)
- `ScoreContext` is read-only during scoring -- safe to share across threads
- Use `std.Thread.Pool` to dispatch work

**Thread safety design:**
- Each thread gets its own `scratch: ArrayListUnmanaged(u32)` buffer
- `model.atoms.items` is shared; writes are disjoint by mover (no locks needed)
- brute-force and greedy cliques remain sequential (they mutate multiple movers' atoms)

**Execution model:**
- Collect all singleton cliques into a work list
- Dispatch to thread pool with `spawnWg` (Zig's WaitGroup pattern)
- After all singletons complete, apply best orientations
- Dispatch fine search the same way
- brute-force/greedy cliques run on the main thread before singletons

**Changes:**
- `optimize()`: restructure clique processing order (brute-force/greedy first, then parallel singletons)
- Add thread pool initialization, per-thread scratch allocation
- Fine search loop dispatched to thread pool

**Expected:** ~1.3s -> ~0.5s (assuming 8 cores, ~3-4x on the parallel portion)

### Phase 3: SIMD Vectorization

**Goal:** Speed up Vec3 math and scoring computations.

**3a: Vec3 SIMD backend**
- Change `Vec3(f32)` internal representation from `{x, y, z}` to `@Vector(4, f32)` with w=0 padding
- `sub()`, `dot()`, `length()` map to hardware SIMD (NEON on ARM64, SSE on x86)
- Public API unchanged -- all callers work without modification
- 4th lane (w) stays 0; `dot` uses `@reduce(.Add, self.v * other.v)` which includes the 0-contributing w lane (no harm)

**3b: Fast exp() approximation for contact scoring**
- `scorePair` uses `@exp(-ratio * ratio)` for contact scores (gap in 0..0.5A range)
- Replace with 2nd-order polynomial or Schraudolph's fast exp approximation
- Contact scoring is relative ranking only; ~1% approximation error is acceptable
- Validate: run full test suite + compare output on benchmark structures

**Changes:**
- `math.zig`: Vec3 internal representation
- `optimizer.zig` or `math.zig`: fast_exp utility function
- `scorePair`: use fast_exp instead of `@exp`

**Expected:** ~0.5s -> ~0.3s

### Phase 4: SoA Scoring Arrays

**Goal:** Improve cache efficiency in scoring inner loops.

**Approach:**
- Scoring only needs `pos`, `vdw_radius`, `flags` per atom (17 bytes)
- Current: fetches full Atom struct (~64+ bytes) per neighbor, wasting cache lines
- Add to `ScoreContext`:
  - `static_radii: []f32` -- VDW radii for static atoms
  - `static_flags: []element.AtomFlags` -- flags for static atoms
- `scorePair` becomes `scorePairSoA(pos_a, r_a, flags_a, pos_b, r_b, flags_b)`
- Static positions already stored; add radii and flags alongside

**Changes:**
- `ScoreContext`: add `static_radii`, `static_flags` fields
- `buildScoreContext`: extract radii/flags from atoms
- `scoreMover`: pass SoA data to scoring instead of `atoms[oi]`
- Atom struct itself is unchanged (no ripple through parser/writer/placer)

**Expected:** ~0.3s -> ~0.25s

### Phase 5: Writer Optimization

**Goal:** Reduce mmCIF output time from 11% to ~3%.

**Approach:**
- Implement fixed-precision float-to-string for common mmCIF fields:
  - Coordinates: 3 decimal places (multiply by 1000, integer division)
  - B-factor: 2 decimal places
  - Occupancy: 2 decimal places
- Bypass `std.fmt.formatFloat` which handles general cases (scientific notation, NaN, etc.)
- Direct write to output buffer with integer arithmetic

**Changes:**
- `writer/mmcif_writer.zig`: add `writeFixedFloat3` and `writeFixedFloat2` helpers
- Replace `std.fmt.formatFloat` calls in atom_site loop

**Expected:** ~0.25s -> ~0.2s

## Benchmark Plan

Measure after each phase on all test structures:

| Structure | Residues | Movers | Baseline |
|-----------|----------|--------|----------|
| AF-C1P619 | 16 | 12 | 0.008s |
| AF-P0A9J6 | 309 | 299 | 0.07s |
| AF-P22523 | 1486 | 1320 | 1.1s |
| AF-P76347 | 2339 | 2434 | 3.1s |

## Implementation Order

Each phase is a separate PR:
1. Phase 1 (early-exit) -- highest ROI, simplest change
2. Phase 2 (multithread) -- second highest ROI
3. Phase 3 (SIMD) -- moderate ROI, low risk
4. Phase 4 (SoA) -- incremental improvement
5. Phase 5 (writer) -- polish

## Risks and Mitigations

| Risk | Mitigation |
|------|-----------|
| Centroid margin too tight -> missed interactions | Use conservative margin (2.0A); compare outputs before/after |
| Thread safety bugs in atom position writes | Disjoint write guarantee via debug assert; run tests with ThreadSanitizer |
| SIMD Vec3 padding changes dot product results | w=0 contributes 0 to dot; validate with existing 233 tests |
| Fast exp approximation changes optimization results | Compare full outputs on all benchmark structures; allow <0.1% score deviation |
| Writer fixed-point loses precision | Use 3 decimal places (same as PDB format standard) |
