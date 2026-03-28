# Issue #16: Optimizer Upgrade — Fine Search + Iterative Greedy

## Problem

The optimizer uses coarse discrete orientations only (12 for single-H, 3 for
methyl/NH3). The greedy fallback for large cliques optimizes each mover
independently, ignoring inter-mover interactions. Both limit hydrogen
placement quality.

## Scope

Phase 1: Fine angular search after coarse optimization, and iterative greedy
for large cliques. Deferred: score caching, vertex-cut decomposition,
delete-mask/conditional suppression (Issue #17 or later).

## Design

### Part A: Fine Search

After coarse optimization (singleton or brute-force), refine each mover's
best orientation with a fine angular search around the coarse best.

**For single_h_rotator** (coarse: 12 orientations at 30° intervals):
- Fine search: ±15° around best, at 5° steps = 6 fine positions
- Generates positions by rotating H around the same axis at fine offsets

**For methyl_rotator / nh3_rotator** (coarse: 3 orientations at 60° intervals):
- Fine search: ±30° around best, at 10° steps = 6 fine positions
- All 3 H atoms rotated together

**Implementation:**
- Add `refineSingleH(allocator, atoms, mover, best_angle)` to `rotator.zig`
- Add `refineGroupRotator(allocator, atoms, mover, best_angle)` to `rotator.zig`
- Both return temporary fine `[]Orientation` for scoring
- `optimizer.zig` adds `fineSearchMover()` that scores fine orientations and
  updates `best_orientation` if improvement found
- Fine search runs after coarse optimization for singletons and after
  brute-force for each mover in the clique

### Part B: Iterative Greedy

Replace `optimizeGreedy` (independent singleton per mover) with iterative:

1. Initialize all movers to orientation 0
2. For each mover in clique: optimize as singleton (others fixed at current best)
3. Repeat until convergence or max 3 iterations
4. Apply fine search to each mover after convergence

### Changed Files

| File | Changes |
|------|---------|
| `src/optimize/rotator.zig` | Add `refineSingleH`, `refineGroupRotator` |
| `src/optimize/optimizer.zig` | Add `fineSearchMover`, update `optimizeGreedy` to iterative |

### Test Plan

1. Fine search finds better score than coarse best for a single-H rotator
   positioned near an H-bond partner
2. Iterative greedy produces better result than independent greedy for a
   coupled 2-mover clique
3. Existing optimizer tests pass unchanged
4. Fine search does not regress quality (score never worse than coarse)

## Out of Scope

- Score caching
- Vertex-cut decomposition
- Delete-mask / conditional atom suppression
- Flip-specific API extensions (Issue #17)
