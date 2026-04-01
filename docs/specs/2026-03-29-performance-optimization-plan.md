# Performance Optimization Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Achieve 10-15x speedup on the optimizer and writer, bringing 2339-residue structures from 3.1s to 0.2-0.3s.

**Architecture:** Five sequential phases — algorithm improvement (mover centroid early-exit), multithreading (thread pool for independent movers), SIMD (Vec3 @Vector backend + fast exp), SoA scoring arrays (cache-friendly layout), and writer optimization (fixed-point float formatting). Each phase is a separate PR with benchmarks.

**Tech Stack:** Zig 0.15, std.Thread, @Vector SIMD intrinsics, ARM64 NEON / x86 SSE

---

## File Map

| File | Role | Phases |
|------|------|--------|
| `src/optimize/optimizer.zig` | Scoring loop, ScoreContext, optimize() | 1, 2, 4 |
| `src/math.zig` | Vec3 type, rotation utilities | 3 |
| `src/writer/mmcif_writer.zig` | mmCIF output, atom_site formatting | 5 |

---

## Task 1: Mover Centroid Early-Exit

**Files:**
- Modify: `src/optimize/optimizer.zig:33-43` (ScoreContext struct)
- Modify: `src/optimize/optimizer.zig:336-343` (mover-vs-mover loop in scoreMover)
- Modify: `src/optimize/optimizer.zig:373-435` (buildScoreContext)

- [ ] **Step 1: Write failing test for centroid computation**

Add test in `src/optimize/optimizer.zig` after the existing tests:

```zig
test "buildScoreContext computes mover centroids from parent atom positions" {
    const allocator = testing.allocator;

    var model = Model.init(allocator);
    defer model.deinit();
    // atom 0: static at origin
    try model.atoms.append(allocator, .{ .pos = .{ .x = 0, .y = 0, .z = 0 } });
    // atom 1: mover H atom (parent is atom 0 conceptually, but centroid uses atom's own pos)
    try model.atoms.append(allocator, .{ .pos = .{ .x = 1, .y = 0, .z = 0 } });
    // atom 2: another mover H atom
    try model.atoms.append(allocator, .{ .pos = .{ .x = 10, .y = 0, .z = 0 } });

    var m0 = try makeTestMover(allocator, 1, &.{
        .{ .x = 1, .y = 0, .z = 0 },
    }, &.{0});
    defer m0.deinit();
    var m1 = try makeTestMover(allocator, 2, &.{
        .{ .x = 10, .y = 0, .z = 0 },
    }, &.{0});
    defer m1.deinit();

    const movers = [_]Mover{ m0, m1 };
    var score_ctx = try buildScoreContext(allocator, &movers, model.atoms.items);
    defer score_ctx.deinit(allocator);

    // Centroids should be at initial atom positions
    try testing.expectApproxEqAbs(@as(f32, 1.0), score_ctx.mover_centroids[0].x, 0.01);
    try testing.expectApproxEqAbs(@as(f32, 10.0), score_ctx.mover_centroids[1].x, 0.01);
}
```

- [ ] **Step 2: Run test to verify it fails**

Run: `zig build test --summary all 2>&1 | grep -A2 "mover centroids"`
Expected: compilation error — `ScoreContext` has no field `mover_centroids`

- [ ] **Step 3: Add mover_centroids to ScoreContext**

In `src/optimize/optimizer.zig`, update ScoreContext:

```zig
const ScoreContext = struct {
    cell_list: ?CellList,
    static_positions: []math_mod.Vec3(f32),
    static_atom_indices: []u32,
    mover_centroids: []math_mod.Vec3(f32),

    fn deinit(self: *ScoreContext, allocator: Allocator) void {
        if (self.cell_list) |*cl| cl.deinit();
        allocator.free(self.static_positions);
        allocator.free(self.static_atom_indices);
        allocator.free(self.mover_centroids);
    }
};
```

- [ ] **Step 4: Compute centroids in buildScoreContext**

At the end of `buildScoreContext`, before the `return` statement, add centroid computation. The centroid for each mover is the mean position of its controlled atoms at construction time:

```zig
    const mover_centroids = try allocator.alloc(math_mod.Vec3(f32), movers.len);
    errdefer allocator.free(mover_centroids);
    for (movers, 0..) |m, mi| {
        var cx: f32 = 0;
        var cy: f32 = 0;
        var cz: f32 = 0;
        for (m.atom_indices) |ai| {
            cx += atoms[ai].pos.x;
            cy += atoms[ai].pos.y;
            cz += atoms[ai].pos.z;
        }
        const n: f32 = @floatFromInt(m.atom_indices.len);
        mover_centroids[mi] = .{ .x = cx / n, .y = cy / n, .z = cz / n };
    }

    return .{
        .cell_list = cell_list,
        .static_positions = static_positions,
        .static_atom_indices = static_atom_indices,
        .mover_centroids = mover_centroids,
    };
```

- [ ] **Step 5: Run test to verify it passes**

Run: `zig build test --summary all`
Expected: all tests pass (including new centroid test)

- [ ] **Step 6: Write failing test for early-exit correctness**

The early-exit must not change optimization results. Add a test that verifies scoring produces the same result with and without centroid filtering — we do this by testing that the existing singleton test still picks the correct orientation:

```zig
test "mover centroid early-exit skips distant movers in scoring" {
    const allocator = testing.allocator;

    var model = Model.init(allocator);
    defer model.deinit();
    // Atom 0: static obstacle at origin
    try model.atoms.append(allocator, .{ .pos = .{ .x = 0, .y = 0, .z = 0 }, .vdw_radius = 1.7 });
    // Atom 1: mover 0's H (near obstacle)
    try model.atoms.append(allocator, .{ .pos = .{ .x = 5, .y = 0, .z = 0 }, .vdw_radius = 1.2 });
    // Atom 2: mover 1's H (far away at x=100, should be skipped)
    try model.atoms.append(allocator, .{ .pos = .{ .x = 100, .y = 0, .z = 0 }, .vdw_radius = 1.2 });

    var m0 = try makeTestMover(allocator, 1, &.{
        .{ .x = 1.5, .y = 0, .z = 0 }, // bumps obstacle
        .{ .x = 10.0, .y = 0, .z = 0 }, // far from obstacle
    }, &.{ 0, 0 });
    defer m0.deinit();
    var m1 = try makeTestMover(allocator, 2, &.{
        .{ .x = 100, .y = 0, .z = 0 },
    }, &.{0});
    defer m1.deinit();

    var movers = [_]Mover{ m0, m1 };
    var score_ctx = try buildScoreContext(allocator, &movers, model.atoms.items);
    defer score_ctx.deinit(allocator);

    var scratch = std.ArrayListUnmanaged(u32).empty;
    defer scratch.deinit(allocator);

    // m1 is 100A away; centroid early-exit should skip it when scoring m0
    optimizeSingleton(&movers, 0, &model, .{}, &score_ctx, allocator, &scratch);
    try testing.expectEqual(@as(u16, 1), movers[0].best_orientation);
}
```

- [ ] **Step 7: Run test — should pass (early-exit not yet implemented, but result is same)**

Run: `zig build test --summary all`
Expected: PASS — the test validates correctness, not performance

- [ ] **Step 8: Implement centroid early-exit in scoreMover**

In `scoreMover`, replace the mover-vs-mover loop (lines 336-343) with a centroid-gated version:

```zig
        // Mover-controlled atoms: skip distant movers via centroid early-exit.
        // H atoms are at most ~1.5A from their initial position, so centroid + margin
        // is a conservative bounding sphere for any mover orientation.
        const centroid_cutoff = search_radius + 2.0; // 5.0 + 2.0 margin for H displacement
        const centroid_cutoff2 = centroid_cutoff * centroid_cutoff;
        const my_centroid = score_ctx.mover_centroids[mover_idx];

        for (movers, 0..) |other_m, other_idx| {
            if (other_idx == mover_idx) continue;
            // Early-exit: skip movers whose centroid is too far from ours
            const cdiff = my_centroid.sub(score_ctx.mover_centroids[other_idx]);
            if (cdiff.dot(cdiff) > centroid_cutoff2) continue;
            for (other_m.atom_indices) |oi| {
                if (oi == ai) continue;
                if (isMoverAtom(m, oi)) continue;
                total += scorePair(a, atoms[oi], config);
            }
        }
```

- [ ] **Step 9: Run all tests**

Run: `zig build test --summary all`
Expected: all tests pass

- [ ] **Step 10: Benchmark**

Run:
```bash
zig build -Doptimize=ReleaseFast
time ./zig-out/bin/zreduce examples/data/AF-P0A9J6-F1-model_v6.cif -o /dev/null
time ./zig-out/bin/zreduce examples/data/AF-P22523-F1-model_v6.cif -o /dev/null
time ./zig-out/bin/zreduce examples/data/AF-P76347-F1-model_v6.cif -o /dev/null
```
Expected: AF-P76347 drops from 3.1s to ~1.0-1.5s

- [ ] **Step 11: Compare output correctness**

```bash
./zig-out/bin/zreduce examples/data/AF-P76347-F1-model_v6.cif -o /tmp/phase1.cif
diff <(wc -l < /tmp/phase1.cif) <(wc -l < examples/result/AF-P76347-F1-zreduce.cif)
```
Expected: identical or ±1 line (minor orientation tie-breaking differences are acceptable)

- [ ] **Step 12: Commit**

```bash
git checkout -b perf/phase1-centroid-early-exit
git add src/optimize/optimizer.zig
git commit -m "perf: add mover centroid early-exit to skip distant mover-vs-mover scoring"
```

---

## Task 2: Multithreaded Singleton and Fine Search

**Files:**
- Modify: `src/optimize/optimizer.zig:46-113` (optimize function)

- [ ] **Step 1: Write test for thread-safe singleton optimization**

```zig
test "parallel singleton optimization produces valid results" {
    // Same setup as "optimize singleton picks best orientation" but verifies
    // that the thread-pool path works correctly.
    const allocator = testing.allocator;

    var model = Model.init(allocator);
    defer model.deinit();
    try model.atoms.append(allocator, .{ .pos = .{ .x = 0, .y = 0, .z = 0 }, .vdw_radius = 1.7 });
    try model.atoms.append(allocator, .{ .pos = .{ .x = 5, .y = 0, .z = 0 }, .vdw_radius = 1.7 });
    try model.atoms.append(allocator, .{ .pos = .{ .x = 0, .y = 5, .z = 0 }, .vdw_radius = 1.7 });

    var m0 = try makeTestMover(allocator, 1, &.{
        .{ .x = 1.5, .y = 0, .z = 0 },
        .{ .x = 10.0, .y = 0, .z = 0 },
    }, &.{ 0, 0 });
    defer m0.deinit();
    var m1 = try makeTestMover(allocator, 2, &.{
        .{ .x = 0, .y = 1.5, .z = 0 },
        .{ .x = 0, .y = 10.0, .z = 0 },
    }, &.{ 0, 0 });
    defer m1.deinit();

    var movers = [_]Mover{ m0, m1 };
    const opt_result = try optimize(allocator, &movers, &model, .{});

    // Both should pick orientation 1 (away from obstacle)
    try testing.expectEqual(@as(u16, 1), movers[0].best_orientation);
    try testing.expectEqual(@as(u16, 1), movers[1].best_orientation);
    try testing.expectEqual(@as(u32, 2), opt_result.n_singletons);
}
```

- [ ] **Step 2: Run test to verify it passes with current sequential code**

Run: `zig build test --summary all`
Expected: PASS (the test uses the public `optimize()` API)

- [ ] **Step 3: Restructure optimize() for parallel dispatch**

Restructure `optimize()` to:
1. Process brute-force/greedy cliques sequentially first
2. Collect singleton clique indices into a list
3. Dispatch singletons to thread pool
4. Apply best orientations
5. Dispatch fine search to thread pool

```zig
pub fn optimize(
    allocator: Allocator,
    movers: []Mover,
    model: *Model,
    config: OptConfig,
) !OptResult {
    var result = OptResult{};

    var score_ctx = try buildScoreContext(allocator, movers, model.atoms.items);
    defer score_ctx.deinit(allocator);

    // Build interaction graph
    var graph = try clique_mod.buildInteractionGraph(
        allocator, movers, model.atoms.items, config.interaction_cutoff,
    );
    defer graph.deinit();

    const cliques = try clique_mod.findCliques(allocator, &graph);
    defer {
        for (cliques) |c| allocator.free(c);
        allocator.free(cliques);
    }
    result.total_cliques = @intCast(cliques.len);

    // Collect singleton indices; process multi-mover cliques sequentially
    var singleton_indices = std.ArrayListUnmanaged(u32).empty;
    defer singleton_indices.deinit(allocator);

    var scratch = std.ArrayListUnmanaged(u32).empty;
    defer scratch.deinit(allocator);

    for (cliques) |clq| {
        if (clq.len == 1) {
            try singleton_indices.append(allocator, clq[0]);
            result.n_singletons += 1;
        } else if (totalStates(movers, clq) <= config.brute_force_limit) {
            optimizeBruteForce(allocator, movers, clq, model, config, &score_ctx, &scratch) catch |err| switch (err) {
                error.OutOfMemory => {
                    std.log.warn("brute-force OOM for clique of {d} movers; falling back to greedy", .{clq.len});
                    for (clq) |mi| {
                        try singleton_indices.append(allocator, mi);
                        result.n_singletons += 1;
                    }
                    result.n_vertex_cut += 1;
                    continue;
                },
            };
            result.n_brute_force += 1;
        } else {
            optimizeIterativeGreedy(movers, clq, model, config, &score_ctx, allocator, &scratch);
            result.n_vertex_cut += 1;
        }
    }

    // Parallel singleton optimization
    const n_threads = @min(std.Thread.getCpuCount() catch 1, 8);
    if (n_threads > 1 and singleton_indices.items.len > 1) {
        var pool: std.Thread.Pool = undefined;
        try pool.init(.{ .allocator = allocator, .n_jobs = n_threads });
        defer pool.deinit();

        var wg = std.Thread.WaitGroup{};
        for (singleton_indices.items) |mi| {
            pool.spawnWg(&wg, parallelSingleton, .{ movers, mi, model, config, &score_ctx, allocator });
        }
        wg.wait();
    } else {
        for (singleton_indices.items) |mi| {
            optimizeSingleton(movers, mi, model, config, &score_ctx, allocator, &scratch);
        }
    }

    // Apply best coarse orientations
    for (movers) |*m| {
        m.applyOrientation(model.atoms.items, m.best_orientation);
        m.current_orientation = m.best_orientation;
    }

    // Parallel fine search
    if (n_threads > 1 and movers.len > 1) {
        var pool: std.Thread.Pool = undefined;
        try pool.init(.{ .allocator = allocator, .n_jobs = n_threads });
        defer pool.deinit();

        var wg = std.Thread.WaitGroup{};
        for (0..movers.len) |mi| {
            pool.spawnWg(&wg, parallelFineSearch, .{ allocator, movers, @as(u32, @intCast(mi)), model, config, &score_ctx });
        }
        wg.wait();
    } else {
        for (0..movers.len) |mi| {
            fineSearchMover(allocator, movers, @intCast(mi), model, config, &score_ctx, &scratch);
        }
    }

    return result;
}
```

- [ ] **Step 4: Implement parallelSingleton and parallelFineSearch wrappers**

These wrappers allocate thread-local scratch buffers:

```zig
fn parallelSingleton(
    movers: []Mover,
    mover_idx: u32,
    model: *Model,
    config: OptConfig,
    score_ctx: *const ScoreContext,
    allocator: Allocator,
) void {
    var scratch = std.ArrayListUnmanaged(u32).empty;
    defer scratch.deinit(allocator);
    optimizeSingleton(movers, mover_idx, model, config, score_ctx, allocator, &scratch);
}

fn parallelFineSearch(
    allocator: Allocator,
    movers: []Mover,
    mover_idx: u32,
    model: *Model,
    config: OptConfig,
    score_ctx: *const ScoreContext,
) void {
    var scratch = std.ArrayListUnmanaged(u32).empty;
    defer scratch.deinit(allocator);
    fineSearchMover(allocator, movers, mover_idx, model, config, score_ctx, &scratch);
}
```

- [ ] **Step 5: Run all tests**

Run: `zig build test --summary all`
Expected: all tests pass

- [ ] **Step 6: Benchmark**

```bash
zig build -Doptimize=ReleaseFast
time ./zig-out/bin/zreduce examples/data/AF-P76347-F1-model_v6.cif -o /dev/null
```
Expected: drops from ~1.3s to ~0.4-0.6s

- [ ] **Step 7: Compare output correctness**

```bash
./zig-out/bin/zreduce examples/data/AF-P76347-F1-model_v6.cif -o /tmp/phase2.cif
wc -l /tmp/phase2.cif
```
Expected: same line count as phase 1 output

- [ ] **Step 8: Commit**

```bash
git add src/optimize/optimizer.zig
git commit -m "perf: parallelize singleton and fine-search optimization with thread pool"
```

---

## Task 3: SIMD Vec3 Backend

**Files:**
- Modify: `src/math.zig:8-74` (Vec3 type)

- [ ] **Step 1: Run existing math tests as baseline**

Run: `zig build test --summary all 2>&1 | grep -c "passed"`
Expected: 233 (or current count) — all pass. These serve as regression tests.

- [ ] **Step 2: Replace Vec3 internals with @Vector(4, f32)**

Replace the Vec3 struct body. The public API stays the same (`x`, `y`, `z` fields, all methods). Use `@Vector` internally:

```zig
pub fn Vec3(comptime T: type) type {
    return struct {
        const Self = @This();
        v: @Vector(4, T),

        pub const zero = Self{ .v = @splat(0) };

        pub fn init(x: T, y: T, z: T) Self {
            return .{ .v = .{ x, y, z, 0 } };
        }

        // Field accessors for compatibility
        pub inline fn getX(self: Self) T { return self.v[0]; }
        pub inline fn getY(self: Self) T { return self.v[1]; }
        pub inline fn getZ(self: Self) T { return self.v[2]; }

        pub fn add(self: Self, other: Self) Self {
            return .{ .v = self.v + other.v };
        }

        pub fn sub(self: Self, other: Self) Self {
            return .{ .v = self.v - other.v };
        }

        pub fn scale(self: Self, s: T) Self {
            return .{ .v = self.v * @as(@Vector(4, T), @splat(s)) };
        }

        pub fn scaleTo(self: Self, len: T) Self {
            const l = self.length();
            if (l < 1e-10) return Self.zero;
            return self.scale(len / l);
        }

        pub fn dot(self: Self, other: Self) T {
            const prod = self.v * other.v;
            return prod[0] + prod[1] + prod[2]; // ignore w lane
        }

        pub fn cross(self: Self, other: Self) Self {
            return .{ .v = .{
                self.v[1] * other.v[2] - self.v[2] * other.v[1],
                self.v[2] * other.v[0] - self.v[0] * other.v[2],
                self.v[0] * other.v[1] - self.v[1] * other.v[0],
                0,
            } };
        }

        pub fn length(self: Self) T {
            return @sqrt(self.dot(self));
        }

        pub fn distance(self: Self, other: Self) T {
            return self.sub(other).length();
        }

        pub fn normalize(self: Self) Self {
            const l = self.length();
            if (l < 1e-10) return Self.zero;
            return self.scale(1.0 / l);
        }

        pub fn negate(self: Self) Self {
            return .{ .v = -self.v };
        }

        pub fn cast(self: Self, comptime U: type) Vec3(U) {
            return .{ .v = .{
                @floatCast(self.v[0]),
                @floatCast(self.v[1]),
                @floatCast(self.v[2]),
                0,
            } };
        }
    };
}
```

**IMPORTANT:** This changes the struct layout from `{x, y, z}` to `{v: @Vector(4, T)}`. All call sites that use `.x`, `.y`, `.z` field access or struct literal `{.x=, .y=, .z=}` must be updated. This is a large ripple.

**Alternative (lower risk):** Keep the `x, y, z` field API, add SIMD only to `dot` and `sub` via inline helper:

```zig
pub fn dot(self: Self, other: Self) T {
    const a: @Vector(4, T) = .{ self.x, self.y, self.z, 0 };
    const b: @Vector(4, T) = .{ other.x, other.y, other.z, 0 };
    const prod = a * b;
    return prod[0] + prod[1] + prod[2];
}

pub fn sub(self: Self, other: Self) Self {
    const a: @Vector(4, T) = .{ self.x, self.y, self.z, 0 };
    const b: @Vector(4, T) = .{ other.x, other.y, other.z, 0 };
    const r = a - b;
    return .{ .x = r[0], .y = r[1], .z = r[2] };
}
```

**Use the lower-risk alternative** — it keeps all existing `.x`, `.y`, `.z` usage working and only changes the hot path methods.

- [ ] **Step 3: Run all tests**

Run: `zig build test --summary all`
Expected: all tests pass (same behavior, SIMD under the hood)

- [ ] **Step 4: Implement fast_exp approximation**

Add to `src/math.zig`:

```zig
/// Fast approximation of exp(x) for x in [-10, 0] range.
/// Uses Schraudolph's method with bias correction.
/// Max relative error ~1% which is sufficient for contact scoring.
pub fn fastExp(x: f32) f32 {
    // Clamp to avoid overflow in integer conversion
    const clamped = @max(x, -87.0);
    // Schraudolph's approximation: reinterpret float bits
    const v: i32 = @intFromFloat(12102203.0 * clamped + 1065353216.0);
    return @bitCast(@as(u32, @intCast(@max(v, 0))));
}
```

- [ ] **Step 5: Write test for fast_exp accuracy**

```zig
test "fastExp approximation within 2% for scoring range" {
    // Test the range used in contact scoring: exp(-r*r) where r = gap/0.25, gap in [0, 0.5]
    // So x ranges from 0 to -(0.5/0.25)^2 = -4
    const test_values = [_]f32{ 0.0, -0.5, -1.0, -2.0, -4.0 };
    for (test_values) |x| {
        const exact = @exp(x);
        const approx = fastExp(x);
        const rel_err = @abs(approx - exact) / exact;
        try std.testing.expect(rel_err < 0.02);
    }
}
```

- [ ] **Step 6: Run tests**

Run: `zig build test --summary all`
Expected: all pass

- [ ] **Step 7: Use fast_exp in scorePair**

In `src/optimize/optimizer.zig`, in `scorePair`, replace:
```zig
return @exp(-ratio * ratio);
```
with:
```zig
return math_mod.fastExp(-ratio * ratio);
```

- [ ] **Step 8: Run all tests + benchmark**

Run: `zig build test --summary all && zig build -Doptimize=ReleaseFast`
Then: `time ./zig-out/bin/zreduce examples/data/AF-P76347-F1-model_v6.cif -o /dev/null`
Expected: all tests pass; benchmark ~0.3s

- [ ] **Step 9: Commit**

```bash
git add src/math.zig src/optimize/optimizer.zig
git commit -m "perf: add SIMD Vec3 dot/sub and fast exp approximation for scoring"
```

---

## Task 4: SoA Scoring Arrays

**Files:**
- Modify: `src/optimize/optimizer.zig:33-43` (ScoreContext)
- Modify: `src/optimize/optimizer.zig:300-345` (scoreMover)
- Modify: `src/optimize/optimizer.zig:348-371` (scorePair -> scorePairSoA)
- Modify: `src/optimize/optimizer.zig:373-435` (buildScoreContext)

- [ ] **Step 1: Add SoA fields to ScoreContext**

```zig
const ScoreContext = struct {
    cell_list: ?CellList,
    static_positions: []math_mod.Vec3(f32),
    static_atom_indices: []u32,
    static_radii: []f32,
    static_flags: []element.AtomFlags,
    mover_centroids: []math_mod.Vec3(f32),

    fn deinit(self: *ScoreContext, allocator: Allocator) void {
        if (self.cell_list) |*cl| cl.deinit();
        allocator.free(self.static_positions);
        allocator.free(self.static_atom_indices);
        allocator.free(self.static_radii);
        allocator.free(self.static_flags);
        allocator.free(self.mover_centroids);
    }
};
```

- [ ] **Step 2: Extract radii and flags in buildScoreContext**

After building `static_positions` and `static_atom_indices`, add:

```zig
    const static_radii = try allocator.alloc(f32, static_count);
    errdefer allocator.free(static_radii);
    const static_flags = try allocator.alloc(element.AtomFlags, static_count);
    errdefer allocator.free(static_flags);

    // Reset out_i and fill all SoA arrays in a single pass
    out_i = 0;
    for (atoms, 0..) |a, i| {
        if (moved_atoms[i]) continue;
        static_positions[out_i] = a.pos;
        static_atom_indices[out_i] = @intCast(i);
        static_radii[out_i] = a.vdw_radius;
        static_flags[out_i] = a.flags;
        out_i += 1;
    }
```

Note: merge the existing position/index fill loop with this one so there's only one pass.

- [ ] **Step 3: Add scorePairSoA function**

```zig
fn scorePairSoA(
    pos_a: math_mod.Vec3(f32),
    radius_a: f32,
    flags_a: element.AtomFlags,
    pos_b: math_mod.Vec3(f32),
    radius_b: f32,
    flags_b: element.AtomFlags,
    config: OptConfig,
) f32 {
    const diff = pos_a.sub(pos_b);
    const dist2 = diff.dot(diff);
    const sum_r = radius_a + radius_b;
    const sum_r2 = sum_r * sum_r;

    if (dist2 < sum_r2) {
        const dist = @sqrt(dist2);
        const gap = dist - sum_r;
        if (scorer_mod.isHBond(flags_a, flags_b, gap, config.scoring_params)) {
            return config.scoring_params.hb_weight * (-0.5 * gap);
        }
        return -config.scoring_params.bump_weight * (-0.5 * gap);
    }

    const threshold = sum_r + 0.5;
    if (dist2 < threshold * threshold) {
        const dist = @sqrt(dist2);
        const gap = dist - sum_r;
        const ratio = gap / config.scoring_params.gap_scale;
        return math_mod.fastExp(-ratio * ratio);
    }
    return 0.0;
}
```

- [ ] **Step 4: Update scoreMover to use SoA for static atoms**

In the static atoms scoring path, replace `scorePair(a, atoms[oi], config)` with:

```zig
            for (scratch.items) |static_idx| {
                total += scorePairSoA(
                    a.pos, a.vdw_radius, a.flags,
                    score_ctx.static_positions[static_idx],
                    score_ctx.static_radii[static_idx],
                    score_ctx.static_flags[static_idx],
                    config,
                );
            }
```

And in the pairwise fallback:

```zig
            for (score_ctx.static_atom_indices, 0..) |_, si| {
                total += scorePairSoA(
                    a.pos, a.vdw_radius, a.flags,
                    score_ctx.static_positions[si],
                    score_ctx.static_radii[si],
                    score_ctx.static_flags[si],
                    config,
                );
            }
```

Keep `scorePair` for the mover-vs-mover path (those still need `atoms[oi]`).

- [ ] **Step 5: Run all tests + benchmark**

Run: `zig build test --summary all && zig build -Doptimize=ReleaseFast`
Then: `time ./zig-out/bin/zreduce examples/data/AF-P76347-F1-model_v6.cif -o /dev/null`
Expected: all pass; benchmark ~0.25s

- [ ] **Step 6: Commit**

```bash
git add src/optimize/optimizer.zig
git commit -m "perf: use SoA arrays for static atom scoring to improve cache efficiency"
```

---

## Task 5: Writer Fixed-Point Float Formatting

**Files:**
- Modify: `src/writer/mmcif_writer.zig`

- [ ] **Step 1: Write test for fixed-point float formatting**

Add test in `src/writer/mmcif_writer.zig`:

```zig
test "writeFixedFloat3 formats coordinates correctly" {
    var buf: [32]u8 = undefined;
    var fbs = std.io.fixedBufferStream(&buf);

    try writeFixedFloat3(fbs.writer(), 12.345);
    try testing.expectEqualStrings("12.345", fbs.getWritten());

    fbs.reset();
    try writeFixedFloat3(fbs.writer(), -0.001);
    try testing.expectEqualStrings("-0.001", fbs.getWritten());

    fbs.reset();
    try writeFixedFloat3(fbs.writer(), 0.0);
    try testing.expectEqualStrings("0.000", fbs.getWritten());

    fbs.reset();
    try writeFixedFloat3(fbs.writer(), 123.4567);
    try testing.expectEqualStrings("123.457", fbs.getWritten());
}

test "writeFixedFloat2 formats b-factor correctly" {
    var buf: [32]u8 = undefined;
    var fbs = std.io.fixedBufferStream(&buf);

    try writeFixedFloat2(fbs.writer(), 45.67);
    try testing.expectEqualStrings("45.67", fbs.getWritten());

    fbs.reset();
    try writeFixedFloat2(fbs.writer(), 1.0);
    try testing.expectEqualStrings("1.00", fbs.getWritten());
}
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `zig build test --summary all 2>&1 | grep "writeFixedFloat"`
Expected: compilation error — function not defined

- [ ] **Step 3: Implement writeFixedFloat3 and writeFixedFloat2**

```zig
/// Write f32 with exactly 3 decimal places using integer arithmetic.
/// Much faster than std.fmt for fixed-precision output.
fn writeFixedFloat3(writer: anytype, val: f32) !void {
    if (val < 0) {
        try writer.writeByte('-');
        return writeFixedFloat3(writer, -val);
    }
    const scaled: u64 = @intFromFloat(@round(val * 1000.0));
    const int_part = scaled / 1000;
    const frac_part = scaled % 1000;
    try writer.print("{d}", .{int_part});
    try writer.writeByte('.');
    if (frac_part < 10) {
        try writer.writeAll("00");
    } else if (frac_part < 100) {
        try writer.writeByte('0');
    }
    try writer.print("{d}", .{frac_part});
}

/// Write f32 with exactly 2 decimal places using integer arithmetic.
fn writeFixedFloat2(writer: anytype, val: f32) !void {
    if (val < 0) {
        try writer.writeByte('-');
        return writeFixedFloat2(writer, -val);
    }
    const scaled: u64 = @intFromFloat(@round(val * 100.0));
    const int_part = scaled / 100;
    const frac_part = scaled % 100;
    try writer.print("{d}", .{int_part});
    try writer.writeByte('.');
    if (frac_part < 10) {
        try writer.writeByte('0');
    }
    try writer.print("{d}", .{frac_part});
}
```

- [ ] **Step 4: Run tests to verify they pass**

Run: `zig build test --summary all`
Expected: all pass

- [ ] **Step 5: Replace std.fmt calls in writeAtomSitePreserving**

In `src/writer/mmcif_writer.zig`, replace the coordinate/b-factor formatting lines:

Replace:
```zig
try writer.print("{d:.3}", .{atom.pos.x});
```
With:
```zig
try writeFixedFloat3(writer, atom.pos.x);
```

Same for `pos.y`, `pos.z`. And for occupancy/b_factor:
```zig
try writeFixedFloat2(writer, atom.occupancy);
try writeFixedFloat2(writer, atom.b_factor);
```

Also update the fallback `writeAtomSite` function if it has similar formatting.

- [ ] **Step 6: Run all tests + benchmark**

Run: `zig build test --summary all && zig build -Doptimize=ReleaseFast`
Then: `time ./zig-out/bin/zreduce examples/data/AF-P76347-F1-model_v6.cif -o /dev/null`
Expected: all pass; benchmark ~0.2s

- [ ] **Step 7: Verify output matches**

```bash
./zig-out/bin/zreduce examples/data/AF-P76347-F1-model_v6.cif -o /tmp/phase5.cif
diff /tmp/phase5.cif /tmp/phase4.cif
```
Expected: identical (same precision)

- [ ] **Step 8: Commit**

```bash
git add src/writer/mmcif_writer.zig
git commit -m "perf: use fixed-point integer arithmetic for mmCIF float formatting"
```

---

## Final: Create PR with All Phases

- [ ] **Step 1: Push branch and create PR**

```bash
git push -u origin perf/phase1-centroid-early-exit
gh pr create --title "perf: optimize scoring, threading, SIMD, and writer" --body "..."
```

Include benchmark table in PR body showing before/after for all test structures.
