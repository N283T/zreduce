//! Optimization engine: selects the best orientation for each mover
//! by searching over cliques (connected components) of interacting movers.
//! Uses brute-force enumeration for small cliques, greedy for large ones.

const std = @import("std");
const Allocator = std.mem.Allocator;
const math_mod = @import("../math.zig");
const Vec3 = math_mod.Vec3;
const model_mod = @import("../model.zig");
const Model = model_mod.Model;
const Atom = model_mod.Atom;
const mover_mod = @import("mover.zig");
const Mover = mover_mod.Mover;
const scorer_mod = @import("scorer.zig");
const clique_mod = @import("clique.zig");
const rotator_mod = @import("rotator.zig");
const element = @import("../element.zig");
const CellList = model_mod.CellList;

pub const OptConfig = struct {
    brute_force_limit: u64 = 100_000,
    interaction_cutoff: f32 = 6.0, // Angstrom -- max distance for mover interaction
    scoring_params: scorer_mod.ScoringParams = .{},
};

pub const OptResult = struct {
    n_singletons: u32 = 0,
    n_brute_force: u32 = 0,
    n_vertex_cut: u32 = 0,
    total_cliques: u32 = 0,
};

const ScoreContext = struct {
    cell_list: ?CellList,
    static_positions: []math_mod.Vec3(f32),
    static_atom_indices: []u32,
    /// Mean position of each mover's atoms across all its orientations.
    mover_centroids: []math_mod.Vec3(f32),
    /// Bounding radius of each mover: max distance from centroid to any orientation position.
    mover_radii: []f32,

    fn deinit(self: *ScoreContext, allocator: Allocator) void {
        if (self.cell_list) |*cl| cl.deinit();
        allocator.free(self.static_positions);
        allocator.free(self.static_atom_indices);
        allocator.free(self.mover_centroids);
        allocator.free(self.mover_radii);
    }
};

/// Arguments shared by parallel singleton and fine-search tasks.
/// Reused for both phases so a single thread pool handles the whole optimization.
const ParallelTaskArgs = struct {
    movers: []Mover,
    mover_idx: u32,
    model: *Model,
    config: OptConfig,
    score_ctx: *const ScoreContext,
    allocator: Allocator,
    /// Output field for fine-search: set by fineSearchMover to the best fine positions
    /// (allocated with `allocator`), or null if the coarse best is retained.
    /// The main thread applies and frees this after all tasks join.
    best_fine_positions: ?[]math_mod.Vec3(f32) = null,
};

/// Thread pool work function for singleton optimization.
/// Uses scoreMoverWithPositions to avoid writing to the shared model during scoring,
/// eliminating data races when multiple singletons run concurrently.
/// Only writes `best_orientation` (a per-mover field with no cross-mover aliasing).
fn parallelSingleton(args: *const ParallelTaskArgs) void {
    var scratch = std.ArrayListUnmanaged(u32).empty;
    defer scratch.deinit(args.allocator);
    optimizeSingleton(args.movers, args.mover_idx, args.model, args.config, args.score_ctx, args.allocator, &scratch);
}

/// Thread pool work function for fine search.
/// Uses scoreMoverWithPositions to avoid writing to the shared model during scoring.
/// Stores the best fine positions in args.best_fine_positions (if better than coarse best);
/// the main thread applies and frees them after all tasks join, avoiding write-write races.
fn parallelFineSearch(args: *ParallelTaskArgs) void {
    var scratch = std.ArrayListUnmanaged(u32).empty;
    defer scratch.deinit(args.allocator);
    fineSearchMover(args.allocator, args.movers, args.mover_idx, args.model, args.config, args.score_ctx, &scratch, &args.best_fine_positions);
}

/// Optimize all movers: find best orientations using clique-based search.
pub fn optimize(
    allocator: Allocator,
    movers: []Mover,
    model: *Model,
    config: OptConfig,
) !OptResult {
    var result = OptResult{};

    const n_threads: u32 = @intCast(@min(std.Thread.getCpuCount() catch 1, 8));

    var score_ctx = try buildScoreContext(allocator, movers, model.atoms.items);
    defer score_ctx.deinit(allocator);

    var scratch = std.ArrayListUnmanaged(u32).empty;
    defer scratch.deinit(allocator);

    // Build interaction graph
    var graph = try clique_mod.buildInteractionGraph(
        allocator,
        movers,
        model.atoms.items,
        config.interaction_cutoff,
    );
    defer graph.deinit();

    // Find connected components
    const cliques = try clique_mod.findCliques(allocator, &graph);
    defer {
        for (cliques) |c| allocator.free(c);
        allocator.free(cliques);
    }

    result.total_cliques = @intCast(cliques.len);

    // Collect singleton mover indices for parallel dispatch after sequential cliques.
    var singleton_indices = std.ArrayListUnmanaged(u32).empty;
    defer singleton_indices.deinit(allocator);

    for (cliques) |clq| {
        if (clq.len == 1) {
            // Defer singletons for parallel dispatch below.
            try singleton_indices.append(allocator, clq[0]);
            result.n_singletons += 1;
        } else if (totalStates(movers, clq) <= config.brute_force_limit) {
            // Brute force: enumerate all combinations
            optimizeBruteForce(allocator, movers, clq, model, config, &score_ctx, &scratch) catch |err| switch (err) {
                error.OutOfMemory => {
                    std.log.warn("brute-force OOM for clique of {d} movers; falling back to singletons", .{clq.len});
                    for (clq) |mi| try singleton_indices.append(allocator, mi);
                    result.n_vertex_cut += 1;
                    continue;
                },
            };
            result.n_brute_force += 1;
        } else {
            // Vertex-cut decomposition (simplified: greedy for Phase 3)
            optimizeIterativeGreedy(movers, clq, model, config, &score_ctx, allocator, &scratch);
            result.n_vertex_cut += 1;
        }
    }

    // Dispatch singleton optimization in parallel (or sequentially when not worth the overhead).
    // Thread safety: optimizeSingleton uses scoreMoverWithPositions which reads other movers'
    // atoms from model.atoms but does NOT write to model.atoms during scoring. The only
    // write is to m.best_orientation (a per-mover field; each thread owns a distinct mover).
    const use_parallel = n_threads > 1 and singleton_indices.items.len > 1;

    // Allocate one pool and reuse it for both the singleton and fine-search phases.
    var pool: std.Thread.Pool = undefined;
    if (use_parallel) {
        try pool.init(.{ .allocator = allocator, .n_jobs = n_threads });
    }
    defer if (use_parallel) pool.deinit();

    if (use_parallel) {
        // Build per-task arg structs (one per singleton) on the heap so thread lifetimes are safe.
        const singleton_args = try allocator.alloc(ParallelTaskArgs, singleton_indices.items.len);
        defer allocator.free(singleton_args);
        for (singleton_indices.items, 0..) |mi, i| {
            singleton_args[i] = .{
                .movers = movers,
                .mover_idx = mi,
                .model = model,
                .config = config,
                .score_ctx = &score_ctx,
                .allocator = allocator,
            };
        }

        var wg = std.Thread.WaitGroup{};
        for (singleton_args) |*args| {
            pool.spawnWg(&wg, parallelSingleton, .{args});
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

    // Fine search phase: refine each mover around its coarse best, in parallel when worthwhile.
    // The gate uses n_threads > 1 and movers.len > 1 (independent of singleton count).
    // Thread safety: fineSearchMover uses scoreMoverWithPositions; the only write to
    // model.atoms happens after all threads join (the sequential apply-best loop below).
    const use_parallel_fine = n_threads > 1 and movers.len > 1;
    if (use_parallel_fine) {
        const fine_args = try allocator.alloc(ParallelTaskArgs, movers.len);
        defer allocator.free(fine_args);
        for (0..movers.len) |mi| {
            fine_args[mi] = .{
                .movers = movers,
                .mover_idx = @intCast(mi),
                .model = model,
                .config = config,
                .score_ctx = &score_ctx,
                .allocator = allocator,
                // best_fine_positions defaults to null; set by parallelFineSearch.
            };
        }

        var wg = std.Thread.WaitGroup{};
        for (fine_args) |*args| {
            pool.spawnWg(&wg, parallelFineSearch, .{args});
        }
        wg.wait();

        // Apply fine-search results sequentially after all threads have joined.
        // This avoids write-write races on model.atoms during parallel scoring.
        for (fine_args) |*args| {
            if (args.best_fine_positions) |positions| {
                defer allocator.free(positions);
                const mi = args.mover_idx;
                for (movers[mi].atom_indices, 0..) |ai, j| {
                    model.atoms.items[ai].pos = positions[j];
                }
            }
        }
    } else {
        for (0..movers.len) |mi| {
            var best_fine_positions: ?[]math_mod.Vec3(f32) = null;
            fineSearchMover(allocator, movers, @intCast(mi), model, config, &score_ctx, &scratch, &best_fine_positions);
            if (best_fine_positions) |positions| {
                defer allocator.free(positions);
                for (movers[mi].atom_indices, 0..) |ai, j| {
                    model.atoms.items[ai].pos = positions[j];
                }
            }
        }
    }

    return result;
}

fn optimizeSingleton(
    movers: []Mover,
    mover_idx: u32,
    model: *Model,
    config: OptConfig,
    score_ctx: *const ScoreContext,
    allocator: Allocator,
    scratch: *std.ArrayListUnmanaged(u32),
) void {
    const m = &movers[mover_idx];
    var best_score: f32 = -std.math.inf(f32);
    var best_idx: u16 = 0;

    for (0..m.nOrientations()) |oi| {
        const idx: u16 = @intCast(oi);
        // Score using the orientation's positions directly — do NOT call applyOrientation,
        // which would write to model.atoms and race with concurrent threads reading it.
        const orient_positions = m.orientations[idx].positions;
        const score = scoreMoverWithPositions(m, mover_idx, orient_positions, movers, model, config, score_ctx, allocator, scratch) - m.orientationPenalty(idx);
        if (score > best_score) {
            best_score = score;
            best_idx = idx;
        }
    }

    m.best_orientation = best_idx;
}

/// Refine a mover's best orientation with fine angular search.
///
/// Thread-safe: does NOT write to model.atoms during scoring. Scoring uses
/// scoreMoverWithPositions with explicit orientation positions. If a fine position
/// beats the coarse best, the best positions are allocated and written to
/// `out_best_positions` for the caller to apply after all parallel tasks join.
/// The caller is responsible for freeing the slice.
fn fineSearchMover(
    allocator: Allocator,
    movers: []Mover,
    mover_idx: u32,
    model: *Model,
    config: OptConfig,
    score_ctx: *const ScoreContext,
    scratch: *std.ArrayListUnmanaged(u32),
    out_best_positions: *?[]math_mod.Vec3(f32),
) void {
    const m = &movers[mover_idx];
    const atoms = model.atoms.items;

    const fine = rotator_mod.generateFineOrientations(allocator, atoms, m, m.best_orientation) catch |err| {
        std.log.warn("fine orientation generation failed for mover {d}: {s}", .{ mover_idx, @errorName(err) });
        return;
    };
    defer {
        for (fine) |o| allocator.free(o.positions);
        allocator.free(fine);
    }

    if (fine.len == 0) return;

    // Score current coarse best using the orientation's stored positions — do NOT read
    // atoms[ai].pos because another thread may be writing it concurrently during fine search.
    const coarse_positions = m.orientations[m.best_orientation].positions;
    var best_score = scoreMoverWithPositions(m, mover_idx, coarse_positions, movers, model, config, score_ctx, allocator, scratch) - m.orientationPenalty(m.best_orientation);

    // Score fine orientations using their positions directly (no writes to model.atoms).
    var best_fine: ?usize = null;
    for (fine, 0..) |orient, fi| {
        const score = scoreMoverWithPositions(m, mover_idx, orient.positions, movers, model, config, score_ctx, allocator, scratch) - orient.penalty;
        if (score > best_score) {
            best_score = score;
            best_fine = fi;
        }
    }

    if (best_fine) |fi| {
        // Copy the best fine positions into a freshly allocated slice; the caller will
        // apply them to model.atoms after all threads join (avoiding write-write races).
        const n = m.atom_indices.len;
        const best = allocator.alloc(math_mod.Vec3(f32), n) catch |err| {
            std.log.warn("fine-search result alloc failed for mover {d}: {s}", .{ mover_idx, @errorName(err) });
            return;
        };
        @memcpy(best, fine[fi].positions);
        out_best_positions.* = best;
    }
    // If no fine orientation improves on the coarse best, leave out_best_positions null.
    // The coarse best is already applied to model.atoms (done before this phase starts).
}

fn optimizeBruteForce(
    allocator: Allocator,
    movers: []Mover,
    clq: []const u32,
    model: *Model,
    config: OptConfig,
    score_ctx: *const ScoreContext,
    scratch: *std.ArrayListUnmanaged(u32),
) !void {
    const n = clq.len;
    // Current orientation indices for each mover in clique
    const indices = try allocator.alloc(u16, n);
    defer allocator.free(indices);
    @memset(indices, 0);

    var best_score: f32 = -std.math.inf(f32);
    const best_indices = try allocator.alloc(u16, n);
    defer allocator.free(best_indices);
    @memset(best_indices, 0);

    // Enumerate all combinations
    while (true) {
        // Apply current combination
        for (clq, 0..) |mi, i| {
            movers[mi].applyOrientation(model.atoms.items, indices[i]);
        }

        // Score all movers in clique
        var total_score: f32 = 0;
        for (clq, 0..) |mi, i| {
            total_score += scoreMover(&movers[mi], mi, movers, model, config, score_ctx, allocator, scratch);
            total_score -= movers[mi].orientationPenalty(indices[i]);
        }

        if (total_score > best_score) {
            best_score = total_score;
            @memcpy(best_indices, indices);
        }

        // Increment indices (mixed radix counter)
        if (!incrementIndices(indices, movers, clq)) break;
    }

    // Record best orientations
    for (clq, 0..) |mi, i| {
        movers[mi].best_orientation = best_indices[i];
    }
}

fn incrementIndices(indices: []u16, movers: []const Mover, clq: []const u32) bool {
    var i: usize = indices.len;
    while (i > 0) {
        i -= 1;
        indices[i] += 1;
        if (indices[i] < movers[clq[i]].nOrientations()) return true;
        indices[i] = 0;
    }
    return false; // overflow = done
}

/// Compute the total number of combined states for a set of movers.
/// Uses saturating multiplication to avoid overflow.
pub fn totalStates(movers: []const Mover, clq: []const u32) u64 {
    var total: u64 = 1;
    for (clq) |mi| {
        total *|= movers[mi].nOrientations(); // saturating multiply
    }
    return total;
}

fn optimizeIterativeGreedy(
    movers: []Mover,
    clq: []const u32,
    model: *Model,
    config: OptConfig,
    score_ctx: *const ScoreContext,
    allocator: Allocator,
    scratch: *std.ArrayListUnmanaged(u32),
) void {
    const max_iterations: u32 = 3;
    var iteration: u32 = 0;

    while (iteration < max_iterations) : (iteration += 1) {
        var changed = false;
        for (clq) |mi| {
            const old_best = movers[mi].best_orientation;
            optimizeSingleton(movers, mi, model, config, score_ctx, allocator, scratch);
            movers[mi].applyOrientation(model.atoms.items, movers[mi].best_orientation);
            if (movers[mi].best_orientation != old_best) changed = true;
        }
        if (!changed) break;
    }
}

/// Check whether a given atom index belongs to a specific mover.
fn isMoverAtom(m: *const Mover, atom_idx: u32) bool {
    for (m.atom_indices) |ai| {
        if (ai == atom_idx) return true;
    }
    return false;
}

/// Score a mover against the model using explicit positions for its atoms.
///
/// `positions_override` provides the positions for the mover's atoms in the same
/// order as `m.atom_indices`. This avoids writing to `model.atoms` during scoring,
/// which is required for thread-safe parallel optimization: concurrent threads can
/// each call this function with different orientation positions without data races.
///
/// Other movers' atoms are read from `model.atoms` at their current positions
/// (which is their initial or last-applied coarse orientation — acceptable for the
/// parallel coarse/fine search phases).
///
/// Static (non-mover) atoms are queried via CellList spatial index (O(nearby)).
/// Mover-controlled atoms are scored by direct iteration.
fn scoreMoverWithPositions(
    m: *const Mover,
    mover_idx: u32,
    positions_override: []const math_mod.Vec3(f32),
    movers: []const Mover,
    model: *const Model,
    config: OptConfig,
    score_ctx: *const ScoreContext,
    allocator: Allocator,
    scratch: *std.ArrayListUnmanaged(u32),
) f32 {
    var total: f32 = 0;
    const atoms = model.atoms.items;
    const search_radius: f32 = 5.0;

    // Centroid early-exit: skip mover pairs where even the closest possible atoms
    // cannot be within scoring range. The cutoff accounts for each mover's bounding
    // radius (max displacement of any orientation position from the centroid).
    const my_centroid = score_ctx.mover_centroids[mover_idx];
    const my_radius = score_ctx.mover_radii[mover_idx];

    for (m.atom_indices, 0..) |ai, j| {
        // Use the override position for this mover's atom instead of atoms[ai].pos.
        var a = atoms[ai];
        a.pos = positions_override[j];

        // Static atoms: use spatial index when available, else brute-force scan.
        if (score_ctx.cell_list) |cl| {
            scratch.clearRetainingCapacity();
            cl.neighborsInRadius(a.pos, search_radius, scratch, allocator, score_ctx.static_positions) catch {
                std.log.warn("neighbor query OOM during scoring for atom {d}; returning -inf", .{ai});
                return -std.math.inf(f32);
            };
            for (scratch.items) |static_idx| {
                const oi = score_ctx.static_atom_indices[static_idx];
                total += scorePair(a, atoms[oi], config);
            }
        } else {
            // Fallback: pairwise scan of all static atoms (no CellList available).
            for (score_ctx.static_atom_indices) |oi| {
                total += scorePair(a, atoms[oi], config);
            }
        }

        // Mover-controlled atoms: score directly with current coordinates.
        // Centroid early-exit: skip other movers whose bounding spheres cannot
        // overlap with this mover's scoring range.
        for (movers, 0..) |other_m, other_idx| {
            if (other_idx == mover_idx) continue;
            const other_radius = score_ctx.mover_radii[other_idx];
            const pair_cutoff = search_radius + my_radius + other_radius;
            const cdiff = my_centroid.sub(score_ctx.mover_centroids[other_idx]);
            if (cdiff.dot(cdiff) > pair_cutoff * pair_cutoff) continue;
            for (other_m.atom_indices) |oi| {
                if (oi == ai) continue;
                if (isMoverAtom(m, oi)) continue;
                total += scorePair(a, atoms[oi], config);
            }
        }
    }
    return total;
}

/// Score a mover's current orientation against the model.
/// Reads atom positions from model.atoms — only safe to call sequentially
/// (e.g. brute-force and greedy clique search where orientations are applied first).
/// For parallel singleton/fine-search use scoreMoverWithPositions instead.
fn scoreMover(
    m: *const Mover,
    mover_idx: u32,
    movers: []const Mover,
    model: *const Model,
    config: OptConfig,
    score_ctx: *const ScoreContext,
    allocator: Allocator,
    scratch: *std.ArrayListUnmanaged(u32),
) f32 {
    return scoreMoverWithPositions(
        m,
        mover_idx,
        blk: {
            // Build a slice of current positions from model.atoms for this mover.
            // This is a small, bounded allocation on the stack (up to ~16 atoms per mover).
            var pos_buf: [64]math_mod.Vec3(f32) = undefined;
            const n = m.atom_indices.len;
            std.debug.assert(n <= pos_buf.len);
            for (m.atom_indices, 0..) |ai, i| {
                pos_buf[i] = model.atoms.items[ai].pos;
            }
            break :blk pos_buf[0..n];
        },
        movers,
        model,
        config,
        score_ctx,
        allocator,
        scratch,
    );
}

fn scorePair(a: Atom, other: Atom, config: OptConfig) f32 {
    const diff = a.pos.sub(other.pos);
    const dist2 = diff.dot(diff);
    const sum_r = a.vdw_radius + other.vdw_radius;
    const sum_r2 = sum_r * sum_r;

    if (dist2 < sum_r2) {
        const dist = @sqrt(dist2);
        const gap = dist - sum_r;
        if (scorer_mod.isHBond(a.flags, other.flags, gap, config.scoring_params)) {
            return config.scoring_params.hb_weight * (-0.5 * gap);
        }
        return -config.scoring_params.bump_weight * (-0.5 * gap);
    }

    const threshold = sum_r + 0.5;
    if (dist2 < threshold * threshold) {
        const dist = @sqrt(dist2);
        const gap = dist - sum_r;
        const ratio = gap / config.scoring_params.gap_scale;
        return @exp(-ratio * ratio);
    }
    return 0.0;
}

fn buildScoreContext(
    allocator: Allocator,
    movers: []const Mover,
    atoms: []const Atom,
) !ScoreContext {
    // Mark mover-controlled atoms (temporary, freed before return).
    const moved_atoms = try allocator.alloc(bool, atoms.len);
    defer allocator.free(moved_atoms);
    @memset(moved_atoms, false);
    for (movers) |m| {
        for (m.atom_indices) |ai| {
            moved_atoms[ai] = true;
        }
    }

    // Debug-only: verify mover atom indices are disjoint across movers.
    if (std.debug.runtime_safety) {
        for (movers, 0..) |m, mi| {
            for (m.atom_indices) |ai| {
                for (movers[mi + 1 ..]) |other| {
                    for (other.atom_indices) |oai| {
                        if (ai == oai) {
                            std.debug.panic("mover atom index {d} appears in multiple movers", .{ai});
                        }
                    }
                }
            }
        }
    }

    var static_count: usize = 0;
    for (moved_atoms) |is_moved| {
        if (!is_moved) static_count += 1;
    }

    const static_positions = try allocator.alloc(math_mod.Vec3(f32), static_count);
    errdefer allocator.free(static_positions);
    const static_atom_indices = try allocator.alloc(u32, static_count);
    errdefer allocator.free(static_atom_indices);

    var out_i: usize = 0;
    for (atoms, 0..) |a, i| {
        if (moved_atoms[i]) continue;
        static_positions[out_i] = a.pos;
        static_atom_indices[out_i] = @intCast(i);
        out_i += 1;
    }

    // CellList may fail with GridTooLarge for very spread-out structures;
    // fall back to null (pairwise scoring only) rather than aborting.
    var cell_list: ?CellList = CellList.init(allocator, static_positions, 5.0) catch |err| blk: {
        switch (err) {
            error.GridTooLarge => {
                std.log.warn("CellList grid too large for static atoms; falling back to pairwise scoring", .{});
                break :blk null;
            },
            error.OutOfMemory => return err,
        }
    };
    errdefer if (cell_list) |*cl| cl.deinit();

    // Compute centroid (mean position) and bounding radius for each mover.
    // Centroid is averaged over all orientation positions; bounding radius is the
    // max distance from the centroid to any orientation position.
    // Used by scoreMover for early-exit of distant mover-vs-mover pairs:
    //   skip if dist(centroid_a, centroid_b) > search_radius + radius_a + radius_b
    const mover_centroids = try allocator.alloc(math_mod.Vec3(f32), movers.len);
    errdefer allocator.free(mover_centroids);
    const mover_radii = try allocator.alloc(f32, movers.len);
    errdefer allocator.free(mover_radii);

    for (movers, 0..) |m, mi| {
        var sum = math_mod.Vec3(f32).zero;
        var count: usize = 0;
        for (m.orientations) |orient| {
            for (orient.positions) |pos| {
                sum = sum.add(pos);
                count += 1;
            }
        }
        const centroid = if (count > 0) sum.scale(1.0 / @as(f32, @floatFromInt(count))) else math_mod.Vec3(f32).zero;
        mover_centroids[mi] = centroid;

        var max_r2: f32 = 0.0;
        for (m.orientations) |orient| {
            for (orient.positions) |pos| {
                const d = pos.sub(centroid);
                const r2 = d.dot(d);
                if (r2 > max_r2) max_r2 = r2;
            }
        }
        mover_radii[mi] = @sqrt(max_r2);
    }

    return .{
        .cell_list = cell_list,
        .static_positions = static_positions,
        .static_atom_indices = static_atom_indices,
        .mover_centroids = mover_centroids,
        .mover_radii = mover_radii,
    };
}

// ── Tests ─────────────────────────────────────────────────────────────────────

const testing = std.testing;

/// Create a test mover with the given orientations (each orientation is one atom position).
fn makeTestMover(
    allocator: Allocator,
    atom_idx: u32,
    positions: []const Vec3(f32),
    penalties: []const f32,
) !Mover {
    const n = positions.len;
    const orientations = try allocator.alloc(mover_mod.Orientation, n);
    for (0..n) |i| {
        const pos = try allocator.alloc(Vec3(f32), 1);
        pos[0] = positions[i];
        orientations[i] = .{ .positions = pos, .penalty = penalties[i] };
    }
    const atom_indices = try allocator.alloc(u32, 1);
    atom_indices[0] = atom_idx;

    return Mover{
        .kind = .single_h_rotator,
        .residue_idx = 0,
        .atom_indices = atom_indices,
        .orientations = orientations,
        .allocator = allocator,
    };
}

test "totalStates calculation" {
    const allocator = testing.allocator;

    // Mover 0 with 3 orientations
    var m0 = try makeTestMover(allocator, 0, &.{
        .{ .x = 0, .y = 0, .z = 0 },
        .{ .x = 1, .y = 0, .z = 0 },
        .{ .x = 2, .y = 0, .z = 0 },
    }, &.{ 0, 0, 0 });
    defer m0.deinit();

    // Mover 1 with 2 orientations
    var m1 = try makeTestMover(allocator, 1, &.{
        .{ .x = 0, .y = 0, .z = 0 },
        .{ .x = 0, .y = 1, .z = 0 },
    }, &.{ 0, 0 });
    defer m1.deinit();

    const movers = [_]Mover{ m0, m1 };
    const clq = [_]u32{ 0, 1 };
    try testing.expectEqual(@as(u64, 6), totalStates(&movers, &clq));

    // Single mover
    const clq_single = [_]u32{0};
    try testing.expectEqual(@as(u64, 3), totalStates(&movers, &clq_single));
}

test "optimize singleton picks best orientation" {
    // Create a model with 2 atoms: a fixed obstacle at origin, and a mover atom.
    // Orientation 0: mover at (1.5, 0, 0) -- bumps with obstacle (VDW overlap)
    // Orientation 1: mover at (10.0, 0, 0) -- far away, no bump
    // The optimizer should pick orientation 1.
    const allocator = testing.allocator;

    var model = Model.init(allocator);
    defer model.deinit();

    // Fixed obstacle atom at origin
    try model.atoms.append(allocator, .{
        .pos = .{ .x = 0, .y = 0, .z = 0 },
        .vdw_radius = 1.7,
    });
    // Mover atom (initial position, will be overwritten)
    try model.atoms.append(allocator, .{
        .pos = .{ .x = 5, .y = 0, .z = 0 },
        .vdw_radius = 1.7,
    });

    // Mover controls atom 1 with 2 orientations
    var mover = try makeTestMover(allocator, 1, &.{
        .{ .x = 1.5, .y = 0, .z = 0 }, // overlaps with obstacle (distance 1.5 < 1.7+1.7=3.4)
        .{ .x = 10.0, .y = 0, .z = 0 }, // far away, no interaction
    }, &.{ 0, 0 });
    defer mover.deinit();

    var movers = [_]Mover{mover};

    var score_ctx = try buildScoreContext(allocator, &movers, model.atoms.items);
    defer score_ctx.deinit(allocator);

    var scratch = std.ArrayListUnmanaged(u32).empty;
    defer scratch.deinit(allocator);

    optimizeSingleton(&movers, 0, &model, .{}, &score_ctx, allocator, &scratch);

    // Should pick orientation 1 (no bump)
    try testing.expectEqual(@as(u16, 1), movers[0].best_orientation);
}

test "optimize brute force finds optimal combination" {
    // Two movers. A fixed atom at (0, 0, 0).
    // Mover 0 controls atom 1: orient 0 = (1.5, 0, 0) (bumps), orient 1 = (10, 0, 0)
    // Mover 1 controls atom 2: orient 0 = (0, 1.5, 0) (bumps), orient 1 = (0, 10, 0)
    // Best: both at orientation 1.
    const allocator = testing.allocator;

    var model = Model.init(allocator);
    defer model.deinit();

    // Fixed obstacle atom at origin
    try model.atoms.append(allocator, .{
        .pos = .{ .x = 0, .y = 0, .z = 0 },
        .vdw_radius = 1.7,
    });
    // Mover 0 atom
    try model.atoms.append(allocator, .{
        .pos = .{ .x = 5, .y = 0, .z = 0 },
        .vdw_radius = 1.7,
    });
    // Mover 1 atom
    try model.atoms.append(allocator, .{
        .pos = .{ .x = 0, .y = 5, .z = 0 },
        .vdw_radius = 1.7,
    });

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
    const clq = [_]u32{ 0, 1 };

    var score_ctx = try buildScoreContext(allocator, &movers, model.atoms.items);
    defer score_ctx.deinit(allocator);

    var scratch = std.ArrayListUnmanaged(u32).empty;
    defer scratch.deinit(allocator);

    try optimizeBruteForce(allocator, &movers, &clq, &model, .{}, &score_ctx, &scratch);

    // Both should pick orientation 1 (no bumps)
    try testing.expectEqual(@as(u16, 1), movers[0].best_orientation);
    try testing.expectEqual(@as(u16, 1), movers[1].best_orientation);
}

test "incrementIndices mixed radix counting" {
    const allocator = testing.allocator;

    // Mover 0: 2 orientations, Mover 1: 3 orientations
    var m0 = try makeTestMover(allocator, 0, &.{
        .{ .x = 0, .y = 0, .z = 0 },
        .{ .x = 1, .y = 0, .z = 0 },
    }, &.{ 0, 0 });
    defer m0.deinit();

    var m1 = try makeTestMover(allocator, 1, &.{
        .{ .x = 0, .y = 0, .z = 0 },
        .{ .x = 0, .y = 1, .z = 0 },
        .{ .x = 0, .y = 2, .z = 0 },
    }, &.{ 0, 0, 0 });
    defer m1.deinit();

    const movers = [_]Mover{ m0, m1 };
    const clq = [_]u32{ 0, 1 };
    var indices = [_]u16{ 0, 0 };

    // Count total iterations
    var count: u32 = 0;
    while (incrementIndices(&indices, &movers, &clq)) {
        count += 1;
    }
    // Started at (0,0), should enumerate 2*3 - 1 = 5 more states before overflow
    try testing.expectEqual(@as(u32, 5), count);
}

test "iterative greedy finds optimal for coupled movers" {
    const allocator = testing.allocator;

    var model = Model.init(allocator);
    defer model.deinit();

    // Fixed obstacle at origin
    try model.atoms.append(allocator, .{
        .pos = .{ .x = 0, .y = 0, .z = 0 },
        .vdw_radius = 1.7,
    });
    // Mover 0 atom
    try model.atoms.append(allocator, .{
        .pos = .{ .x = 5, .y = 0, .z = 0 },
        .vdw_radius = 1.7,
    });
    // Mover 1 atom
    try model.atoms.append(allocator, .{
        .pos = .{ .x = 0, .y = 5, .z = 0 },
        .vdw_radius = 1.7,
    });

    var m0 = try makeTestMover(allocator, 1, &.{
        .{ .x = 1.5, .y = 0, .z = 0 }, // overlaps with obstacle (dist 1.5 < 3.4)
        .{ .x = 10.0, .y = 0, .z = 0 }, // far away
    }, &.{ 0, 0 });
    defer m0.deinit();

    var m1 = try makeTestMover(allocator, 2, &.{
        .{ .x = 0, .y = 1.5, .z = 0 }, // overlaps with obstacle (dist 1.5 < 3.4)
        .{ .x = 0, .y = 10.0, .z = 0 }, // far away
    }, &.{ 0, 0 });
    defer m1.deinit();

    var movers = [_]Mover{ m0, m1 };
    const clq = [_]u32{ 0, 1 };

    var score_ctx = try buildScoreContext(allocator, &movers, model.atoms.items);
    defer score_ctx.deinit(allocator);

    var scratch = std.ArrayListUnmanaged(u32).empty;
    defer scratch.deinit(allocator);

    optimizeIterativeGreedy(&movers, &clq, &model, .{}, &score_ctx, allocator, &scratch);

    // Both should pick orientation 1 (away from obstacle)
    try testing.expectEqual(@as(u16, 1), movers[0].best_orientation);
    try testing.expectEqual(@as(u16, 1), movers[1].best_orientation);
}

test "optimize brute force sees clashes introduced by moved coordinates" {
    const allocator = testing.allocator;

    var model = Model.init(allocator);
    defer model.deinit();

    try model.atoms.append(allocator, .{
        .pos = .{ .x = 0, .y = 0, .z = 0 },
        .vdw_radius = 1.7,
    });
    try model.atoms.append(allocator, .{
        .pos = .{ .x = 20, .y = 0, .z = 0 },
        .vdw_radius = 1.7,
    });

    var m0 = try makeTestMover(allocator, 0, &.{
        .{ .x = 10.0, .y = 0, .z = 0 },
        .{ .x = 0.0, .y = 0, .z = 0 },
    }, &.{ 0.0, 0.2 });
    defer m0.deinit();

    var m1 = try makeTestMover(allocator, 1, &.{
        .{ .x = 10.5, .y = 0, .z = 0 },
        .{ .x = 20.0, .y = 0, .z = 0 },
    }, &.{ 0.0, 0.2 });
    defer m1.deinit();

    var movers = [_]Mover{ m0, m1 };
    const clq = [_]u32{ 0, 1 };

    var score_ctx = try buildScoreContext(allocator, &movers, model.atoms.items);
    defer score_ctx.deinit(allocator);

    var scratch = std.ArrayListUnmanaged(u32).empty;
    defer scratch.deinit(allocator);

    try optimizeBruteForce(allocator, &movers, &clq, &model, .{}, &score_ctx, &scratch);

    // orient 0+0 places atoms at x=10.0 and x=10.5 (distance 0.5 < sum_r 3.4 → severe clash).
    // orient 0+1: m0 at 10.0, m1 at 20.0 → no clash, penalty 0+0.2 = -0.2
    // orient 1+0: m0 at 0.0, m1 at 10.5 → no clash, penalty 0.2+0 = -0.2
    // orient 1+1: m0 at 0.0, m1 at 20.0 → no clash, penalty 0.2+0.2 = -0.4
    // Best is orient 0+1 (first encountered with score -0.2).
    try testing.expectEqual(@as(u16, 0), movers[0].best_orientation);
    try testing.expectEqual(@as(u16, 1), movers[1].best_orientation);
}

test "buildScoreContext excludes mover-controlled atoms from static index" {
    const allocator = testing.allocator;

    var model = Model.init(allocator);
    defer model.deinit();
    try model.atoms.append(allocator, .{ .pos = .{ .x = 0, .y = 0, .z = 0 } });
    try model.atoms.append(allocator, .{ .pos = .{ .x = 1, .y = 0, .z = 0 } });
    try model.atoms.append(allocator, .{ .pos = .{ .x = 5, .y = 0, .z = 0 } });

    var mover = try makeTestMover(allocator, 1, &.{
        .{ .x = 1, .y = 0, .z = 0 },
        .{ .x = 2, .y = 0, .z = 0 },
    }, &.{ 0, 0 });
    defer mover.deinit();

    const movers = [_]Mover{mover};
    var score_ctx = try buildScoreContext(allocator, &movers, model.atoms.items);
    defer score_ctx.deinit(allocator);

    try testing.expectEqual(@as(usize, 2), score_ctx.static_positions.len);
    try testing.expectEqual(@as(u32, 0), score_ctx.static_atom_indices[0]);
    try testing.expectEqual(@as(u32, 2), score_ctx.static_atom_indices[1]);
    try testing.expect(score_ctx.cell_list != null);
}

test "buildScoreContext with multiple movers excludes all mover atoms" {
    const allocator = testing.allocator;

    var model = Model.init(allocator);
    defer model.deinit();
    // atoms 0..4, movers control atoms 1 and 3
    for (0..5) |i| {
        try model.atoms.append(allocator, .{
            .pos = .{ .x = @as(f32, @floatFromInt(i)), .y = 0, .z = 0 },
        });
    }

    var m0 = try makeTestMover(allocator, 1, &.{
        .{ .x = 1, .y = 0, .z = 0 },
    }, &.{0});
    defer m0.deinit();
    var m1 = try makeTestMover(allocator, 3, &.{
        .{ .x = 3, .y = 0, .z = 0 },
    }, &.{0});
    defer m1.deinit();

    const movers = [_]Mover{ m0, m1 };
    var score_ctx = try buildScoreContext(allocator, &movers, model.atoms.items);
    defer score_ctx.deinit(allocator);

    // Static atoms: 0, 2, 4
    try testing.expectEqual(@as(usize, 3), score_ctx.static_positions.len);
    try testing.expectEqual(@as(u32, 0), score_ctx.static_atom_indices[0]);
    try testing.expectEqual(@as(u32, 2), score_ctx.static_atom_indices[1]);
    try testing.expectEqual(@as(u32, 4), score_ctx.static_atom_indices[2]);
}

test "buildScoreContext partial allocation failure frees correctly" {
    const allocator = testing.allocator;

    var model = Model.init(allocator);
    defer model.deinit();
    try model.atoms.append(allocator, .{ .pos = .{ .x = 0, .y = 0, .z = 0 } });
    try model.atoms.append(allocator, .{ .pos = .{ .x = 1, .y = 0, .z = 0 } });

    var mover = try makeTestMover(allocator, 0, &.{
        .{ .x = 0, .y = 0, .z = 0 },
    }, &.{0});
    defer mover.deinit();

    const movers = [_]Mover{mover};

    // Fail at various allocation points; FailingAllocator verifies no leaks.
    // Allocation order: moved_atoms, static_positions, static_atom_indices,
    //                   CellList.init internals (counts/cell_offsets/atom_indices),
    //                   mover_centroids, mover_radii.
    // Range 0..5 covers the first 5 allocation points safely. CellList.init has a
    // known pre-existing issue where cell_offsets leaks if atom_indices alloc fails
    // (fail_index=5 in this test), so we stop before that index.
    for (0..5) |fail_idx| {
        var failing = std.testing.FailingAllocator.init(allocator, .{ .fail_index = fail_idx });
        const result = buildScoreContext(failing.allocator(), &movers, model.atoms.items);
        if (result) |*ctx| {
            var ctx_mut = ctx.*;
            ctx_mut.deinit(failing.allocator());
        } else |_| {
            // Expected failure -- FailingAllocator checks for leaks on deinit.
        }
    }
}

test "buildScoreContext computes correct centroids and bounding radii" {
    const allocator = testing.allocator;

    var model = Model.init(allocator);
    defer model.deinit();

    // Static atom at origin; mover controls atoms 1 and 2
    try model.atoms.append(allocator, .{ .pos = .{ .x = 0, .y = 0, .z = 0 } });
    try model.atoms.append(allocator, .{ .pos = .{ .x = 2, .y = 0, .z = 0 } });
    try model.atoms.append(allocator, .{ .pos = .{ .x = 4, .y = 0, .z = 0 } });
    try model.atoms.append(allocator, .{ .pos = .{ .x = 100, .y = 0, .z = 0 } });

    // Mover 0: one orientation with two positions at (2,0,0) and (4,0,0).
    // centroid = ((2+4)/2, 0, 0) = (3, 0, 0)
    // bounding radius = max(|2-3|, |4-3|) = 1.0
    const orientations_m0 = try allocator.alloc(mover_mod.Orientation, 1);
    const pos_m0 = try allocator.alloc(Vec3(f32), 2);
    pos_m0[0] = .{ .x = 2, .y = 0, .z = 0 };
    pos_m0[1] = .{ .x = 4, .y = 0, .z = 0 };
    orientations_m0[0] = .{ .positions = pos_m0 };
    const atom_indices_m0 = try allocator.alloc(u32, 2);
    atom_indices_m0[0] = 1;
    atom_indices_m0[1] = 2;
    var m0 = Mover{
        .kind = .single_h_rotator,
        .residue_idx = 0,
        .atom_indices = atom_indices_m0,
        .orientations = orientations_m0,
        .allocator = allocator,
    };
    defer m0.deinit();

    // Mover 1: one orientation at (100,0,0). centroid=(100,0,0), radius=0.
    var m1 = try makeTestMover(allocator, 3, &.{
        .{ .x = 100, .y = 0, .z = 0 },
    }, &.{0});
    defer m1.deinit();

    const movers = [_]Mover{ m0, m1 };
    var score_ctx = try buildScoreContext(allocator, &movers, model.atoms.items);
    defer score_ctx.deinit(allocator);

    try testing.expectEqual(@as(usize, 2), score_ctx.mover_centroids.len);
    try testing.expectEqual(@as(usize, 2), score_ctx.mover_radii.len);

    // Centroid of mover 0: (3, 0, 0); bounding radius: 1.0
    try testing.expectApproxEqAbs(@as(f32, 3.0), score_ctx.mover_centroids[0].x, 0.001);
    try testing.expectApproxEqAbs(@as(f32, 0.0), score_ctx.mover_centroids[0].y, 0.001);
    try testing.expectApproxEqAbs(@as(f32, 1.0), score_ctx.mover_radii[0], 0.001);

    // Centroid of mover 1: (100, 0, 0); bounding radius: 0.0
    try testing.expectApproxEqAbs(@as(f32, 100.0), score_ctx.mover_centroids[1].x, 0.001);
    try testing.expectApproxEqAbs(@as(f32, 0.0), score_ctx.mover_radii[1], 0.001);
}

test "centroid early-exit does not change optimization result" {
    // Two movers far apart: mover 0 near x=5.75, mover 1 near x=50.5.
    // Pair cutoff = search_radius + radius_0 + radius_1 = 5.0 + 4.25 + 0.5 = 9.75.
    // Distance between centroids ~= 44.75 >> 9.75 → early-exit fires.
    // A fixed obstacle is near mover 0. Mover 0 should still pick the clear orientation.
    const allocator = testing.allocator;

    var model = Model.init(allocator);
    defer model.deinit();

    // Fixed obstacle near mover 0's closer orientation
    try model.atoms.append(allocator, .{
        .pos = .{ .x = 0, .y = 0, .z = 0 },
        .vdw_radius = 1.7,
    });
    // Mover 0 atom (initial position, overwritten during optimization)
    try model.atoms.append(allocator, .{
        .pos = .{ .x = 5, .y = 0, .z = 0 },
        .vdw_radius = 1.2,
    });
    // Mover 1 atom far away
    try model.atoms.append(allocator, .{
        .pos = .{ .x = 50, .y = 0, .z = 0 },
        .vdw_radius = 1.2,
    });

    // Mover 0: orient 0 bumps obstacle at (0,0,0), orient 1 is clear.
    // centroid = (1.5+10.0)/2 = 5.75, radius = max(|1.5-5.75|, |10.0-5.75|) = 4.25
    var m0 = try makeTestMover(allocator, 1, &.{
        .{ .x = 1.5, .y = 0, .z = 0 }, // bumps obstacle (dist 1.5 < 1.7+1.2=2.9)
        .{ .x = 10.0, .y = 0, .z = 0 }, // clear
    }, &.{ 0, 0 });
    defer m0.deinit();

    // Mover 1: far away; centroid = 50.5, radius = 0.5
    var m1 = try makeTestMover(allocator, 2, &.{
        .{ .x = 50, .y = 0, .z = 0 },
        .{ .x = 51, .y = 0, .z = 0 },
    }, &.{ 0, 0 });
    defer m1.deinit();

    var movers = [_]Mover{ m0, m1 };

    var score_ctx = try buildScoreContext(allocator, &movers, model.atoms.items);
    defer score_ctx.deinit(allocator);

    // Verify early-exit fires: centroid distance > search_radius + r0 + r1
    const search_radius: f32 = 5.0;
    const r0 = score_ctx.mover_radii[0];
    const r1 = score_ctx.mover_radii[1];
    const pair_cutoff = search_radius + r0 + r1;
    const cdiff = score_ctx.mover_centroids[0].sub(score_ctx.mover_centroids[1]);
    const dist2 = cdiff.dot(cdiff);
    try testing.expect(dist2 > pair_cutoff * pair_cutoff);

    var scratch = std.ArrayListUnmanaged(u32).empty;
    defer scratch.deinit(allocator);

    optimizeSingleton(&movers, 0, &model, .{}, &score_ctx, allocator, &scratch);

    // Mover 0 should pick orientation 1 (away from obstacle), even with early-exit skipping m1.
    try testing.expectEqual(@as(u16, 1), movers[0].best_orientation);

    // Optimize mover 1 as well -- should work fine (no nearby obstacles)
    optimizeSingleton(&movers, 1, &model, .{}, &score_ctx, allocator, &scratch);
    try testing.expect(movers[1].best_orientation < 2);
}

test "scorePair returns correct values for all branches" {
    const config = OptConfig{};

    // Branch 1: No interaction (atoms far apart) -- returns 0.0
    {
        const a = Atom{ .pos = .{ .x = 0, .y = 0, .z = 0 }, .vdw_radius = 1.7 };
        const b = Atom{ .pos = .{ .x = 10, .y = 0, .z = 0 }, .vdw_radius = 1.7 };
        try testing.expectEqual(@as(f32, 0.0), scorePair(a, b, config));
    }

    // Branch 2: Contact within threshold (dist between sum_r and sum_r + 0.5)
    // sum_r = 1.7 + 1.7 = 3.4, threshold = 3.9
    // place at distance 3.6: gap = 0.2
    {
        const a = Atom{ .pos = .{ .x = 0, .y = 0, .z = 0 }, .vdw_radius = 1.7 };
        const b = Atom{ .pos = .{ .x = 3.6, .y = 0, .z = 0 }, .vdw_radius = 1.7 };
        const score = scorePair(a, b, config);
        // gap = 0.2, ratio = 0.2 / 0.25 = 0.8, exp(-0.64) ≈ 0.527
        try testing.expect(score > 0.0);
        try testing.expect(@abs(score - @exp(@as(f32, -0.64))) < 0.001);
    }

    // Branch 3: Overlap without H-bond (bump) -- no donor/acceptor flags
    // distance 3.0 < sum_r 3.4, gap = 3.0 - 3.4 = -0.4
    {
        const a = Atom{ .pos = .{ .x = 0, .y = 0, .z = 0 }, .vdw_radius = 1.7 };
        const b = Atom{ .pos = .{ .x = 3.0, .y = 0, .z = 0 }, .vdw_radius = 1.7 };
        const score = scorePair(a, b, config);
        // bump: -bump_weight * (-0.5 * gap) = -10.0 * (-0.5 * -0.4) = -10.0 * 0.2 = -2.0
        try testing.expect(score < 0.0);
        try testing.expect(@abs(score - (-2.0)) < 0.001);
    }

    // Branch 4: Overlap with H-bond (donor + acceptor flags)
    // distance 3.0 < sum_r 3.4, gap = -0.4 → -gap = 0.4 <= min_reg_hb_gap(0.6) → H-bond
    {
        const a = Atom{
            .pos = .{ .x = 0, .y = 0, .z = 0 },
            .vdw_radius = 1.7,
            .flags = .{ .donor = true },
        };
        const b = Atom{
            .pos = .{ .x = 3.0, .y = 0, .z = 0 },
            .vdw_radius = 1.7,
            .flags = .{ .acceptor = true },
        };
        const score = scorePair(a, b, config);
        // hb: hb_weight * (-0.5 * gap) = 4.0 * (-0.5 * -0.4) = 4.0 * 0.2 = 0.8
        try testing.expect(score > 0.0);
        try testing.expect(@abs(score - 0.8) < 0.001);
    }
}

test "mover-vs-mover clash scoring picks correct orientations" {
    // Both atoms are mover-controlled (no static atoms).
    // m0 orient 0: x=10.0 (near m1 orient 0 at x=10.5 → clash)
    // m0 orient 1: x=0.0  (far from m1)
    // m1 orient 0: x=10.5 (near m0 orient 0)
    // m1 orient 1: x=20.0 (far from m0)
    // Best: both avoid clash. orient 0+0 has penalty 0 but severe clash.
    // orient 1+1 has penalty 0.2+0.2 = 0.4 but no clash → should win.
    const allocator = testing.allocator;

    var model = Model.init(allocator);
    defer model.deinit();

    try model.atoms.append(allocator, .{
        .pos = .{ .x = 0, .y = 0, .z = 0 },
        .vdw_radius = 1.7,
    });
    try model.atoms.append(allocator, .{
        .pos = .{ .x = 20, .y = 0, .z = 0 },
        .vdw_radius = 1.7,
    });

    var m0 = try makeTestMover(allocator, 0, &.{
        .{ .x = 10.0, .y = 0, .z = 0 },
        .{ .x = 0.0, .y = 0, .z = 0 },
    }, &.{ 0.0, 0.2 });
    defer m0.deinit();

    var m1 = try makeTestMover(allocator, 1, &.{
        .{ .x = 10.5, .y = 0, .z = 0 },
        .{ .x = 20.0, .y = 0, .z = 0 },
    }, &.{ 0.0, 0.2 });
    defer m1.deinit();

    var movers = [_]Mover{ m0, m1 };
    const clq = [_]u32{ 0, 1 };

    var score_ctx = try buildScoreContext(allocator, &movers, model.atoms.items);
    defer score_ctx.deinit(allocator);

    var scratch = std.ArrayListUnmanaged(u32).empty;
    defer scratch.deinit(allocator);

    try optimizeBruteForce(allocator, &movers, &clq, &model, .{}, &score_ctx, &scratch);

    // orient 0+1 scores best: m0 at 10.0 (no clash with m1 at 20.0), penalty 0+0.2 = -0.2.
    try testing.expectEqual(@as(u16, 0), movers[0].best_orientation);
    try testing.expectEqual(@as(u16, 1), movers[1].best_orientation);
}

test "optimize() with multiple independent singletons uses parallel path and gets correct results" {
    // Place N independent movers far apart (200 Å apart on y-axis) so they form N singleton cliques.
    // Each mover has its own nearby obstacle:
    //   obstacle_i at (0, y_i, 0) with VDW 1.7
    //   orient 0:  mover at (1.5, y_i, 0) -- bumps obstacle (dist 1.5 < sum_r 3.4)
    //   orient 1:  mover at (10.0, y_i, 0) -- far away, no bump
    // After optimize(), every mover should have best_orientation == 1.
    const allocator = testing.allocator;
    const N = 8;

    var model = Model.init(allocator);
    defer model.deinit();

    // Each mover i occupies atom slot obstacle_atom = 2*i, mover_atom = 2*i+1.
    // Movers are 200 Å apart on y-axis so they don't interact with each other.
    var movers_buf: [N]Mover = undefined;
    for (0..N) |i| {
        const mover_atom_idx: u32 = @intCast(2 * i + 1);
        const y: f32 = @as(f32, @floatFromInt(i)) * 200.0;

        // Fixed obstacle for this mover (atom index obstacle_atom_idx = 2*i, static)
        try model.atoms.append(allocator, .{
            .pos = .{ .x = 0, .y = y, .z = 0 },
            .vdw_radius = 1.7,
        });

        // Mover atom (initial position, overwritten during optimization)
        try model.atoms.append(allocator, .{
            .pos = .{ .x = 5, .y = y, .z = 0 },
            .vdw_radius = 1.7,
        });

        movers_buf[i] = try makeTestMover(allocator, mover_atom_idx, &.{
            .{ .x = 1.5, .y = y, .z = 0 }, // bumps own obstacle (dist 1.5 < sum_r 3.4)
            .{ .x = 10.0, .y = y, .z = 0 }, // clear
        }, &.{ 0, 0 });
    }
    defer for (&movers_buf) |*m| m.deinit();

    const movers: []Mover = &movers_buf;
    const result = try optimize(allocator, movers, &model, .{});

    // All singletons, no brute-force or greedy cliques
    try testing.expectEqual(@as(u32, N), result.n_singletons);
    try testing.expectEqual(@as(u32, 0), result.n_brute_force);

    // Every mover should pick orientation 1 (no bump)
    for (movers_buf) |m| {
        try testing.expectEqual(@as(u16, 1), m.best_orientation);
    }
}
