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

/// Optimize all movers: find best orientations using clique-based search.
pub fn optimize(
    allocator: Allocator,
    movers: []Mover,
    model: *Model,
    config: OptConfig,
) !OptResult {
    var result = OptResult{};

    const positions = try allocator.alloc(math_mod.Vec3(f32), model.atoms.items.len);
    defer allocator.free(positions);
    syncPositions(positions, model.atoms.items);

    var cell_list = try CellList.init(allocator, positions, 5.0);
    defer cell_list.deinit();

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

    for (cliques) |clq| {
        if (clq.len == 1) {
            // Singleton: score all orientations, pick best
            optimizeSingleton(movers, clq[0], model, config, &cell_list, positions, allocator, &scratch);
            result.n_singletons += 1;
        } else if (totalStates(movers, clq) <= config.brute_force_limit) {
            // Brute force: enumerate all combinations
            optimizeBruteForce(allocator, movers, clq, model, config, &cell_list, positions, &scratch) catch |err| switch (err) {
                error.OutOfMemory, error.GridTooLarge => {
                    // Fallback to greedy on allocation or grid-sizing failure
                    std.debug.print("Warning: brute-force optimization failed ({s}), falling back to greedy for clique of {d} movers\n", .{ @errorName(err), clq.len });
                    for (clq) |mi| optimizeSingleton(movers, mi, model, config, &cell_list, positions, allocator, &scratch);
                    result.n_vertex_cut += 1;
                    continue;
                },
            };
            result.n_brute_force += 1;
        } else {
            // Vertex-cut decomposition (simplified: greedy for Phase 3)
            optimizeIterativeGreedy(movers, clq, model, config, &cell_list, positions, allocator, &scratch);
            result.n_vertex_cut += 1;
        }
    }

    // Apply best coarse orientations
    for (movers) |*m| {
        m.applyOrientation(model.atoms.items, m.best_orientation);
        m.current_orientation = m.best_orientation;
    }

    // Fine search phase: refine each mover around its coarse best
    for (0..movers.len) |mi| {
        fineSearchMover(allocator, movers, @intCast(mi), model, config, &cell_list, positions, &scratch);
    }

    return result;
}

fn optimizeSingleton(
    movers: []Mover,
    mover_idx: u32,
    model: *Model,
    config: OptConfig,
    cell_list: *CellList,
    positions: []math_mod.Vec3(f32),
    allocator: Allocator,
    scratch: *std.ArrayListUnmanaged(u32),
) void {
    const m = &movers[mover_idx];
    var best_score: f32 = -std.math.inf(f32);
    var best_idx: u16 = 0;

    for (0..m.nOrientations()) |oi| {
        const idx: u16 = @intCast(oi);
        m.applyOrientation(model.atoms.items, idx);
        rebuildCellList(allocator, cell_list, positions, model.atoms.items) catch |err| {
            std.debug.print("Warning: CellList rebuild failed during singleton optimization: {s}\n", .{@errorName(err)});
            break;
        };
        const score = scoreMover(m, mover_idx, movers, model, config, cell_list, positions, allocator, scratch) - m.orientationPenalty(idx);
        if (score > best_score) {
            best_score = score;
            best_idx = idx;
        }
    }

    m.best_orientation = best_idx;
}

/// Refine a mover's best orientation with fine angular search.
/// Directly updates atom positions if a fine position improves the score.
fn fineSearchMover(
    allocator: Allocator,
    movers: []Mover,
    mover_idx: u32,
    model: *Model,
    config: OptConfig,
    cell_list: *CellList,
    positions: []math_mod.Vec3(f32),
    scratch: *std.ArrayListUnmanaged(u32),
) void {
    const m = &movers[mover_idx];
    const atoms = model.atoms.items;

    const fine = rotator_mod.generateFineOrientations(allocator, atoms, m, m.best_orientation) catch return;
    defer {
        for (fine) |o| allocator.free(o.positions);
        allocator.free(fine);
    }

    if (fine.len == 0) return;

    // Rebuild CellList to reflect current coarse-best positions before baseline scoring
    rebuildCellList(allocator, cell_list, positions, atoms) catch |err| {
        std.debug.print("Warning: CellList rebuild failed during fine search: {s}\n", .{@errorName(err)});
        return;
    };

    // Score current coarse best (already applied)
    var best_score = scoreMover(m, mover_idx, movers, model, config, cell_list, positions, allocator, scratch) - m.orientationPenalty(m.best_orientation);

    // Score fine orientations
    var best_fine: ?usize = null;
    for (fine, 0..) |orient, fi| {
        for (m.atom_indices, 0..) |ai, j| {
            atoms[ai].pos = orient.positions[j];
        }
        rebuildCellList(allocator, cell_list, positions, atoms) catch |err| {
            std.debug.print("Warning: CellList rebuild failed during fine search: {s}\n", .{@errorName(err)});
            // Restore coarse-best positions before returning
            m.applyOrientation(atoms, m.best_orientation);
            return;
        };
        const score = scoreMover(m, mover_idx, movers, model, config, cell_list, positions, allocator, scratch) - orient.penalty;
        if (score > best_score) {
            best_score = score;
            best_fine = fi;
        }
    }

    if (best_fine) |fi| {
        // Apply the best fine orientation
        for (m.atom_indices, 0..) |ai, j| {
            atoms[ai].pos = fine[fi].positions[j];
        }
    } else {
        // Restore coarse best
        m.applyOrientation(atoms, m.best_orientation);
    }
}

fn optimizeBruteForce(
    allocator: Allocator,
    movers: []Mover,
    clq: []const u32,
    model: *Model,
    config: OptConfig,
    cell_list: *CellList,
    positions: []math_mod.Vec3(f32),
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
        try rebuildCellList(allocator, cell_list, positions, model.atoms.items);

        // Score all movers in clique
        var total_score: f32 = 0;
        for (clq, 0..) |mi, i| {
            total_score += scoreMover(&movers[mi], mi, movers, model, config, cell_list, positions, allocator, scratch);
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
    cell_list: *CellList,
    positions: []math_mod.Vec3(f32),
    allocator: Allocator,
    scratch: *std.ArrayListUnmanaged(u32),
) void {
    const max_iterations: u32 = 3;
    var iteration: u32 = 0;

    while (iteration < max_iterations) : (iteration += 1) {
        var changed = false;
        for (clq) |mi| {
            const old_best = movers[mi].best_orientation;
            optimizeSingleton(movers, mi, model, config, cell_list, positions, allocator, scratch);
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

/// Score a mover's current orientation against the model.
/// Uses a lightweight distance-based VDW overlap scoring.
/// Skips atoms belonging to the same mover.
/// Uses CellList spatial index for fast neighbor queries (O(N) instead of O(N²)).
fn scoreMover(
    m: *const Mover,
    mover_idx: u32,
    movers: []const Mover,
    model: *const Model,
    config: OptConfig,
    cell_list: *const CellList,
    positions: []const math_mod.Vec3(f32),
    allocator: Allocator,
    scratch: *std.ArrayListUnmanaged(u32),
) f32 {
    _ = movers;
    _ = mover_idx;
    var total: f32 = 0;
    const atoms = model.atoms.items;
    // Safe for H-mover atoms: max H(1.22) + largest neighbor(2.75) + 0.5 contact margin = 4.47
    const search_radius: f32 = 5.0;

    for (m.atom_indices) |ai| {
        const a = atoms[ai];

        // Query CellList for nearby atoms using current mover position
        scratch.clearRetainingCapacity();
        cell_list.neighborsInRadius(a.pos, search_radius, scratch, allocator, positions) catch
            return -std.math.inf(f32);

        for (scratch.items) |oi| {
            if (oi == ai) continue;
            if (isMoverAtom(m, oi)) continue;

            const other = atoms[oi];
            const diff = a.pos.sub(other.pos);
            const dist2 = diff.dot(diff);
            const sum_r = a.vdw_radius + other.vdw_radius;
            const sum_r2 = sum_r * sum_r;

            if (dist2 < sum_r2) {
                const dist = @sqrt(dist2);
                const gap = dist - sum_r;
                if (scorer_mod.isHBond(a.flags, other.flags, gap, config.scoring_params)) {
                    total += config.scoring_params.hb_weight * (-0.5 * gap);
                } else {
                    total -= config.scoring_params.bump_weight * (-0.5 * gap);
                }
            } else {
                const threshold = sum_r + 0.5;
                if (dist2 < threshold * threshold) {
                    const dist = @sqrt(dist2);
                    const gap = dist - sum_r;
                    const ratio = gap / config.scoring_params.gap_scale;
                    total += @exp(-ratio * ratio);
                }
            }
        }
    }
    return total;
}

fn syncPositions(positions: []math_mod.Vec3(f32), atoms: []const Atom) void {
    std.debug.assert(positions.len == atoms.len);
    for (atoms, 0..) |a, i| {
        positions[i] = a.pos;
    }
}

fn rebuildCellList(
    allocator: Allocator,
    cell_list: *CellList,
    positions: []math_mod.Vec3(f32),
    atoms: []const Atom,
) !void {
    syncPositions(positions, atoms);
    const new_cell_list = try CellList.init(allocator, positions, 5.0);
    cell_list.deinit();
    cell_list.* = new_cell_list;
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

    // Build spatial index from current atom positions
    const pos = try allocator.alloc(math_mod.Vec3(f32), model.atoms.items.len);
    defer allocator.free(pos);
    syncPositions(pos, model.atoms.items);
    var cl = try CellList.init(allocator, pos, 5.0);
    defer cl.deinit();

    var scratch = std.ArrayListUnmanaged(u32).empty;
    defer scratch.deinit(allocator);

    optimizeSingleton(&movers, 0, &model, .{}, &cl, pos, allocator, &scratch);

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

    // Build spatial index from current atom positions
    const pos = try allocator.alloc(math_mod.Vec3(f32), model.atoms.items.len);
    defer allocator.free(pos);
    syncPositions(pos, model.atoms.items);
    var cl = try CellList.init(allocator, pos, 5.0);
    defer cl.deinit();

    var scratch = std.ArrayListUnmanaged(u32).empty;
    defer scratch.deinit(allocator);

    try optimizeBruteForce(allocator, &movers, &clq, &model, .{}, &cl, pos, &scratch);

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

    // Build spatial index from current atom positions
    const pos = try allocator.alloc(math_mod.Vec3(f32), model.atoms.items.len);
    defer allocator.free(pos);
    syncPositions(pos, model.atoms.items);
    var cl = try CellList.init(allocator, pos, 5.0);
    defer cl.deinit();

    var scratch = std.ArrayListUnmanaged(u32).empty;
    defer scratch.deinit(allocator);

    optimizeIterativeGreedy(&movers, &clq, &model, .{}, &cl, pos, allocator, &scratch);

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

    const pos = try allocator.alloc(math_mod.Vec3(f32), model.atoms.items.len);
    defer allocator.free(pos);
    syncPositions(pos, model.atoms.items);

    var cl = try CellList.init(allocator, pos, 5.0);
    defer cl.deinit();

    var scratch = std.ArrayListUnmanaged(u32).empty;
    defer scratch.deinit(allocator);

    try optimizeBruteForce(allocator, &movers, &clq, &model, .{}, &cl, pos, &scratch);

    try testing.expect(!(movers[0].best_orientation == 0 and movers[1].best_orientation == 0));
}

test "rebuildCellList leaves previous index usable on allocation failure" {
    const allocator = testing.allocator;

    var model = Model.init(allocator);
    defer model.deinit();
    try model.atoms.append(allocator, .{ .pos = .{ .x = 0, .y = 0, .z = 0 } });
    try model.atoms.append(allocator, .{ .pos = .{ .x = 1, .y = 0, .z = 0 } });

    const positions = try allocator.alloc(math_mod.Vec3(f32), model.atoms.items.len);
    defer allocator.free(positions);
    syncPositions(positions, model.atoms.items);

    var cl = try CellList.init(allocator, positions, 5.0);
    defer cl.deinit();

    var failing = std.testing.FailingAllocator.init(allocator, .{ .fail_index = 0 });
    try testing.expectError(error.OutOfMemory, rebuildCellList(failing.allocator(), &cl, positions, model.atoms.items));

    var result = std.ArrayListUnmanaged(u32).empty;
    defer result.deinit(allocator);
    try cl.neighborsInRadius(model.atoms.items[0].pos, 2.0, &result, allocator, positions);
    try testing.expect(result.items.len >= 2);
}
