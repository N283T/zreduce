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
const element = @import("../element.zig");

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
            optimizeSingleton(movers, clq[0], model, config);
            result.n_singletons += 1;
        } else if (totalStates(movers, clq) <= config.brute_force_limit) {
            // Brute force: enumerate all combinations
            optimizeBruteForce(allocator, movers, clq, model, config) catch {
                // Fallback to greedy if allocation fails
                for (clq) |mi| optimizeSingleton(movers, mi, model, config);
            };
            result.n_brute_force += 1;
        } else {
            // Vertex-cut decomposition (simplified: greedy for Phase 3)
            optimizeGreedy(movers, clq, model, config);
            result.n_vertex_cut += 1;
        }
    }

    // Apply best orientations
    for (movers) |*m| {
        m.applyOrientation(model.atoms.items, m.best_orientation);
        m.current_orientation = m.best_orientation;
    }

    return result;
}

fn optimizeSingleton(movers: []Mover, mover_idx: u32, model: *Model, config: OptConfig) void {
    const m = &movers[mover_idx];
    var best_score: f32 = -std.math.inf(f32);
    var best_idx: u16 = 0;

    for (0..m.nOrientations()) |oi| {
        const idx: u16 = @intCast(oi);
        m.applyOrientation(model.atoms.items, idx);
        const score = scoreMover(m, mover_idx, movers, model, config) - m.orientationPenalty(idx);
        if (score > best_score) {
            best_score = score;
            best_idx = idx;
        }
    }

    m.best_orientation = best_idx;
}

fn optimizeBruteForce(
    allocator: Allocator,
    movers: []Mover,
    clq: []const u32,
    model: *Model,
    config: OptConfig,
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
            total_score += scoreMover(&movers[mi], mi, movers, model, config);
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

fn optimizeGreedy(movers: []Mover, clq: []const u32, model: *Model, config: OptConfig) void {
    // Simple greedy: optimize each mover in the clique independently
    for (clq) |mi| {
        optimizeSingleton(movers, mi, model, config);
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
fn scoreMover(
    m: *const Mover,
    mover_idx: u32,
    movers: []const Mover,
    model: *const Model,
    config: OptConfig,
) f32 {
    _ = mover_idx;
    var total: f32 = 0;
    const atoms = model.atoms.items;
    for (m.atom_indices) |ai| {
        const a = atoms[ai];
        for (atoms, 0..) |other, oi| {
            if (oi == ai) continue;
            if (isMoverAtom(m, @intCast(oi))) continue;

            const diff = a.pos.sub(other.pos);
            const dist2 = diff.dot(diff);
            const sum_r = a.vdw_radius + other.vdw_radius;
            const sum_r2 = sum_r * sum_r;

            if (dist2 < sum_r2) {
                const dist = @sqrt(dist2);
                const gap = dist - sum_r;
                // Check H-bond
                if (scorer_mod.isHBond(a.flags, other.flags, gap, config.scoring_params)) {
                    total += config.scoring_params.hb_weight * (-0.5 * gap);
                } else {
                    total -= config.scoring_params.bump_weight * (-0.5 * gap);
                }
            } else {
                const threshold = sum_r + 0.5;
                if (dist2 < threshold * threshold) {
                    // Contact region
                    const dist = @sqrt(dist2);
                    const gap = dist - sum_r;
                    const ratio = gap / config.scoring_params.gap_scale;
                    total += @exp(-ratio * ratio);
                }
            }
        }
    }
    // Also check interactions with other movers' atoms
    for (movers) |*other_m| {
        if (@intFromPtr(other_m) == @intFromPtr(m)) continue;
        for (m.atom_indices) |ai| {
            const a = atoms[ai];
            for (other_m.atom_indices) |bi| {
                if (ai == bi) continue;
                const b = atoms[bi];
                const diff = a.pos.sub(b.pos);
                const dist2 = diff.dot(diff);
                const sum_r = a.vdw_radius + b.vdw_radius;
                const sum_r2 = sum_r * sum_r;

                if (dist2 < sum_r2) {
                    const dist = @sqrt(dist2);
                    const gap = dist - sum_r;
                    if (scorer_mod.isHBond(a.flags, b.flags, gap, config.scoring_params)) {
                        total += config.scoring_params.hb_weight * (-0.5 * gap);
                    } else {
                        total -= config.scoring_params.bump_weight * (-0.5 * gap);
                    }
                }
            }
        }
    }
    return total;
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
    optimizeSingleton(&movers, 0, &model, .{});

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
    try optimizeBruteForce(allocator, &movers, &clq, &model, .{});

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
