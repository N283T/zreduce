//! Interaction graph and clique (connected component) detection for movers.
//! Two movers interact if any of their controlled atoms are within a cutoff distance.

const std = @import("std");
const Allocator = std.mem.Allocator;
const math_mod = @import("../math.zig");
const Vec3 = math_mod.Vec3;
const model_mod = @import("../model.zig");
const Atom = model_mod.Atom;
const mover_mod = @import("mover.zig");
const Mover = mover_mod.Mover;

pub const InteractionGraph = struct {
    /// Adjacency list: adj[i] = list of mover indices that interact with mover i
    adjacency: [][]u32,
    n_movers: u32,
    allocator: Allocator,

    pub fn deinit(self: *InteractionGraph) void {
        for (self.adjacency) |neighbors| {
            self.allocator.free(neighbors);
        }
        self.allocator.free(self.adjacency);
    }
};

/// Build interaction graph: edge between movers whose atoms can interact.
/// cutoff: maximum distance between any pair of atoms for two movers to be considered interacting.
/// Typically VDW_radius_sum + probe_radius + overlap_margin (~6-8 Å).
pub fn buildInteractionGraph(
    allocator: Allocator,
    movers: []const Mover,
    atoms: []const Atom,
    cutoff: f32,
) !InteractionGraph {
    const n = movers.len;
    var adj_lists = try allocator.alloc(std.ArrayListUnmanaged(u32), n);
    for (adj_lists) |*al| al.* = .empty;
    errdefer {
        for (adj_lists) |*al| al.deinit(allocator);
        allocator.free(adj_lists);
    }

    // For each pair of movers, check if any atoms are within cutoff
    for (0..n) |i| {
        for (i + 1..n) |j| {
            if (moversInteract(movers[i], movers[j], atoms, cutoff)) {
                try adj_lists[i].append(allocator, @intCast(j));
                try adj_lists[j].append(allocator, @intCast(i));
            }
        }
    }

    // Convert to owned slices
    const adjacency = try allocator.alloc([]u32, n);
    for (0..n) |i| {
        adjacency[i] = try adj_lists[i].toOwnedSlice(allocator);
    }
    allocator.free(adj_lists);

    return InteractionGraph{
        .adjacency = adjacency,
        .n_movers = @intCast(n),
        .allocator = allocator,
    };
}

fn moversInteract(a: Mover, b: Mover, atoms: []const Atom, cutoff: f32) bool {
    const cutoff2 = cutoff * cutoff;
    for (a.atom_indices) |ai| {
        for (b.atom_indices) |bi| {
            const diff = atoms[ai].pos.sub(atoms[bi].pos);
            if (diff.dot(diff) <= cutoff2) return true;
        }
    }
    return false;
}

/// Find connected components in the interaction graph.
/// Returns a list of cliques, where each clique is a list of mover indices.
pub fn findCliques(allocator: Allocator, graph: *const InteractionGraph) ![][]u32 {
    const n = graph.n_movers;
    const visited = try allocator.alloc(bool, n);
    defer allocator.free(visited);
    @memset(visited, false);

    var cliques = std.ArrayListUnmanaged([]u32).empty;
    errdefer {
        for (cliques.items) |c| allocator.free(c);
        cliques.deinit(allocator);
    }

    for (0..n) |i| {
        if (visited[i]) continue;
        // BFS from i
        var component = std.ArrayListUnmanaged(u32).empty;
        var queue = std.ArrayListUnmanaged(u32).empty;
        defer queue.deinit(allocator);

        try queue.append(allocator, @intCast(i));
        visited[i] = true;

        while (queue.items.len > 0) {
            const node = queue.orderedRemove(0);
            try component.append(allocator, node);
            for (graph.adjacency[node]) |neighbor| {
                if (!visited[neighbor]) {
                    visited[neighbor] = true;
                    try queue.append(allocator, neighbor);
                }
            }
        }

        try cliques.append(allocator, try component.toOwnedSlice(allocator));
    }

    return try cliques.toOwnedSlice(allocator);
}

// ── Tests ─────────────────────────────────────────────────────────────────────

const testing = std.testing;

/// Create a minimal Mover for testing. Only atom_indices matters for interaction checks.
fn makeMover(allocator: Allocator, atom_indices: []const u32) !Mover {
    const indices = try allocator.dupe(u32, atom_indices);
    const orientations = try allocator.alloc(mover_mod.Orientation, 0);
    return Mover{
        .kind = .single_h_rotator,
        .residue_idx = 0,
        .atom_indices = indices,
        .orientations = orientations,
        .allocator = allocator,
    };
}

test "single isolated mover forms own clique" {
    const allocator = testing.allocator;

    var atoms = [_]Atom{
        .{ .pos = .{ .x = 0.0, .y = 0.0, .z = 0.0 } },
    };

    var m0 = try makeMover(allocator, &.{0});
    defer m0.deinit();

    const movers = [_]Mover{m0};
    var graph = try buildInteractionGraph(allocator, &movers, &atoms, 6.0);
    defer graph.deinit();

    const cliques = try findCliques(allocator, &graph);
    defer {
        for (cliques) |c| allocator.free(c);
        allocator.free(cliques);
    }

    try testing.expectEqual(@as(usize, 1), cliques.len);
    try testing.expectEqual(@as(usize, 1), cliques[0].len);
    try testing.expectEqual(@as(u32, 0), cliques[0][0]);
}

test "two interacting movers form one clique" {
    const allocator = testing.allocator;

    // Two atoms close together (distance = 2.0, well within 6.0 cutoff)
    var atoms = [_]Atom{
        .{ .pos = .{ .x = 0.0, .y = 0.0, .z = 0.0 } },
        .{ .pos = .{ .x = 2.0, .y = 0.0, .z = 0.0 } },
    };

    var m0 = try makeMover(allocator, &.{0});
    defer m0.deinit();
    var m1 = try makeMover(allocator, &.{1});
    defer m1.deinit();

    const movers = [_]Mover{ m0, m1 };
    var graph = try buildInteractionGraph(allocator, &movers, &atoms, 6.0);
    defer graph.deinit();

    const cliques = try findCliques(allocator, &graph);
    defer {
        for (cliques) |c| allocator.free(c);
        allocator.free(cliques);
    }

    try testing.expectEqual(@as(usize, 1), cliques.len);
    try testing.expectEqual(@as(usize, 2), cliques[0].len);
}

test "two distant movers form two cliques" {
    const allocator = testing.allocator;

    // Two atoms far apart (distance = 100.0, beyond 6.0 cutoff)
    var atoms = [_]Atom{
        .{ .pos = .{ .x = 0.0, .y = 0.0, .z = 0.0 } },
        .{ .pos = .{ .x = 100.0, .y = 0.0, .z = 0.0 } },
    };

    var m0 = try makeMover(allocator, &.{0});
    defer m0.deinit();
    var m1 = try makeMover(allocator, &.{1});
    defer m1.deinit();

    const movers = [_]Mover{ m0, m1 };
    var graph = try buildInteractionGraph(allocator, &movers, &atoms, 6.0);
    defer graph.deinit();

    const cliques = try findCliques(allocator, &graph);
    defer {
        for (cliques) |c| allocator.free(c);
        allocator.free(cliques);
    }

    try testing.expectEqual(@as(usize, 2), cliques.len);
    try testing.expectEqual(@as(usize, 1), cliques[0].len);
    try testing.expectEqual(@as(usize, 1), cliques[1].len);
}

test "chain of interactions" {
    const allocator = testing.allocator;

    // 3 movers: A at x=0, B at x=4 (within 6 of A), C at x=8 (within 6 of B, not of A directly)
    // A-B distance = 4, B-C distance = 4, A-C distance = 8
    var atoms = [_]Atom{
        .{ .pos = .{ .x = 0.0, .y = 0.0, .z = 0.0 } }, // mover A
        .{ .pos = .{ .x = 4.0, .y = 0.0, .z = 0.0 } }, // mover B
        .{ .pos = .{ .x = 8.0, .y = 0.0, .z = 0.0 } }, // mover C
    };

    var m0 = try makeMover(allocator, &.{0});
    defer m0.deinit();
    var m1 = try makeMover(allocator, &.{1});
    defer m1.deinit();
    var m2 = try makeMover(allocator, &.{2});
    defer m2.deinit();

    const movers = [_]Mover{ m0, m1, m2 };
    var graph = try buildInteractionGraph(allocator, &movers, &atoms, 6.0);
    defer graph.deinit();

    // Verify A-C do not directly interact (distance 8 > 6 cutoff)
    try testing.expectEqual(@as(usize, 1), graph.adjacency[0].len); // A only connects to B
    try testing.expectEqual(@as(usize, 2), graph.adjacency[1].len); // B connects to A and C

    const cliques = try findCliques(allocator, &graph);
    defer {
        for (cliques) |c| allocator.free(c);
        allocator.free(cliques);
    }

    // All three form one clique through transitivity
    try testing.expectEqual(@as(usize, 1), cliques.len);
    try testing.expectEqual(@as(usize, 3), cliques[0].len);
}
