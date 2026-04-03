//! Spatial cell list for efficient neighbor lookups.

const std = @import("std");
const Allocator = std.mem.Allocator;
const math = @import("../math.zig");
const Vec3 = math.Vec3;

pub const CellList = struct {
    /// Atom indices sorted by cell.
    atom_indices: []u32,
    /// Offset into atom_indices for each cell (length = nx*ny*nz + 1, inclusive end).
    cell_offsets: []u32,
    nx: u32,
    ny: u32,
    nz: u32,
    cell_size: f32,
    x_min: f32,
    y_min: f32,
    z_min: f32,
    allocator: Allocator,

    pub fn init(allocator: Allocator, positions: []const Vec3(f32), cell_size: f32) !CellList {
        std.debug.assert(cell_size > 0);
        if (positions.len == 0) {
            const empty_indices = try allocator.alloc(u32, 0);
            errdefer allocator.free(empty_indices);
            const empty_offsets = try allocator.alloc(u32, 1);
            return CellList{
                .atom_indices = empty_indices,
                .cell_offsets = empty_offsets,
                .nx = 0,
                .ny = 0,
                .nz = 0,
                .cell_size = cell_size,
                .x_min = 0,
                .y_min = 0,
                .z_min = 0,
                .allocator = allocator,
            };
        }

        // Reject non-finite coordinates before building the grid.
        for (positions) |p| {
            if (!std.math.isFinite(p.x) or !std.math.isFinite(p.y) or !std.math.isFinite(p.z)) {
                return error.NonFiniteCoordinate;
            }
        }

        // Compute bounding box.
        var x_min = positions[0].x;
        var y_min = positions[0].y;
        var z_min = positions[0].z;
        var x_max = positions[0].x;
        var y_max = positions[0].y;
        var z_max = positions[0].z;

        for (positions[1..]) |p| {
            if (p.x < x_min) x_min = p.x;
            if (p.y < y_min) y_min = p.y;
            if (p.z < z_min) z_min = p.z;
            if (p.x > x_max) x_max = p.x;
            if (p.y > y_max) y_max = p.y;
            if (p.z > z_max) z_max = p.z;
        }

        // Add a small margin so boundary atoms don't fall outside.
        const margin = cell_size;
        x_min -= margin;
        y_min -= margin;
        z_min -= margin;
        x_max += margin;
        y_max += margin;
        z_max += margin;

        const nx_f = @ceil((x_max - x_min) / cell_size);
        const ny_f = @ceil((y_max - y_min) / cell_size);
        const nz_f = @ceil((z_max - z_min) / cell_size);
        const max_dim: f32 = @floatFromInt(std.math.maxInt(u32));
        if (nx_f > max_dim or ny_f > max_dim or nz_f > max_dim) return error.GridTooLarge;
        const nx: u32 = @max(1, @as(u32, @intFromFloat(nx_f)));
        const ny: u32 = @max(1, @as(u32, @intFromFloat(ny_f)));
        const nz: u32 = @max(1, @as(u32, @intFromFloat(nz_f)));

        const total_cells = std.math.mul(u32, nx, ny) catch return error.GridTooLarge;
        const total_cells_3d = std.math.mul(u32, total_cells, nz) catch return error.GridTooLarge;

        // Count atoms per cell (counting sort phase 1).
        const counts = try allocator.alloc(u32, total_cells_3d);
        defer allocator.free(counts);
        @memset(counts, 0);

        for (positions) |p| {
            const cell = cellIndex(p, x_min, y_min, z_min, cell_size, nx, ny, nz);
            counts[cell] += 1;
        }

        // Build prefix sum (offsets array has total_cells_3d+1 entries).
        const cell_offsets = try allocator.alloc(u32, total_cells_3d + 1);
        errdefer allocator.free(cell_offsets);
        cell_offsets[0] = 0;
        for (0..total_cells_3d) |i| {
            cell_offsets[i + 1] = cell_offsets[i] + counts[i];
        }

        // Fill atom_indices (counting sort phase 2).
        const atom_indices = try allocator.alloc(u32, positions.len);
        @memset(counts, 0); // reuse as cursor

        for (positions, 0..) |p, idx| {
            const cell = cellIndex(p, x_min, y_min, z_min, cell_size, nx, ny, nz);
            const pos = cell_offsets[cell] + counts[cell];
            atom_indices[pos] = @intCast(idx);
            counts[cell] += 1;
        }

        return CellList{
            .atom_indices = atom_indices,
            .cell_offsets = cell_offsets,
            .nx = nx,
            .ny = ny,
            .nz = nz,
            .cell_size = cell_size,
            .x_min = x_min,
            .y_min = y_min,
            .z_min = z_min,
            .allocator = allocator,
        };
    }

    pub fn deinit(self: *CellList) void {
        self.allocator.free(self.atom_indices);
        self.allocator.free(self.cell_offsets);
    }

    /// Find all atom indices within radius of query point.
    /// Appends to result (does not clear it first).
    pub fn neighborsInRadius(
        self: *const CellList,
        query: Vec3(f32),
        radius: f32,
        result: *std.ArrayListUnmanaged(u32),
        allocator: Allocator,
        positions: []const Vec3(f32),
    ) !void {
        std.debug.assert(!std.math.isNan(query.x) and !std.math.isNan(query.y) and !std.math.isNan(query.z));
        if (self.nx == 0) return;

        const r2 = radius * radius;
        const cell_steps: i32 = @intFromFloat(@ceil(radius / self.cell_size));

        // Find the cell of the query point.
        const qcx: i32 = @intFromFloat(@floor((query.x - self.x_min) / self.cell_size));
        const qcy: i32 = @intFromFloat(@floor((query.y - self.y_min) / self.cell_size));
        const qcz: i32 = @intFromFloat(@floor((query.z - self.z_min) / self.cell_size));

        const nx: i32 = @intCast(self.nx);
        const ny: i32 = @intCast(self.ny);
        const nz: i32 = @intCast(self.nz);

        var dx: i32 = -cell_steps;
        while (dx <= cell_steps) : (dx += 1) {
            var dy: i32 = -cell_steps;
            while (dy <= cell_steps) : (dy += 1) {
                var dz: i32 = -cell_steps;
                while (dz <= cell_steps) : (dz += 1) {
                    const cx = qcx + dx;
                    const cy = qcy + dy;
                    const cz = qcz + dz;

                    if (cx < 0 or cy < 0 or cz < 0) continue;
                    if (cx >= nx or cy >= ny or cz >= nz) continue;

                    const cell: u32 = @intCast(
                        @as(i32, @intCast(cx)) +
                            @as(i32, @intCast(cy)) * nx +
                            @as(i32, @intCast(cz)) * nx * ny,
                    );

                    const start = self.cell_offsets[cell];
                    const end = self.cell_offsets[cell + 1];

                    for (self.atom_indices[start..end]) |idx| {
                        const diff = query.sub(positions[idx]);
                        if (diff.dot(diff) <= r2) {
                            try result.append(allocator, idx);
                        }
                    }
                }
            }
        }
    }
};

fn cellIndex(
    p: Vec3(f32),
    x_min: f32,
    y_min: f32,
    z_min: f32,
    cell_size: f32,
    nx: u32,
    ny: u32,
    nz: u32,
) u32 {
    const cx: u32 = @intCast(@min(
        nx - 1,
        @as(u32, @intFromFloat(@floor((p.x - x_min) / cell_size))),
    ));
    const cy: u32 = @intCast(@min(
        ny - 1,
        @as(u32, @intFromFloat(@floor((p.y - y_min) / cell_size))),
    ));
    const cz: u32 = @intCast(@min(
        nz - 1,
        @as(u32, @intFromFloat(@floor((p.z - z_min) / cell_size))),
    ));
    return cx + cy * nx + cz * nx * ny;
}

// ── Tests ─────────────────────────────────────────────────────────────────────

test "CellList basic neighbor search" {
    const allocator = std.testing.allocator;

    const positions = [_]Vec3(f32){
        .{ .x = 0, .y = 0, .z = 0 },
        .{ .x = 1, .y = 0, .z = 0 },
        .{ .x = 10, .y = 10, .z = 10 },
    };

    var cl = try CellList.init(allocator, &positions, 3.0);
    defer cl.deinit();

    var result = std.ArrayListUnmanaged(u32).empty;
    defer result.deinit(allocator);

    const query = Vec3(f32){ .x = 0, .y = 0, .z = 0 };
    try cl.neighborsInRadius(query, 2.0, &result, allocator, &positions);

    // Should find index 0 (distance 0) and index 1 (distance 1), but NOT index 2 (distance ~17)
    try std.testing.expectEqual(@as(usize, 2), result.items.len);

    // Verify the found indices are 0 and 1 (order may vary)
    var found0 = false;
    var found1 = false;
    for (result.items) |idx| {
        if (idx == 0) found0 = true;
        if (idx == 1) found1 = true;
    }
    try std.testing.expect(found0);
    try std.testing.expect(found1);
}

test "CellList empty positions" {
    const allocator = std.testing.allocator;
    var cl = try CellList.init(allocator, &[_]Vec3(f32){}, 3.0);
    defer cl.deinit();

    var result = std.ArrayListUnmanaged(u32).empty;
    defer result.deinit(allocator);

    const query = Vec3(f32){ .x = 0, .y = 0, .z = 0 };
    try cl.neighborsInRadius(query, 2.0, &result, allocator, &[_]Vec3(f32){});
    try std.testing.expectEqual(@as(usize, 0), result.items.len);
}
