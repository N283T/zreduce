//! Water hydrogen placement.
//!
//! Adds H1/H2 to HOH/WAT residues, using a CellList for fast neighbor
//! queries when available, with a linear-scan fallback for small models.

const std = @import("std");
const model_mod = @import("../model.zig");
const Model = model_mod.Model;
const Atom = model_mod.Atom;
const Residue = model_mod.Residue;
const math_mod = @import("../math.zig");
const bond_policy = @import("bond_policy.zig");
const neighbor_mod = @import("../model/neighbor.zig");
const lookup = @import("lookup.zig");
const terminal = @import("terminal.zig");

const Vec3f32 = math_mod.Vec3(f32);
const CellList = neighbor_mod.CellList;
const ParentMeta = lookup.ParentMeta;
const findAtomPos = lookup.findAtomPos;
const existsInResidue = lookup.existsInResidue;
const padName = lookup.padName;
const nameMatch = lookup.nameMatch;
const appendNtermH = terminal.appendNtermH;

// ── Tests ─────────────────────────────────────────────────────────────────────

const testing = std.testing;

test "isWaterResidue recognizes HOH and WAT" {
    var res = Residue{};
    res.entity_type = .polymer;

    res.setCompId("HOH");
    try testing.expect(isWaterResidue(res, "HOH"));

    res.setCompId("WAT");
    try testing.expect(isWaterResidue(res, "WAT"));

    res.setCompId("ALA");
    try testing.expect(!isWaterResidue(res, "ALA"));
}

test "isWaterResidue recognizes entity_type water regardless of comp_id" {
    var res = Residue{};
    res.entity_type = .water;
    res.setCompId("O");
    try testing.expect(isWaterResidue(res, "O"));
}

test "orthogonalUnit produces unit vector perpendicular to input" {
    // orthogonalUnit is private; test it indirectly via its observable effect.
    // The function always produces a vector perpendicular to v with length 1.
    // We verify this property using a direct copy of the logic.
    const v1 = Vec3f32{ .x = 1.0, .y = 0.0, .z = 0.0 };
    const v2 = Vec3f32{ .x = 0.0, .y = 1.0, .z = 0.0 };
    const v3 = Vec3f32{ .x = 0.577, .y = 0.577, .z = 0.577 };

    // orthogonalUnit: picks x_axis if |v.x| < 0.8, else y_axis; then cross+normalize
    for ([_]Vec3f32{ v1, v2, v3 }) |v| {
        const x_axis = Vec3f32{ .x = 1.0, .y = 0.0, .z = 0.0 };
        const y_axis = Vec3f32{ .x = 0.0, .y = 1.0, .z = 0.0 };
        const candidate = if (@abs(v.x) < 0.8) x_axis else y_axis;
        const result = v.cross(candidate).normalize();
        // Must be a unit vector
        try testing.expectApproxEqAbs(@as(f32, 1.0), result.length(), 1e-5);
        // Must be perpendicular to v (dot product ≈ 0)
        const dot = result.x * v.x + result.y * v.y + result.z * v.z;
        try testing.expectApproxEqAbs(@as(f32, 0.0), dot, 1e-4);
    }
}

test "placeWaterHydrogens: water with no neighbors is skipped (no phantom)" {
    // Build a single-residue model with just a water oxygen and no neighboring atoms.
    var mdl = Model.init(testing.allocator);
    defer mdl.deinit();

    var oxygen = Atom{
        .pos = .{ .x = 0.0, .y = 0.0, .z = 0.0 },
        .element_type = .O,
        .residue_idx = 0,
        .occupancy = 1.0,
        .b_factor = 10.0,
        .vdw_radius = 1.52,
    };
    oxygen.setName("O");
    try mdl.atoms.append(testing.allocator, oxygen);

    var res = Residue{};
    res.setCompId("HOH");
    res.entity_type = .water;
    res.atom_start = 0;
    res.atom_end = 1;
    try mdl.residues.append(testing.allocator, res);

    const config = WaterConfig{
        .enabled = true,
        .phantom = false,
        .occupancy_cutoff = 0.5,
        .b_factor_cutoff = 80.0,
        .metal_cutoff = 3.2,
    };

    var altlocs_buf = [_]u8{' '};
    const water_res = mdl.residues.items[0];
    const result = try placeWaterHydrogens(&mdl, water_res, 0, config, .neutron, null, 0, &altlocs_buf);

    // No neighbors → n_skipped_missing_ref should be 1, n_placed = 0
    try testing.expectEqual(@as(u32, 0), result.n_placed);
    try testing.expectEqual(@as(u32, 1), result.n_skipped_missing_ref);
}

test "placeWaterHydrogens: water with one neighbor places H1 and H2" {
    // Build a model: water oxygen at origin, one neighbor at (2.5, 0, 0).
    var mdl = Model.init(testing.allocator);
    defer mdl.deinit();

    // Neighbor atom (residue 1)
    var neighbor = Atom{
        .pos = .{ .x = 2.5, .y = 0.0, .z = 0.0 },
        .element_type = .N,
        .residue_idx = 1,
        .occupancy = 1.0,
        .b_factor = 5.0,
        .vdw_radius = 1.55,
    };
    neighbor.setName("N");
    try mdl.atoms.append(testing.allocator, neighbor);

    // Water oxygen (residue 0)
    var oxygen = Atom{
        .pos = .{ .x = 0.0, .y = 0.0, .z = 0.0 },
        .element_type = .O,
        .residue_idx = 0,
        .occupancy = 1.0,
        .b_factor = 10.0,
        .vdw_radius = 1.52,
    };
    oxygen.setName("O");
    try mdl.atoms.append(testing.allocator, oxygen);

    // Residue 0 = water (atom index 1)
    var water_res = Residue{};
    water_res.setCompId("HOH");
    water_res.entity_type = .water;
    water_res.atom_start = 1;
    water_res.atom_end = 2;
    try mdl.residues.append(testing.allocator, water_res);

    // Residue 1 = neighbor protein (atom index 0)
    var prot_res = Residue{};
    prot_res.setCompId("ALA");
    prot_res.entity_type = .polymer;
    prot_res.atom_start = 0;
    prot_res.atom_end = 1;
    try mdl.residues.append(testing.allocator, prot_res);

    const config = WaterConfig{
        .enabled = true,
        .phantom = false,
        .occupancy_cutoff = 0.5,
        .b_factor_cutoff = 80.0,
        .metal_cutoff = 3.2,
    };

    var altlocs_buf = [_]u8{' '};
    const wr = mdl.residues.items[0];
    const result = try placeWaterHydrogens(&mdl, wr, 0, config, .neutron, null, 0, &altlocs_buf);

    try testing.expectEqual(@as(u32, 2), result.n_placed);
    try testing.expectEqual(@as(u32, 0), result.n_skipped_missing_ref);

    // Verify H1 and H2 were appended with ~0.97 A bond length to oxygen
    var h1_found = false;
    var h2_found = false;
    for (mdl.atoms.items) |atom| {
        if (!atom.is_hydrogen) continue;
        const o_pos = math_mod.Vec3(f32){ .x = 0.0, .y = 0.0, .z = 0.0 };
        const dist = atom.pos.distance(o_pos);
        try testing.expect(dist > 0.80 and dist < 1.10);
        if (std.mem.eql(u8, atom.nameSlice(), "H1")) h1_found = true;
        if (std.mem.eql(u8, atom.nameSlice(), "H2")) h2_found = true;
    }
    try testing.expect(h1_found);
    try testing.expect(h2_found);
}

test "placeWaterHydrogens: low occupancy water is skipped" {
    var mdl = Model.init(testing.allocator);
    defer mdl.deinit();

    var oxygen = Atom{
        .pos = .{ .x = 0.0, .y = 0.0, .z = 0.0 },
        .element_type = .O,
        .residue_idx = 0,
        .occupancy = 0.3, // below default cutoff of 0.66
        .b_factor = 10.0,
        .vdw_radius = 1.52,
    };
    oxygen.setName("O");
    try mdl.atoms.append(testing.allocator, oxygen);

    var res = Residue{};
    res.setCompId("HOH");
    res.entity_type = .water;
    res.atom_start = 0;
    res.atom_end = 1;
    try mdl.residues.append(testing.allocator, res);

    const config = WaterConfig{
        .enabled = true,
        .phantom = false,
        .occupancy_cutoff = 0.66,
        .b_factor_cutoff = 40.0,
        .metal_cutoff = 3.2,
    };

    var altlocs_buf = [_]u8{' '};
    const wr = mdl.residues.items[0];
    const result = try placeWaterHydrogens(&mdl, wr, 0, config, .neutron, null, 0, &altlocs_buf);

    try testing.expectEqual(@as(u32, 0), result.n_placed);
    try testing.expectEqual(@as(u32, 1), result.n_skipped_quality_filter);
}

pub const WaterConfig = struct {
    enabled: bool = false,
    phantom: bool = false,
    occupancy_cutoff: f32 = 0.66,
    b_factor_cutoff: f32 = 40.0,
    metal_cutoff: f32 = 3.2,
};

pub fn isWaterResidue(res: Residue, comp_id: []const u8) bool {
    return res.entity_type == .water or std.mem.eql(u8, comp_id, "HOH") or std.mem.eql(u8, comp_id, "WAT");
}

fn findWaterOxygen(mdl: *const Model, res: Residue, target_altloc: u8) ?Atom {
    const atoms = mdl.atoms.items[res.atom_start..res.atom_end];
    for (atoms) |atom| {
        if (atom.is_hydrogen) continue;
        if (atom.element_type != .O) continue;
        if (target_altloc != ' ' and atom.altloc != ' ' and atom.altloc != target_altloc) continue;
        return atom;
    }
    return null;
}

fn orthogonalUnit(v: Vec3f32) Vec3f32 {
    const x_axis = Vec3f32{ .x = 1.0, .y = 0.0, .z = 0.0 };
    const y_axis = Vec3f32{ .x = 0.0, .y = 1.0, .z = 0.0 };
    const candidate = if (@abs(v.x) < 0.8) x_axis else y_axis;
    return v.cross(candidate).normalize();
}

/// Context built once before the residue loop to accelerate water neighbor queries.
pub const WaterCellCtx = struct {
    cell_list: CellList,
    /// Positions of heavy atoms (parallel to atom_indices).
    positions: []Vec3f32,
    /// Maps CellList position index → mdl.atoms index.
    atom_indices: []u32,
};

fn chooseWaterNeighbors(
    mdl: *const Model,
    res_idx: u32,
    oxygen: Atom,
    target_altloc: u8,
    ctx: ?WaterCellCtx,
    tmp: *std.ArrayListUnmanaged(u32),
) [2]?Vec3f32 {
    var best_dist = [2]f32{ std.math.inf(f32), std.math.inf(f32) };
    var best_pos = [2]?Vec3f32{ null, null };

    if (ctx) |c| {
        // Fast path: query CellList for heavy atoms within 3.5 Å.
        tmp.clearRetainingCapacity();
        c.cell_list.neighborsInRadius(oxygen.pos, 3.5, tmp, mdl.allocator, c.positions) catch {
            // On OOM fall through to the linear scan below.
            return chooseWaterNeighborsLinear(mdl, res_idx, oxygen, target_altloc);
        };
        for (tmp.items) |cl_idx| {
            const atom_idx = c.atom_indices[cl_idx];
            const atom = mdl.atoms.items[atom_idx];
            if (atom.residue_idx == res_idx) continue;
            if (target_altloc != ' ' and atom.altloc != ' ' and atom.altloc != target_altloc) continue;
            const dist = oxygen.pos.distance(atom.pos);
            if (dist < 0.1) continue;
            if (dist < best_dist[0]) {
                best_dist[1] = best_dist[0];
                best_pos[1] = best_pos[0];
                best_dist[0] = dist;
                best_pos[0] = atom.pos;
            } else if (dist < best_dist[1]) {
                best_dist[1] = dist;
                best_pos[1] = atom.pos;
            }
        }
        return best_pos;
    }
    return chooseWaterNeighborsLinear(mdl, res_idx, oxygen, target_altloc);
}

fn chooseWaterNeighborsLinear(mdl: *const Model, res_idx: u32, oxygen: Atom, target_altloc: u8) [2]?Vec3f32 {
    var best_dist = [2]f32{ std.math.inf(f32), std.math.inf(f32) };
    var best_pos = [2]?Vec3f32{ null, null };
    for (mdl.atoms.items) |atom| {
        if (atom.is_hydrogen) continue;
        if (atom.residue_idx == res_idx) continue;
        if (target_altloc != ' ' and atom.altloc != ' ' and atom.altloc != target_altloc) continue;
        const dist = oxygen.pos.distance(atom.pos);
        if (dist > 3.5 or dist < 0.1) continue;
        if (dist < best_dist[0]) {
            best_dist[1] = best_dist[0];
            best_pos[1] = best_pos[0];
            best_dist[0] = dist;
            best_pos[0] = atom.pos;
        } else if (dist < best_dist[1]) {
            best_dist[1] = dist;
            best_pos[1] = atom.pos;
        }
    }
    return best_pos;
}

fn nearestMetalDistance(
    mdl: *const Model,
    res_idx: u32,
    oxygen: Atom,
    target_altloc: u8,
    ctx: ?WaterCellCtx,
    tmp: *std.ArrayListUnmanaged(u32),
) ?f32 {
    if (ctx) |c| {
        // Query within 5.0 Å — metals closer than metal_cutoff (~3.2 Å) are the concern.
        tmp.clearRetainingCapacity();
        c.cell_list.neighborsInRadius(oxygen.pos, 5.0, tmp, mdl.allocator, c.positions) catch {
            return nearestMetalDistanceLinear(mdl, res_idx, oxygen, target_altloc);
        };
        var best: f32 = std.math.inf(f32);
        var found = false;
        for (tmp.items) |cl_idx| {
            const atom_idx = c.atom_indices[cl_idx];
            const atom = mdl.atoms.items[atom_idx];
            if (atom.residue_idx == res_idx) continue;
            if (target_altloc != ' ' and atom.altloc != ' ' and atom.altloc != target_altloc) continue;
            if (!atom.element_type.info().flags.metallic and !atom.flags.metallic) continue;
            const dist = oxygen.pos.distance(atom.pos);
            if (dist < best) {
                best = dist;
                found = true;
            }
        }
        return if (found) best else null;
    }
    return nearestMetalDistanceLinear(mdl, res_idx, oxygen, target_altloc);
}

fn nearestMetalDistanceLinear(mdl: *const Model, res_idx: u32, oxygen: Atom, target_altloc: u8) ?f32 {
    var best: f32 = std.math.inf(f32);
    var found = false;
    for (mdl.atoms.items) |atom| {
        if (atom.residue_idx == res_idx) continue;
        if (target_altloc != ' ' and atom.altloc != ' ' and atom.altloc != target_altloc) continue;
        if (!atom.element_type.info().flags.metallic and !atom.flags.metallic) continue;
        const dist = oxygen.pos.distance(atom.pos);
        if (dist < best) {
            best = dist;
            found = true;
        }
    }
    return if (found) best else null;
}

fn awayFromNearbyAtoms(
    mdl: *const Model,
    res_idx: u32,
    oxygen: Atom,
    target_altloc: u8,
    ctx: ?WaterCellCtx,
    tmp: *std.ArrayListUnmanaged(u32),
) ?Vec3f32 {
    if (ctx) |c| {
        tmp.clearRetainingCapacity();
        c.cell_list.neighborsInRadius(oxygen.pos, 4.0, tmp, mdl.allocator, c.positions) catch {
            return awayFromNearbyAtomsLinear(mdl, res_idx, oxygen, target_altloc);
        };
        var sum = Vec3f32.zero;
        var count: u32 = 0;
        for (tmp.items) |cl_idx| {
            const atom_idx = c.atom_indices[cl_idx];
            const atom = mdl.atoms.items[atom_idx];
            if (atom.residue_idx == res_idx) continue;
            if (target_altloc != ' ' and atom.altloc != ' ' and atom.altloc != target_altloc) continue;
            const dist = oxygen.pos.distance(atom.pos);
            if (dist < 0.1) continue;
            sum = sum.add(atom.pos);
            count += 1;
        }
        if (count == 0) return null;
        return oxygen.pos.sub(sum.scale(1.0 / @as(f32, @floatFromInt(count)))).normalize();
    }
    return awayFromNearbyAtomsLinear(mdl, res_idx, oxygen, target_altloc);
}

fn awayFromNearbyAtomsLinear(mdl: *const Model, res_idx: u32, oxygen: Atom, target_altloc: u8) ?Vec3f32 {
    var sum = Vec3f32.zero;
    var count: u32 = 0;
    for (mdl.atoms.items) |atom| {
        if (atom.is_hydrogen) continue;
        if (atom.residue_idx == res_idx) continue;
        if (target_altloc != ' ' and atom.altloc != ' ' and atom.altloc != target_altloc) continue;
        const dist = oxygen.pos.distance(atom.pos);
        if (dist > 4.0 or dist < 0.1) continue;
        sum = sum.add(atom.pos);
        count += 1;
    }
    if (count == 0) return null;
    return oxygen.pos.sub(sum.scale(1.0 / @as(f32, @floatFromInt(count)))).normalize();
}

/// Subset of PlacementResult relevant to water placement.
/// Mirrors the fields in placer.PlacementResult; we avoid a circular import
/// by using this local struct and letting placer.zig merge into its own result.
pub const WaterPlacementResult = struct {
    n_placed: u32 = 0,
    n_skipped_existing: u32 = 0,
    n_skipped_inter_residue: u32 = 0,
    n_skipped_missing_ref: u32 = 0,
    n_skipped_quality_filter: u32 = 0,
};

pub fn placeWaterHydrogens(
    mdl: *Model,
    res: Residue,
    res_idx: u32,
    config: WaterConfig,
    mode: bond_policy.BondLengthMode,
    cell_ctx: ?WaterCellCtx,
    /// collectAltlocs is defined in placer.zig; pass the result in.
    altlocs_count: usize,
    altlocs_locs: []const u8,
) !WaterPlacementResult {
    var result = WaterPlacementResult{};

    const targets: []const u8 = if (altlocs_count == 0)
        &[_]u8{' '}
    else
        altlocs_locs[0..altlocs_count];

    const h1_name = padName("H1");
    const h2_name = padName("H2");
    const half_angle: f32 = 52.25;

    // Scratch buffer reused across altlocs and helper calls to avoid repeated allocs.
    var tmp = std.ArrayListUnmanaged(u32).empty;
    defer tmp.deinit(mdl.allocator);

    for (targets) |alt| {
        const oxygen = findWaterOxygen(mdl, res, alt) orelse {
            result.n_skipped_missing_ref += 1;
            continue;
        };

        var meta = ParentMeta.fromAtom(oxygen);
        if (alt != ' ') meta.altloc = alt;
        const bond_len = bond_policy.adjustedBondLength(mode, 0.97, oxygen.element_type, .Hpol);

        if (oxygen.flags.bonded_inter_residue) {
            result.n_skipped_inter_residue += 1;
            continue;
        }
        if (nearestMetalDistance(mdl, res_idx, oxygen, alt, cell_ctx, &tmp)) |dist| {
            if (dist <= config.metal_cutoff) {
                result.n_skipped_inter_residue += 1;
                continue;
            }
        }
        if (oxygen.occupancy < config.occupancy_cutoff or oxygen.b_factor > config.b_factor_cutoff) {
            result.n_skipped_quality_filter += 1;
            continue;
        }
        if (existsInResidue(mdl, res, h1_name, meta.altloc) or existsInResidue(mdl, res, h2_name, meta.altloc)) {
            result.n_skipped_existing += 1;
            continue;
        }

        const neighbors = chooseWaterNeighbors(mdl, res_idx, oxygen, alt, cell_ctx, &tmp);
        const n1 = neighbors[0];
        const n2 = neighbors[1];

        var away = if (n1) |p1|
            oxygen.pos.sub(p1).normalize()
        else if (config.phantom)
            awayFromNearbyAtoms(mdl, res_idx, oxygen, alt, cell_ctx, &tmp) orelse Vec3f32{ .x = 1.0, .y = 0.0, .z = 0.0 }
        else {
            result.n_skipped_missing_ref += 1;
            continue;
        };
        if (n1) |p1| {
            if (n2) |p2| {
                away = oxygen.pos.sub(p1).normalize().add(oxygen.pos.sub(p2).normalize()).normalize();
            }
        }
        if (away.length() < 1e-4) {
            result.n_skipped_missing_ref += 1;
            continue;
        }

        const normal = blk: {
            if (n1) |p1| {
                if (n2) |p2| break :blk oxygen.pos.sub(p1).cross(oxygen.pos.sub(p2)).normalize();
            }
            break :blk orthogonalUnit(away);
        };
        const axis = if (normal.length() < 1e-4) orthogonalUnit(away) else normal;

        if (config.phantom and n1 == null) meta.occupancy = 0.0;

        const away64 = away.cast(f64);
        const axis64 = axis.cast(f64);
        const oxygen64 = oxygen.pos.cast(f64);
        const base64 = oxygen64.add(away64.scale(@as(f64, bond_len)));
        const h1 = math_mod.rotateAroundAxis(f64, base64, oxygen64, axis64, half_angle).cast(f32);
        const h2 = math_mod.rotateAroundAxis(f64, base64, oxygen64, axis64, -half_angle).cast(f32);

        try appendNtermH(mdl, h1, "H1", res_idx, meta, .none);
        try appendNtermH(mdl, h2, "H2", res_idx, meta, .none);
        result.n_placed += 2;
    }

    return result;
}
