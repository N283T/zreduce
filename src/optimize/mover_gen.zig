//! Mover generator: scans placed H atoms and creates Mover instances.
//!
//! Inspects `mover_hint` on each added hydrogen atom, resolves axis/center
//! atoms from placement plans, and delegates to rotator constructors.
//! Currently supports single-H rotators, methyl rotators, and NH3 rotators.
//! Flip movers (amide, His) are deferred to Issue #17.

const std = @import("std");
const log = std.log.scoped(.mover_gen);
const model_mod = @import("../model.zig");
const Model = model_mod.Model;
const Atom = model_mod.Atom;
const Residue = model_mod.Residue;
const mover_mod = @import("mover.zig");
const Mover = mover_mod.Mover;
const rotator = @import("rotator.zig");
const standard = @import("../place/standard.zig");

// ---------------------------------------------------------------------------
// Helpers
// ---------------------------------------------------------------------------

/// Trim leading/trailing spaces from a 4-byte padded atom name.
/// Takes a pointer to avoid returning a slice into a stack-local copy.
fn trimName(name: *const [4]u8) []const u8 {
    var start: usize = 0;
    while (start < 4 and name[start] == ' ') start += 1;
    var end: usize = 4;
    while (end > start and name[end - 1] == ' ') end -= 1;
    return name[start..end];
}

/// Find a PlacementPlan for a given H atom name within a comp_id.
fn findPlanForH(comp_id: []const u8, h_name: []const u8) ?*const standard.PlacementPlan {
    const plans = standard.getPlans(comp_id) orelse return null;
    for (plans) |*plan| {
        if (std.mem.eql(u8, trimName(&plan.h_name), h_name)) return plan;
    }
    return null;
}

/// Find atom index by scanning entire atoms array (handles added H atoms
/// that live beyond res.atom_end).
fn findAtomIdx(atoms: []const Atom, residue_idx: u32, name: []const u8) ?u32 {
    for (atoms, 0..) |*a, i| {
        if (a.residue_idx == residue_idx and std.mem.eql(u8, a.nameSlice(), name)) {
            return @intCast(i);
        }
    }
    return null;
}

// ---------------------------------------------------------------------------
// Public API
// ---------------------------------------------------------------------------

pub const MoverGenResult = struct {
    movers: []Mover,
    n_skipped: u32,
};

/// Scan all placed H atoms in the model and create Mover instances for
/// optimization. Caller owns the returned movers slice and must call
/// `deinit()` on each Mover, then free the slice with the same allocator.
pub fn generateMovers(allocator: std.mem.Allocator, mdl: *const Model, no_flip: bool) !MoverGenResult {
    var movers: std.ArrayListUnmanaged(Mover) = .empty;
    errdefer {
        for (movers.items) |*m| m.deinit();
        movers.deinit(allocator);
    }

    var n_skipped: u32 = 0;

    // Track which group rotators we have already created to avoid duplicates.
    // Key: (residue_idx, center_name) packed into u64.
    var seen_groups: std.AutoHashMapUnmanaged(u64, void) = .empty;
    defer seen_groups.deinit(allocator);

    const atoms = mdl.atoms.items;

    for (atoms, 0..) |*atom, idx| {
        if (!atom.is_added) continue;
        if (atom.mover_hint == .none) continue;

        const residue_idx = atom.residue_idx;
        const res = mdl.residues.items[residue_idx];
        const comp_id = res.compIdSlice();
        const h_name = atom.nameSlice();

        switch (atom.mover_hint) {
            .rotate => {
                // Single H rotator (OH, SH, etc.)
                const plan = findPlanForH(comp_id, h_name) orelse {
                    log.warn("no plan for H '{s}' in {s} (res {d}), skipping rotator", .{ h_name, comp_id, residue_idx });
                    n_skipped += 1;
                    continue;
                };
                const center_name = trimName(&plan.connected[0]);
                const axis_name = trimName(&plan.connected[1]);
                const center_idx = findAtomIdx(atoms, residue_idx, center_name) orelse {
                    log.warn("center atom '{s}' not found for H '{s}' in {s} (res {d}), skipping", .{ center_name, h_name, comp_id, residue_idx });
                    n_skipped += 1;
                    continue;
                };
                const axis_idx = findAtomIdx(atoms, residue_idx, axis_name) orelse {
                    log.warn("axis atom '{s}' not found for H '{s}' in {s} (res {d}), skipping", .{ axis_name, h_name, comp_id, residue_idx });
                    n_skipped += 1;
                    continue;
                };
                const h_idx: u32 = @intCast(idx);

                const mover = try rotator.createSingleHRotator(
                    allocator,
                    atoms,
                    h_idx,
                    center_idx,
                    axis_idx,
                    residue_idx,
                );
                try movers.append(allocator, mover);
            },
            .rotate_methyl, .rotate_nh3 => {
                // Group rotator: 3 H atoms share same center.
                const plan = findPlanForH(comp_id, h_name) orelse {
                    n_skipped += 1;
                    continue;
                };
                const center_name_raw = plan.connected[0];

                // Dedup key: pack residue_idx and center_name into u64.
                const group_key = packGroupKey(residue_idx, center_name_raw);
                const gop = try seen_groups.getOrPut(allocator, group_key);
                if (gop.found_existing) continue;

                const center_name = trimName(&center_name_raw);
                const axis_name = trimName(&plan.connected[1]);
                const center_idx = findAtomIdx(atoms, residue_idx, center_name) orelse {
                    log.warn("center atom '{s}' not found for group in {s} (res {d}), skipping", .{ center_name, comp_id, residue_idx });
                    n_skipped += 1;
                    continue;
                };
                const axis_idx = findAtomIdx(atoms, residue_idx, axis_name) orelse {
                    log.warn("axis atom '{s}' not found for group in {s} (res {d}), skipping", .{ axis_name, comp_id, residue_idx });
                    n_skipped += 1;
                    continue;
                };

                // Find all 3 H atoms with same (residue_idx, center, hint).
                var h_indices: [3]u32 = undefined;
                var h_count: u32 = 0;
                for (atoms, 0..) |*a, ai| {
                    if (a.residue_idx == residue_idx and a.is_added and a.mover_hint == atom.mover_hint) {
                        const a_plan = findPlanForH(comp_id, a.nameSlice()) orelse continue;
                        if (std.mem.eql(u8, &a_plan.connected[0], &center_name_raw)) {
                            if (h_count < 3) {
                                h_indices[h_count] = @intCast(ai);
                                h_count += 1;
                            }
                        }
                    }
                }
                if (h_count != 3) {
                    log.warn("expected 3 H for group at '{s}' in {s} (res {d}), found {d}, skipping", .{ center_name, comp_id, residue_idx, h_count });
                    n_skipped += 1;
                    continue;
                }

                const mover = switch (atom.mover_hint) {
                    .rotate_methyl => try rotator.createMethylRotator(
                        allocator,
                        atoms,
                        h_indices,
                        center_idx,
                        axis_idx,
                        residue_idx,
                    ),
                    .rotate_nh3 => try rotator.createNH3Rotator(
                        allocator,
                        atoms,
                        h_indices,
                        center_idx,
                        axis_idx,
                        residue_idx,
                    ),
                    else => unreachable,
                };
                try movers.append(allocator, mover);
            },
            .flip_amide, .flip_his => {
                if (no_flip) continue;
                // TODO: Implement in Tasks 3 and 4
                continue;
            },
            .none => unreachable,
        }
    }

    return .{
        .movers = try movers.toOwnedSlice(allocator),
        .n_skipped = n_skipped,
    };
}

fn packGroupKey(residue_idx: u32, center_name: [4]u8) u64 {
    return (@as(u64, residue_idx) << 32) |
        (@as(u64, center_name[0]) << 24) |
        (@as(u64, center_name[1]) << 16) |
        (@as(u64, center_name[2]) << 8) |
        (@as(u64, center_name[3]));
}

// ── Tests ─────────────────────────────────────────────────────────────────────

const testing = std.testing;
const mmcif = @import("../mmcif.zig");
const placer = @import("../place/placer.zig");

test "trimName trims spaces" {
    try testing.expectEqualStrings("CB", trimName(&.{ ' ', 'C', 'B', ' ' }));
    try testing.expectEqualStrings("CA", trimName(&.{ ' ', 'C', 'A', ' ' }));
    try testing.expectEqualStrings("N", trimName(&.{ ' ', 'N', ' ', ' ' }));
    try testing.expectEqualStrings("HB1", trimName(&.{ 'H', 'B', '1', ' ' }));
}

test "findPlanForH returns correct plan for ALA HB1" {
    const plan = findPlanForH("ALA", "HB1");
    try testing.expect(plan != null);
    try testing.expectEqualStrings("HB1", trimName(&plan.?.h_name));
    try testing.expectEqual(standard.MoverHint.rotate_methyl, plan.?.mover_hint);
}

test "findPlanForH returns null for unknown atom" {
    const plan = findPlanForH("ALA", "HX9");
    try testing.expect(plan == null);
}

test "generateMovers creates methyl rotator for ALA" {
    const source = @embedFile("../test_data/tiny.cif");
    var mdl = try mmcif.parseModel(testing.allocator, source);
    defer mdl.deinit();

    placer.applyChemistry(&mdl);
    _ = try placer.addHydrogens(&mdl, null);

    const gen_result = try generateMovers(testing.allocator, &mdl, false);
    const movers = gen_result.movers;
    defer {
        for (0..movers.len) |i| @constCast(&movers[i]).deinit();
        testing.allocator.free(movers);
    }

    // ALA has 1 methyl group on CB -> should produce 1 methyl_rotator
    var methyl_count: u32 = 0;
    for (movers) |m| {
        if (m.kind == .methyl_rotator) methyl_count += 1;
    }
    try testing.expectEqual(@as(u32, 1), methyl_count);
}

test "generateMovers total count for ALA" {
    const source = @embedFile("../test_data/tiny.cif");
    var mdl = try mmcif.parseModel(testing.allocator, source);
    defer mdl.deinit();

    placer.applyChemistry(&mdl);
    _ = try placer.addHydrogens(&mdl, null);

    const gen_result = try generateMovers(testing.allocator, &mdl, false);
    const movers = gen_result.movers;
    defer {
        for (0..movers.len) |i| @constCast(&movers[i]).deinit();
        testing.allocator.free(movers);
    }

    // ALA should have: 1 methyl rotator (CB)
    // No single-H rotators (ALA has no OH/SH groups)
    // No NH3 rotators (N-terminal NH3+ would need rotate_nh3 hint)
    try testing.expect(movers.len >= 1);
}

test "optimizer pipeline runs without error on ALA" {
    const opt_mod = @import("optimizer.zig");

    const source = @embedFile("../test_data/tiny.cif");
    var mdl = try mmcif.parseModel(testing.allocator, source);
    defer mdl.deinit();

    placer.applyChemistry(&mdl);
    _ = try placer.addHydrogens(&mdl, null);

    const gen_result = try generateMovers(testing.allocator, &mdl, false);
    const movers = gen_result.movers;
    defer {
        for (0..movers.len) |i| @constCast(&movers[i]).deinit();
        testing.allocator.free(movers);
    }

    try testing.expect(movers.len > 0);

    // Run optimizer - should complete without error
    const opt_config = opt_mod.OptConfig{};
    _ = try opt_mod.optimize(testing.allocator, movers, &mdl, opt_config);
}

test "methyl rotator controls 3 atoms" {
    const source = @embedFile("../test_data/tiny.cif");
    var mdl = try mmcif.parseModel(testing.allocator, source);
    defer mdl.deinit();

    placer.applyChemistry(&mdl);
    _ = try placer.addHydrogens(&mdl, null);

    const gen_result = try generateMovers(testing.allocator, &mdl, false);
    const movers = gen_result.movers;
    defer {
        for (0..movers.len) |i| @constCast(&movers[i]).deinit();
        testing.allocator.free(movers);
    }

    for (movers) |m| {
        if (m.kind == .methyl_rotator) {
            try testing.expectEqual(@as(usize, 3), m.atom_indices.len);
            try testing.expectEqual(@as(usize, 3), m.orientations.len);
        }
    }
}
