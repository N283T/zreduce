//! Mover generator: scans placed H atoms and creates Mover instances.
//!
//! Inspects `mover_hint` on each added hydrogen atom, resolves axis/center
//! atoms from placement plans, and delegates to rotator constructors.
//! Supports single-H rotators, methyl rotators, NH3 rotators,
//! amide flippers (ASN/GLN), and histidine ring flippers.

const std = @import("std");
const log = std.log.scoped(.mover_gen);
const model_mod = @import("../model.zig");
const Model = model_mod.Model;
const Atom = model_mod.Atom;
const Residue = model_mod.Residue;
const mover_mod = @import("mover.zig");
const Mover = mover_mod.Mover;
const rotator = @import("rotator.zig");
const flipper = @import("flipper.zig");
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
/// Prefers an exact altloc match and falls back to blank altloc for shared atoms.
fn findAtomIdx(atoms: []const Atom, residue_idx: u32, name: []const u8, target_altloc: u8) ?u32 {
    var blank_match: ?u32 = null;
    for (atoms, 0..) |*a, i| {
        if (a.residue_idx != residue_idx) continue;
        if (!std.mem.eql(u8, a.nameSlice(), name)) continue;
        if (a.altloc == target_altloc) return @intCast(i);
        if (a.altloc == ' ' and blank_match == null) blank_match = @intCast(i);
    }
    return blank_match;
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
    var seen_groups: std.AutoHashMapUnmanaged(GroupKey, void) = .empty;
    defer seen_groups.deinit(allocator);

    const atoms = mdl.atoms.items;

    for (atoms, 0..) |*atom, idx| {
        if (!atom.is_added) continue;
        if (atom.mover_hint == .none) continue;

        const residue_idx = atom.residue_idx;
        const res = mdl.residues.items[residue_idx];
        const comp_id = res.compIdSlice();
        const h_name = atom.nameSlice();
        const target_altloc = atom.altloc;

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
                const center_idx = findAtomIdx(atoms, residue_idx, center_name, target_altloc) orelse {
                    log.warn("center atom '{s}' not found for H '{s}' in {s} (res {d}), skipping", .{ center_name, h_name, comp_id, residue_idx });
                    n_skipped += 1;
                    continue;
                };
                const axis_idx = findAtomIdx(atoms, residue_idx, axis_name, target_altloc) orelse {
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
                const group_key = GroupKey{ .residue_idx = residue_idx, .center_name = center_name_raw, .altloc = target_altloc };
                const gop = try seen_groups.getOrPut(allocator, group_key);
                if (gop.found_existing) continue;

                const center_name = trimName(&center_name_raw);
                const axis_name = trimName(&plan.connected[1]);
                const center_idx = findAtomIdx(atoms, residue_idx, center_name, target_altloc) orelse {
                    log.warn("center atom '{s}' not found for group in {s} (res {d}), skipping", .{ center_name, comp_id, residue_idx });
                    n_skipped += 1;
                    continue;
                };
                const axis_idx = findAtomIdx(atoms, residue_idx, axis_name, target_altloc) orelse {
                    log.warn("axis atom '{s}' not found for group in {s} (res {d}), skipping", .{ axis_name, comp_id, residue_idx });
                    n_skipped += 1;
                    continue;
                };

                // Find all 3 H atoms with same (residue_idx, center, hint).
                var h_indices: [3]u32 = undefined;
                var h_count: u32 = 0;
                for (atoms, 0..) |*a, ai| {
                    if (a.residue_idx == residue_idx and a.is_added and a.mover_hint == atom.mover_hint and a.altloc == target_altloc) {
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
            .flip_amide => {
                if (no_flip) continue;

                // Deduplicate: one amide flipper per residue
                const group_key = GroupKey{ .residue_idx = residue_idx, .center_name = .{ 'A', 'M', 'D', ' ' }, .altloc = target_altloc };
                const gop = try seen_groups.getOrPut(allocator, group_key);
                if (gop.found_existing) continue;

                const plan = findPlanForH(comp_id, h_name) orelse {
                    log.warn("no plan for flip H '{s}' in {s} (res {d}), skipping", .{ h_name, comp_id, residue_idx });
                    n_skipped += 1;
                    continue;
                };

                // N atom = connected[0] (ND2 for ASN, NE2 for GLN)
                const n_name = trimName(&plan.connected[0]);
                // C atom = connected[1] (CG for ASN, CD for GLN)
                const c_name = trimName(&plan.connected[1]);
                // O atom: OD1 for ASN, OE1 for GLN
                const o_name: []const u8 = if (std.mem.eql(u8, comp_id, "ASN")) "OD1" else "OE1";

                const n_at_idx = findAtomIdx(atoms, residue_idx, n_name, target_altloc) orelse {
                    log.warn("N atom '{s}' not found for amide flip in {s} (res {d})", .{ n_name, comp_id, residue_idx });
                    n_skipped += 1;
                    continue;
                };
                const o_at_idx = findAtomIdx(atoms, residue_idx, o_name, target_altloc) orelse {
                    log.warn("O atom '{s}' not found for amide flip in {s} (res {d})", .{ o_name, comp_id, residue_idx });
                    n_skipped += 1;
                    continue;
                };
                const c_at_idx = findAtomIdx(atoms, residue_idx, c_name, target_altloc) orelse {
                    log.warn("C atom '{s}' not found for amide flip in {s} (res {d})", .{ c_name, comp_id, residue_idx });
                    n_skipped += 1;
                    continue;
                };

                // Find the two H atoms for this amide group
                var h_indices: [2]u32 = undefined;
                var h_count: u32 = 0;
                const all_plans = standard.getPlans(comp_id) orelse continue;
                for (all_plans) |*p| {
                    if (p.mover_hint != .flip_amide) continue;
                    const p_h_name = trimName(&p.h_name);
                    if (findAtomIdx(atoms, residue_idx, p_h_name, target_altloc)) |hi| {
                        if (h_count < 2) {
                            h_indices[h_count] = hi;
                            h_count += 1;
                        }
                    }
                }
                if (h_count != 2) {
                    log.warn("expected 2 H for amide flip in {s} (res {d}), found {d}", .{ comp_id, residue_idx, h_count });
                    n_skipped += 1;
                    continue;
                }

                const m = try flipper.createAmideFlipper(
                    allocator,
                    atoms,
                    o_at_idx,
                    n_at_idx,
                    h_indices[0],
                    h_indices[1],
                    c_at_idx,
                    residue_idx,
                );
                try movers.append(allocator, m);
            },
            .flip_his => {
                if (no_flip) continue;

                // Deduplicate: one His flipper per residue
                const group_key = GroupKey{ .residue_idx = residue_idx, .center_name = .{ 'H', 'I', 'S', ' ' }, .altloc = target_altloc };
                const gop = try seen_groups.getOrPut(allocator, group_key);
                if (gop.found_existing) continue;

                // Find ring heavy atoms by name
                const nd1_idx = findAtomIdx(atoms, residue_idx, "ND1", target_altloc) orelse {
                    log.warn("ND1 not found for His flip (res {d})", .{residue_idx});
                    n_skipped += 1;
                    continue;
                };
                const cd2_idx = findAtomIdx(atoms, residue_idx, "CD2", target_altloc) orelse {
                    log.warn("CD2 not found for His flip (res {d})", .{residue_idx});
                    n_skipped += 1;
                    continue;
                };
                const ce1_idx = findAtomIdx(atoms, residue_idx, "CE1", target_altloc) orelse {
                    log.warn("CE1 not found for His flip (res {d})", .{residue_idx});
                    n_skipped += 1;
                    continue;
                };
                const ne2_idx = findAtomIdx(atoms, residue_idx, "NE2", target_altloc) orelse {
                    log.warn("NE2 not found for His flip (res {d})", .{residue_idx});
                    n_skipped += 1;
                    continue;
                };

                // H atoms are optional
                const hd1_idx = findAtomIdx(atoms, residue_idx, "HD1", target_altloc);
                const he2_idx = findAtomIdx(atoms, residue_idx, "HE2", target_altloc);

                const m = try flipper.createHisFlipper(
                    allocator,
                    atoms,
                    nd1_idx,
                    cd2_idx,
                    ce1_idx,
                    ne2_idx,
                    hd1_idx,
                    he2_idx,
                    residue_idx,
                );
                try movers.append(allocator, m);
            },
            .none => unreachable,
        }
    }

    return .{
        .movers = try movers.toOwnedSlice(allocator),
        .n_skipped = n_skipped,
    };
}

const GroupKey = struct {
    residue_idx: u32,
    center_name: [4]u8,
    altloc: u8,
};

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

test "findAtomIdx prefers exact altloc and falls back to blank" {
    const source = @embedFile("../test_data/ala_altloc.cif");
    var mdl = try mmcif.parseModel(testing.allocator, source);
    defer mdl.deinit();

    try testing.expectEqual(findAtomIdx(mdl.atoms.items, 0, "CA", 'A').?, 1);
    try testing.expectEqual(findAtomIdx(mdl.atoms.items, 0, "CA", 'B').?, 2);
    try testing.expectEqual(findAtomIdx(mdl.atoms.items, 0, "N", 'A').?, 0);
    try testing.expectEqual(findAtomIdx(mdl.atoms.items, 0, "N", 'B').?, 0);
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

test "generateMovers creates separate methyl rotators per altloc conformer" {
    const source = @embedFile("../test_data/ala_altloc.cif");
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

    var methyl_count: u32 = 0;
    var center_a: u32 = 0;
    var center_b: u32 = 0;
    for (movers) |m| {
        if (m.kind != .methyl_rotator) continue;
        methyl_count += 1;
        const center_idx = m.center_idx.?;
        if (mdl.atoms.items[center_idx].altloc == 'A') center_a += 1;
        if (mdl.atoms.items[center_idx].altloc == 'B') center_b += 1;
    }

    try testing.expectEqual(@as(u32, 2), methyl_count);
    try testing.expectEqual(@as(u32, 1), center_a);
    try testing.expectEqual(@as(u32, 1), center_b);
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

test "generateMovers creates amide flipper for ASN" {
    const source = @embedFile("../test_data/asn.cif");
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

    var amide_count: u32 = 0;
    for (movers) |m| {
        if (m.kind == .amide_flip) amide_count += 1;
    }
    try testing.expectEqual(@as(u32, 1), amide_count);
}

test "no_flip suppresses flip movers" {
    const source = @embedFile("../test_data/asn.cif");
    var mdl = try mmcif.parseModel(testing.allocator, source);
    defer mdl.deinit();

    placer.applyChemistry(&mdl);
    _ = try placer.addHydrogens(&mdl, null);

    const gen_result = try generateMovers(testing.allocator, &mdl, true);
    const movers = gen_result.movers;
    defer {
        for (0..movers.len) |i| @constCast(&movers[i]).deinit();
        testing.allocator.free(movers);
    }

    for (movers) |m| {
        try testing.expect(m.kind != .amide_flip);
        try testing.expect(m.kind != .his_flip);
    }
}

test "generateMovers creates His flipper" {
    const source = @embedFile("../test_data/his.cif");
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

    var his_count: u32 = 0;
    for (movers) |m| {
        if (m.kind == .his_flip) {
            his_count += 1;
            try testing.expectEqual(@as(usize, 6), m.orientations.len);
        }
    }
    try testing.expectEqual(@as(u32, 1), his_count);
}
