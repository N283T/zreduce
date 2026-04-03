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
const nucleotide = @import("../place/nucleotide.zig");
const modified = @import("../place/modified.zig");
const protonation = @import("../place/protonation.zig");
const bond_policy = @import("../place/bond_policy.zig");
const ccd_mod = @import("../ccd.zig");
const ComponentDict = ccd_mod.ComponentDict;

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
    const plans = standard.getPlans(comp_id) orelse
        nucleotide.getPlans(comp_id) orelse
        modified.getPlans(comp_id) orelse
        return null;
    for (plans) |*plan| {
        if (std.mem.eql(u8, trimName(&plan.h_name), h_name)) return plan;
    }
    return null;
}

/// Find atom index for a given atom name and altloc within a residue.
///
/// Performs two bounded scans instead of a full O(N) scan over all atoms:
///   - Pass 1: atoms[res.atom_start..res.atom_end] for original atoms (no
///             residue_idx check needed since the range is per-residue).
///   - Pass 2: atoms[original_atom_count..] for added H atoms, filtered by
///             residue_idx (added atoms live after all original atoms and carry
///             a residue_idx back-pointer).
///
/// Prefers an exact altloc match and falls back to blank altloc for shared atoms.
fn findAtomIdx(
    atoms: []const Atom,
    res: Residue,
    residue_idx: u32,
    name: []const u8,
    target_altloc: u8,
    original_atom_count: u32,
) ?u32 {
    var blank_match: ?u32 = null;

    // Pass 1: original atoms for this residue (tight range, no residue_idx filter needed)
    for (atoms[res.atom_start..res.atom_end], res.atom_start..) |*a, i| {
        if (!std.mem.eql(u8, a.nameSlice(), name)) continue;
        if (a.altloc == target_altloc) return @intCast(i);
        if (a.altloc == ' ' and blank_match == null) blank_match = @intCast(i);
    }

    // Pass 2: added H atoms appended after all original atoms
    for (atoms[original_atom_count..], original_atom_count..) |*a, i| {
        if (a.residue_idx != residue_idx) continue;
        if (!std.mem.eql(u8, a.nameSlice(), name)) continue;
        if (a.altloc == target_altloc) return @intCast(i);
        if (a.altloc == ' ' and blank_match == null) blank_match = @intCast(i);
    }

    return blank_match;
}

/// Resolve center and axis atom names for a rotatable H using CCD topology.
/// For a single-H rotator: center = heavy atom bonded to H, axis = heavy atom bonded to center.
fn resolveCcdRotatorAtoms(
    ccd_dict: ?*const ComponentDict,
    comp_id: []const u8,
    h_name: []const u8,
    center_out: *[4]u8,
    axis_out: *[4]u8,
) bool {
    const dict = ccd_dict orelse return false;
    const component = dict.get(comp_id) orelse return false;

    // Find H atom index in CCD component
    var h_ccd_idx: ?u16 = null;
    for (component.atoms, 0..) |a, i| {
        const a_name = trimName(&a.name);
        if (std.mem.eql(u8, a_name, h_name)) {
            h_ccd_idx = @intCast(i);
            break;
        }
    }
    const h_idx = h_ccd_idx orelse return false;

    // Find the heavy atom bonded to this H (center)
    var center_idx: ?u16 = null;
    for (component.bonds) |bond| {
        const other: u16 = if (bond.atom_idx_1 == h_idx)
            bond.atom_idx_2
        else if (bond.atom_idx_2 == h_idx)
            bond.atom_idx_1
        else
            continue;
        if (component.atoms[other].element_symbol[0] != 'H') {
            center_idx = other;
            break;
        }
    }
    const ci = center_idx orelse return false;
    center_out.* = component.atoms[ci].name;

    // Find the first heavy atom bonded to center (axis), excluding H atoms
    for (component.bonds) |bond| {
        const other: u16 = if (bond.atom_idx_1 == ci)
            bond.atom_idx_2
        else if (bond.atom_idx_2 == ci)
            bond.atom_idx_1
        else
            continue;
        if (other == h_idx) continue;
        if (component.atoms[other].element_symbol[0] != 'H') {
            axis_out.* = component.atoms[other].name;
            return true;
        }
    }
    return false;
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
/// inline_dict is checked first per component, then ccd_dict as fallback.
pub fn generateMovers(
    allocator: std.mem.Allocator,
    mdl: *const Model,
    no_flip: bool,
    ccd_dict: ?*const ComponentDict,
    inline_dict: ?*const ComponentDict,
    overrides: ?*const protonation.ProtonationOverrides,
    mode: bond_policy.BondLengthMode,
) !MoverGenResult {
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
        const residue_state = if (overrides) |ov| ov.find(mdl, residue_idx) else null;

        // Per-component dict fallback: inline_dict first, then ccd_dict.
        const effective_dict: ?*const ComponentDict = blk: {
            if (inline_dict) |d| {
                if (d.get(comp_id) != null) break :blk inline_dict;
            }
            break :blk ccd_dict;
        };

        switch (atom.mover_hint) {
            .rotate => {
                // Single H rotator (OH, SH, etc.)
                // Try standard plan first, fall back to CCD topology
                var center_name_buf: [4]u8 = undefined;
                var axis_name_buf: [4]u8 = undefined;
                var center_name: []const u8 = undefined;
                var axis_name: []const u8 = undefined;

                if (findPlanForH(comp_id, h_name)) |plan| {
                    center_name_buf = plan.connected[0];
                    axis_name_buf = plan.connected[1];
                    center_name = trimName(&center_name_buf);
                    axis_name = trimName(&axis_name_buf);
                } else if (resolveCcdRotatorAtoms(effective_dict, comp_id, h_name, &center_name_buf, &axis_name_buf)) {
                    center_name = trimName(&center_name_buf);
                    axis_name = trimName(&axis_name_buf);
                } else {
                    log.warn("no plan for H '{s}' in {s} (res {d}), skipping rotator", .{ h_name, comp_id, residue_idx });
                    n_skipped += 1;
                    continue;
                }
                const center_idx = findAtomIdx(atoms, res, residue_idx, center_name, target_altloc, mdl.original_atom_count) orelse {
                    log.warn("center atom '{s}' not found for H '{s}' in {s} (res {d}), skipping", .{ center_name, h_name, comp_id, residue_idx });
                    n_skipped += 1;
                    continue;
                };
                const axis_idx = findAtomIdx(atoms, res, residue_idx, axis_name, target_altloc, mdl.original_atom_count) orelse {
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
                if (residue_state) |state| {
                    if (state == .lys and state.lys == .neutral and atom.mover_hint == .rotate_nh3) continue;
                }
                // Group rotator: 3 H atoms share same center.
                // Try standard plan first, fall back to CCD topology.
                var center_name_raw: [4]u8 = undefined;
                var axis_name_buf: [4]u8 = undefined;
                var center_name: []const u8 = undefined;
                var axis_name: []const u8 = undefined;

                if (findPlanForH(comp_id, h_name)) |plan| {
                    center_name_raw = plan.connected[0];
                    axis_name_buf = plan.connected[1];
                    center_name = trimName(&center_name_raw);
                    axis_name = trimName(&axis_name_buf);
                } else if (resolveCcdRotatorAtoms(effective_dict, comp_id, h_name, &center_name_raw, &axis_name_buf)) {
                    center_name = trimName(&center_name_raw);
                    axis_name = trimName(&axis_name_buf);
                } else {
                    n_skipped += 1;
                    continue;
                }

                // Dedup key
                const group_key = GroupKey{ .residue_idx = residue_idx, .center_name = center_name_raw, .altloc = target_altloc };
                const gop = try seen_groups.getOrPut(allocator, group_key);
                if (gop.found_existing) continue;

                const center_idx = findAtomIdx(atoms, res, residue_idx, center_name, target_altloc, mdl.original_atom_count) orelse {
                    log.warn("center atom '{s}' not found for group in {s} (res {d}), skipping", .{ center_name, comp_id, residue_idx });
                    n_skipped += 1;
                    continue;
                };
                const axis_idx = findAtomIdx(atoms, res, residue_idx, axis_name, target_altloc, mdl.original_atom_count) orelse {
                    log.warn("axis atom '{s}' not found for group in {s} (res {d}), skipping", .{ axis_name, comp_id, residue_idx });
                    n_skipped += 1;
                    continue;
                };

                // Find all 3 H atoms with same (residue_idx, center, hint).
                // Use both standard plan lookup and CCD fallback to match center.
                // Scan only added atoms (atoms[original_atom_count..]) for efficiency.
                var h_indices: [3]u32 = undefined;
                var h_count: u32 = 0;
                for (atoms[mdl.original_atom_count..], mdl.original_atom_count..) |*a, ai| {
                    if (a.residue_idx == residue_idx and a.is_added and a.mover_hint == atom.mover_hint and a.altloc == target_altloc) {
                        // Resolve this H's center via standard plan or CCD
                        var a_center: [4]u8 = undefined;
                        var a_axis: [4]u8 = undefined;
                        const has_center = if (findPlanForH(comp_id, a.nameSlice())) |ap| blk: {
                            a_center = ap.connected[0];
                            break :blk true;
                        } else resolveCcdRotatorAtoms(effective_dict, comp_id, a.nameSlice(), &a_center, &a_axis);

                        if (has_center and std.mem.eql(u8, &a_center, &center_name_raw)) {
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

                const n_at_idx = findAtomIdx(atoms, res, residue_idx, n_name, target_altloc, mdl.original_atom_count) orelse {
                    log.warn("N atom '{s}' not found for amide flip in {s} (res {d})", .{ n_name, comp_id, residue_idx });
                    n_skipped += 1;
                    continue;
                };
                const o_at_idx = findAtomIdx(atoms, res, residue_idx, o_name, target_altloc, mdl.original_atom_count) orelse {
                    log.warn("O atom '{s}' not found for amide flip in {s} (res {d})", .{ o_name, comp_id, residue_idx });
                    n_skipped += 1;
                    continue;
                };
                const c_at_idx = findAtomIdx(atoms, res, residue_idx, c_name, target_altloc, mdl.original_atom_count) orelse {
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
                    if (findAtomIdx(atoms, res, residue_idx, p_h_name, target_altloc, mdl.original_atom_count)) |hi| {
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
                    mode,
                );
                try movers.append(allocator, m);
            },
            .flip_his => {
                if (residue_state) |state| {
                    if (state == .his and state.his != .auto) continue;
                }
                if (no_flip) continue;

                // Deduplicate: one His flipper per residue
                const group_key = GroupKey{ .residue_idx = residue_idx, .center_name = .{ 'H', 'I', 'S', ' ' }, .altloc = target_altloc };
                const gop = try seen_groups.getOrPut(allocator, group_key);
                if (gop.found_existing) continue;

                // Find ring heavy atoms by name
                const nd1_idx = findAtomIdx(atoms, res, residue_idx, "ND1", target_altloc, mdl.original_atom_count) orelse {
                    log.warn("ND1 not found for His flip (res {d})", .{residue_idx});
                    n_skipped += 1;
                    continue;
                };
                const cd2_idx = findAtomIdx(atoms, res, residue_idx, "CD2", target_altloc, mdl.original_atom_count) orelse {
                    log.warn("CD2 not found for His flip (res {d})", .{residue_idx});
                    n_skipped += 1;
                    continue;
                };
                const ce1_idx = findAtomIdx(atoms, res, residue_idx, "CE1", target_altloc, mdl.original_atom_count) orelse {
                    log.warn("CE1 not found for His flip (res {d})", .{residue_idx});
                    n_skipped += 1;
                    continue;
                };
                const ne2_idx = findAtomIdx(atoms, res, residue_idx, "NE2", target_altloc, mdl.original_atom_count) orelse {
                    log.warn("NE2 not found for His flip (res {d})", .{residue_idx});
                    n_skipped += 1;
                    continue;
                };

                // H atoms are optional
                const hd1_idx = findAtomIdx(atoms, res, residue_idx, "HD1", target_altloc, mdl.original_atom_count);
                const he2_idx = findAtomIdx(atoms, res, residue_idx, "HE2", target_altloc, mdl.original_atom_count);

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
                    mode,
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

test "findPlanForH finds nucleotide plans" {
    // Aromatic C-H on adenine
    const h8 = findPlanForH("DA", "H8");
    try testing.expect(h8 != null);
    try testing.expectEqualStrings("H8", trimName(&h8.?.h_name));

    // OH rotator on RNA ribose
    const ho2 = findPlanForH("A", "HO2'");
    try testing.expect(ho2 != null);
    try testing.expectEqual(standard.MoverHint.rotate, ho2.?.mover_hint);

    // Thymine methyl
    const h71 = findPlanForH("DT", "H71");
    try testing.expect(h71 != null);
    try testing.expectEqual(standard.MoverHint.rotate_methyl, h71.?.mover_hint);

    // Standard AA still works
    try testing.expect(findPlanForH("ALA", "HB1") != null);

    // Unknown returns null
    try testing.expect(findPlanForH("DA", "HX9") == null);
}

test "findAtomIdx prefers exact altloc and falls back to blank" {
    const source = @embedFile("../test_data/ala_altloc.cif");
    var mdl = try mmcif.parseModel(testing.allocator, source);
    defer mdl.deinit();

    const res0 = mdl.residues.items[0];
    const orig = mdl.original_atom_count;

    try testing.expectEqual(findAtomIdx(mdl.atoms.items, res0, 0, "CA", 'A', orig).?, 1);
    try testing.expectEqual(findAtomIdx(mdl.atoms.items, res0, 0, "CA", 'B', orig).?, 2);
    try testing.expectEqual(findAtomIdx(mdl.atoms.items, res0, 0, "N", 'A', orig).?, 0);
    try testing.expectEqual(findAtomIdx(mdl.atoms.items, res0, 0, "N", 'B', orig).?, 0);
}

test "generateMovers creates methyl rotator for ALA" {
    const source = @embedFile("../test_data/tiny.cif");
    var mdl = try mmcif.parseModel(testing.allocator, source);
    defer mdl.deinit();

    placer.applyChemistry(&mdl);
    _ = try placer.addHydrogens(&mdl, null, null);

    const gen_result = try generateMovers(testing.allocator, &mdl, false, null, null, null, .neutron);
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
    _ = try placer.addHydrogens(&mdl, null, null);

    const gen_result = try generateMovers(testing.allocator, &mdl, false, null, null, null, .neutron);
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
    _ = try placer.addHydrogens(&mdl, null, null);

    const gen_result = try generateMovers(testing.allocator, &mdl, false, null, null, null, .neutron);
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
    _ = try placer.addHydrogens(&mdl, null, null);

    const gen_result = try generateMovers(testing.allocator, &mdl, false, null, null, null, .neutron);
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
    _ = try placer.addHydrogens(&mdl, null, null);

    const gen_result = try generateMovers(testing.allocator, &mdl, false, null, null, null, .neutron);
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
    _ = try placer.addHydrogens(&mdl, null, null);

    const gen_result = try generateMovers(testing.allocator, &mdl, false, null, null, null, .neutron);
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
    _ = try placer.addHydrogens(&mdl, null, null);

    const gen_result = try generateMovers(testing.allocator, &mdl, true, null, null, null, .neutron);
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
    _ = try placer.addHydrogens(&mdl, null, null);

    const gen_result = try generateMovers(testing.allocator, &mdl, false, null, null, null, .neutron);
    const movers = gen_result.movers;
    defer {
        for (0..movers.len) |i| @constCast(&movers[i]).deinit();
        testing.allocator.free(movers);
    }

    var his_count: u32 = 0;
    for (movers) |m| {
        if (m.kind == .his_flip) {
            his_count += 1;
            try testing.expectEqual(@as(usize, 4), m.orientations.len);
        }
    }
    try testing.expectEqual(@as(u32, 1), his_count);
}

test "generateMovers skips fixed His override" {
    const source = @embedFile("../test_data/his.cif");
    var mdl = try mmcif.parseModel(testing.allocator, source);
    defer mdl.deinit();
    mdl.residues.items[0].auth_seq_id = 1;

    var overrides = try protonation.parseString(testing.allocator,
        \\A:1 HIS HID
    );
    defer overrides.deinit();

    placer.applyChemistryWithConfig(&mdl, .{ .protonation = &overrides });
    _ = try placer.addHydrogensWithConfig(&mdl, null, null, .{ .protonation = &overrides });

    const gen_result = try generateMovers(testing.allocator, &mdl, false, null, null, &overrides, .neutron);
    const movers = gen_result.movers;
    defer {
        for (0..movers.len) |i| @constCast(&movers[i]).deinit();
        testing.allocator.free(movers);
    }

    for (movers) |m| {
        try testing.expect(m.kind != .his_flip);
    }
}
