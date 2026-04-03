//! Unified hydrogen placement entry point.
//!
//! Adds hydrogens to a Model using standard plans for known residues (20 AA)
//! and CCD-derived plans for non-standard residues.

const std = @import("std");
const model_mod = @import("../model.zig");
const Model = model_mod.Model;
const Atom = model_mod.Atom;
const Residue = model_mod.Residue;
const ccd_mod = @import("../ccd.zig");
const ComponentDict = ccd_mod.ComponentDict;
const standard = @import("standard.zig");
const nucleotide = @import("nucleotide.zig");
const modified = @import("modified.zig");
const ccd_derive = @import("ccd_derive.zig");
const geometry = @import("geometry.zig");
const bond_policy = @import("bond_policy.zig");
const math_mod = @import("../math.zig");
const element = @import("../element.zig");

const Vec3f32 = math_mod.Vec3(f32);
const topology = @import("topology.zig");
const chemistry = @import("chemistry.zig");
const protonation = @import("protonation.zig");
const neighbor_mod = @import("../model/neighbor.zig");
const CellList = neighbor_mod.CellList;

pub const lookup = @import("lookup.zig");
pub const terminal = @import("terminal.zig");
pub const execute = @import("execute.zig");
pub const water_mod = @import("water.zig");

const findAtom = lookup.findAtom;
const findAtomPos = lookup.findAtomPos;
const existsInResidue = lookup.existsInResidue;
const padName = lookup.padName;

/// Result of a single executePlan call — why a hydrogen was or wasn't placed.
/// Re-exported from execute.zig; defined there to avoid circular imports.
pub const PlaceResult = execute.PlaceResult;

pub const PlacementResult = struct {
    n_placed: u32 = 0,
    n_skipped_existing: u32 = 0,
    n_skipped_inter_residue: u32 = 0,
    n_skipped_missing_ref: u32 = 0,
    n_skipped_quality_filter: u32 = 0,
    n_residues: u32 = 0,

    pub fn totalSkipped(self: PlacementResult) u32 {
        return self.n_skipped_existing + self.n_skipped_inter_residue + self.n_skipped_missing_ref + self.n_skipped_quality_filter;
    }

    fn tally(self: *PlacementResult, r: PlaceResult) void {
        switch (r) {
            .placed => self.n_placed += 1,
            .existing_h => self.n_skipped_existing += 1,
            .inter_residue => self.n_skipped_inter_residue += 1,
            .missing_parent, .missing_ref => self.n_skipped_missing_ref += 1,
        }
    }
};

pub const WaterConfig = water_mod.WaterConfig;

pub const PlacementConfig = struct {
    water: WaterConfig = .{},
    protonation: ?*const protonation.ProtonationOverrides = null,
    bond_policy: bond_policy.BondPolicy = .{},
};

const AltlocSet = struct { locs: [10]u8, count: usize };

/// Collect distinct non-blank altloc values from a residue's atoms.
fn collectAltlocs(mdl: *const Model, res: Residue) AltlocSet {
    var result: AltlocSet = .{ .locs = .{ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, .count = 0 };
    const atoms = mdl.atoms.items[res.atom_start..res.atom_end];
    for (atoms) |a| {
        if (a.altloc == ' ') continue;
        var found = false;
        for (result.locs[0..result.count]) |existing| {
            if (existing == a.altloc) {
                found = true;
                break;
            }
        }
        if (!found) {
            if (result.count >= result.locs.len) {
                std.log.warn("residue has more than {d} altlocs, excess ignored", .{result.locs.len});
                break;
            }
            result.locs[result.count] = a.altloc;
            result.count += 1;
        }
    }
    return result;
}

/// Add hydrogens to the model.
/// Uses hardcoded plans (standard AA, nucleotides, modified), CCD-derived plans as fallback.
/// For CCD-derived placement, inline_dict is checked first, then ccd_dict (per-component fallback).
/// New atoms are appended to the end of model.atoms. Each new atom carries its residue_idx.
pub fn addHydrogens(
    mdl: *Model,
    ccd_dict: ?*const ComponentDict,
    inline_dict: ?*const ComponentDict,
) !PlacementResult {
    return addHydrogensWithConfig(mdl, ccd_dict, inline_dict, .{});
}

/// Add hydrogens to the model with placement options.
pub fn addHydrogensWithConfig(
    mdl: *Model,
    ccd_dict: ?*const ComponentDict,
    inline_dict: ?*const ComponentDict,
    config: PlacementConfig,
) !PlacementResult {
    var result = PlacementResult{};

    // Build a CellList over all non-hydrogen atoms so that water neighbor
    // queries (max radius 4.0 Å) can run in O(1) instead of O(N).
    // We collect positions and a mapping from CellList index → atom index.
    var water_cell_ctx: ?water_mod.WaterCellCtx = null;
    defer if (water_cell_ctx) |*ctx| {
        ctx.cell_list.deinit();
        mdl.allocator.free(ctx.positions);
        mdl.allocator.free(ctx.atom_indices);
    };

    if (config.water.enabled) blk: {
        var positions = std.ArrayListUnmanaged(Vec3f32){};
        var atom_idx_map = std.ArrayListUnmanaged(u32){};
        errdefer positions.deinit(mdl.allocator);
        errdefer atom_idx_map.deinit(mdl.allocator);

        for (mdl.atoms.items, 0..) |atom, i| {
            if (atom.is_hydrogen) continue;
            try positions.append(mdl.allocator, atom.pos);
            try atom_idx_map.append(mdl.allocator, @intCast(i));
        }

        const pos_slice = try positions.toOwnedSlice(mdl.allocator);
        errdefer mdl.allocator.free(pos_slice);
        const idx_slice = try atom_idx_map.toOwnedSlice(mdl.allocator);
        errdefer mdl.allocator.free(idx_slice);

        const cl = CellList.init(mdl.allocator, pos_slice, 5.0) catch |err| switch (err) {
            error.GridTooLarge => {
                // Fall back to linear scan: free arrays and leave water_cell_ctx null.
                mdl.allocator.free(pos_slice);
                mdl.allocator.free(idx_slice);
                break :blk;
            },
            else => return err,
        };
        water_cell_ctx = water_mod.WaterCellCtx{
            .cell_list = cl,
            .positions = pos_slice,
            .atom_indices = idx_slice,
        };
    }

    const n_residues = mdl.residues.items.len;
    for (0..n_residues) |res_idx| {
        const res = mdl.residues.items[res_idx];
        const comp_id = res.compIdSlice();
        const protonation_state = if (config.protonation) |overrides| overrides.find(mdl, res_idx) else null;

        if (config.water.enabled and water_mod.isWaterResidue(res, comp_id)) {
            const altlocs = collectAltlocs(mdl, res);
            const wr = try water_mod.placeWaterHydrogens(
                mdl,
                res,
                @intCast(res_idx),
                config.water,
                config.bond_policy.mode,
                water_cell_ctx,
                altlocs.count,
                &altlocs.locs,
            );
            result.n_placed += wr.n_placed;
            result.n_skipped_existing += wr.n_skipped_existing;
            result.n_skipped_inter_residue += wr.n_skipped_inter_residue;
            result.n_skipped_missing_ref += wr.n_skipped_missing_ref;
            result.n_skipped_quality_filter += wr.n_skipped_quality_filter;
            result.n_residues += 1;
            continue;
        }

        // Detect terminal residues.
        // The first residue in a chain is always a real N-terminus (gets NH3+/NH2+),
        // even if is_chain_break_before is set due to unobserved leading residues.
        // Mid-chain residues after a gap are terminal for bond topology but keep
        // a single backbone amide H (no NH3+/NH2+ protonation).
        const chain = mdl.chains.items[res.chain_idx];
        const is_real_nterm = (res_idx == chain.residue_start);
        const is_break_nterm = res.is_chain_break_before and (res_idx != chain.residue_start);
        const is_nterm_for_bonding = is_real_nterm or is_break_nterm;
        const is_cterm = (res_idx == chain.residue_end - 1) or
            (res_idx + 1 < n_residues and mdl.residues.items[res_idx + 1].is_chain_break_before);

        if (standard.getPlans(comp_id) orelse
            nucleotide.getPlans(comp_id) orelse
            modified.getPlans(comp_id)) |plans|
        {
            const bonds = topology.getBonds(comp_id); // null for non-standard
            const altlocs = collectAltlocs(mdl, res);

            const targets: []const u8 = if (altlocs.count == 0)
                &[_]u8{' '}
            else
                altlocs.locs[0..altlocs.count];

            // Detect whether plans include a backbone amide H (peptide residues).
            // Nucleotides and residues like PCA (cyclized N) have no backbone H.
            const has_backbone = for (plans) |plan| {
                if (terminal.isBackboneH(&plan)) break true;
            } else false;

            for (targets) |alt| {
                for (plans) |plan| {
                    if (shouldSkipPlanForProtonation(comp_id, &plan, protonation_state)) continue;
                    // Skip backbone amide H only on real N-termini (NH3+/NH2+, not NH).
                    // After a chain break, keep the single amide H.
                    if (is_real_nterm and terminal.isBackboneH(&plan)) continue;

                    result.tally(try execute.executePlan(mdl, res, @intCast(res_idx), &plan, bonds, alt, config.bond_policy.mode));
                }

                // Real N-terminal peptide residues: place NH3+/NH2+ instead of single backbone H.
                // Residues after a chain break keep a single break-amide H.
                if (is_real_nterm and has_backbone) {
                    const nterm = if (std.mem.eql(u8, comp_id, "PRO"))
                        try terminal.placeNtermNH2Pro(mdl, res, @intCast(res_idx), alt, config.bond_policy.mode)
                    else
                        try terminal.placeNtermNH3(mdl, res, @intCast(res_idx), alt, config.bond_policy.mode);
                    result.n_placed += nterm.placed;
                    result.n_skipped_existing += nterm.skipped;
                }

                // 3' terminal nucleotide: place HO3' (leaving atom, absent mid-chain).
                if (is_cterm and nucleotide.getPlans(comp_id) != null) {
                    const oh = try terminal.place3primeOH(mdl, res, @intCast(res_idx), alt, config.bond_policy.mode);
                    result.n_placed += oh.placed;
                    result.n_skipped_existing += oh.skipped;
                }

                if (try execute.placeOverrideHydrogen(mdl, res, @intCast(res_idx), alt, protonation_state, config.bond_policy.mode)) |place_result| {
                    result.tally(place_result);
                }
            }

            result.n_residues += 1;
        } else {
            // Try inline dict first, then fall back to external CCD (per-component fallback)
            const component = if (inline_dict) |d| d.get(comp_id) else null;
            const effective_component = component orelse if (ccd_dict) |d| d.get(comp_id) else null;
            if (effective_component) |comp| {
                const existing = try collectAtomNames(mdl.allocator, mdl, res);
                defer mdl.allocator.free(existing);
                const explicit_ir = try collectInterResidueAtomNames(mdl.allocator, mdl, res);
                defer mdl.allocator.free(explicit_ir);

                // For amino-acid-like residues (have N, CA, C), add implicit peptide
                // bond atoms to inter-residue list. Peptide bonds are not in struct_conn
                // but backbone N (non-N-term) and C (non-C-term) participate in them.
                // This ensures leaving H on backbone N (e.g. H2 on 2MR) are skipped.
                const ir_atoms = try addImplicitPeptideBondAtoms(
                    mdl.allocator,
                    mdl,
                    res,
                    explicit_ir,
                    is_nterm_for_bonding,
                    is_cterm,
                );
                defer mdl.allocator.free(ir_atoms);
                const plans = try ccd_derive.derivePlans(mdl.allocator, &comp, existing, ir_atoms);
                defer mdl.allocator.free(plans);

                // For non-5'-terminal nucleotides, skip H on phosphate O (HOP2/HOP3).
                // Phosphodiester bonds are implicit polymer bonds (not in struct_conn),
                // so phosphate H should only be placed at the 5' terminus.
                const skip_phosphate_h = !is_nterm_for_bonding and terminal.hasPhosphorusAtom(mdl, res);

                for (plans) |plan| {
                    if (skip_phosphate_h and terminal.isPhosphateH(plan.h_name)) {
                        result.n_skipped_inter_residue += 1;
                        continue;
                    }
                    result.tally(try execute.executePlan(mdl, res, @intCast(res_idx), &plan, null, ' ', config.bond_policy.mode));
                }

                // 3' terminal nucleotide (CCD path): place HO3' only if CCD
                // plans did not already include it (avoids duplicate placement).
                if (is_cterm and terminal.isNucleotideResidue(mdl, res)) {
                    const ho3_name = padName("HO3'");
                    const already_in_plans = for (plans) |plan| {
                        if (std.mem.eql(u8, &plan.h_name, &ho3_name)) break true;
                    } else false;
                    if (!already_in_plans) {
                        const oh = try terminal.place3primeOH(mdl, res, @intCast(res_idx), ' ', config.bond_policy.mode);
                        result.n_placed += oh.placed;
                        result.n_skipped_existing += oh.skipped;
                    }
                }

                result.n_residues += 1;
            }
        }
    }

    return result;
}

/// Apply residue/atom-specific chemical annotations to heavy atoms.
/// Updates element_type, flags, and vdw_radius on standard-residue heavy atoms.
/// Also applies terminal charge annotations (N-terminal positive, C-terminal negative).
/// Should be called once after parsing, before hydrogen placement.
pub fn applyChemistry(mdl: *Model) void {
    applyChemistryWithConfig(mdl, .{});
}

pub fn applyChemistryWithConfig(mdl: *Model, config: PlacementConfig) void {
    const n_residues = mdl.residues.items.len;
    for (0..n_residues) |res_idx| {
        const res = mdl.residues.items[res_idx];
        const comp_id = res.compIdSlice();
        const protonation_state = if (config.protonation) |overrides| overrides.find(mdl, res_idx) else null;

        // Detect terminal residues.
        // N-terminal positive charge applies only at real chain starts (first residue).
        // C-terminal negative charge applies at chain ends and before chain breaks.
        const chain = mdl.chains.items[res.chain_idx];
        const is_real_nterm = (res_idx == chain.residue_start);
        const is_cterm = (res_idx == chain.residue_end - 1) or
            (res_idx + 1 < n_residues and mdl.residues.items[res_idx + 1].is_chain_break_before);

        const atoms = mdl.atoms.items[res.atom_start..res.atom_end];
        for (atoms) |*atom| {
            if (atom.is_hydrogen) continue;

            // Apply standard residue annotations (replace flags, but preserve bonded_inter_residue)
            const has_std_ann = if (chemistry.getAnnotationWithOverride(comp_id, atom.nameRaw(), protonation_state)) |ann| blk: {
                atom.element_type = ann.atom_type;
                const keep_bonded = atom.flags.bonded_inter_residue;
                atom.flags = ann.flags;
                atom.flags.bonded_inter_residue = keep_bonded;
                atom.vdw_radius = ann.atom_type.info().explicit_radius;
                break :blk true;
            } else false;

            // Apply terminal annotations (merge flags via OR)
            // Only set element_type/vdw_radius for atoms without standard annotation (e.g. OXT)
            if (chemistry.getTerminalAnnotation(atom.nameRaw(), is_real_nterm, is_cterm)) |term_ann| {
                atom.flags = element.mergeFlags(atom.flags, term_ann.flags);
                if (!has_std_ann) {
                    atom.element_type = term_ann.atom_type;
                    atom.vdw_radius = term_ann.atom_type.info().explicit_radius;
                }
            }
        }
    }
}

fn shouldSkipPlanForProtonation(comp_id: []const u8, plan: *const standard.PlacementPlan, state: ?protonation.ResidueState) bool {
    const s = state orelse return false;

    if (std.mem.eql(u8, comp_id, "HIS") and s == .his) {
        if (s.his == .hid and std.mem.eql(u8, &plan.h_name, " HE2")) return true;
        if (s.his == .hie and std.mem.eql(u8, &plan.h_name, " HD1")) return true;
        return false;
    }
    if (std.mem.eql(u8, comp_id, "LYS") and s == .lys and s.lys == .neutral) {
        return std.mem.eql(u8, &plan.h_name, " HZ3");
    }
    if (std.mem.eql(u8, comp_id, "CYS") and s == .cys and s.cys == .thiolate) {
        return std.mem.eql(u8, &plan.h_name, " HG ");
    }
    return false;
}

/// Find the raw [4]u8 name of an atom in the model (using nameMatch for lookup).
/// Returns the actual stored name bytes, which may differ from PDB-padded format.
fn findAtomName(mdl: *const Model, res: Residue, name: [4]u8, target_altloc: u8) ?[4]u8 {
    const a = findAtom(mdl, res, name, target_altloc) orelse return null;
    return a.nameRaw();
}

/// Extend inter-residue atom names with implicit peptide bond atoms.
/// For amino-acid-like residues (having N, CA, C), backbone N participates in
/// a peptide bond to the previous residue (unless N-terminal), and backbone C
/// participates in a peptide bond to the next residue (unless C-terminal).
/// These implicit bonds are not recorded in struct_conn but affect leaving atom logic.
fn addImplicitPeptideBondAtoms(
    allocator: std.mem.Allocator,
    mdl: *const Model,
    res: Residue,
    explicit_ir: []const [4]u8,
    is_nterm: bool,
    is_cterm: bool,
) ![]const [4]u8 {
    // Check if residue has backbone atoms (amino-acid-like)
    const has_n = findAtomPos(mdl, res, .{ ' ', 'N', ' ', ' ' }, ' ') != null;
    const has_ca = findAtomPos(mdl, res, .{ ' ', 'C', 'A', ' ' }, ' ') != null;
    const has_c = findAtomPos(mdl, res, .{ ' ', 'C', ' ', ' ' }, ' ') != null;

    if (!(has_n and has_ca and has_c)) return try allocator.dupe([4]u8, explicit_ir);

    var list = std.ArrayListUnmanaged([4]u8){};
    defer list.deinit(allocator);
    try list.appendSlice(allocator, explicit_ir);

    // Non-N-terminal: backbone N has implicit peptide bond to previous C.
    // Use model atom name format (setName style, not PDB-padded).
    if (!is_nterm) try list.append(allocator, findAtomName(mdl, res, .{ ' ', 'N', ' ', ' ' }, ' ') orelse .{ 'N', ' ', ' ', ' ' });
    // Non-C-terminal: backbone C has implicit peptide bond to next N
    if (!is_cterm) try list.append(allocator, findAtomName(mdl, res, .{ ' ', 'C', ' ', ' ' }, ' ') orelse .{ 'C', ' ', ' ', ' ' });

    return try allocator.dupe([4]u8, list.items);
}

/// Collect existing atom names in a residue as [4]u8 arrays.
/// Caller must free the returned slice with the provided allocator.
fn collectAtomNames(allocator: std.mem.Allocator, mdl: *const Model, res: Residue) ![][4]u8 {
    const atoms = mdl.atoms.items[res.atom_start..res.atom_end];
    const result = try allocator.alloc([4]u8, atoms.len);
    for (atoms, 0..) |a, i| {
        result[i] = a.nameRaw();
    }
    return result;
}

/// Collect names of heavy atoms in this residue that carry an inter-residue bond.
/// Iterates mdl.bonds so that an atom involved in multiple inter-residue bonds
/// produces multiple entries (one per bond), giving correct valence accounting.
fn collectInterResidueAtomNames(allocator: std.mem.Allocator, mdl: *const Model, res: Residue) ![]const [4]u8 {
    const bond_mod = @import("../model/bond.zig");
    var names = std.ArrayListUnmanaged([4]u8){};
    defer names.deinit(allocator);
    for (mdl.bonds.items) |bond| {
        if (bond.source == bond_mod.BondSource.struct_conn or bond.source == bond_mod.BondSource.branch_link) {
            if (bond.atom_1 >= res.atom_start and bond.atom_1 < res.atom_end) {
                const atom = mdl.atoms.items[bond.atom_1];
                if (!atom.is_hydrogen) try names.append(allocator, atom.nameRaw());
            }
            if (bond.atom_2 >= res.atom_start and bond.atom_2 < res.atom_end) {
                const atom = mdl.atoms.items[bond.atom_2];
                if (!atom.is_hydrogen) try names.append(allocator, atom.nameRaw());
            }
        }
    }
    return try allocator.dupe([4]u8, names.items);
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

const testing = std.testing;
const mmcif = @import("../mmcif.zig");

test "place hydrogens on ALA" {
    const source = @embedFile("../test_data/tiny.cif");
    var mdl = try mmcif.parseModel(testing.allocator, source);
    defer mdl.deinit();

    const initial_count = mdl.atoms.items.len;
    try testing.expectEqual(@as(usize, 5), initial_count);

    const result = try addHydrogens(&mdl, null, null);

    // ALA should get H atoms added
    try testing.expect(result.n_placed > 0);
    try testing.expect(mdl.atoms.items.len > initial_count);
    try testing.expectEqual(@as(u32, 1), result.n_residues);

    // Find HA atom and check bond length to CA (~1.10 A)
    const ca_pos = mdl.atoms.items[1].pos; // CA is index 1
    var found_ha = false;
    for (mdl.atoms.items) |atom| {
        if (std.mem.eql(u8, atom.nameSlice(), "HA")) {
            found_ha = true;
            const dist = atom.pos.distance(ca_pos);
            try testing.expect(dist > 0.8);
            try testing.expect(dist < 1.4);
            break;
        }
    }
    try testing.expect(found_ha);
}

test "xray mode produces shorter C-H bond than neutron" {
    const source = @embedFile("../test_data/tiny.cif");

    // Place in neutron mode (default)
    var mdl_n = try mmcif.parseModel(testing.allocator, source);
    defer mdl_n.deinit();
    _ = try addHydrogensWithConfig(&mdl_n, null, null, .{ .bond_policy = .{ .mode = .neutron } });

    // Place in xray mode
    var mdl_x = try mmcif.parseModel(testing.allocator, source);
    defer mdl_x.deinit();
    _ = try addHydrogensWithConfig(&mdl_x, null, null, .{ .bond_policy = .{ .mode = .xray } });

    // Find HA in both and compare distance to CA
    const ca_pos_n = mdl_n.atoms.items[1].pos;
    const ca_pos_x = mdl_x.atoms.items[1].pos;
    var dist_n: f32 = 0;
    var dist_x: f32 = 0;
    for (mdl_n.atoms.items) |atom| {
        if (std.mem.eql(u8, atom.nameSlice(), "HA")) {
            dist_n = atom.pos.distance(ca_pos_n);
            break;
        }
    }
    for (mdl_x.atoms.items) |atom| {
        if (std.mem.eql(u8, atom.nameSlice(), "HA")) {
            dist_x = atom.pos.distance(ca_pos_x);
            break;
        }
    }
    // Neutron C-H ~1.10, xray C-H ~0.98
    try testing.expect(dist_n > 1.0);
    try testing.expectApproxEqAbs(@as(f32, 0.98), dist_x, 0.02);
    try testing.expect(dist_x < dist_n);
}

test "protonation override fixes HIS tautomer during placement" {
    const source = @embedFile("../test_data/his.cif");
    var mdl = try mmcif.parseModel(testing.allocator, source);
    defer mdl.deinit();
    mdl.residues.items[0].auth_seq_id = 1;

    var overrides = try protonation.parseString(testing.allocator,
        \\A:1 HIS HIE
    );
    defer overrides.deinit();

    applyChemistryWithConfig(&mdl, .{ .protonation = &overrides });
    _ = try addHydrogensWithConfig(&mdl, null, null, .{ .protonation = &overrides });

    try testing.expect(findAddedAtomIdx(&mdl, 0, "HE2") != null);
    try testing.expect(findAddedAtomIdx(&mdl, 0, "HD1") == null);
}

test "protonation override adds ASP sidechain proton" {
    const source =
        \\data_ASP
        \\#
        \\loop_
        \\_atom_site.group_PDB
        \\_atom_site.id
        \\_atom_site.type_symbol
        \\_atom_site.label_atom_id
        \\_atom_site.label_comp_id
        \\_atom_site.label_asym_id
        \\_atom_site.label_seq_id
        \\_atom_site.Cartn_x
        \\_atom_site.Cartn_y
        \\_atom_site.Cartn_z
        \\_atom_site.occupancy
        \\_atom_site.B_iso_or_equiv
        \\_atom_site.label_alt_id
        \\_atom_site.auth_asym_id
        \\_atom_site.auth_seq_id
        \\ATOM 1 N N ASP A 1 0.0 0.0 0.0 1.00 10.0 . A 1
        \\ATOM 2 C CA ASP A 1 1.5 0.0 0.0 1.00 10.0 . A 1
        \\ATOM 3 C C ASP A 1 2.1 1.4 0.0 1.00 10.0 . A 1
        \\ATOM 4 O O ASP A 1 3.3 1.6 0.0 1.00 10.0 . A 1
        \\ATOM 5 C CB ASP A 1 2.0 -0.8 1.2 1.00 10.0 . A 1
        \\ATOM 6 C CG ASP A 1 3.4 -0.4 1.4 1.00 10.0 . A 1
        \\ATOM 7 O OD1 ASP A 1 4.2 0.3 0.7 1.00 10.0 . A 1
        \\ATOM 8 O OD2 ASP A 1 3.8 -0.8 2.6 1.00 10.0 . A 1
        \\#
    ;
    var mdl = try mmcif.parseModel(testing.allocator, source);
    defer mdl.deinit();

    var overrides = try protonation.parseString(testing.allocator,
        \\A:1 ASP OD2
    );
    defer overrides.deinit();

    applyChemistryWithConfig(&mdl, .{ .protonation = &overrides });
    _ = try addHydrogensWithConfig(&mdl, null, null, .{ .protonation = &overrides });

    try testing.expect(findAddedAtomIdx(&mdl, 0, "HD2") != null);
}

test "protonation override adds GLU sidechain proton" {
    const source =
        \\data_GLU
        \\#
        \\loop_
        \\_atom_site.group_PDB
        \\_atom_site.id
        \\_atom_site.type_symbol
        \\_atom_site.label_atom_id
        \\_atom_site.label_comp_id
        \\_atom_site.label_asym_id
        \\_atom_site.label_seq_id
        \\_atom_site.Cartn_x
        \\_atom_site.Cartn_y
        \\_atom_site.Cartn_z
        \\_atom_site.occupancy
        \\_atom_site.B_iso_or_equiv
        \\_atom_site.label_alt_id
        \\_atom_site.auth_asym_id
        \\_atom_site.auth_seq_id
        \\ATOM 1 N N GLU A 1 0.0 0.0 0.0 1.00 10.0 . A 1
        \\ATOM 2 C CA GLU A 1 1.5 0.0 0.0 1.00 10.0 . A 1
        \\ATOM 3 C C GLU A 1 2.1 1.4 0.0 1.00 10.0 . A 1
        \\ATOM 4 O O GLU A 1 3.3 1.6 0.0 1.00 10.0 . A 1
        \\ATOM 5 C CB GLU A 1 2.0 -0.8 1.2 1.00 10.0 . A 1
        \\ATOM 6 C CG GLU A 1 3.4 -0.4 1.4 1.00 10.0 . A 1
        \\ATOM 7 C CD GLU A 1 4.3 -1.3 0.6 1.00 10.0 . A 1
        \\ATOM 8 O OE1 GLU A 1 5.5 -0.9 0.4 1.00 10.0 . A 1
        \\ATOM 9 O OE2 GLU A 1 3.8 -2.4 0.2 1.00 10.0 . A 1
        \\#
    ;
    var mdl = try mmcif.parseModel(testing.allocator, source);
    defer mdl.deinit();

    var overrides = try protonation.parseString(testing.allocator,
        \\A:1 GLU OE2
    );
    defer overrides.deinit();

    applyChemistryWithConfig(&mdl, .{ .protonation = &overrides });
    _ = try addHydrogensWithConfig(&mdl, null, null, .{ .protonation = &overrides });

    try testing.expect(findAddedAtomIdx(&mdl, 0, "HE2") != null);
    try testing.expect(findAddedAtomIdx(&mdl, 0, "HE1") == null);
}

test "protonation override LYS neutral skips HZ3" {
    const source =
        \\data_LYS
        \\#
        \\loop_
        \\_atom_site.group_PDB
        \\_atom_site.id
        \\_atom_site.type_symbol
        \\_atom_site.label_atom_id
        \\_atom_site.label_comp_id
        \\_atom_site.label_asym_id
        \\_atom_site.label_seq_id
        \\_atom_site.Cartn_x
        \\_atom_site.Cartn_y
        \\_atom_site.Cartn_z
        \\_atom_site.occupancy
        \\_atom_site.B_iso_or_equiv
        \\_atom_site.label_alt_id
        \\_atom_site.auth_asym_id
        \\_atom_site.auth_seq_id
        \\ATOM  1 N N   LYS A 1  0.000  0.000  0.000 1.00 10.0 . A 1
        \\ATOM  2 C CA  LYS A 1  1.458  0.000  0.000 1.00 10.0 . A 1
        \\ATOM  3 C C   LYS A 1  2.009  1.420  0.000 1.00 10.0 . A 1
        \\ATOM  4 O O   LYS A 1  3.200  1.600  0.000 1.00 10.0 . A 1
        \\ATOM  5 C CB  LYS A 1  1.986 -0.760  1.220 1.00 10.0 . A 1
        \\ATOM  6 C CG  LYS A 1  3.500 -0.800  1.220 1.00 10.0 . A 1
        \\ATOM  7 C CD  LYS A 1  4.028 -1.560  2.440 1.00 10.0 . A 1
        \\ATOM  8 C CE  LYS A 1  5.542 -1.600  2.440 1.00 10.0 . A 1
        \\ATOM  9 N NZ  LYS A 1  6.070 -2.360  3.660 1.00 10.0 . A 1
        \\#
    ;
    var mdl = try mmcif.parseModel(testing.allocator, source);
    defer mdl.deinit();

    var overrides = try protonation.parseString(testing.allocator,
        \\A:1 LYS NEUTRAL
    );
    defer overrides.deinit();

    applyChemistryWithConfig(&mdl, .{ .protonation = &overrides });
    _ = try addHydrogensWithConfig(&mdl, null, null, .{ .protonation = &overrides });

    try testing.expect(findAddedAtomIdx(&mdl, 0, "HZ1") != null);
    try testing.expect(findAddedAtomIdx(&mdl, 0, "HZ2") != null);
    try testing.expect(findAddedAtomIdx(&mdl, 0, "HZ3") == null);
}

test "protonation override CYS thiolate skips HG" {
    const source =
        \\data_CYS
        \\#
        \\loop_
        \\_atom_site.group_PDB
        \\_atom_site.id
        \\_atom_site.type_symbol
        \\_atom_site.label_atom_id
        \\_atom_site.label_comp_id
        \\_atom_site.label_asym_id
        \\_atom_site.label_seq_id
        \\_atom_site.Cartn_x
        \\_atom_site.Cartn_y
        \\_atom_site.Cartn_z
        \\_atom_site.occupancy
        \\_atom_site.B_iso_or_equiv
        \\_atom_site.label_alt_id
        \\_atom_site.auth_asym_id
        \\_atom_site.auth_seq_id
        \\ATOM 1 N N   CYS A 1  0.000  0.000  0.000 1.00 10.0 . A 1
        \\ATOM 2 C CA  CYS A 1  1.458  0.000  0.000 1.00 10.0 . A 1
        \\ATOM 3 C C   CYS A 1  2.009  1.420  0.000 1.00 10.0 . A 1
        \\ATOM 4 O O   CYS A 1  3.200  1.600  0.000 1.00 10.0 . A 1
        \\ATOM 5 C CB  CYS A 1  1.986 -0.760  1.220 1.00 10.0 . A 1
        \\ATOM 6 S SG  CYS A 1  1.300 -0.200  2.800 1.00 10.0 . A 1
        \\#
    ;
    var mdl = try mmcif.parseModel(testing.allocator, source);
    defer mdl.deinit();

    var overrides = try protonation.parseString(testing.allocator,
        \\A:1 CYS THIOLATE
    );
    defer overrides.deinit();

    applyChemistryWithConfig(&mdl, .{ .protonation = &overrides });
    _ = try addHydrogensWithConfig(&mdl, null, null, .{ .protonation = &overrides });

    try testing.expect(findAddedAtomIdx(&mdl, 0, "HG") == null);
}

fn findAddedAtomIdx(mdl: *const Model, residue_idx: u32, name: []const u8) ?u32 {
    for (mdl.atoms.items, 0..) |atom, idx| {
        if (atom.residue_idx != residue_idx) continue;
        if (!atom.is_added) continue;
        if (std.mem.eql(u8, atom.nameSlice(), name)) return @intCast(idx);
    }
    return null;
}

test "placed atoms have correct metadata" {
    const source = @embedFile("../test_data/tiny.cif");
    var mdl = try mmcif.parseModel(testing.allocator, source);
    defer mdl.deinit();

    _ = try addHydrogens(&mdl, null, null);

    // Check that newly added atoms have correct flags
    for (mdl.atoms.items[5..]) |atom| {
        try testing.expect(atom.is_hydrogen);
        try testing.expect(atom.is_added);
        try testing.expectEqual(@as(u32, 0), atom.residue_idx);
    }
}

test "PlacementResult tracks counts" {
    const source = @embedFile("../test_data/tiny.cif");
    var mdl = try mmcif.parseModel(testing.allocator, source);
    defer mdl.deinit();

    const result = try addHydrogens(&mdl, null, null);

    // ALA has 5 plans; backbone H skipped on N-term but NH3+ (H1,H2,H3) added = 4+3=7
    try testing.expectEqual(@as(u32, 7), result.n_placed + result.totalSkipped());
}

test "placed H inherits parent atom metadata" {
    const source = @embedFile("../test_data/tiny.cif");
    var mdl = try mmcif.parseModel(testing.allocator, source);
    defer mdl.deinit();

    // tiny.cif atoms have occupancy=1.0, b_factor=10.0, altloc=' '
    _ = try addHydrogens(&mdl, null, null);

    // All placed H atoms should inherit b_factor=10.0 from parent
    for (mdl.atoms.items) |atom| {
        if (atom.is_added and atom.is_hydrogen) {
            try testing.expectEqual(@as(f32, 10.0), atom.b_factor);
            try testing.expectEqual(@as(f32, 1.0), atom.occupancy);
        }
    }
}

test "duplicate H atoms are not placed" {
    const source = @embedFile("../test_data/ala_with_h.cif");
    var mdl = try mmcif.parseModel(testing.allocator, source);
    defer mdl.deinit();

    _ = try addHydrogens(&mdl, null, null);

    // Count HA atoms — should be exactly 1 (the pre-existing one)
    var ha_count: u32 = 0;
    for (mdl.atoms.items) |atom| {
        if (std.mem.eql(u8, atom.nameSlice(), "HA")) ha_count += 1;
    }
    try testing.expectEqual(@as(u32, 1), ha_count);

    // The original HA should NOT be overwritten (b_factor should remain 12.0)
    for (mdl.atoms.items) |atom| {
        if (std.mem.eql(u8, atom.nameSlice(), "HA")) {
            try testing.expectEqual(@as(f32, 12.0), atom.b_factor);
            break;
        }
    }
}

test "PlacementResult counts duplicates as skipped" {
    const source = @embedFile("../test_data/ala_with_h.cif");
    var mdl = try mmcif.parseModel(testing.allocator, source);
    defer mdl.deinit();

    const result = try addHydrogens(&mdl, null, null);

    // HA was pre-existing so should be counted as skipped (existing_h)
    // Total plans attempted should still be the same as clean ALA
    try testing.expect(result.n_skipped_existing >= 1);
    try testing.expectEqual(@as(u32, 1), result.n_residues);
}

test "placement succeeds on stretched geometry with bond topology" {
    const source = @embedFile("../test_data/ala_stretched.cif");
    var mdl = try mmcif.parseModel(testing.allocator, source);
    defer mdl.deinit();

    _ = try addHydrogens(&mdl, null, null);

    // With bond topology, HA should be placed even though CB is >1.9A from CA
    // (HA placement type is hxr3 which needs the 3rd neighbor = CB)
    var found_ha = false;
    for (mdl.atoms.items) |atom| {
        if (std.mem.eql(u8, atom.nameSlice(), "HA")) {
            found_ha = true;
            break;
        }
    }
    try testing.expect(found_ha);
}

test "applyChemistry sets backbone C to C_eq_O" {
    const source = @embedFile("../test_data/tiny.cif");
    var mdl = try mmcif.parseModel(testing.allocator, source);
    defer mdl.deinit();

    // Before: C has generic element type
    // tiny.cif atoms: N(0), CA(1), C(2), O(3), CB(4)
    try testing.expectEqual(element.AtomType.C, mdl.atoms.items[2].element_type);

    applyChemistry(&mdl);

    // After: C has carbonyl type
    try testing.expectEqual(element.AtomType.C_eq_O, mdl.atoms.items[2].element_type);
    try testing.expectApproxEqAbs(@as(f32, 1.65), mdl.atoms.items[2].vdw_radius, 1e-6);
}

test "applyChemistry sets backbone O acceptor flag" {
    const source = @embedFile("../test_data/tiny.cif");
    var mdl = try mmcif.parseModel(testing.allocator, source);
    defer mdl.deinit();

    try testing.expect(!mdl.atoms.items[3].flags.acceptor);
    applyChemistry(&mdl);
    try testing.expect(mdl.atoms.items[3].flags.acceptor);
}

test "applyChemistry sets backbone N donor flag" {
    const source = @embedFile("../test_data/tiny.cif");
    var mdl = try mmcif.parseModel(testing.allocator, source);
    defer mdl.deinit();

    try testing.expect(!mdl.atoms.items[0].flags.donor);
    applyChemistry(&mdl);
    try testing.expect(mdl.atoms.items[0].flags.donor);
}

test "placed H atoms have correct flags from element table" {
    const source = @embedFile("../test_data/tiny.cif");
    var mdl = try mmcif.parseModel(testing.allocator, source);
    defer mdl.deinit();

    applyChemistry(&mdl);
    _ = try addHydrogens(&mdl, null, null);

    for (mdl.atoms.items) |atom| {
        if (atom.is_added and atom.is_hydrogen) {
            const type_flags = atom.element_type.info().flags;
            // Hpol atoms should have donor flag
            if (atom.element_type == .Hpol) {
                try testing.expect(atom.flags.donor);
            }
            // All placed H flags should match their element_type flags
            try testing.expectEqual(type_flags.donor, atom.flags.donor);
            try testing.expectEqual(type_flags.aromatic, atom.flags.aromatic);
        }
    }
}

test "applyChemistry adds positive flag to N-terminal N" {
    const source = @embedFile("../test_data/tiny.cif");
    var mdl = try mmcif.parseModel(testing.allocator, source);
    defer mdl.deinit();

    applyChemistry(&mdl);

    // tiny.cif has single ALA — both N-term and C-term
    // N atom (index 0) should have donor (standard) + positive (terminal)
    try testing.expect(mdl.atoms.items[0].flags.donor);
    try testing.expect(mdl.atoms.items[0].flags.positive);
}

test "applyChemistry adds negative flag to C-terminal O" {
    const source = @embedFile("../test_data/tiny.cif");
    var mdl = try mmcif.parseModel(testing.allocator, source);
    defer mdl.deinit();

    applyChemistry(&mdl);

    // O atom (index 3) should have acceptor (standard) + negative (terminal)
    try testing.expect(mdl.atoms.items[3].flags.acceptor);
    try testing.expect(mdl.atoms.items[3].flags.negative);
}

test "applyChemistry annotates OXT as negative acceptor" {
    const source = @embedFile("../test_data/ala_cterm.cif");
    var mdl = try mmcif.parseModel(testing.allocator, source);
    defer mdl.deinit();

    applyChemistry(&mdl);

    // Find OXT atom
    var oxt_found = false;
    for (mdl.atoms.items) |atom| {
        if (std.mem.eql(u8, atom.nameSlice(), "OXT")) {
            oxt_found = true;
            try testing.expectEqual(element.AtomType.O, atom.element_type);
            try testing.expect(atom.flags.negative);
            try testing.expect(atom.flags.acceptor);
            break;
        }
    }
    try testing.expect(oxt_found);
}

test "multi-chain terminal detection is correct" {
    const source = @embedFile("../test_data/multi_chain.cif");
    var mdl = try mmcif.parseModel(testing.allocator, source);
    defer mdl.deinit();

    applyChemistry(&mdl);

    // Chain A: ALA(N-term, res 0) + GLY(C-term, res 1)
    // Chain B: VAL(both N-term and C-term, res 2)

    // ALA N (atom 0): N-terminal → positive + donor
    try testing.expect(mdl.atoms.items[0].flags.positive);
    try testing.expect(mdl.atoms.items[0].flags.donor);

    // GLY N (atom 4): internal-ish (C-terminal residue, but N is not annotated for C-term)
    try testing.expect(!mdl.atoms.items[4].flags.positive);

    // GLY O (atom 7): C-terminal → negative + acceptor
    try testing.expect(mdl.atoms.items[7].flags.negative);
    try testing.expect(mdl.atoms.items[7].flags.acceptor);

    // VAL N (atom 8): N-terminal of chain B → positive + donor
    try testing.expect(mdl.atoms.items[8].flags.positive);
    try testing.expect(mdl.atoms.items[8].flags.donor);
}

test "OXT does not receive hydrogen atoms" {
    const source = @embedFile("../test_data/ala_cterm.cif");
    var mdl = try mmcif.parseModel(testing.allocator, source);
    defer mdl.deinit();

    applyChemistry(&mdl);
    _ = try addHydrogens(&mdl, null, null);

    // No hydrogen should be bonded to OXT
    for (mdl.atoms.items) |atom| {
        if (atom.is_added and atom.is_hydrogen) {
            const name = atom.nameSlice();
            try testing.expect(!std.mem.eql(u8, name, "HOXT"));
        }
    }
}

test "multi-conformer residue places H per conformer" {
    const source = @embedFile("../test_data/ala_altloc.cif");
    var mdl = try mmcif.parseModel(testing.allocator, source);
    defer mdl.deinit();

    applyChemistry(&mdl);
    const result = try addHydrogens(&mdl, null, null);

    // Should place H for both conformers A and B
    var ha_a_count: u32 = 0;
    var ha_b_count: u32 = 0;
    for (mdl.atoms.items) |atom| {
        if (atom.is_added and std.mem.eql(u8, atom.nameSlice(), "HA")) {
            if (atom.altloc == 'A') ha_a_count += 1;
            if (atom.altloc == 'B') ha_b_count += 1;
        }
    }
    try testing.expectEqual(@as(u32, 1), ha_a_count);
    try testing.expectEqual(@as(u32, 1), ha_b_count);
    try testing.expect(result.n_placed > 0);
}

test "conformer A and B H atoms have different positions" {
    const source = @embedFile("../test_data/ala_altloc.cif");
    var mdl = try mmcif.parseModel(testing.allocator, source);
    defer mdl.deinit();

    applyChemistry(&mdl);
    _ = try addHydrogens(&mdl, null, null);

    var ha_a_pos: ?Vec3f32 = null;
    var ha_b_pos: ?Vec3f32 = null;
    for (mdl.atoms.items) |atom| {
        if (atom.is_added and std.mem.eql(u8, atom.nameSlice(), "HA")) {
            if (atom.altloc == 'A') ha_a_pos = atom.pos;
            if (atom.altloc == 'B') ha_b_pos = atom.pos;
        }
    }
    try testing.expect(ha_a_pos != null);
    try testing.expect(ha_b_pos != null);
    const diff = ha_a_pos.?.distance(ha_b_pos.?);
    try testing.expect(diff > 0.01);
}

test "placed H atoms have correct mover_hint" {
    const source = @embedFile("../test_data/tiny.cif");
    var mdl = try mmcif.parseModel(testing.allocator, source);
    defer mdl.deinit();

    applyChemistry(&mdl);
    _ = try addHydrogens(&mdl, null, null);

    // ALA methyl H (HB1/HB2/HB3) should have rotate_methyl hint
    var methyl_count: u32 = 0;
    for (mdl.atoms.items) |atom| {
        if (atom.is_added and atom.mover_hint == .rotate_methyl) {
            methyl_count += 1;
        }
    }
    // ALA has 3 methyl H on CB
    try testing.expectEqual(@as(u32, 3), methyl_count);
}

test "chain break residue keeps single backbone H" {
    const source = @embedFile("../test_data/gap_chain.cif");
    var mdl = try mmcif.parseModel(testing.allocator, source);
    defer mdl.deinit();

    applyChemistry(&mdl);
    _ = try addHydrogens(&mdl, null, null);

    // Second residue (seq_id 3, index 1) is after a chain break.
    // It should keep the single backbone amide H, not gain NH3+.
    var h_count: u32 = 0;
    var h1_count: u32 = 0;
    var h2_count: u32 = 0;
    var h3_count: u32 = 0;
    for (mdl.atoms.items) |atom| {
        if (!atom.is_added or atom.residue_idx != 1) continue;
        if (std.mem.eql(u8, atom.nameSlice(), "H")) h_count += 1;
        if (std.mem.eql(u8, atom.nameSlice(), "H1")) h1_count += 1;
        if (std.mem.eql(u8, atom.nameSlice(), "H2")) h2_count += 1;
        if (std.mem.eql(u8, atom.nameSlice(), "H3")) h3_count += 1;
    }
    try testing.expectEqual(@as(u32, 1), h_count);
    try testing.expectEqual(@as(u32, 0), h1_count);
    try testing.expectEqual(@as(u32, 0), h2_count);
    try testing.expectEqual(@as(u32, 0), h3_count);
}

test "chain break residue is not annotated as positively charged N-terminus" {
    const source = @embedFile("../test_data/gap_chain.cif");
    var mdl = try mmcif.parseModel(testing.allocator, source);
    defer mdl.deinit();

    applyChemistry(&mdl);

    const res1 = mdl.residues.items[1];
    const atoms1 = mdl.atoms.items[res1.atom_start..res1.atom_end];
    var n_positive = false;
    for (atoms1) |atom| {
        if (std.mem.eql(u8, atom.nameSlice(), "N")) {
            n_positive = atom.flags.positive;
        }
    }
    try testing.expect(!n_positive);
}

test "first observed residue with N-terminal disorder gets NH3+" {
    const source = @embedFile("../test_data/nterm_disorder.cif");
    var mdl = try mmcif.parseModel(testing.allocator, source);
    defer mdl.deinit();

    // First observed residue (seq_id 2, index 0) has is_chain_break_before
    // because seq_id 1 is unobserved. It is still the physical N-terminus
    // and should receive NH3+ (H1, H2, H3), not a single amide H.
    try testing.expect(mdl.residues.items[0].is_chain_break_before);

    applyChemistry(&mdl);
    _ = try addHydrogens(&mdl, null, null);

    var h_count: u32 = 0;
    var h1_count: u32 = 0;
    var h2_count: u32 = 0;
    var h3_count: u32 = 0;
    for (mdl.atoms.items) |atom| {
        if (!atom.is_added or atom.residue_idx != 0) continue;
        if (std.mem.eql(u8, atom.nameSlice(), "H")) h_count += 1;
        if (std.mem.eql(u8, atom.nameSlice(), "H1")) h1_count += 1;
        if (std.mem.eql(u8, atom.nameSlice(), "H2")) h2_count += 1;
        if (std.mem.eql(u8, atom.nameSlice(), "H3")) h3_count += 1;
    }
    // NH3+: no single backbone H, but H1/H2/H3 present
    try testing.expectEqual(@as(u32, 0), h_count);
    try testing.expectEqual(@as(u32, 1), h1_count);
    try testing.expectEqual(@as(u32, 1), h2_count);
    try testing.expectEqual(@as(u32, 1), h3_count);
}

test "first observed residue with N-terminal disorder gets positive charge" {
    const source = @embedFile("../test_data/nterm_disorder.cif");
    var mdl = try mmcif.parseModel(testing.allocator, source);
    defer mdl.deinit();

    applyChemistry(&mdl);

    const res0 = mdl.residues.items[0];
    const atoms0 = mdl.atoms.items[res0.atom_start..res0.atom_end];
    var n_positive = false;
    for (atoms0) |atom| {
        if (std.mem.eql(u8, atom.nameSlice(), "N")) {
            n_positive = atom.flags.positive;
        }
    }
    try testing.expect(n_positive);
}

test "residue before chain break gets C-terminal charge" {
    const source = @embedFile("../test_data/gap_chain.cif");
    var mdl = try mmcif.parseModel(testing.allocator, source);
    defer mdl.deinit();

    applyChemistry(&mdl);

    // First residue (seq_id 1, index 0): before gap -> C-terminal -> O gets negative
    const res0 = mdl.residues.items[0];
    const atoms0 = mdl.atoms.items[res0.atom_start..res0.atom_end];
    var o_negative = false;
    for (atoms0) |atom| {
        if (std.mem.eql(u8, atom.nameSlice(), "O")) {
            o_negative = atom.flags.negative;
        }
    }
    try testing.expect(o_negative);
}

test "addHydrogens skips H on bonded_inter_residue atom" {
    const mmcif_mod = @import("../mmcif.zig");
    const source = @embedFile("../test_data/disulfide.cif");
    var mdl = try mmcif_mod.parseModel(testing.allocator, source);
    defer mdl.deinit();

    // Apply chemistry first (it overwrites flags), then set bonded_inter_residue
    applyChemistry(&mdl);

    // Manually set bonded_inter_residue on SG atoms (index 5 and 11)
    mdl.atoms.items[5].flags.bonded_inter_residue = true;
    mdl.atoms.items[11].flags.bonded_inter_residue = true;

    const result = try addHydrogens(&mdl, null, null);

    // Verify no HG was placed on either CYS SG
    for (mdl.atoms.items) |atom| {
        if (atom.is_added) {
            const name = atom.nameSlice();
            // SG should not have HG placed
            try std.testing.expect(!std.mem.eql(u8, name, "HG"));
        }
    }
    _ = result;
}

test "water placement adds two hydrogens when enabled" {
    const water_cif =
        \\data_WATER
        \\#
        \\loop_
        \\_entity.id
        \\_entity.type
        \\1 polymer
        \\2 water
        \\#
        \\loop_
        \\_atom_site.group_PDB
        \\_atom_site.id
        \\_atom_site.type_symbol
        \\_atom_site.label_atom_id
        \\_atom_site.label_comp_id
        \\_atom_site.label_asym_id
        \\_atom_site.label_entity_id
        \\_atom_site.label_seq_id
        \\_atom_site.Cartn_x
        \\_atom_site.Cartn_y
        \\_atom_site.Cartn_z
        \\_atom_site.occupancy
        \\_atom_site.B_iso_or_equiv
        \\_atom_site.label_alt_id
        \\ATOM   1  N  N   GLY A 1 1  0.000 0.000 0.000 1.00 10.0 .
        \\ATOM   2  O  O   GLY A 1 1  1.200 1.100 0.000 1.00 10.0 .
        \\ATOM   3  O  O   HOH B 2 1  2.700 0.200 0.000 1.00 12.0 .
        \\#
    ;

    var mdl = try mmcif.parseModel(testing.allocator, water_cif);
    defer mdl.deinit();

    const result = try addHydrogensWithConfig(&mdl, null, null, .{
        .water = .{ .enabled = true },
    });

    try testing.expectEqual(@as(u32, 2), result.n_placed);

    const oxygen_pos = lookup.findAtomPos(&mdl, mdl.residues.items[1], padName("O"), ' ') orelse unreachable;
    var h_positions: [2]Vec3f32 = undefined;
    var h_count: u32 = 0;
    for (mdl.atoms.items) |atom| {
        if (atom.is_added and atom.residue_idx == 1) {
            try testing.expect(atom.is_hydrogen);
            try testing.expectApproxEqAbs(@as(f32, 1.0), atom.occupancy, 1e-6);
            try testing.expectApproxEqAbs(@as(f32, 12.0), atom.b_factor, 1e-6);
            if (h_count < 2) h_positions[h_count] = atom.pos;
            h_count += 1;
        }
    }
    try testing.expectEqual(@as(u32, 2), h_count);

    // Verify O-H bond length (~0.97 A)
    try testing.expectApproxEqAbs(@as(f32, 0.97), oxygen_pos.distance(h_positions[0]), 0.02);
    try testing.expectApproxEqAbs(@as(f32, 0.97), oxygen_pos.distance(h_positions[1]), 0.02);

    // Verify H-O-H angle (~104.5 degrees)
    const hoh_angle = math_mod.angle(f32, h_positions[0], oxygen_pos, h_positions[1]);
    try testing.expectApproxEqAbs(@as(f32, 104.5), hoh_angle, 1.0);
}

test "water placement respects occupancy cutoff" {
    const water_cif =
        \\data_WATER_OCC
        \\#
        \\loop_
        \\_entity.id
        \\_entity.type
        \\1 polymer
        \\2 water
        \\#
        \\loop_
        \\_atom_site.group_PDB
        \\_atom_site.id
        \\_atom_site.type_symbol
        \\_atom_site.label_atom_id
        \\_atom_site.label_comp_id
        \\_atom_site.label_asym_id
        \\_atom_site.label_entity_id
        \\_atom_site.label_seq_id
        \\_atom_site.Cartn_x
        \\_atom_site.Cartn_y
        \\_atom_site.Cartn_z
        \\_atom_site.occupancy
        \\_atom_site.B_iso_or_equiv
        \\_atom_site.label_alt_id
        \\ATOM   1  N  N   GLY A 1 1  0.000 0.000 0.000 1.00 10.0 .
        \\ATOM   2  O  O   GLY A 1 1  1.200 1.100 0.000 1.00 10.0 .
        \\ATOM   3  O  O   HOH B 2 1  2.700 0.200 0.000 0.50 12.0 .
        \\#
    ;

    var mdl = try mmcif.parseModel(testing.allocator, water_cif);
    defer mdl.deinit();

    const result = try addHydrogensWithConfig(&mdl, null, null, .{
        .water = .{
            .enabled = true,
            .occupancy_cutoff = 0.66,
        },
    });

    try testing.expectEqual(@as(u32, 0), result.n_placed);
    try testing.expectEqual(@as(u32, 1), result.n_skipped_quality_filter);
}

test "water placement skips coordinated water oxygen" {
    const water_cif =
        \\data_WATER_METAL
        \\#
        \\loop_
        \\_entity.id
        \\_entity.type
        \\1 polymer
        \\2 water
        \\#
        \\loop_
        \\_atom_site.group_PDB
        \\_atom_site.id
        \\_atom_site.type_symbol
        \\_atom_site.label_atom_id
        \\_atom_site.label_comp_id
        \\_atom_site.label_asym_id
        \\_atom_site.label_entity_id
        \\_atom_site.label_seq_id
        \\_atom_site.Cartn_x
        \\_atom_site.Cartn_y
        \\_atom_site.Cartn_z
        \\_atom_site.occupancy
        \\_atom_site.B_iso_or_equiv
        \\_atom_site.label_alt_id
        \\ATOM   1  N  N   GLY A 1 1  0.000 0.000 0.000 1.00 10.0 .
        \\ATOM   2  O  O   GLY A 1 1  1.200 1.100 0.000 1.00 10.0 .
        \\ATOM   3  O  O   HOH B 2 1  2.700 0.200 0.000 1.00 12.0 .
        \\#
    ;

    var mdl = try mmcif.parseModel(testing.allocator, water_cif);
    defer mdl.deinit();
    mdl.atoms.items[2].flags.bonded_inter_residue = true;

    const result = try addHydrogensWithConfig(&mdl, null, null, .{
        .water = .{ .enabled = true },
    });

    try testing.expectEqual(@as(u32, 0), result.n_placed);
    try testing.expect(result.n_skipped_inter_residue > 0);
}

test "water phantom mode places zero-occupancy hydrogens for isolated water" {
    const water_cif =
        \\data_WATER_PHANTOM
        \\#
        \\loop_
        \\_entity.id
        \\_entity.type
        \\1 water
        \\#
        \\loop_
        \\_atom_site.group_PDB
        \\_atom_site.id
        \\_atom_site.type_symbol
        \\_atom_site.label_atom_id
        \\_atom_site.label_comp_id
        \\_atom_site.label_asym_id
        \\_atom_site.label_entity_id
        \\_atom_site.label_seq_id
        \\_atom_site.Cartn_x
        \\_atom_site.Cartn_y
        \\_atom_site.Cartn_z
        \\_atom_site.occupancy
        \\_atom_site.B_iso_or_equiv
        \\_atom_site.label_alt_id
        \\HETATM 1  O  O   HOH A 1 1  0.000 0.000 0.000 1.00 12.0 .
        \\#
    ;

    var mdl = try mmcif.parseModel(testing.allocator, water_cif);
    defer mdl.deinit();

    const result = try addHydrogensWithConfig(&mdl, null, null, .{
        .water = .{
            .enabled = true,
            .phantom = true,
        },
    });

    try testing.expectEqual(@as(u32, 2), result.n_placed);
    var zero_occ: u32 = 0;
    for (mdl.atoms.items) |atom| {
        if (atom.is_added) {
            try testing.expectApproxEqAbs(@as(f32, 0.0), atom.occupancy, 1e-6);
            zero_occ += 1;
        }
    }
    try testing.expectEqual(@as(u32, 2), zero_occ);
}

test "water placement skips water near metal by distance" {
    const water_cif =
        \\data_WATER_METAL_DIST
        \\#
        \\loop_
        \\_entity.id
        \\_entity.type
        \\1 water
        \\2 non-polymer
        \\#
        \\loop_
        \\_atom_site.group_PDB
        \\_atom_site.id
        \\_atom_site.type_symbol
        \\_atom_site.label_atom_id
        \\_atom_site.label_comp_id
        \\_atom_site.label_asym_id
        \\_atom_site.label_entity_id
        \\_atom_site.label_seq_id
        \\_atom_site.Cartn_x
        \\_atom_site.Cartn_y
        \\_atom_site.Cartn_z
        \\_atom_site.occupancy
        \\_atom_site.B_iso_or_equiv
        \\_atom_site.label_alt_id
        \\HETATM 1  O  O   HOH A 1 1  0.000 0.000 0.000 1.00 12.0 .
        \\HETATM 2  ZN ZN   ZN B 2 1  2.500 0.000 0.000 1.00 10.0 .
        \\#
    ;

    var mdl = try mmcif.parseModel(testing.allocator, water_cif);
    defer mdl.deinit();

    const result = try addHydrogensWithConfig(&mdl, null, null, .{
        .water = .{ .enabled = true },
    });

    try testing.expectEqual(@as(u32, 0), result.n_placed);
    try testing.expect(result.n_skipped_inter_residue > 0);
}

test "water placement respects B-factor cutoff" {
    const water_cif =
        \\data_WATER_BFAC
        \\#
        \\loop_
        \\_entity.id
        \\_entity.type
        \\1 polymer
        \\2 water
        \\#
        \\loop_
        \\_atom_site.group_PDB
        \\_atom_site.id
        \\_atom_site.type_symbol
        \\_atom_site.label_atom_id
        \\_atom_site.label_comp_id
        \\_atom_site.label_asym_id
        \\_atom_site.label_entity_id
        \\_atom_site.label_seq_id
        \\_atom_site.Cartn_x
        \\_atom_site.Cartn_y
        \\_atom_site.Cartn_z
        \\_atom_site.occupancy
        \\_atom_site.B_iso_or_equiv
        \\_atom_site.label_alt_id
        \\ATOM   1  N  N   GLY A 1 1  0.000 0.000 0.000 1.00 10.0 .
        \\ATOM   2  O  O   GLY A 1 1  1.200 1.100 0.000 1.00 10.0 .
        \\ATOM   3  O  O   HOH B 2 1  2.700 0.200 0.000 1.00 50.0 .
        \\#
    ;

    var mdl = try mmcif.parseModel(testing.allocator, water_cif);
    defer mdl.deinit();

    const result = try addHydrogensWithConfig(&mdl, null, null, .{
        .water = .{
            .enabled = true,
            .b_factor_cutoff = 40.0,
        },
    });

    try testing.expectEqual(@as(u32, 0), result.n_placed);
    try testing.expectEqual(@as(u32, 1), result.n_skipped_quality_filter);
}

test "water placement skips water with existing H atoms" {
    const water_cif =
        \\data_WATER_EXIST_H
        \\#
        \\loop_
        \\_entity.id
        \\_entity.type
        \\1 polymer
        \\2 water
        \\#
        \\loop_
        \\_atom_site.group_PDB
        \\_atom_site.id
        \\_atom_site.type_symbol
        \\_atom_site.label_atom_id
        \\_atom_site.label_comp_id
        \\_atom_site.label_asym_id
        \\_atom_site.label_entity_id
        \\_atom_site.label_seq_id
        \\_atom_site.Cartn_x
        \\_atom_site.Cartn_y
        \\_atom_site.Cartn_z
        \\_atom_site.occupancy
        \\_atom_site.B_iso_or_equiv
        \\_atom_site.label_alt_id
        \\ATOM   1  N  N   GLY A 1 1  0.000 0.000 0.000 1.00 10.0 .
        \\ATOM   2  O  O   GLY A 1 1  1.200 1.100 0.000 1.00 10.0 .
        \\ATOM   3  O  O   HOH B 2 1  2.700 0.200 0.000 1.00 12.0 .
        \\ATOM   4  H  H1  HOH B 2 1  3.100 0.800 0.500 1.00 12.0 .
        \\#
    ;

    var mdl = try mmcif.parseModel(testing.allocator, water_cif);
    defer mdl.deinit();

    const result = try addHydrogensWithConfig(&mdl, null, null, .{
        .water = .{ .enabled = true },
    });

    try testing.expectEqual(@as(u32, 0), result.n_placed);
    try testing.expectEqual(@as(u32, 1), result.n_skipped_existing);
}

test "backbone NH placed in peptide plane using C(i-1)" {
    // Two-residue ALA-ALA with realistic geometry.
    // ALA 1: N-term (gets NH3+, no backbone H)
    // ALA 2: should get backbone H in the C(1)-N(2)-CA(2) peptide plane.
    const two_ala_cif =
        \\data_TWO_ALA
        \\#
        \\loop_
        \\_atom_site.group_PDB
        \\_atom_site.id
        \\_atom_site.type_symbol
        \\_atom_site.label_atom_id
        \\_atom_site.label_comp_id
        \\_atom_site.label_asym_id
        \\_atom_site.label_seq_id
        \\_atom_site.Cartn_x
        \\_atom_site.Cartn_y
        \\_atom_site.Cartn_z
        \\_atom_site.occupancy
        \\_atom_site.B_iso_or_equiv
        \\_atom_site.label_alt_id
        \\ATOM 1  N  N   ALA A 1  -1.200  0.000  0.000 1.00 10.0 .
        \\ATOM 2  C  CA  ALA A 1   0.000  0.000  0.000 1.00 10.0 .
        \\ATOM 3  C  C   ALA A 1   0.550  1.420  0.000 1.00 10.0 .
        \\ATOM 4  O  O   ALA A 1   1.720  1.600  0.000 1.00 10.0 .
        \\ATOM 5  C  CB  ALA A 1   0.550 -0.760  1.200 1.00 10.0 .
        \\ATOM 6  N  N   ALA A 2   -0.100  2.500  0.000 1.00 10.0 .
        \\ATOM 7  C  CA  ALA A 2   0.400  3.870  0.000 1.00 10.0 .
        \\ATOM 8  C  C   ALA A 2   1.920  3.900  0.000 1.00 10.0 .
        \\ATOM 9  O  O   ALA A 2   2.500  4.980  0.000 1.00 10.0 .
        \\ATOM 10 C  CB  ALA A 2  -0.100  4.600  1.200 1.00 10.0 .
        \\#
    ;

    const mmcif_mod = @import("../mmcif.zig");
    var mdl = try mmcif_mod.parseModel(testing.allocator, two_ala_cif);
    defer mdl.deinit();

    applyChemistry(&mdl);
    _ = try addHydrogens(&mdl, null, null);

    // Find the backbone H on residue 2 (ALA 2)
    var backbone_h_pos: ?Vec3f32 = null;
    for (mdl.atoms.items) |atom| {
        if (atom.is_added and atom.residue_idx == 1 and
            std.mem.eql(u8, atom.nameSlice(), "H"))
        {
            backbone_h_pos = atom.pos;
        }
    }
    // Backbone H must exist on residue 2
    const h_pos = backbone_h_pos orelse return error.TestUnexpectedResult;

    // Get reference atoms
    const n2 = lookup.findAtomPos(&mdl, mdl.residues.items[1], .{ ' ', 'N', ' ', ' ' }, ' ') orelse unreachable;
    const ca2 = lookup.findAtomPos(&mdl, mdl.residues.items[1], .{ ' ', 'C', 'A', ' ' }, ' ') orelse unreachable;
    const c1 = lookup.findAtomPos(&mdl, mdl.residues.items[0], .{ ' ', 'C', ' ', ' ' }, ' ') orelse unreachable;

    // Check H-N bond length (~1.02 A)
    try testing.expectApproxEqAbs(h_pos.distance(n2), 1.02, 0.05);

    // Check C(i-1)-N-H and CA-N-H angles are approximately equal (bisector)
    const cn_h = math_mod.angle(f32, c1, n2, h_pos);
    const ca_n_h = math_mod.angle(f32, ca2, n2, h_pos);
    try testing.expectApproxEqAbs(cn_h, ca_n_h, 5.0);

    // Both angles should be roughly 119° (peptide plane bisector)
    try testing.expect(cn_h > 110.0 and cn_h < 130.0);

    // H should lie approximately in the C(i-1)-N-CA plane (z ≈ 0 for this fixture)
    try testing.expectApproxEqAbs(h_pos.z, 0.0, 0.1);
}
