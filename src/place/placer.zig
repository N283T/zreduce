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

/// Result of a single executePlan call — why a hydrogen was or wasn't placed.
pub const PlaceResult = enum {
    placed,
    existing_h, // hydrogen already present in residue
    inter_residue, // parent atom bonded_inter_residue (disulfide, glycosidic)
    missing_parent, // parent heavy atom (connected[0]) not found
    missing_ref, // reference neighbor or geometric lookup failed
};

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

pub const WaterConfig = struct {
    enabled: bool = false,
    phantom: bool = false,
    occupancy_cutoff: f32 = 0.66,
    b_factor_cutoff: f32 = 40.0,
    metal_cutoff: f32 = 3.2,
};

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
        if (!found and result.count < 10) {
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

    const n_residues = mdl.residues.items.len;
    for (0..n_residues) |res_idx| {
        const res = mdl.residues.items[res_idx];
        const comp_id = res.compIdSlice();
        const protonation_state = if (config.protonation) |overrides| overrides.find(mdl, res_idx) else null;

        if (config.water.enabled and isWaterResidue(res, comp_id)) {
            const water_result = try placeWaterHydrogens(mdl, res, @intCast(res_idx), config.water, config.bond_policy.mode);
            result.n_placed += water_result.n_placed;
            result.n_skipped_existing += water_result.n_skipped_existing;
            result.n_skipped_inter_residue += water_result.n_skipped_inter_residue;
            result.n_skipped_missing_ref += water_result.n_skipped_missing_ref;
            result.n_skipped_quality_filter += water_result.n_skipped_quality_filter;
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
                if (isBackboneH(&plan)) break true;
            } else false;

            for (targets) |alt| {
                for (plans) |plan| {
                    if (shouldSkipPlanForProtonation(comp_id, &plan, protonation_state)) continue;
                    // Skip backbone amide H only on real N-termini (NH3+/NH2+, not NH).
                    // After a chain break, keep the single amide H.
                    if (is_real_nterm and isBackboneH(&plan)) continue;

                    result.tally(try executePlan(mdl, res, @intCast(res_idx), &plan, bonds, alt, config.bond_policy.mode));
                }

                // Real N-terminal peptide residues: place NH3+/NH2+ instead of single backbone H.
                // Residues after a chain break keep a single break-amide H.
                if (is_real_nterm and has_backbone) {
                    const nterm = if (std.mem.eql(u8, comp_id, "PRO"))
                        try placeNtermNH2Pro(mdl, res, @intCast(res_idx), alt, config.bond_policy.mode)
                    else
                        try placeNtermNH3(mdl, res, @intCast(res_idx), alt, config.bond_policy.mode);
                    result.n_placed += nterm.placed;
                    result.n_skipped_existing += nterm.skipped;
                }

                // 3' terminal nucleotide: place HO3' (leaving atom, absent mid-chain).
                if (is_cterm and nucleotide.getPlans(comp_id) != null) {
                    const oh = try place3primeOH(mdl, res, @intCast(res_idx), alt, config.bond_policy.mode);
                    result.n_placed += oh.placed;
                    result.n_skipped_existing += oh.skipped;
                }

                if (try placeOverrideHydrogen(mdl, res, @intCast(res_idx), alt, protonation_state, config.bond_policy.mode)) |place_result| {
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
                const skip_phosphate_h = !is_nterm_for_bonding and hasPhosphorusAtom(mdl, res);

                for (plans) |plan| {
                    if (skip_phosphate_h and isPhosphateH(plan.h_name)) {
                        result.n_skipped_inter_residue += 1;
                        continue;
                    }
                    result.tally(try executePlan(mdl, res, @intCast(res_idx), &plan, null, ' ', config.bond_policy.mode));
                }

                // 3' terminal nucleotide (CCD path): place HO3' only if CCD
                // plans did not already include it (avoids duplicate placement).
                if (is_cterm and isNucleotideResidue(mdl, res)) {
                    const ho3_name = padName("HO3'");
                    const already_in_plans = for (plans) |plan| {
                        if (std.mem.eql(u8, &plan.h_name, &ho3_name)) break true;
                    } else false;
                    if (!already_in_plans) {
                        const oh = try place3primeOH(mdl, res, @intCast(res_idx), ' ', config.bond_policy.mode);
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
            const has_std_ann = if (chemistry.getAnnotationWithOverride(comp_id, atom.name, protonation_state)) |ann| blk: {
                atom.element_type = ann.atom_type;
                const keep_bonded = atom.flags.bonded_inter_residue;
                atom.flags = ann.flags;
                atom.flags.bonded_inter_residue = keep_bonded;
                atom.vdw_radius = ann.atom_type.info().explicit_radius;
                break :blk true;
            } else false;

            // Apply terminal annotations (merge flags via OR)
            // Only set element_type/vdw_radius for atoms without standard annotation (e.g. OXT)
            if (chemistry.getTerminalAnnotation(atom.name, is_real_nterm, is_cterm)) |term_ann| {
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

fn placeOverrideHydrogen(mdl: *Model, res: Residue, res_idx: u32, target_altloc: u8, state: ?protonation.ResidueState, mode: bond_policy.BondLengthMode) !?PlaceResult {
    const s = state orelse return null;
    const comp_id = res.compIdSlice();

    if (std.mem.eql(u8, comp_id, "ASP") and s == .asp) {
        return switch (s.asp) {
            .deprotonated => null,
            .atom1 => try placeCarboxylHydrogen(mdl, res, res_idx, target_altloc, mode, "OD1", "OD2", "CG", "HD1"),
            .atom2 => try placeCarboxylHydrogen(mdl, res, res_idx, target_altloc, mode, "OD2", "OD1", "CG", "HD2"),
        };
    }
    if (std.mem.eql(u8, comp_id, "GLU") and s == .glu) {
        return switch (s.glu) {
            .deprotonated => null,
            .atom1 => try placeCarboxylHydrogen(mdl, res, res_idx, target_altloc, mode, "OE1", "OE2", "CD", "HE1"),
            .atom2 => try placeCarboxylHydrogen(mdl, res, res_idx, target_altloc, mode, "OE2", "OE1", "CD", "HE2"),
        };
    }
    return null;
}

fn placeCarboxylHydrogen(
    mdl: *Model,
    res: Residue,
    res_idx: u32,
    target_altloc: u8,
    mode: bond_policy.BondLengthMode,
    protonated_o_name: []const u8,
    partner_o_name: []const u8,
    carbon_name: []const u8,
    h_name: []const u8,
) !PlaceResult {
    const protonated_o = findAtom(mdl, res, padName(protonated_o_name), target_altloc) orelse return .missing_parent;
    if (protonated_o.flags.bonded_inter_residue) return .inter_residue;

    var meta = ParentMeta.fromAtom(protonated_o);
    if (target_altloc != ' ') meta.altloc = target_altloc;
    const h_padded = padName(h_name);
    if (existsInResidue(mdl, res, h_padded, meta.altloc)) return .existing_h;

    const o_pos = protonated_o.pos;
    const partner_o_pos = findAtomPos(mdl, res, padName(partner_o_name), target_altloc) orelse return .missing_ref;
    const carbon_pos = findAtomPos(mdl, res, padName(carbon_name), target_altloc) orelse return .missing_ref;
    const h_pos = geometry.placeHXR2Planar(
        o_pos.cast(f64),
        carbon_pos.cast(f64),
        partner_o_pos.cast(f64),
        bond_policy.adjustedBondLength(mode, 0.97, protonated_o.element_type, .Hpol),
        0.0,
    );
    try appendHydrogenNamed(mdl, h_pos.cast(f32), h_name, .Hpol, .none, res_idx, meta);
    return .placed;
}

// ---------------------------------------------------------------------------
// Plan execution
// ---------------------------------------------------------------------------

/// Execute a single placement plan: find reference atoms, compute H position, add to model.
/// Returns true if placed, false if skipped (duplicate H or missing reference atoms).
fn executePlan(mdl: *Model, res: Residue, res_idx: u32, plan: *const standard.PlacementPlan, bonds: ?[]const topology.BondEntry, target_altloc: u8, mode: bond_policy.BondLengthMode) !PlaceResult {
    // Resolve parent heavy atom (connected[0]) for metadata and position
    const base_atom = findAtom(mdl, res, plan.connected[0], target_altloc) orelse return .missing_parent;

    // Skip H placement if parent atom is involved in an inter-residue bond
    // (e.g. disulfide SG, glycosidic leaving O) — the bond already satisfies valence.
    // Only applies to standard/hardcoded plans (bonds != null). CCD-derived plans
    // already account for inter-residue bonds via analyzeBonds extra_bonds, and
    // atoms like glycan C1 may still have free valence for H despite the flag.
    if (bonds != null and base_atom.flags.bonded_inter_residue) return .inter_residue;

    var meta = ParentMeta.fromAtom(base_atom);
    // Override altloc when iterating conformers: if parent has blank altloc
    // (shared backbone) but we're targeting a specific conformer, the placed H
    // should inherit the target conformer's altloc, not blank.
    if (target_altloc != ' ') meta.altloc = target_altloc;

    // Skip if this hydrogen already exists in the residue
    if (existsInResidue(mdl, res, plan.h_name, meta.altloc)) return .existing_h;

    const effective_bond_len = bond_policy.adjustedBondLength(mode, plan.bond_len, base_atom.element_type, plan.atom_type);

    switch (plan.placement_type) {
        .hxr3 => {
            // connected[0]=center, connected[1..2]=two known neighbors
            // Need to find 3rd heavy-atom neighbor of center
            const center_pos = base_atom.pos;
            const n1_pos = findAtomPos(mdl, res, plan.connected[1], target_altloc) orelse return .missing_ref;
            const n2_pos = findAtomPos(mdl, res, plan.connected[2], target_altloc) orelse return .missing_ref;
            const n3_pos = (if (bonds) |b|
                findThirdBondedNeighbor(mdl, res, b, plan.connected[0], plan.connected[1], plan.connected[2], target_altloc)
            else
                findThirdNeighbor(mdl, res, plan.connected[0], plan.connected[1], plan.connected[2], target_altloc)) orelse return .missing_ref;

            const h_pos = geometry.placeHXR3(
                center_pos.cast(f64),
                n1_pos.cast(f64),
                n2_pos.cast(f64),
                n3_pos.cast(f64),
                effective_bond_len,
            );
            try appendHydrogen(mdl, h_pos.cast(f32), plan, res_idx, meta);
            return .placed;
        },

        .h2xr2 => {
            // connected[0]=center, connected[1]=reference neighbor
            const center_pos = base_atom.pos;
            const n1_pos = findAtomPos(mdl, res, plan.connected[1], target_altloc) orelse return .missing_ref;
            const n2_pos = (if (bonds) |b|
                findBondedNeighbor(mdl, res, b, plan.connected[0], plan.connected[1], target_altloc)
            else
                findOtherNeighbor(mdl, res, plan.connected[0], plan.connected[1], target_altloc)) orelse return .missing_ref;

            const h_pos = geometry.placeH2XR2(
                center_pos.cast(f64),
                n1_pos.cast(f64),
                n2_pos.cast(f64),
                effective_bond_len,
                plan.angle,
                plan.dihedral,
            );
            try appendHydrogen(mdl, h_pos.cast(f32), plan, res_idx, meta);
            return .placed;
        },

        .h3xr => {
            // connected[0]=a1 (center), connected[1]=a2, connected[2]=a3
            const a1_pos = base_atom.pos;
            const a2_pos = findAtomPos(mdl, res, plan.connected[1], target_altloc) orelse return .missing_ref;

            // Backbone amide NH: place in peptide plane using C(i-1) from previous residue.
            // H bisects the C(i-1)-N-CA exterior angle, lying in the peptide plane.
            // Only matches " H  " (single amide H), not " H1 "/" H2 "/" H3 " (N-terminal NH3+).
            if (isBackboneAmideH(plan)) {
                if (findPrevResAtomPos(mdl, res_idx, .{ ' ', 'C', ' ', ' ' }, target_altloc)) |prev_c_pos| {
                    const h_pos = geometry.placeHXR2Planar(
                        a1_pos.cast(f64),
                        prev_c_pos.cast(f64),
                        a2_pos.cast(f64),
                        effective_bond_len,
                        0.0,
                    );
                    try appendHydrogen(mdl, h_pos.cast(f32), plan, res_idx, meta);
                    return .placed;
                }
                // No previous C available: first residue, chain break, or missing atom.
                // Fall through to dihedral-based placement as degraded fallback.
            }

            const a3_pos = if (!isBlank(plan.connected[2]))
                findAtomPos(mdl, res, plan.connected[2], target_altloc) orelse return .missing_ref
            else if (bonds) |b|
                findBondedNeighbor(mdl, res, b, plan.connected[1], plan.connected[0], target_altloc) orelse return .missing_ref
            else
                findOtherNeighbor(mdl, res, plan.connected[1], plan.connected[0], target_altloc) orelse return .missing_ref;

            const h_pos = geometry.placeH3XR(
                a1_pos.cast(f64),
                a2_pos.cast(f64),
                a3_pos.cast(f64),
                effective_bond_len,
                plan.angle,
                plan.dihedral,
            );
            try appendHydrogen(mdl, h_pos.cast(f32), plan, res_idx, meta);
            return .placed;
        },

        .hxr2_planar => {
            // connected[0] and connected[1] are neighbors; center is the atom between them
            const n1_pos = base_atom.pos;
            const n2_pos = findAtomPos(mdl, res, plan.connected[1], target_altloc) orelse return .missing_ref;
            const center_pos = (if (bonds) |b|
                findBondedAtomBetween(mdl, res, b, plan.connected[0], plan.connected[1], target_altloc)
            else
                findAtomBetween(mdl, res, plan.connected[0], plan.connected[1], target_altloc)) orelse return .missing_ref;

            const h_pos = geometry.placeHXR2Planar(
                center_pos.cast(f64),
                n1_pos.cast(f64),
                n2_pos.cast(f64),
                effective_bond_len,
                plan.fudge,
            );
            try appendHydrogen(mdl, h_pos.cast(f32), plan, res_idx, meta);
            return .placed;
        },

        .hxr2_frac => {
            const a1_pos = base_atom.pos;
            const a2_pos = findAtomPos(mdl, res, plan.connected[1], target_altloc) orelse return .missing_ref;
            const a3_pos = findAtomPos(mdl, res, plan.connected[2], target_altloc) orelse return .missing_ref;

            const h_pos = geometry.placeHXR2Frac(
                a1_pos.cast(f64),
                a2_pos.cast(f64),
                a3_pos.cast(f64),
                effective_bond_len,
                plan.fudge,
            );
            try appendHydrogen(mdl, h_pos.cast(f32), plan, res_idx, meta);
            return .placed;
        },

        .hxy => {
            const center_pos = base_atom.pos;
            const neighbor_pos = findAtomPos(mdl, res, plan.connected[1], target_altloc) orelse return .missing_ref;

            const h_pos = geometry.placeHXY(
                center_pos.cast(f64),
                neighbor_pos.cast(f64),
                effective_bond_len,
            );
            try appendHydrogen(mdl, h_pos.cast(f32), plan, res_idx, meta);
            return .placed;
        },
    }
}

// ---------------------------------------------------------------------------
// Name comparison helpers
// ---------------------------------------------------------------------------

/// Compare a 4-char PDB-padded name with a model atom's nameSlice().
/// Trims leading/trailing spaces from the PDB name for comparison.
fn nameMatch(pdb_name: [4]u8, atom_name_slice: []const u8) bool {
    // Trim leading spaces
    var start: usize = 0;
    while (start < 4 and pdb_name[start] == ' ') start += 1;
    // Trim trailing spaces
    var end: usize = 4;
    while (end > start and pdb_name[end - 1] == ' ') end -= 1;
    const trimmed_len = end - start;
    if (trimmed_len != atom_name_slice.len) return false;
    for (start..end, 0..) |i, j| {
        if (pdb_name[i] != atom_name_slice[j]) return false;
    }
    return true;
}

/// Check if a 4-char name is blank (all spaces).
fn isBlank(name: [4]u8) bool {
    return name[0] == ' ' and name[1] == ' ' and name[2] == ' ' and name[3] == ' ';
}

// ---------------------------------------------------------------------------
// Atom lookup helpers
// ---------------------------------------------------------------------------

/// Find an atom by name in the previous residue within the same chain.
/// Returns null if there is no previous residue (first in chain or chain break).
fn findPrevResAtomPos(mdl: *const Model, res_idx: u32, name: [4]u8, target_altloc: u8) ?Vec3f32 {
    if (res_idx == 0) return null;
    const cur_res = mdl.residues.items[res_idx];
    if (cur_res.is_chain_break_before) return null;
    const prev_res = mdl.residues.items[res_idx - 1];
    if (prev_res.chain_idx != cur_res.chain_idx) return null;
    return findAtomPos(mdl, prev_res, name, target_altloc);
}

/// Find an atom by 4-char PDB name within a residue. Returns its position.
/// When target_altloc is specified, prefers an atom with matching altloc,
/// but falls back to an atom with blank (' ') altloc if no exact match.
fn findAtomPos(mdl: *const Model, res: Residue, name: [4]u8, target_altloc: u8) ?Vec3f32 {
    if (isBlank(name)) return null;
    const atoms = mdl.atoms.items[res.atom_start..res.atom_end];
    var blank_match: ?Vec3f32 = null;
    for (atoms) |a| {
        if (nameMatch(name, a.nameSlice())) {
            if (a.altloc == target_altloc) return a.pos;
            if (a.altloc == ' ' and blank_match == null) blank_match = a.pos;
        }
    }
    // Note: we only search the original residue atom range. Newly appended
    // H atoms are not referenced by plans (plans only reference heavy atoms).
    return blank_match;
}

/// Find an atom by 4-char PDB name within a residue. Returns the full Atom.
/// When target_altloc is specified, prefers an atom with matching altloc,
/// but falls back to an atom with blank (' ') altloc if no exact match.
fn findAtom(mdl: *const Model, res: Residue, name: [4]u8, target_altloc: u8) ?Atom {
    if (isBlank(name)) return null;
    const atoms = mdl.atoms.items[res.atom_start..res.atom_end];
    var blank_match: ?Atom = null;
    for (atoms) |a| {
        if (nameMatch(name, a.nameSlice())) {
            if (a.altloc == target_altloc) return a;
            if (a.altloc == ' ' and blank_match == null) blank_match = a;
        }
    }
    return blank_match;
}

/// Check if an atom with the given name and altloc already exists in a residue.
/// Used to prevent duplicate hydrogen placement.
/// Only searches the original residue atom range, not newly appended H atoms.
/// altloc ' ' (blank) matches only blank; 'A' matches only 'A'.
fn existsInResidue(mdl: *const Model, res: Residue, name: [4]u8, altloc: u8) bool {
    const atoms = mdl.atoms.items[res.atom_start..res.atom_end];
    for (atoms) |a| {
        if (nameMatch(name, a.nameSlice()) and a.altloc == altloc) {
            return true;
        }
    }
    return false;
}

/// Find a heavy-atom neighbor of `center_name` that is NOT `exclude_name`.
/// Uses distance-based bonding (within 1.9 A of center).
fn findOtherNeighbor(mdl: *const Model, res: Residue, center_name: [4]u8, exclude_name: [4]u8, target_altloc: u8) ?Vec3f32 {
    const center_pos = findAtomPos(mdl, res, center_name, target_altloc) orelse return null;
    const bond_cutoff: f32 = 1.9;

    const atoms = mdl.atoms.items[res.atom_start..res.atom_end];
    for (atoms) |a| {
        if (a.is_hydrogen) continue;
        if (a.altloc != target_altloc and a.altloc != ' ') continue;
        const aname = a.nameSlice();
        if (nameMatch(center_name, aname)) continue;
        if (nameMatch(exclude_name, aname)) continue;
        if (a.pos.distance(center_pos) < bond_cutoff) {
            return a.pos;
        }
    }
    return null;
}

/// Find the 3rd heavy-atom neighbor of center (for HXR3 placement).
/// Excludes the two already-known neighbors.
fn findThirdNeighbor(mdl: *const Model, res: Residue, center_name: [4]u8, n1_name: [4]u8, n2_name: [4]u8, target_altloc: u8) ?Vec3f32 {
    const center_pos = findAtomPos(mdl, res, center_name, target_altloc) orelse return null;
    const bond_cutoff: f32 = 1.9;

    const atoms = mdl.atoms.items[res.atom_start..res.atom_end];
    for (atoms) |a| {
        if (a.is_hydrogen) continue;
        if (a.altloc != target_altloc and a.altloc != ' ') continue;
        const aname = a.nameSlice();
        if (nameMatch(center_name, aname)) continue;
        if (nameMatch(n1_name, aname)) continue;
        if (nameMatch(n2_name, aname)) continue;
        if (a.pos.distance(center_pos) < bond_cutoff) {
            return a.pos;
        }
    }
    return null;
}

/// Find an atom that is bonded to BOTH named atoms (the atom "between" them).
/// Used for planar placement where the center atom is implicit.
fn findAtomBetween(mdl: *const Model, res: Residue, name1: [4]u8, name2: [4]u8, target_altloc: u8) ?Vec3f32 {
    const pos1 = findAtomPos(mdl, res, name1, target_altloc) orelse return null;
    const pos2 = findAtomPos(mdl, res, name2, target_altloc) orelse return null;
    const bond_cutoff: f32 = 1.9;

    const atoms = mdl.atoms.items[res.atom_start..res.atom_end];
    for (atoms) |a| {
        if (a.is_hydrogen) continue;
        if (a.altloc != target_altloc and a.altloc != ' ') continue;
        const aname = a.nameSlice();
        if (nameMatch(name1, aname)) continue;
        if (nameMatch(name2, aname)) continue;
        if (a.pos.distance(pos1) < bond_cutoff and a.pos.distance(pos2) < bond_cutoff) {
            return a.pos;
        }
    }
    return null;
}

// ---------------------------------------------------------------------------
// Bond-aware neighbor queries
// ---------------------------------------------------------------------------

/// Given a bond entry, return the partner of `atom_name`, or null if not involved.
fn bondPartner(bond: topology.BondEntry, atom_name: [4]u8) ?[4]u8 {
    if (std.mem.eql(u8, &bond.a1, &atom_name)) return bond.a2;
    if (std.mem.eql(u8, &bond.a2, &atom_name)) return bond.a1;
    return null;
}

/// Find a bonded neighbor of `center_name` that is NOT `exclude_name`, using topology.
/// Returns the first match; result depends on bond table ordering.
fn findBondedNeighbor(
    mdl: *const Model,
    res: Residue,
    bonds: []const topology.BondEntry,
    center_name: [4]u8,
    exclude_name: [4]u8,
    target_altloc: u8,
) ?Vec3f32 {
    for (bonds) |bond| {
        const partner = bondPartner(bond, center_name) orelse continue;
        if (std.mem.eql(u8, &partner, &exclude_name)) continue;
        if (findAtomPos(mdl, res, partner, target_altloc)) |pos| return pos;
    }
    return null;
}

/// Find the 3rd bonded neighbor of `center_name`, excluding `n1_name` and `n2_name`.
fn findThirdBondedNeighbor(
    mdl: *const Model,
    res: Residue,
    bonds: []const topology.BondEntry,
    center_name: [4]u8,
    n1_name: [4]u8,
    n2_name: [4]u8,
    target_altloc: u8,
) ?Vec3f32 {
    for (bonds) |bond| {
        const partner = bondPartner(bond, center_name) orelse continue;
        if (std.mem.eql(u8, &partner, &n1_name)) continue;
        if (std.mem.eql(u8, &partner, &n2_name)) continue;
        if (findAtomPos(mdl, res, partner, target_altloc)) |pos| return pos;
    }
    return null;
}

/// Find an atom bonded to BOTH `name1` and `name2` using topology.
fn findBondedAtomBetween(
    mdl: *const Model,
    res: Residue,
    bonds: []const topology.BondEntry,
    name1: [4]u8,
    name2: [4]u8,
    target_altloc: u8,
) ?Vec3f32 {
    // Collect atoms bonded to name1 (max 8 — no standard AA atom has more)
    var bonded_to_1: [8][4]u8 = undefined;
    var count_1: usize = 0;
    for (bonds) |bond| {
        const partner = bondPartner(bond, name1) orelse continue;
        if (count_1 < bonded_to_1.len) {
            bonded_to_1[count_1] = partner;
            count_1 += 1;
        }
    }

    // Check which are also bonded to name2
    for (bonds) |bond| {
        const partner = bondPartner(bond, name2) orelse continue;
        for (bonded_to_1[0..count_1]) |candidate| {
            if (std.mem.eql(u8, &partner, &candidate)) {
                if (findAtomPos(mdl, res, partner, target_altloc)) |pos| return pos;
            }
        }
    }
    return null;
}

/// Metadata inherited from the parent heavy atom.
const ParentMeta = struct {
    altloc: u8 = ' ',
    occupancy: f32 = 1.0,
    b_factor: f32 = 0.0,

    /// Extract metadata from an atom.
    fn fromAtom(a: Atom) ParentMeta {
        return .{
            .altloc = a.altloc,
            .occupancy = a.occupancy,
            .b_factor = a.b_factor,
        };
    }
};

/// Pad a short atom name (e.g. "H1") to a 4-char PDB-padded name.
/// name must be at most 4 bytes.
fn padName(name: []const u8) [4]u8 {
    std.debug.assert(name.len <= 4);
    var padded: [4]u8 = .{ ' ', ' ', ' ', ' ' };
    for (name, 0..) |c, i| padded[i] = c;
    return padded;
}

/// Check if a placement plan is for a backbone amide H (single NH or N-terminal H1/H2/H3).
/// Used by N-terminal skip logic — matches " H  ", " H1 ", " H2 ", " H3 " but only
/// when the parent atom (connected[0]) is backbone nitrogen " N  ".
/// This prevents false matches on nucleotide ring H (e.g. guanine H1 on N1 via C6).
fn isBackboneH(plan: *const standard.PlacementPlan) bool {
    const parent = plan.connected[0];
    if (!(parent[0] == ' ' and parent[1] == 'N' and parent[2] == ' ' and parent[3] == ' ')) return false;
    const h = plan.h_name;
    if (h[0] == ' ' and h[1] == 'H' and h[2] == ' ' and h[3] == ' ') return true;
    if (h[0] == ' ' and h[1] == 'H' and h[2] == '1' and h[3] == ' ') return true;
    if (h[0] == ' ' and h[1] == 'H' and h[2] == '2' and h[3] == ' ') return true;
    if (h[0] == ' ' and h[1] == 'H' and h[2] == '3' and h[3] == ' ') return true;
    return false;
}

/// Check if a plan is specifically the single backbone amide H (" H  ").
/// Used by peptide-plane placement — must NOT match N-terminal H1/H2/H3.
fn isBackboneAmideH(plan: *const standard.PlacementPlan) bool {
    const h = plan.h_name;
    return h[0] == ' ' and h[1] == 'H' and h[2] == ' ' and h[3] == ' ';
}

/// Append an N-terminal H atom to the model.
fn appendNtermH(mdl: *Model, h_pos: Vec3f32, name: []const u8, res_idx: u32, meta: ParentMeta) !void {
    const hpol_info = element.AtomType.Hpol.info();
    var atom = Atom{
        .pos = h_pos,
        .element_type = .Hpol,
        .residue_idx = res_idx,
        .is_hydrogen = true,
        .is_added = true,
        .vdw_radius = hpol_info.explicit_radius,
        .altloc = meta.altloc,
        .occupancy = meta.occupancy,
        .b_factor = meta.b_factor,
        .flags = .{ .donor = true },
    };
    atom.setName(name);
    try mdl.atoms.append(mdl.allocator, atom);
}

const NtermResult = struct { placed: u32, skipped: u32 };

/// Place NH3+ hydrogens (H1, H2, H3) on the N-terminal residue.
/// Uses h3xr (dihedral-controlled) placement around the N-CA bond.
fn placeNtermNH3(mdl: *Model, res: Residue, res_idx: u32, target_altloc: u8, mode: bond_policy.BondLengthMode) !NtermResult {
    const n_atom = findAtom(mdl, res, .{ ' ', 'N', ' ', ' ' }, target_altloc) orelse return .{ .placed = 0, .skipped = 0 };
    const ca_pos = findAtomPos(mdl, res, .{ ' ', 'C', 'A', ' ' }, target_altloc) orelse return .{ .placed = 0, .skipped = 0 };
    const c_pos = findAtomPos(mdl, res, .{ ' ', 'C', ' ', ' ' }, target_altloc) orelse return .{ .placed = 0, .skipped = 0 };

    var meta = ParentMeta.fromAtom(n_atom);
    if (target_altloc != ' ') meta.altloc = target_altloc;

    const n64 = n_atom.pos.cast(f64);
    const ca64 = math_mod.Vec3(f64){ .x = ca_pos.x, .y = ca_pos.y, .z = ca_pos.z };
    const c64 = math_mod.Vec3(f64){ .x = c_pos.x, .y = c_pos.y, .z = c_pos.z };

    const bond_len: f64 = bond_policy.adjustedBondLength(mode, 1.00, n_atom.element_type, .Hpol);
    const angle_deg: f64 = 109.5;
    const dihedrals = [3]f64{ 180.0, 60.0, -60.0 };
    const names = [3][]const u8{ "H1", "H2", "H3" };

    var placed: u32 = 0;
    var skipped: u32 = 0;
    for (names, dihedrals) |name, dihedral| {
        if (existsInResidue(mdl, res, padName(name), meta.altloc)) {
            skipped += 1;
            continue;
        }

        const h64 = geometry.placeH3XR(n64, ca64, c64, bond_len, angle_deg, dihedral);
        const h_pos = Vec3f32{ .x = @floatCast(h64.x), .y = @floatCast(h64.y), .z = @floatCast(h64.z) };
        try appendNtermH(mdl, h_pos, name, res_idx, meta);
        placed += 1;
    }
    return .{ .placed = placed, .skipped = skipped };
}

/// Place NH2+ hydrogens (H2, H3) on N-terminal PRO.
/// PRO has CD bonded to N (secondary amine), so only 2 H positions.
fn placeNtermNH2Pro(mdl: *Model, res: Residue, res_idx: u32, target_altloc: u8, mode: bond_policy.BondLengthMode) !NtermResult {
    const n_atom = findAtom(mdl, res, .{ ' ', 'N', ' ', ' ' }, target_altloc) orelse return .{ .placed = 0, .skipped = 0 };
    const ca_pos = findAtomPos(mdl, res, .{ ' ', 'C', 'A', ' ' }, target_altloc) orelse return .{ .placed = 0, .skipped = 0 };
    const cd_pos = findAtomPos(mdl, res, .{ ' ', 'C', 'D', ' ' }, target_altloc) orelse return .{ .placed = 0, .skipped = 0 };

    var meta = ParentMeta.fromAtom(n_atom);
    if (target_altloc != ' ') meta.altloc = target_altloc;

    // PRO N is sp3 with 2 heavy-atom neighbors (CA, CD) and 2 H.
    // Use h2xr2 (two H on atom with 2 heavy neighbors).
    const n64 = n_atom.pos.cast(f64);
    const ca64 = math_mod.Vec3(f64){ .x = ca_pos.x, .y = ca_pos.y, .z = ca_pos.z };
    const cd64 = math_mod.Vec3(f64){ .x = cd_pos.x, .y = cd_pos.y, .z = cd_pos.z };

    const bond_len: f64 = bond_policy.adjustedBondLength(mode, 1.00, n_atom.element_type, .Hpol);
    const angle_deg: f64 = 109.5;
    const names = [2][]const u8{ "H2", "H3" };
    const dihedrals = [2]f64{ 120.0, -120.0 };

    var placed: u32 = 0;
    var skipped: u32 = 0;
    for (names, dihedrals) |name, dihedral| {
        if (existsInResidue(mdl, res, padName(name), meta.altloc)) {
            skipped += 1;
            continue;
        }

        const h64 = geometry.placeH2XR2(n64, ca64, cd64, bond_len, angle_deg, dihedral);
        const h_pos = Vec3f32{ .x = @floatCast(h64.x), .y = @floatCast(h64.y), .z = @floatCast(h64.z) };
        try appendNtermH(mdl, h_pos, name, res_idx, meta);
        placed += 1;
    }
    return .{ .placed = placed, .skipped = skipped };
}

/// Place HO3' on 3' terminal nucleotide residues.
/// O3' is sp3 with 2 heavy-atom neighbors (C3', and normally the next P — absent at 3' terminus).
/// Uses h3xr geometry: O3' center, C3' and C4' as references, tetrahedral angle.
fn place3primeOH(mdl: *Model, res: Residue, res_idx: u32, target_altloc: u8, mode: bond_policy.BondLengthMode) !NtermResult {
    const o3_atom = findAtom(mdl, res, .{ ' ', 'O', '3', '\'' }, target_altloc) orelse return .{ .placed = 0, .skipped = 0 };
    const c3_pos = findAtomPos(mdl, res, .{ ' ', 'C', '3', '\'' }, target_altloc) orelse return .{ .placed = 0, .skipped = 0 };
    const c4_pos = findAtomPos(mdl, res, .{ ' ', 'C', '4', '\'' }, target_altloc) orelse return .{ .placed = 0, .skipped = 0 };

    var meta = ParentMeta.fromAtom(o3_atom);
    if (target_altloc != ' ') meta.altloc = target_altloc;

    const ho3_name = padName("HO3'");
    if (existsInResidue(mdl, res, ho3_name, meta.altloc)) {
        return .{ .placed = 0, .skipped = 1 };
    }

    const o3_64 = o3_atom.pos.cast(f64);
    const c3_64 = math_mod.Vec3(f64){ .x = c3_pos.x, .y = c3_pos.y, .z = c3_pos.z };
    const c4_64 = math_mod.Vec3(f64){ .x = c4_pos.x, .y = c4_pos.y, .z = c4_pos.z };

    const bond_len: f64 = bond_policy.adjustedBondLength(mode, 0.97, o3_atom.element_type, .Hpol);
    const h64 = geometry.placeH3XR(o3_64, c3_64, c4_64, bond_len, 109.5, 180.0);
    const h_pos = Vec3f32{ .x = @floatCast(h64.x), .y = @floatCast(h64.y), .z = @floatCast(h64.z) };

    const hpol_info = element.AtomType.Hpol.info();
    var atom = Atom{
        .pos = h_pos,
        .element_type = .Hpol,
        .residue_idx = res_idx,
        .is_hydrogen = true,
        .is_added = true,
        .vdw_radius = hpol_info.explicit_radius,
        .altloc = meta.altloc,
        .occupancy = meta.occupancy,
        .b_factor = meta.b_factor,
        .flags = .{ .donor = true },
        .mover_hint = .rotate,
    };
    atom.setName("HO3'");
    try mdl.atoms.append(mdl.allocator, atom);
    return .{ .placed = 1, .skipped = 0 };
}

/// Check if a residue has a nucleotide sugar-phosphate backbone (has O3' atom).
fn isNucleotideResidue(mdl: *const Model, res: Residue) bool {
    const atoms = mdl.atoms.items[res.atom_start..res.atom_end];
    for (atoms) |a| {
        if (nameMatch(.{ ' ', 'O', '3', '\'' }, a.nameSlice())) return true;
    }
    return false;
}

/// Check if a hydrogen name is a phosphate H (HOP2, HOP3).
/// In CCD, these are H atoms on phosphate oxygens OP2/OP3.
fn isPhosphateH(h_name: [4]u8) bool {
    return (h_name[0] == 'H' and h_name[1] == 'O' and h_name[2] == 'P');
}

/// Check if a residue has atom named "P" (phosphorus — nucleotide backbone).
/// Matches both setName("P") format {'P',' ',' ',' '} and PDB-padded {' ','P',' ',' '}.
fn hasPhosphorusAtom(mdl: *const Model, res: Residue) bool {
    const atoms = mdl.atoms.items[res.atom_start..res.atom_end];
    for (atoms) |a| {
        if (nameMatch(.{ ' ', 'P', ' ', ' ' }, a.nameSlice())) return true;
    }
    return false;
}

fn isWaterResidue(res: Residue, comp_id: []const u8) bool {
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

fn chooseWaterNeighbors(mdl: *const Model, res_idx: u32, oxygen: Atom, target_altloc: u8) [2]?Vec3f32 {
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

fn nearestMetalDistance(mdl: *const Model, res_idx: u32, oxygen: Atom, target_altloc: u8) ?f32 {
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

fn awayFromNearbyAtoms(mdl: *const Model, res_idx: u32, oxygen: Atom, target_altloc: u8) ?Vec3f32 {
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

fn placeWaterHydrogens(mdl: *Model, res: Residue, res_idx: u32, config: WaterConfig, mode: bond_policy.BondLengthMode) !PlacementResult {
    var result = PlacementResult{};
    const altlocs = collectAltlocs(mdl, res);
    const targets: []const u8 = if (altlocs.count == 0)
        &[_]u8{' '}
    else
        altlocs.locs[0..altlocs.count];

    const h1_name = padName("H1");
    const h2_name = padName("H2");
    const half_angle: f32 = 52.25;

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
        if (nearestMetalDistance(mdl, res_idx, oxygen, alt)) |dist| {
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

        const neighbors = chooseWaterNeighbors(mdl, res_idx, oxygen, alt);
        const n1 = neighbors[0];
        const n2 = neighbors[1];

        var away = if (n1) |p1|
            oxygen.pos.sub(p1).normalize()
        else if (config.phantom)
            awayFromNearbyAtoms(mdl, res_idx, oxygen, alt) orelse Vec3f32{ .x = 1.0, .y = 0.0, .z = 0.0 }
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

        try appendNtermH(mdl, h1, "H1", res_idx, meta);
        try appendNtermH(mdl, h2, "H2", res_idx, meta);
        result.n_placed += 2;
    }

    return result;
}

/// Find the raw [4]u8 name of an atom in the model (using nameMatch for lookup).
/// Returns the actual stored name bytes, which may differ from PDB-padded format.
fn findAtomName(mdl: *const Model, res: Residue, name: [4]u8, target_altloc: u8) ?[4]u8 {
    const a = findAtom(mdl, res, name, target_altloc) orelse return null;
    return a.name;
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
        result[i] = a.name;
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
                if (!atom.is_hydrogen) try names.append(allocator, atom.name);
            }
            if (bond.atom_2 >= res.atom_start and bond.atom_2 < res.atom_end) {
                const atom = mdl.atoms.items[bond.atom_2];
                if (!atom.is_hydrogen) try names.append(allocator, atom.name);
            }
        }
    }
    return try allocator.dupe([4]u8, names.items);
}

/// Append a new hydrogen atom to the model.
fn appendHydrogen(mdl: *Model, pos: Vec3f32, plan: *const standard.PlacementPlan, res_idx: u32, meta: ParentMeta) !void {
    try appendHydrogenNamed(mdl, pos, trimPlanName(&plan.h_name), plan.atom_type, plan.mover_hint, res_idx, meta);
}

fn appendHydrogenNamed(
    mdl: *Model,
    pos: Vec3f32,
    name: []const u8,
    atom_type: element.AtomType,
    mover_hint: standard.MoverHint,
    res_idx: u32,
    meta: ParentMeta,
) !void {
    var atom = Atom{
        .pos = pos,
        .element_type = atom_type,
        .residue_idx = res_idx,
        .is_hydrogen = true,
        .is_added = true,
        .vdw_radius = atom_type.info().explicit_radius,
        .flags = atom_type.info().flags,
        .altloc = meta.altloc,
        .occupancy = meta.occupancy,
        .b_factor = meta.b_factor,
        .mover_hint = mover_hint,
    };
    atom.setName(name);
    try mdl.atoms.append(mdl.allocator, atom);
}

fn trimPlanName(name: *const [4]u8) []const u8 {
    var start: usize = 0;
    while (start < 4 and name[start] == ' ') start += 1;
    var end: usize = 4;
    while (end > start and name[end - 1] == ' ') end -= 1;
    return name[start..end];
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

const testing = std.testing;
const mmcif = @import("../mmcif.zig");

test "nameMatch trims PDB-style names correctly" {
    try testing.expect(nameMatch(.{ ' ', 'N', ' ', ' ' }, "N"));
    try testing.expect(nameMatch(.{ ' ', 'C', 'A', ' ' }, "CA"));
    try testing.expect(nameMatch(.{ 'H', 'G', '1', '1' }, "HG11"));
    try testing.expect(!nameMatch(.{ ' ', 'N', ' ', ' ' }, "CA"));
}

test "isBlank detects blank names" {
    try testing.expect(isBlank(.{ ' ', ' ', ' ', ' ' }));
    try testing.expect(!isBlank(.{ 'N', ' ', ' ', ' ' }));
}

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

test "findAtom returns full atom with metadata" {
    const source = @embedFile("../test_data/tiny.cif");
    var mdl = try mmcif.parseModel(testing.allocator, source);
    defer mdl.deinit();

    const res = mdl.residues.items[0];

    // CA should be found with correct metadata
    const ca = findAtom(&mdl, res, .{ ' ', 'C', 'A', ' ' }, ' ');
    try testing.expect(ca != null);
    try testing.expectEqual(@as(f32, 1.0), ca.?.occupancy);
    try testing.expectEqual(@as(f32, 10.0), ca.?.b_factor);
    try testing.expectEqual(@as(u8, ' '), ca.?.altloc);

    // Non-existent atom returns null
    const xx = findAtom(&mdl, res, .{ ' ', 'X', 'X', ' ' }, ' ');
    try testing.expect(xx == null);
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

test "existsInResidue checks name and altloc" {
    const source = @embedFile("../test_data/tiny.cif");
    var mdl = try mmcif.parseModel(testing.allocator, source);
    defer mdl.deinit();

    const res = mdl.residues.items[0];

    // "N" with altloc ' ' exists in tiny.cif
    try testing.expect(existsInResidue(&mdl, res, .{ ' ', 'N', ' ', ' ' }, ' '));
    // "CA" with altloc ' ' exists
    try testing.expect(existsInResidue(&mdl, res, .{ ' ', 'C', 'A', ' ' }, ' '));
    // "HA" does not exist
    try testing.expect(!existsInResidue(&mdl, res, .{ ' ', 'H', 'A', ' ' }, ' '));
    // "N" with altloc 'A' does not exist (tiny.cif atoms have altloc=' ')
    try testing.expect(!existsInResidue(&mdl, res, .{ ' ', 'N', ' ', ' ' }, 'A'));
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

test "findBondedNeighbor returns correct neighbor from topology" {
    const source = @embedFile("../test_data/tiny.cif");
    var mdl = try mmcif.parseModel(testing.allocator, source);
    defer mdl.deinit();

    const res = mdl.residues.items[0];
    const bonds = topology.getBonds("ALA").?;

    // CA is bonded to N, C, CB. Excluding N, should find C or CB.
    const result = findBondedNeighbor(&mdl, res, bonds, .{ ' ', 'C', 'A', ' ' }, .{ ' ', 'N', ' ', ' ' }, ' ');
    try testing.expect(result != null);

    // Verify it's C or CB (both bonded to CA, excluding N)
    const c_pos = findAtomPos(&mdl, res, .{ ' ', 'C', ' ', ' ' }, ' ');
    const cb_pos = findAtomPos(&mdl, res, .{ ' ', 'C', 'B', ' ' }, ' ');
    const pos = result.?;
    const is_c = c_pos != null and pos.x == c_pos.?.x and pos.y == c_pos.?.y and pos.z == c_pos.?.z;
    const is_cb = cb_pos != null and pos.x == cb_pos.?.x and pos.y == cb_pos.?.y and pos.z == cb_pos.?.z;
    try testing.expect(is_c or is_cb);
}

test "findThirdBondedNeighbor finds the third bonded atom" {
    const source = @embedFile("../test_data/tiny.cif");
    var mdl = try mmcif.parseModel(testing.allocator, source);
    defer mdl.deinit();

    const res = mdl.residues.items[0];
    const bonds = topology.getBonds("ALA").?;

    // CA bonded to N, C, CB. Excluding N and C → CB.
    const result = findThirdBondedNeighbor(&mdl, res, bonds, .{ ' ', 'C', 'A', ' ' }, .{ ' ', 'N', ' ', ' ' }, .{ ' ', 'C', ' ', ' ' }, ' ');
    try testing.expect(result != null);
    const cb_pos = findAtomPos(&mdl, res, .{ ' ', 'C', 'B', ' ' }, ' ').?;
    try testing.expectEqual(cb_pos.x, result.?.x);
    try testing.expectEqual(cb_pos.y, result.?.y);
    try testing.expectEqual(cb_pos.z, result.?.z);
}

test "findBondedAtomBetween finds atom bonded to both" {
    const source = @embedFile("../test_data/tiny.cif");
    var mdl = try mmcif.parseModel(testing.allocator, source);
    defer mdl.deinit();

    const res = mdl.residues.items[0];
    const bonds = topology.getBonds("ALA").?;

    // N and C are both bonded to CA
    const result = findBondedAtomBetween(&mdl, res, bonds, .{ ' ', 'N', ' ', ' ' }, .{ ' ', 'C', ' ', ' ' }, ' ');
    try testing.expect(result != null);
    const ca_pos = findAtomPos(&mdl, res, .{ ' ', 'C', 'A', ' ' }, ' ').?;
    try testing.expectEqual(ca_pos.x, result.?.x);
    try testing.expectEqual(ca_pos.y, result.?.y);
    try testing.expectEqual(ca_pos.z, result.?.z);
}

test "bond-based query finds stretched CB that distance-based misses" {
    const source = @embedFile("../test_data/ala_stretched.cif");
    var mdl = try mmcif.parseModel(testing.allocator, source);
    defer mdl.deinit();

    const res = mdl.residues.items[0];
    const bonds = topology.getBonds("ALA").?;

    // Distance-based: CB is >1.9A from CA, should NOT be found
    const dist_result = findThirdNeighbor(&mdl, res, .{ ' ', 'C', 'A', ' ' }, .{ ' ', 'N', ' ', ' ' }, .{ ' ', 'C', ' ', ' ' }, ' ');
    try testing.expect(dist_result == null);

    // Bond-based: CB is bonded to CA in topology, SHOULD be found
    const bond_result = findThirdBondedNeighbor(&mdl, res, bonds, .{ ' ', 'C', 'A', ' ' }, .{ ' ', 'N', ' ', ' ' }, .{ ' ', 'C', ' ', ' ' }, ' ');
    try testing.expect(bond_result != null);
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

test "findAtomPos with altloc prefers matching conformer" {
    const source = @embedFile("../test_data/ala_altloc.cif");
    var mdl = try mmcif.parseModel(testing.allocator, source);
    defer mdl.deinit();

    const res = mdl.residues.items[0];

    // CA with altloc 'A' should return conformer A position (0, 0, 0)
    const ca_a = findAtomPos(&mdl, res, .{ ' ', 'C', 'A', ' ' }, 'A');
    try testing.expect(ca_a != null);
    try testing.expectApproxEqAbs(@as(f32, 0.0), ca_a.?.x, 1e-3);

    // CA with altloc 'B' should return conformer B position (0.1, 0.1, 0.1)
    const ca_b = findAtomPos(&mdl, res, .{ ' ', 'C', 'A', ' ' }, 'B');
    try testing.expect(ca_b != null);
    try testing.expectApproxEqAbs(@as(f32, 0.1), ca_b.?.x, 1e-3);

    // N has blank altloc — should be found by any target_altloc (fallback)
    const n_a = findAtomPos(&mdl, res, .{ ' ', 'N', ' ', ' ' }, 'A');
    try testing.expect(n_a != null);
    const n_b = findAtomPos(&mdl, res, .{ ' ', 'N', ' ', ' ' }, 'B');
    try testing.expect(n_b != null);
    try testing.expectApproxEqAbs(n_a.?.x, n_b.?.x, 1e-6);
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

    const oxygen_pos = findAtomPos(&mdl, mdl.residues.items[1], padName("O"), ' ') orelse unreachable;
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
    const n2 = findAtomPos(&mdl, mdl.residues.items[1], .{ ' ', 'N', ' ', ' ' }, ' ') orelse unreachable;
    const ca2 = findAtomPos(&mdl, mdl.residues.items[1], .{ ' ', 'C', 'A', ' ' }, ' ') orelse unreachable;
    const c1 = findAtomPos(&mdl, mdl.residues.items[0], .{ ' ', 'C', ' ', ' ' }, ' ') orelse unreachable;

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
