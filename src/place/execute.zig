//! Plan execution — core geometry dispatch for hydrogen placement.
//!
//! Translates placement plans into atom positions and appends
//! new hydrogen atoms to the model.

const std = @import("std");
const model_mod = @import("../model.zig");
const Model = model_mod.Model;
const Atom = model_mod.Atom;
const Residue = model_mod.Residue;
const math_mod = @import("../math.zig");
const element = @import("../element.zig");
const bond_policy = @import("bond_policy.zig");
const geometry = @import("geometry.zig");
const standard = @import("standard.zig");
const topology = @import("topology.zig");
const protonation = @import("protonation.zig");
const lookup = @import("lookup.zig");
const terminal = @import("terminal.zig");

const Vec3f32 = math_mod.Vec3(f32);
const ParentMeta = lookup.ParentMeta;
const findAtom = lookup.findAtom;
const findAtomPos = lookup.findAtomPos;
const existsInResidue = lookup.existsInResidue;
const findOtherNeighbor = lookup.findOtherNeighbor;
const findThirdNeighbor = lookup.findThirdNeighbor;
const findAtomBetween = lookup.findAtomBetween;
const findBondedNeighbor = lookup.findBondedNeighbor;
const findThirdBondedNeighbor = lookup.findThirdBondedNeighbor;
const findBondedAtomBetween = lookup.findBondedAtomBetween;
const padName = lookup.padName;
const trimPlanName = lookup.trimPlanName;
const isBackboneAmideH = terminal.isBackboneAmideH;
const findPrevResAtomPos = lookup.findPrevResAtomPos;

/// Result of a single executePlan call — why a hydrogen was or wasn't placed.
/// Defined here to avoid circular dependency with placer.zig.
pub const PlaceResult = enum {
    placed,
    existing_h, // hydrogen already present in residue
    inter_residue, // parent atom bonded_inter_residue (disulfide, glycosidic)
    missing_parent, // parent heavy atom (connected[0]) not found
    missing_ref, // reference neighbor or geometric lookup failed
};

/// Append a new hydrogen atom to the model.
pub fn appendHydrogen(mdl: *Model, pos: Vec3f32, plan: *const standard.PlacementPlan, res_idx: u32, meta: ParentMeta) !void {
    try appendHydrogenNamed(mdl, pos, trimPlanName(&plan.h_name), plan.atom_type, plan.mover_hint, res_idx, meta);
}

pub fn appendHydrogenNamed(
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

pub fn placeOverrideHydrogen(mdl: *Model, res: Residue, res_idx: u32, target_altloc: u8, state: ?protonation.ResidueState, mode: bond_policy.BondLengthMode) !?PlaceResult {
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

/// Execute a single placement plan: find reference atoms, compute H position, add to model.
/// Returns the placement result (placed / skipped / missing ref).
pub fn executePlan(mdl: *Model, res: Residue, res_idx: u32, plan: *const standard.PlacementPlan, bonds: ?[]const topology.BondEntry, target_altloc: u8, mode: bond_policy.BondLengthMode) !PlaceResult {
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

            const a3_pos = if (!lookup.isBlank(plan.connected[2]))
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
