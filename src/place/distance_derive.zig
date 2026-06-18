//! Distance-based bond inferrer for residues with no CCD dictionary entry.
//!
//! Produces a synthetic `ccd.Component` from model atom coordinates so that the
//! existing `ccd_derive.derivePlans()` pipeline can still generate hydrogen
//! placement plans for unknown ligands.
//!
//! Algorithm:
//!   1. Collect non-hydrogen, non-metal heavy atoms from the residue.
//!   2. Detect bonds by covalent-radius sum with a 1.3× tolerance and a 2.5 Å cap.
//!   3. Promote single bonds to double where both atoms are under-valent.
//!   4. Return a synthetic Component with ideal coordinates taken from actual positions.

const std = @import("std");
const ccd = @import("../ccd.zig");
const model_mod = @import("../model.zig");
const element_mod = @import("../element.zig");
const math_mod = @import("../math.zig");

const Model = model_mod.Model;
const Atom = model_mod.Atom;
const Residue = model_mod.Residue;
const AtomType = element_mod.AtomType;
const Component = ccd.Component;
const CompAtom = ccd.CompAtom;
const CompBond = ccd.CompBond;
const BondOrder = ccd.BondOrder;

// ---------------------------------------------------------------------------
// Constants
// ---------------------------------------------------------------------------

/// Tolerance multiplier on the sum of covalent radii.
const COV_TOLERANCE: f32 = 1.3;

/// Hard upper limit on bond distance regardless of radii.
const MAX_BOND_DIST: f32 = 2.5;

// ---------------------------------------------------------------------------
// Public API
// ---------------------------------------------------------------------------

/// Derive a synthetic Component from model atom coordinates for a residue with
/// no CCD dictionary entry. Uses distance-based bond detection and valence-rule
/// bond order inference.
///
/// Returns null if the residue has fewer than 2 non-hydrogen, non-metal heavy
/// atoms (nothing useful to infer).
///
/// The returned Component's `atoms`, `bonds`, `comp_id`, and `comp_type` slices
/// are all allocated with `allocator`. The caller owns the memory and must free
/// each slice independently.
pub fn deriveComponentFromCoordinates(
    allocator: std.mem.Allocator,
    mdl: *const Model,
    res: Residue,
) !?Component {
    // ── Step 1: collect eligible heavy atoms ────────────────────────────────
    var comp_atoms = std.ArrayListUnmanaged(CompAtom).empty;
    defer comp_atoms.deinit(allocator);

    // Map from local CompAtom index → model atom index within the residue slice.
    // We only need this during construction; not stored in the Component.
    var local_to_model = std.ArrayListUnmanaged(u32).empty;
    defer local_to_model.deinit(allocator);

    const res_atoms = mdl.residueAtoms(res);
    for (res_atoms, 0..) |*atom, i| {
        if (atom.is_hydrogen) continue;
        if (atom.element_type.info().flags.metallic) continue;

        var ca = CompAtom{};
        ca.name = atom.name.buf;
        ca.name_len = @intCast(atom.name.len);
        ca.element_symbol = elementSymbolFromType(atom.element_type);
        ca.ideal_x = atom.pos.x;
        ca.ideal_y = atom.pos.y;
        ca.ideal_z = atom.pos.z;
        // charge and leaving remain at zero/false — unknown from coordinates alone.

        try comp_atoms.append(allocator, ca);
        try local_to_model.append(allocator, @intCast(i));
    }

    if (comp_atoms.items.len < 2) return null;

    // ── Step 2: detect bonds by covalent radius sum ─────────────────────────
    var comp_bonds = std.ArrayListUnmanaged(CompBond).empty;
    defer comp_bonds.deinit(allocator);

    const n = comp_atoms.items.len;
    for (0..n) |i| {
        const ai_idx = local_to_model.items[i];
        const ai = &res_atoms[ai_idx];
        const ri = ai.element_type.info().covalent_radius;

        for (i + 1..n) |j| {
            const aj_idx = local_to_model.items[j];
            const aj = &res_atoms[aj_idx];
            const rj = aj.element_type.info().covalent_radius;

            const dist = ai.pos.distance(aj.pos);
            const threshold = @min((ri + rj) * COV_TOLERANCE, MAX_BOND_DIST);

            if (dist <= threshold) {
                try comp_bonds.append(allocator, .{
                    .atom_idx_1 = @intCast(i),
                    .atom_idx_2 = @intCast(j),
                    .order = .single,
                    .aromatic = false,
                });
            }
        }
    }

    // ── Step 3: bond order promotion ────────────────────────────────────────
    promoteBondOrders(comp_atoms.items, comp_bonds.items);

    // ── Step 3b: add synthetic H atoms ────────────────────────────────────
    // ccd_derive.derivePlans() iterates H atoms in the Component to generate
    // PlacementPlans. We must add explicit H entries (with bonds to their
    // parent heavy atom) so that derivePlans can see them.
    try addSyntheticHydrogens(allocator, &comp_atoms, &comp_bonds);

    // ── Step 4: build the Component ─────────────────────────────────────────
    const comp_id = try allocator.dupe(u8, res.compIdSlice());
    errdefer allocator.free(comp_id);

    const comp_type = try allocator.dupe(u8, "non-polymer");
    errdefer allocator.free(comp_type);

    const owned_atoms = try comp_atoms.toOwnedSlice(allocator);
    errdefer allocator.free(owned_atoms);

    const owned_bonds = try comp_bonds.toOwnedSlice(allocator);

    return Component{
        .comp_id = comp_id,
        .comp_type = comp_type,
        .atoms = owned_atoms,
        .bonds = owned_bonds,
    };
}

// ---------------------------------------------------------------------------
// Element symbol helper
// ---------------------------------------------------------------------------

/// Return a 2-byte element symbol for the given AtomType in the format expected
/// by CompAtom.element_symbol: uppercase first letter, space or lowercase second.
fn elementSymbolFromType(at: AtomType) [2]u8 {
    return switch (at) {
        .H, .Har, .Hpol, .Ha_p, .HOd => .{ 'H', ' ' },
        .C, .Car, .C_eq_O => .{ 'C', ' ' },
        .N, .Nacc => .{ 'N', ' ' },
        .O => .{ 'O', ' ' },
        .P => .{ 'P', ' ' },
        .S => .{ 'S', ' ' },
        .Se => .{ 'S', 'e' },
        .F => .{ 'F', ' ' },
        .Cl => .{ 'C', 'l' },
        .Br => .{ 'B', 'r' },
        .I => .{ 'I', ' ' },
        .Li => .{ 'L', 'i' },
        .Na => .{ 'N', 'a' },
        .Mg => .{ 'M', 'g' },
        .K => .{ 'K', ' ' },
        .Ca => .{ 'C', 'a' },
        .Mn => .{ 'M', 'n' },
        .Fe => .{ 'F', 'e' },
        .Co => .{ 'C', 'o' },
        .Ni => .{ 'N', 'i' },
        .Cu => .{ 'C', 'u' },
        .Zn => .{ 'Z', 'n' },
        .As => .{ 'A', 's' },
        .Rb => .{ 'R', 'b' },
        .Sr => .{ 'S', 'r' },
        .Mo => .{ 'M', 'o' },
        .Ag => .{ 'A', 'g' },
        .Cd => .{ 'C', 'd' },
        .Sn => .{ 'S', 'n' },
        .Cs => .{ 'C', 's' },
        .Ba => .{ 'B', 'a' },
        .W => .{ 'W', ' ' },
        .Pt => .{ 'P', 't' },
        .Au => .{ 'A', 'u' },
        .Hg => .{ 'H', 'g' },
        .Pb => .{ 'P', 'b' },
        .U => .{ 'U', ' ' },
        .unknown => .{ ' ', ' ' },
    };
}

// ---------------------------------------------------------------------------
// Synthetic hydrogen generation
// ---------------------------------------------------------------------------

/// Add synthetic H atoms to the component for each heavy atom that has
/// unfilled valence. These H entries are required by ccd_derive.derivePlans().
fn addSyntheticHydrogens(
    allocator: std.mem.Allocator,
    comp_atoms: *std.ArrayListUnmanaged(CompAtom),
    comp_bonds: *std.ArrayListUnmanaged(CompBond),
) !void {
    const n_heavy = comp_atoms.items.len;
    var h_counter: u16 = 0;

    for (0..n_heavy) |i| {
        const atom = comp_atoms.items[i];
        const mv = maxValence(atom.element_symbol);
        if (mv == 0) continue; // unknown element, skip

        const cv = currentValence(comp_bonds.items, i);
        if (cv >= mv) continue; // fully saturated

        const n_h = mv - cv;
        for (0..n_h) |_| {
            h_counter += 1;

            // Generate H atom name: H1, H2, H3, ...
            var name_buf: [8]u8 = undefined;
            const name_str = std.fmt.bufPrint(&name_buf, "H{d}", .{h_counter}) catch unreachable;
            const copy_len: u4 = @intCast(@min(name_str.len, 4));
            var h_name: [4]u8 = .{ ' ', ' ', ' ', ' ' };
            @memcpy(h_name[0..copy_len], name_str[0..copy_len]);

            try comp_atoms.append(allocator, CompAtom{
                .name = h_name,
                .name_len = copy_len,
                .element_symbol = .{ 'H', ' ' },
            });

            try comp_bonds.append(allocator, CompBond{
                .atom_idx_1 = @intCast(i),
                .atom_idx_2 = @intCast(comp_atoms.items.len - 1),
                .order = .single,
            });
        }
    }
}

// ---------------------------------------------------------------------------
// Bond order promotion
// ---------------------------------------------------------------------------

/// Maximum covalent bond valence for each heavy element.
/// Only covers elements that may appear in organic/bio ligands; others stay at 1.
fn maxValence(symbol: [2]u8) u8 {
    if (symbol[0] == 'C') {
        if (symbol[1] == ' ') return 4; // C
        if (symbol[1] == 'l') return 1; // Cl
        if (symbol[1] == 'o') return 2; // Co (metal, won't be here but safe)
        if (symbol[1] == 'u') return 2; // Cu
        if (symbol[1] == 'a') return 2; // Ca
        if (symbol[1] == 'd') return 2; // Cd
        if (symbol[1] == 's') return 1; // Cs
        return 4;
    }
    if (symbol[0] == 'N') {
        if (symbol[1] == ' ') return 3; // N (formal; +1 charge → 4, not attempted here)
        if (symbol[1] == 'i') return 2; // Ni
        if (symbol[1] == 'a') return 1; // Na
        return 3;
    }
    if (symbol[0] == 'O') return 2; // O
    if (symbol[0] == 'S') {
        if (symbol[1] == ' ') return 2; // S (not trying to detect S(IV)/S(VI))
        if (symbol[1] == 'e') return 2; // Se
        if (symbol[1] == 'r') return 2; // Sr
        if (symbol[1] == 'n') return 4; // Sn
        return 2;
    }
    if (symbol[0] == 'P') {
        if (symbol[1] == ' ') return 5; // P (can be 3 or 5; use 5 for safety)
        if (symbol[1] == 't') return 4; // Pt
        if (symbol[1] == 'b') return 4; // Pb
        return 5;
    }
    if (symbol[0] == 'B') {
        if (symbol[1] == 'r') return 1; // Br
        if (symbol[1] == 'a') return 2; // Ba
        return 3; // B
    }
    if (symbol[0] == 'F') return 1; // F
    if (symbol[0] == 'I') return 1; // I
    return 1; // conservative default for all other elements
}

/// Bond-order value as an integer for valence arithmetic.
fn bondOrderValue(order: BondOrder) u8 {
    return switch (order) {
        .single => 1,
        .double => 2,
        .triple => 3,
        .aromatic, .delocalized => 2, // approximate
        .unknown => 1,
    };
}

/// Compute current valence of atom at index `idx` given the current bond list.
fn currentValence(bonds: []const CompBond, idx: usize) u8 {
    var val: u8 = 0;
    for (bonds) |b| {
        if (b.atom_idx_1 == idx or b.atom_idx_2 == idx) {
            val += bondOrderValue(b.order);
        }
    }
    return val;
}

/// Priority for promotion: prefer C–O, then C–N, then C–C, then anything else.
/// Lower number = higher priority.
fn promotionPriority(sym1: [2]u8, sym2: [2]u8) u8 {
    const c1_is_c = sym1[0] == 'C' and sym1[1] == ' ';
    const c2_is_c = sym2[0] == 'C' and sym2[1] == ' ';
    const c1_is_o = sym1[0] == 'O' and sym1[1] == ' ';
    const c2_is_o = sym2[0] == 'O' and sym2[1] == ' ';
    const c1_is_n = sym1[0] == 'N' and sym1[1] == ' ';
    const c2_is_n = sym2[0] == 'N' and sym2[1] == ' ';

    const c_o = (c1_is_c and c2_is_o) or (c2_is_c and c1_is_o);
    const c_n = (c1_is_c and c2_is_n) or (c2_is_c and c1_is_n);
    const c_c = c1_is_c and c2_is_c;

    if (c_o) return 0;
    if (c_n) return 1;
    if (c_c) return 2;
    return 3;
}

/// A bond candidate for promotion: index into bonds slice plus priority.
const PromotionCandidate = struct {
    bond_idx: usize,
    priority: u8,
};

/// Iteratively promote single bonds to double where both endpoints are
/// under-valent. Promotes highest-priority bonds first (C=O > C=N > C=C).
///
/// We do at most one promotion pass so the algorithm stays O(n²) without
/// complex fixpoint iteration. In practice, ligands rarely need more than
/// one pass to assign all obvious double bonds.
fn promoteBondOrders(atoms: []CompAtom, bonds: []CompBond) void {
    // Gather candidates sorted by priority.
    var candidates: [256]PromotionCandidate = undefined;
    var ncand: usize = 0;

    for (bonds, 0..) |b, bi| {
        if (b.order != .single) continue;
        if (ncand >= candidates.len) break;
        const sym1 = atoms[b.atom_idx_1].element_symbol;
        const sym2 = atoms[b.atom_idx_2].element_symbol;
        candidates[ncand] = .{
            .bond_idx = bi,
            .priority = promotionPriority(sym1, sym2),
        };
        ncand += 1;
    }

    // Sort by priority ascending (lower = preferred).
    const cands = candidates[0..ncand];
    std.sort.pdq(PromotionCandidate, cands, {}, struct {
        fn lt(_: void, a: PromotionCandidate, b: PromotionCandidate) bool {
            return a.priority < b.priority;
        }
    }.lt);

    for (cands) |cand| {
        const b = &bonds[cand.bond_idx];
        const idx1: usize = b.atom_idx_1;
        const idx2: usize = b.atom_idx_2;

        const mv1 = maxValence(atoms[idx1].element_symbol);
        const mv2 = maxValence(atoms[idx2].element_symbol);
        const cv1 = currentValence(bonds, idx1);
        const cv2 = currentValence(bonds, idx2);

        // Both atoms must have room for an extra bond order.
        if (cv1 < mv1 and cv2 < mv2) {
            b.order = .double;
        }
    }
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

const testing = std.testing;

/// Count non-hydrogen atoms in a Component (for test assertions on heavy atom topology).
fn countHeavyAtoms(comp: Component) usize {
    var n: usize = 0;
    for (comp.atoms) |a| {
        if (a.element_symbol[0] != 'H') n += 1;
    }
    return n;
}

/// Count bonds between non-hydrogen atoms (for test assertions on heavy atom topology).
fn countHeavyBonds(comp: Component) usize {
    var n: usize = 0;
    for (comp.bonds) |b| {
        if (comp.atoms[b.atom_idx_1].element_symbol[0] != 'H' and
            comp.atoms[b.atom_idx_2].element_symbol[0] != 'H')
        {
            n += 1;
        }
    }
    return n;
}

/// Count synthetic hydrogen atoms in a Component.
fn countSyntheticH(comp: Component) usize {
    var n: usize = 0;
    for (comp.atoms) |a| {
        if (a.element_symbol[0] == 'H') n += 1;
    }
    return n;
}

/// Descriptor for a single test atom.
const TestAtomDesc = struct {
    name: []const u8,
    elem: AtomType,
    x: f32,
    y: f32,
    z: f32,
    is_h: bool = false,
};

/// Build a minimal Model from an array of TestAtomDesc entries.
fn buildModel(allocator: std.mem.Allocator, entries: []const TestAtomDesc) !Model {
    var mdl = Model.init(allocator);
    errdefer mdl.deinit();

    for (entries) |e| {
        var atom = Atom{
            .pos = .{ .x = e.x, .y = e.y, .z = e.z },
            .element_type = e.elem,
            .is_hydrogen = e.is_h,
        };
        atom.setName(e.name);
        try mdl.atoms.append(mdl.allocator, atom);
    }

    return mdl;
}

test "null return: zero heavy atoms" {
    var mdl = Model.init(testing.allocator);
    defer mdl.deinit();

    // No atoms at all — residue window is empty.
    const res = Residue{ .atom_start = 0, .atom_end = 0 };

    const result = try deriveComponentFromCoordinates(testing.allocator, &mdl, res);
    try testing.expectEqual(@as(?Component, null), result);
}

test "null return: one heavy atom only" {
    const entries = [_]TestAtomDesc{
        .{ .name = "C1", .elem = .C, .x = 0, .y = 0, .z = 0 },
    };

    var mdl = try buildModel(testing.allocator, &entries);
    defer mdl.deinit();

    var res = Residue{ .atom_start = 0, .atom_end = 1 };
    res.setCompId("LIG");

    const result = try deriveComponentFromCoordinates(testing.allocator, &mdl, res);
    try testing.expectEqual(@as(?Component, null), result);
}

test "null return: only hydrogen atoms" {
    const entries = [_]TestAtomDesc{
        .{ .name = "H1", .elem = .H, .x = 0, .y = 0, .z = 0, .is_h = true },
        .{ .name = "H2", .elem = .H, .x = 1, .y = 0, .z = 0, .is_h = true },
    };

    var mdl = try buildModel(testing.allocator, &entries);
    defer mdl.deinit();

    var res = Residue{ .atom_start = 0, .atom_end = 2 };
    res.setCompId("LIG");

    const result = try deriveComponentFromCoordinates(testing.allocator, &mdl, res);
    try testing.expectEqual(@as(?Component, null), result);
}

test "bond detection threshold: atoms within range bond" {
    // C at origin, O at 1.2 Å — well within (0.77+0.66)*1.3 = 1.859 Å
    const entries = [_]TestAtomDesc{
        .{ .name = "C1", .elem = .C, .x = 0.0, .y = 0.0, .z = 0.0 },
        .{ .name = "O1", .elem = .O, .x = 1.2, .y = 0.0, .z = 0.0 },
    };

    var mdl = try buildModel(testing.allocator, &entries);
    defer mdl.deinit();

    var res = Residue{ .atom_start = 0, .atom_end = 2 };
    res.setCompId("LIG");

    const comp = (try deriveComponentFromCoordinates(testing.allocator, &mdl, res)).?;
    defer {
        testing.allocator.free(comp.atoms);
        testing.allocator.free(comp.bonds);
        testing.allocator.free(comp.comp_id);
        testing.allocator.free(comp.comp_type);
    }

    try testing.expectEqual(@as(usize, 2), countHeavyAtoms(comp));
    try testing.expectEqual(@as(usize, 1), countHeavyBonds(comp));
    // C=O (double) → C has valence 2, max 4 → 2 syn-H; O full → 0. Total: 2.
    try testing.expectEqual(@as(usize, 2), countSyntheticH(comp));
}

test "bond detection threshold: atoms outside range do not bond" {
    // C at origin, O at 3.0 Å — beyond MAX_BOND_DIST cap of 2.5 Å.
    const entries = [_]TestAtomDesc{
        .{ .name = "C1", .elem = .C, .x = 0.0, .y = 0.0, .z = 0.0 },
        .{ .name = "O1", .elem = .O, .x = 3.0, .y = 0.0, .z = 0.0 },
    };

    var mdl = try buildModel(testing.allocator, &entries);
    defer mdl.deinit();

    var res = Residue{ .atom_start = 0, .atom_end = 2 };
    res.setCompId("LIG");

    const comp = (try deriveComponentFromCoordinates(testing.allocator, &mdl, res)).?;
    defer {
        testing.allocator.free(comp.atoms);
        testing.allocator.free(comp.bonds);
        testing.allocator.free(comp.comp_id);
        testing.allocator.free(comp.comp_type);
    }

    try testing.expectEqual(@as(usize, 2), countHeavyAtoms(comp));
    try testing.expectEqual(@as(usize, 0), countHeavyBonds(comp));
    // No bonds → C gets 4 syn-H, O gets 2 syn-H = 6 total.
    try testing.expectEqual(@as(usize, 6), countSyntheticH(comp));
}

test "metal exclusion: Fe atom not bonded and not included" {
    // Fe at origin, C nearby — Fe is metallic so excluded entirely.
    // Only C should appear; with one atom we get null.
    const entries = [_]TestAtomDesc{
        .{ .name = "FE", .elem = .Fe, .x = 0.0, .y = 0.0, .z = 0.0 },
        .{ .name = "C1", .elem = .C, .x = 1.0, .y = 0.0, .z = 0.0 },
    };

    var mdl = try buildModel(testing.allocator, &entries);
    defer mdl.deinit();

    var res = Residue{ .atom_start = 0, .atom_end = 2 };
    res.setCompId("HEM");

    // Fe is excluded: only C remains → 1 heavy atom → null.
    const result = try deriveComponentFromCoordinates(testing.allocator, &mdl, res);
    try testing.expectEqual(@as(?Component, null), result);
}

test "metal exclusion: two carbons bonded, Fe ignored" {
    // C1 at (0,0,0), C2 at (1.5,0,0), Fe at (0.5,0,0) within bonding distance of both.
    // Fe must be excluded — only C1–C2 bond should appear.
    const entries = [_]TestAtomDesc{
        .{ .name = "C1", .elem = .C, .x = 0.0, .y = 0.0, .z = 0.0 },
        .{ .name = "FE", .elem = .Fe, .x = 0.5, .y = 0.0, .z = 0.0 },
        .{ .name = "C2", .elem = .C, .x = 1.5, .y = 0.0, .z = 0.0 },
    };

    var mdl = try buildModel(testing.allocator, &entries);
    defer mdl.deinit();

    var res = Residue{ .atom_start = 0, .atom_end = 3 };
    res.setCompId("HEM");

    const comp = (try deriveComponentFromCoordinates(testing.allocator, &mdl, res)).?;
    defer {
        testing.allocator.free(comp.atoms);
        testing.allocator.free(comp.bonds);
        testing.allocator.free(comp.comp_id);
        testing.allocator.free(comp.comp_type);
    }

    // Only 2 atoms (both carbons); Fe excluded.
    try testing.expectEqual(@as(usize, 2), countHeavyAtoms(comp));
    // C1–C2 at 1.5 Å: (0.77+0.77)*1.3 = 2.002 Å threshold. 1.5 < 2.002 → bonded.
    try testing.expectEqual(@as(usize, 1), countHeavyBonds(comp));
    // C1: 1 bond (max 4) → C=C double promoted → valence 2, needs 2 syn-H.
    // C2: same → 2 syn-H. Total: 4.
    try testing.expectEqual(@as(usize, 4), countSyntheticH(comp));
}

test "simple sp3 carbon: tetrahedral methane-like with 3 heavy neighbors" {
    // C1 at origin bonded to N1, O1, S1 in a 120° star topology at 1.4 Å.
    // Neighbour-to-neighbour distances = 2*1.4*sin(60°) = 2.425 Å, which exceeds
    // every pairwise covalent-radius threshold:
    //   N–O: (0.70+0.66)*1.3 = 1.768 Å < 2.425 → no bond
    //   N–S: (0.70+1.04)*1.3 = 2.262 Å < 2.425 → no bond
    //   O–S: (0.66+1.04)*1.3 = 2.210 Å < 2.425 → no bond
    // So exactly 3 bonds are expected: C1–N1, C1–O1, C1–S1.
    const entries = [_]TestAtomDesc{
        .{ .name = "C1", .elem = .C, .x = 0.0, .y = 0.0, .z = 0.0 },
        .{ .name = "N1", .elem = .N, .x = 1.4, .y = 0.0, .z = 0.0 },
        .{ .name = "O1", .elem = .O, .x = -0.7, .y = 1.212, .z = 0.0 },
        .{ .name = "S1", .elem = .S, .x = -0.7, .y = -1.212, .z = 0.0 },
    };

    var mdl = try buildModel(testing.allocator, &entries);
    defer mdl.deinit();

    var res = Residue{ .atom_start = 0, .atom_end = 4 };
    res.setCompId("LIG");

    const comp = (try deriveComponentFromCoordinates(testing.allocator, &mdl, res)).?;
    defer {
        testing.allocator.free(comp.atoms);
        testing.allocator.free(comp.bonds);
        testing.allocator.free(comp.comp_id);
        testing.allocator.free(comp.comp_type);
    }

    try testing.expectEqual(@as(usize, 4), countHeavyAtoms(comp));
    // Expect exactly 3 heavy bonds: C1–N1, C1–O1, C1–S1.
    try testing.expectEqual(@as(usize, 3), countHeavyBonds(comp));
    // After promotion: C1–O1 promoted to double (C-O priority 0).
    // C1: double(O)+single(N)+single(S) = valence 4, max 4 → 0 syn-H.
    // N1: 1 single, max 3 → 2 syn-H.
    // O1: 1 double, valence 2, max 2 → 0 syn-H.
    // S1: 1 single, max 2 → 1 syn-H.
    // Total: 0+2+0+1 = 3.
    try testing.expectEqual(@as(usize, 3), countSyntheticH(comp));
}

test "bond order promotion: C–O becomes double (formaldehyde-like)" {
    // C at (0,0,0), O at (1.2,0,0). C has 1 heavy bond (max 4), O has 1 (max 2).
    // Both under-valent → promote C–O to double.
    const entries = [_]TestAtomDesc{
        .{ .name = "C1", .elem = .C, .x = 0.0, .y = 0.0, .z = 0.0 },
        .{ .name = "O1", .elem = .O, .x = 1.2, .y = 0.0, .z = 0.0 },
    };

    var mdl = try buildModel(testing.allocator, &entries);
    defer mdl.deinit();

    var res = Residue{ .atom_start = 0, .atom_end = 2 };
    res.setCompId("LIG");

    const comp = (try deriveComponentFromCoordinates(testing.allocator, &mdl, res)).?;
    defer {
        testing.allocator.free(comp.atoms);
        testing.allocator.free(comp.bonds);
        testing.allocator.free(comp.comp_id);
        testing.allocator.free(comp.comp_type);
    }

    try testing.expectEqual(@as(usize, 1), countHeavyBonds(comp));
    // Find the heavy C-O bond and verify it's double.
    for (comp.bonds) |b| {
        if (comp.atoms[b.atom_idx_1].element_symbol[0] != 'H' and
            comp.atoms[b.atom_idx_2].element_symbol[0] != 'H')
        {
            try testing.expectEqual(BondOrder.double, b.order);
        }
    }
    // C=O double → C valence 2 max 4, needs 2 syn-H; O valence 2 max 2, full → 0.
    try testing.expectEqual(@as(usize, 2), countSyntheticH(comp));
}

test "bond order promotion: CO2-like, both C=O bonds promoted" {
    // C bonded to two O atoms (linear CO2-like fragment without the other bonds).
    // C at origin, O1 at +1.2 Å, O2 at -1.2 Å. O1–O2 distance 2.4 Å > threshold 1.716 Å.
    // After initial singles: C valence=2, O1 valence=1, O2 valence=1.
    // Both C–O bonds have priority 0. In sorted order, the first promotion:
    //   C valence=2 < 4, O1 valence=1 < 2 → promote C–O1 to double. C valence=3.
    // Second C–O: C valence=3 < 4, O2 valence=1 < 2 → also promoted. C valence=4.
    // Result: two C=O double bonds.
    const entries = [_]TestAtomDesc{
        .{ .name = "C1", .elem = .C, .x = 0.0, .y = 0.0, .z = 0.0 },
        .{ .name = "O1", .elem = .O, .x = 1.2, .y = 0.0, .z = 0.0 },
        .{ .name = "O2", .elem = .O, .x = -1.2, .y = 0.0, .z = 0.0 },
    };

    var mdl = try buildModel(testing.allocator, &entries);
    defer mdl.deinit();

    var res = Residue{ .atom_start = 0, .atom_end = 3 };
    res.setCompId("LIG");

    const comp = (try deriveComponentFromCoordinates(testing.allocator, &mdl, res)).?;
    defer {
        testing.allocator.free(comp.atoms);
        testing.allocator.free(comp.bonds);
        testing.allocator.free(comp.comp_id);
        testing.allocator.free(comp.comp_type);
    }

    // Two heavy bonds: C–O1 and C–O2. O1–O2 at 2.4 Å; threshold (0.66+0.66)*1.3=1.716 < 2.4 → no O–O bond.
    try testing.expectEqual(@as(usize, 2), countHeavyBonds(comp));
    // Count double bonds among heavy bonds.
    var n_double: usize = 0;
    for (comp.bonds) |b| {
        if (b.order == .double and
            comp.atoms[b.atom_idx_1].element_symbol[0] != 'H' and
            comp.atoms[b.atom_idx_2].element_symbol[0] != 'H')
        {
            n_double += 1;
        }
    }
    // Both C=O bonds should be promoted since C has room for both.
    try testing.expectEqual(@as(usize, 2), n_double);
    // C: 2 double = valence 4, full → 0 syn-H. O1: valence 2, full → 0. O2: same → 0.
    try testing.expectEqual(@as(usize, 0), countSyntheticH(comp));
}

test "bond order promotion: C–N promoted when both under-valent" {
    const entries = [_]TestAtomDesc{
        .{ .name = "C1", .elem = .C, .x = 0.0, .y = 0.0, .z = 0.0 },
        .{ .name = "N1", .elem = .N, .x = 1.3, .y = 0.0, .z = 0.0 },
    };

    var mdl = try buildModel(testing.allocator, &entries);
    defer mdl.deinit();

    var res = Residue{ .atom_start = 0, .atom_end = 2 };
    res.setCompId("LIG");

    const comp = (try deriveComponentFromCoordinates(testing.allocator, &mdl, res)).?;
    defer {
        testing.allocator.free(comp.atoms);
        testing.allocator.free(comp.bonds);
        testing.allocator.free(comp.comp_id);
        testing.allocator.free(comp.comp_type);
    }

    try testing.expectEqual(@as(usize, 1), countHeavyBonds(comp));
    // C has 1 bond (max 4), N has 1 bond (max 3) → both under-valent → double.
    for (comp.bonds) |b| {
        if (comp.atoms[b.atom_idx_1].element_symbol[0] != 'H' and
            comp.atoms[b.atom_idx_2].element_symbol[0] != 'H')
        {
            try testing.expectEqual(BondOrder.double, b.order);
        }
    }
    // C=N double: C valence 2 max 4 → 2 syn-H; N valence 2 max 3 → 1 syn-H. Total: 3.
    try testing.expectEqual(@as(usize, 3), countSyntheticH(comp));
}

test "bond order promotion: C–O preferred over C–N (priority ordering)" {
    // Three atoms: C bonded to O and N. Both bonds promoted since C (max 4) has room.
    // C at (0,0,0), O at (1.2,0,0), N at (-1.3,0,0).
    // O–N distance = 2.5 Å; threshold = min((0.66+0.70)*1.3, 2.5) = min(1.768, 2.5) = 1.768 Å.
    // 2.5 > 1.768 → no O–N bond. Exactly 2 bonds: C–O and C–N.
    const entries = [_]TestAtomDesc{
        .{ .name = "C1", .elem = .C, .x = 0.0, .y = 0.0, .z = 0.0 },
        .{ .name = "O1", .elem = .O, .x = 1.2, .y = 0.0, .z = 0.0 },
        .{ .name = "N1", .elem = .N, .x = -1.3, .y = 0.0, .z = 0.0 },
    };

    var mdl = try buildModel(testing.allocator, &entries);
    defer mdl.deinit();

    var res = Residue{ .atom_start = 0, .atom_end = 3 };
    res.setCompId("GLX");

    const comp = (try deriveComponentFromCoordinates(testing.allocator, &mdl, res)).?;
    defer {
        testing.allocator.free(comp.atoms);
        testing.allocator.free(comp.bonds);
        testing.allocator.free(comp.comp_id);
        testing.allocator.free(comp.comp_type);
    }

    try testing.expectEqualStrings("GLX", comp.comp_id);
    try testing.expectEqualStrings("non-polymer", comp.comp_type);

    // O–N distance: 2.5 Å exactly equals MAX_BOND_DIST → ≤ threshold → bonded.
    // But (0.66+0.70)*1.3=1.768 Å threshold, MAX_BOND_DIST=2.5 Å cap: min(1.768, 2.5)=1.768.
    // 2.5 > 1.768 → NO O–N bond. Good.
    try testing.expectEqual(@as(usize, 2), countHeavyBonds(comp));

    // Both C–O and C–N should be promoted.
    var n_double: usize = 0;
    for (comp.bonds) |b| {
        if (b.order == .double and
            comp.atoms[b.atom_idx_1].element_symbol[0] != 'H' and
            comp.atoms[b.atom_idx_2].element_symbol[0] != 'H')
        {
            n_double += 1;
        }
    }
    try testing.expectEqual(@as(usize, 2), n_double);
}

test "elementSymbolFromType covers all non-metal elements" {
    // Spot-check a few element symbols.
    try testing.expectEqual([2]u8{ 'H', ' ' }, elementSymbolFromType(.H));
    try testing.expectEqual([2]u8{ 'H', ' ' }, elementSymbolFromType(.Hpol));
    try testing.expectEqual([2]u8{ 'C', ' ' }, elementSymbolFromType(.C));
    try testing.expectEqual([2]u8{ 'C', ' ' }, elementSymbolFromType(.Car));
    try testing.expectEqual([2]u8{ 'N', ' ' }, elementSymbolFromType(.N));
    try testing.expectEqual([2]u8{ 'O', ' ' }, elementSymbolFromType(.O));
    try testing.expectEqual([2]u8{ 'S', ' ' }, elementSymbolFromType(.S));
    try testing.expectEqual([2]u8{ 'S', 'e' }, elementSymbolFromType(.Se));
    try testing.expectEqual([2]u8{ 'C', 'l' }, elementSymbolFromType(.Cl));
    try testing.expectEqual([2]u8{ 'B', 'r' }, elementSymbolFromType(.Br));
    try testing.expectEqual([2]u8{ 'F', 'e' }, elementSymbolFromType(.Fe));
    try testing.expectEqual([2]u8{ ' ', ' ' }, elementSymbolFromType(.unknown));
}

test "comp_atom fields correctly populated from model atom" {
    const entries = [_]TestAtomDesc{
        .{ .name = "CA", .elem = .C, .x = 1.5, .y = 2.5, .z = 3.5 },
        .{ .name = "N", .elem = .N, .x = 0.0, .y = 0.0, .z = 0.0 },
    };

    var mdl = try buildModel(testing.allocator, &entries);
    defer mdl.deinit();

    var res = Residue{ .atom_start = 0, .atom_end = 2 };
    res.setCompId("TST");

    const comp = (try deriveComponentFromCoordinates(testing.allocator, &mdl, res)).?;
    defer {
        testing.allocator.free(comp.atoms);
        testing.allocator.free(comp.bonds);
        testing.allocator.free(comp.comp_id);
        testing.allocator.free(comp.comp_type);
    }

    try testing.expectEqual(@as(usize, 2), countHeavyAtoms(comp));

    // First atom should be CA.
    try testing.expectEqualSlices(u8, "CA  ", &comp.atoms[0].name);
    try testing.expectEqual(@as(u4, 2), comp.atoms[0].name_len);
    try testing.expectEqual([2]u8{ 'C', ' ' }, comp.atoms[0].element_symbol);
    try testing.expectApproxEqAbs(@as(f32, 1.5), comp.atoms[0].ideal_x, 1e-5);
    try testing.expectApproxEqAbs(@as(f32, 2.5), comp.atoms[0].ideal_y, 1e-5);
    try testing.expectApproxEqAbs(@as(f32, 3.5), comp.atoms[0].ideal_z, 1e-5);
    try testing.expectEqual(@as(i8, 0), comp.atoms[0].charge);
    try testing.expectEqual(false, comp.atoms[0].leaving);

    // Second atom: N.
    try testing.expectEqual([2]u8{ 'N', ' ' }, comp.atoms[1].element_symbol);
    try testing.expectEqual(@as(u4, 1), comp.atoms[1].name_len);
}

test "hydrogen atoms excluded from CompAtom list" {
    const entries = [_]TestAtomDesc{
        .{ .name = "C1", .elem = .C, .x = 0.0, .y = 0.0, .z = 0.0 },
        .{ .name = "H1", .elem = .H, .x = 1.1, .y = 0.0, .z = 0.0, .is_h = true },
        .{ .name = "O1", .elem = .O, .x = -1.2, .y = 0.0, .z = 0.0 },
    };

    var mdl = try buildModel(testing.allocator, &entries);
    defer mdl.deinit();

    var res = Residue{ .atom_start = 0, .atom_end = 3 };
    res.setCompId("LIG");

    const comp = (try deriveComponentFromCoordinates(testing.allocator, &mdl, res)).?;
    defer {
        testing.allocator.free(comp.atoms);
        testing.allocator.free(comp.bonds);
        testing.allocator.free(comp.comp_id);
        testing.allocator.free(comp.comp_type);
    }

    // Only C and O as heavy atoms; model H1 excluded from heavy atom list.
    try testing.expectEqual(@as(usize, 2), countHeavyAtoms(comp));
    // Verify no heavy atom is named H1 (model H was excluded).
    for (comp.atoms) |a| {
        if (a.element_symbol[0] != 'H') {
            try testing.expect(!std.mem.eql(u8, a.nameSlice(), "H1"));
        }
    }
    // Synthetic H atoms should still be present (derived from valence).
    try testing.expect(countSyntheticH(comp) > 0);
}
