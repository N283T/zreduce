//! mmCIF extraction: parses _atom_site loop data into a Model.

const std = @import("std");
const Allocator = std.mem.Allocator;

const cif = @import("cif.zig");
const model_mod = @import("model.zig");
const element = @import("element.zig");
const bond_mod = @import("model/bond.zig");
const ccd_mod = @import("ccd.zig");

const Model = model_mod.Model;
const Atom = model_mod.Atom;
const Residue = model_mod.Residue;
const Chain = model_mod.Chain;

pub const MmcifError = error{
    NoAtomSiteLoop,
    MissingCoordinateField,
    InvalidCoordinateValue,
    CifParseError,
    OutOfMemory,
};

/// Column indices for _atom_site loop fields.
const AtomSiteColumns = struct {
    cartn_x: ?usize = null,
    cartn_y: ?usize = null,
    cartn_z: ?usize = null,
    type_symbol: ?usize = null,
    label_atom_id: ?usize = null,
    label_comp_id: ?usize = null,
    label_asym_id: ?usize = null,
    label_seq_id: ?usize = null,
    auth_asym_id: ?usize = null,
    label_alt_id: ?usize = null,
    occupancy: ?usize = null,
    b_iso_or_equiv: ?usize = null,
    group_pdb: ?usize = null,
    id: ?usize = null,
    pdbx_PDB_ins_code: ?usize = null,
    label_entity_id: ?usize = null,
    auth_seq_id: ?usize = null,
};

/// Compare a Chain's label_asym_id with a string from the CIF data.
fn chainAsymMatches(chain: Chain, asym: []const u8) bool {
    return std.mem.eql(u8, chain.labelSlice(), asym);
}

/// Scan _pdbx_poly_seq_scheme to detect sequence gaps and mark residues
/// that follow an unobserved residue (auth_seq_num = '?' or '.').
fn detectChainBreaks(
    mdl: *Model,
    pss: *const cif.Loop,
    col_asym: usize,
    col_seq: usize,
    col_auth_seq: usize,
) void {
    const nrows = pss.length();
    var prev_asym: []const u8 = "";
    var gap_pending = false;

    for (0..nrows) |row| {
        const asym = cif.asString(pss.val(row, col_asym) orelse continue);
        const seq_str = pss.val(row, col_seq) orelse continue;
        const auth_seq = cif.asString(pss.val(row, col_auth_seq) orelse "?");

        // Reset on chain change
        if (!std.mem.eql(u8, asym, prev_asym)) {
            gap_pending = false;
            prev_asym = asym;
        }

        // Unobserved residue (asString converts '?' and '.' to empty)
        if (auth_seq.len == 0) {
            gap_pending = true;
            continue;
        }

        // Observed residue after gap -> find and mark
        if (gap_pending) {
            const seq_id = cif.value.asIntOr(i32, seq_str, 0);
            for (mdl.residues.items) |*res| {
                const chain = mdl.chains.items[res.chain_idx];
                if (chainAsymMatches(chain, asym) and res.seq_id == seq_id) {
                    res.is_chain_break_before = true;
                    break;
                }
            }
            gap_pending = false;
        }
    }
}

const EntityType = model_mod.residue.EntityType;

fn entityTypeFromString(s: []const u8) EntityType {
    if (std.ascii.eqlIgnoreCase(s, "polymer")) return .polymer;
    if (std.ascii.eqlIgnoreCase(s, "non-polymer")) return .non_polymer;
    if (std.ascii.eqlIgnoreCase(s, "branched")) return .branched;
    if (std.ascii.eqlIgnoreCase(s, "water")) return .water;
    return .unknown;
}

fn applyEntityTypes(mdl: *Model, block: *const cif.Block) void {
    const ent = block.findLoop("_entity.id") orelse return;
    const col_id = ent.findTag("_entity.id") orelse return;
    const col_type = ent.findTag("_entity.type") orelse return;

    for (mdl.residues.items) |*res| {
        const chain = mdl.chains.items[res.chain_idx];
        const entity_id = chain.entityIdSlice();
        if (entity_id.len == 0) continue;

        for (0..ent.length()) |row| {
            const eid = cif.asString(ent.val(row, col_id) orelse continue);
            if (std.mem.eql(u8, eid, entity_id)) {
                const etype = cif.asString(ent.val(row, col_type) orelse continue);
                res.entity_type = entityTypeFromString(etype);
                break;
            }
        }
    }
}

/// Parse an mmCIF source string and extract all _atom_site records into a Model.
pub fn parseModel(allocator: Allocator, source: []const u8) MmcifError!Model {
    var doc = cif.readString(allocator, source) catch |err| switch (err) {
        error.OutOfMemory => return MmcifError.OutOfMemory,
        else => return MmcifError.CifParseError,
    };
    defer doc.deinit();

    if (doc.blocks.items.len == 0) return MmcifError.NoAtomSiteLoop;

    const block = &doc.blocks.items[0];
    const loop = block.findLoop("_atom_site.Cartn_x") orelse return MmcifError.NoAtomSiteLoop;

    // Map column indices
    var cols = AtomSiteColumns{};
    cols.cartn_x = loop.findTag("_atom_site.Cartn_x");
    cols.cartn_y = loop.findTag("_atom_site.Cartn_y");
    cols.cartn_z = loop.findTag("_atom_site.Cartn_z");
    cols.type_symbol = loop.findTag("_atom_site.type_symbol");
    cols.label_atom_id = loop.findTag("_atom_site.label_atom_id");
    cols.label_comp_id = loop.findTag("_atom_site.label_comp_id");
    cols.label_asym_id = loop.findTag("_atom_site.label_asym_id");
    cols.label_seq_id = loop.findTag("_atom_site.label_seq_id");
    cols.auth_asym_id = loop.findTag("_atom_site.auth_asym_id");
    cols.label_alt_id = loop.findTag("_atom_site.label_alt_id");
    cols.occupancy = loop.findTag("_atom_site.occupancy");
    cols.b_iso_or_equiv = loop.findTag("_atom_site.B_iso_or_equiv");
    cols.group_pdb = loop.findTag("_atom_site.group_PDB");
    cols.id = loop.findTag("_atom_site.id");
    cols.pdbx_PDB_ins_code = loop.findTag("_atom_site.pdbx_PDB_ins_code");
    cols.label_entity_id = loop.findTag("_atom_site.label_entity_id");
    cols.auth_seq_id = loop.findTag("_atom_site.auth_seq_id");

    // Require x, y, z
    if (cols.cartn_x == null or cols.cartn_y == null or cols.cartn_z == null) {
        return MmcifError.MissingCoordinateField;
    }

    var mdl = Model.init(allocator);
    errdefer mdl.deinit();

    const nrows = loop.length();

    // State tracking for chain/residue boundaries
    var cur_label_asym_id: []const u8 = "";
    var cur_seq_id: []const u8 = "";
    var cur_comp_id: []const u8 = "";
    var cur_ins_code: []const u8 = "";
    var cur_auth_seq: []const u8 = "";
    var res_atom_start: u32 = 0;
    var in_chain = false;
    var in_residue = false;

    for (0..nrows) |row| {
        const x_str = loop.val(row, cols.cartn_x.?) orelse return MmcifError.InvalidCoordinateValue;
        const y_str = loop.val(row, cols.cartn_y.?) orelse return MmcifError.InvalidCoordinateValue;
        const z_str = loop.val(row, cols.cartn_z.?) orelse return MmcifError.InvalidCoordinateValue;

        const x = cif.asFloat(x_str) orelse return MmcifError.InvalidCoordinateValue;
        const y = cif.asFloat(y_str) orelse return MmcifError.InvalidCoordinateValue;
        const z = cif.asFloat(z_str) orelse return MmcifError.InvalidCoordinateValue;

        // Get string fields
        const type_sym = if (cols.type_symbol) |c| cif.asString(loop.val(row, c) orelse ".") else "";
        const atom_name = if (cols.label_atom_id) |c| cif.asString(loop.val(row, c) orelse ".") else "";
        const comp_id = if (cols.label_comp_id) |c| cif.asString(loop.val(row, c) orelse ".") else "";
        const label_asym = if (cols.label_asym_id) |c| cif.asString(loop.val(row, c) orelse ".") else "";
        const seq_id = if (cols.label_seq_id) |c| cif.asString(loop.val(row, c) orelse ".") else "";
        const auth_seq = if (cols.auth_seq_id) |c| cif.asString(loop.val(row, c) orelse ".") else "";
        const auth_asym = if (cols.auth_asym_id) |c| cif.asString(loop.val(row, c) orelse ".") else "";
        const alt_loc_str = if (cols.label_alt_id) |c| loop.val(row, c) orelse "." else ".";
        const occ_str = if (cols.occupancy) |c| loop.val(row, c) orelse "1.0" else "1.0";
        const bfac_str = if (cols.b_iso_or_equiv) |c| loop.val(row, c) orelse "0.0" else "0.0";
        const serial_str = if (cols.id) |c| loop.val(row, c) orelse "0" else "0";
        const ins_code_str = if (cols.pdbx_PDB_ins_code) |c| cif.asString(loop.val(row, c) orelse ".") else "";

        // Chain boundary detection
        if (!in_chain or !std.mem.eql(u8, label_asym, cur_label_asym_id)) {
            // Close previous residue
            if (in_residue) {
                const atom_end: u32 = @intCast(mdl.atoms.items.len);
                const res_idx = mdl.residues.items.len - 1;
                mdl.residues.items[res_idx].atom_end = atom_end;
                in_residue = false;
            }
            // Close previous chain
            if (in_chain) {
                const res_end: u32 = @intCast(mdl.residues.items.len);
                const chain_idx = mdl.chains.items.len - 1;
                mdl.chains.items[chain_idx].residue_end = res_end;
            }
            // Start new chain
            var new_chain = Chain{};
            new_chain.setLabelAsymId(label_asym);
            new_chain.setAuthAsymId(auth_asym);
            const entity_id_str = if (cols.label_entity_id) |c| cif.asString(loop.val(row, c) orelse ".") else "";
            new_chain.setEntityId(entity_id_str);
            new_chain.residue_start = @intCast(mdl.residues.items.len);
            try mdl.chains.append(mdl.allocator, new_chain);
            cur_label_asym_id = label_asym;
            in_chain = true;
            // Force new residue
            cur_seq_id = "";
            cur_comp_id = "";
            cur_ins_code = "";
            cur_auth_seq = "";
        }

        // Residue boundary detection (includes insertion code and auth_seq for non-polymer)
        const new_residue = !in_residue or
            !std.mem.eql(u8, seq_id, cur_seq_id) or
            !std.mem.eql(u8, comp_id, cur_comp_id) or
            !std.mem.eql(u8, ins_code_str, cur_ins_code) or
            !std.mem.eql(u8, auth_seq, cur_auth_seq);

        if (new_residue) {
            // Close previous residue
            if (in_residue) {
                const atom_end: u32 = @intCast(mdl.atoms.items.len);
                const res_idx = mdl.residues.items.len - 1;
                mdl.residues.items[res_idx].atom_end = atom_end;
            }
            // Start new residue
            res_atom_start = @intCast(mdl.atoms.items.len);
            var new_res = Residue{};
            new_res.setCompId(comp_id);
            const label_seq_int = cif.value.asIntOr(i32, seq_id, 0);
            const auth_seq_int = cif.value.asIntOr(i32, auth_seq, 0);
            new_res.seq_id = if (seq_id.len == 0) auth_seq_int else label_seq_int;
            new_res.auth_seq_id = auth_seq_int;
            new_res.ins_code = if (ins_code_str.len > 0) ins_code_str[0] else ' ';
            new_res.chain_idx = @intCast(mdl.chains.items.len - 1);
            new_res.atom_start = res_atom_start;
            new_res.atom_end = res_atom_start; // will be updated
            try mdl.residues.append(mdl.allocator, new_res);
            cur_seq_id = seq_id;
            cur_comp_id = comp_id;
            cur_ins_code = ins_code_str;
            cur_auth_seq = auth_seq;
            in_residue = true;
        }

        // Build atom
        var atom = Atom{
            .pos = .{ .x = x, .y = y, .z = z },
        };
        atom.setName(atom_name);
        atom.element_type = element.elementFromSymbol(type_sym);
        atom.is_hydrogen = switch (atom.element_type) {
            .H, .Har, .Hpol, .Ha_p, .HOd => true,
            else => false,
        };
        atom.occupancy = cif.asFloatOr(occ_str, 1.0);
        atom.b_factor = cif.asFloatOr(bfac_str, 0.0);
        atom.serial = @intCast(cif.value.asIntOr(u32, serial_str, 0));
        atom.residue_idx = @intCast(mdl.residues.items.len - 1);

        // Alt location: first char if not '.'
        const alt_raw = cif.asString(alt_loc_str);
        atom.altloc = if (alt_raw.len > 0) alt_raw[0] else ' ';

        try mdl.atoms.append(mdl.allocator, atom);
    }

    // Close final residue
    if (in_residue) {
        const atom_end: u32 = @intCast(mdl.atoms.items.len);
        const res_idx = mdl.residues.items.len - 1;
        mdl.residues.items[res_idx].atom_end = atom_end;
    }

    // Close final chain
    if (in_chain) {
        const res_end: u32 = @intCast(mdl.residues.items.len);
        const chain_idx = mdl.chains.items.len - 1;
        mdl.chains.items[chain_idx].residue_end = res_end;
    }

    // Parse _pdbx_poly_seq_scheme for chain-break detection (optional)
    if (block.findLoop("_pdbx_poly_seq_scheme.seq_id")) |pss| {
        const col_asym = pss.findTag("_pdbx_poly_seq_scheme.asym_id");
        const col_seq = pss.findTag("_pdbx_poly_seq_scheme.seq_id");
        const col_auth_seq = pss.findTag("_pdbx_poly_seq_scheme.auth_seq_num");

        if (col_asym != null and col_seq != null and col_auth_seq != null) {
            detectChainBreaks(&mdl, pss, col_asym.?, col_seq.?, col_auth_seq.?);
        }
    }

    // Parse _pdbx_unobs_or_zero_occ_atoms count (optional, diagnostic)
    if (block.findLoop("_pdbx_unobs_or_zero_occ_atoms.label_atom_id")) |unobs| {
        mdl.n_unobs_atoms = @intCast(unobs.length());
    }

    // Parse _entity loop for entity types (optional)
    applyEntityTypes(&mdl, block);

    return mdl;
}

// ── Atom Lookup ───────────────────────────────────────────────────────────────

/// Key for looking up a Model atom index from CIF identifiers.
pub const AtomLookupKey = struct {
    label_asym_id: []const u8,
    seq_id: []const u8,
    atom_name: []const u8,
};

/// Hash context for AtomLookupKey using Wyhash with null-byte separators.
pub const AtomLookupContext = struct {
    pub fn hash(_: AtomLookupContext, key: AtomLookupKey) u64 {
        var hasher = std.hash.Wyhash.init(0);
        hasher.update(key.label_asym_id);
        hasher.update(&[_]u8{0});
        hasher.update(key.seq_id);
        hasher.update(&[_]u8{0});
        hasher.update(key.atom_name);
        return hasher.final();
    }

    pub fn eql(_: AtomLookupContext, a: AtomLookupKey, b: AtomLookupKey) bool {
        return std.mem.eql(u8, a.label_asym_id, b.label_asym_id) and
            std.mem.eql(u8, a.seq_id, b.seq_id) and
            std.mem.eql(u8, a.atom_name, b.atom_name);
    }
};

/// HashMap from (label_asym_id, label_seq_id, atom_name) to Model atom index.
pub const AtomLookup = std.HashMap(AtomLookupKey, u32, AtomLookupContext, 80);

/// Build an AtomLookup from a pre-parsed CIF block.
/// Row index in the _atom_site loop == Model atom index (same order as parseModel).
/// For altloc atoms, only the first occurrence is indexed.
/// Also registers auth_seq_id entries for branched entities (label_seq_id == ".").
pub fn buildAtomLookup(allocator: Allocator, block: *const cif.Block) !AtomLookup {
    const loop = block.findLoop("_atom_site.Cartn_x") orelse return MmcifError.NoAtomSiteLoop;

    const col_asym = loop.findTag("_atom_site.label_asym_id") orelse return MmcifError.MissingCoordinateField;
    const col_seq = loop.findTag("_atom_site.label_seq_id") orelse return MmcifError.MissingCoordinateField;
    const col_atom = loop.findTag("_atom_site.label_atom_id") orelse return MmcifError.MissingCoordinateField;
    const col_auth_seq = loop.findTag("_atom_site.auth_seq_id");

    var lookup = AtomLookup.initContext(allocator, AtomLookupContext{});
    errdefer lookup.deinit();

    const nrows = loop.length();
    try lookup.ensureTotalCapacity(@intCast(nrows));

    for (0..nrows) |row| {
        const asym = cif.asString(loop.val(row, col_asym) orelse continue);
        const seq = cif.asString(loop.val(row, col_seq) orelse continue);
        const atom = cif.asString(loop.val(row, col_atom) orelse continue);

        const idx: u32 = @intCast(row);

        // Primary key: label_asym_id + seq_id + atom_name
        // getOrPut skips duplicates (altloc), only first occurrence is indexed.
        const key = AtomLookupKey{
            .label_asym_id = asym,
            .seq_id = seq,
            .atom_name = atom,
        };
        const gop = try lookup.getOrPut(key);
        if (!gop.found_existing) {
            gop.value_ptr.* = idx;
        }

        // Also register auth_seq_id entry for branched entities (label_seq_id == ".")
        if (seq.len == 0) {
            if (col_auth_seq) |c| {
                const auth_seq = cif.asString(loop.val(row, c) orelse ".");
                if (auth_seq.len > 0) {
                    const auth_key = AtomLookupKey{
                        .label_asym_id = asym,
                        .seq_id = auth_seq,
                        .atom_name = atom,
                    };
                    const auth_gop = try lookup.getOrPut(auth_key);
                    if (!auth_gop.found_existing) {
                        auth_gop.value_ptr.* = idx;
                    }
                }
            }
        }
    }

    return lookup;
}

// ── Struct Conn ───────────────────────────────────────────────────────────────

/// Returns true if the connection type is covalent (covale* or disulf).
pub fn isCovalentConnType(conn_type: []const u8) bool {
    var buf: [32]u8 = undefined;
    const len = @min(conn_type.len, buf.len);
    for (0..len) |i| {
        buf[i] = std.ascii.toLower(conn_type[i]);
    }
    const lower = buf[0..len];
    return std.mem.startsWith(u8, lower, "covale") or std.mem.eql(u8, lower, "disulf");
}

/// Parse _struct_conn loop and add inter-residue bonds to Model.bonds.
/// Sets bonded_inter_residue = true on both partner atoms.
pub fn parseStructConn(mdl: *Model, block: *const cif.Block, lookup: *const AtomLookup) !void {
    const sc = block.findLoop("_struct_conn.conn_type_id") orelse return;

    const col_type = sc.findTag("_struct_conn.conn_type_id") orelse return error.MissingCoordinateField;
    const col_asym1 = sc.findTag("_struct_conn.ptnr1_label_asym_id") orelse return error.MissingCoordinateField;
    const col_seq1 = sc.findTag("_struct_conn.ptnr1_label_seq_id") orelse return error.MissingCoordinateField;
    const col_atom1 = sc.findTag("_struct_conn.ptnr1_label_atom_id") orelse return error.MissingCoordinateField;
    const col_asym2 = sc.findTag("_struct_conn.ptnr2_label_asym_id") orelse return error.MissingCoordinateField;
    const col_seq2 = sc.findTag("_struct_conn.ptnr2_label_seq_id") orelse return error.MissingCoordinateField;
    const col_atom2 = sc.findTag("_struct_conn.ptnr2_label_atom_id") orelse return error.MissingCoordinateField;
    const col_sym1 = sc.findTag("_struct_conn.ptnr1_symmetry");
    const col_sym2 = sc.findTag("_struct_conn.ptnr2_symmetry");
    const col_order = sc.findTag("_struct_conn.pdbx_value_order");

    const nrows = sc.length();
    for (0..nrows) |row| {
        const conn_type = cif.asString(sc.val(row, col_type) orelse continue);
        if (!isCovalentConnType(conn_type)) continue;

        // Skip inter-symmetry bonds
        if (col_sym1 != null and col_sym2 != null) {
            const sym1 = cif.asString(sc.val(row, col_sym1.?) orelse ".");
            const sym2 = cif.asString(sc.val(row, col_sym2.?) orelse ".");
            if (sym1.len > 0 and sym2.len > 0 and !std.mem.eql(u8, sym1, sym2)) continue;
        }

        const asym1 = cif.asString(sc.val(row, col_asym1) orelse continue);
        const seq1 = cif.asString(sc.val(row, col_seq1) orelse continue);
        const atom1 = cif.asString(sc.val(row, col_atom1) orelse continue);
        const asym2 = cif.asString(sc.val(row, col_asym2) orelse continue);
        const seq2 = cif.asString(sc.val(row, col_seq2) orelse continue);
        const atom2 = cif.asString(sc.val(row, col_atom2) orelse continue);

        const idx1 = lookup.get(.{ .label_asym_id = asym1, .seq_id = seq1, .atom_name = atom1 }) orelse continue;
        const idx2 = lookup.get(.{ .label_asym_id = asym2, .seq_id = seq2, .atom_name = atom2 }) orelse continue;

        const order_str = if (col_order) |c| cif.asString(sc.val(row, c) orelse ".") else "";
        const order = bond_mod.BondOrder.fromString(order_str);

        try mdl.bonds.append(mdl.allocator, bond_mod.Bond{
            .atom_1 = idx1,
            .atom_2 = idx2,
            .order = order,
            .source = .struct_conn,
        });

        mdl.atoms.items[idx1].flags.bonded_inter_residue = true;
        mdl.atoms.items[idx2].flags.bonded_inter_residue = true;
    }
}

// ── Branch Links ─────────────────────────────────────────────────────────────

/// Parse _pdbx_entity_branch_link loop and add branch bonds to Model.bonds.
/// Sets bonded_inter_residue = true on the leaving atoms of each bond partner.
/// If the leaving atom is not present in the model (e.g. the hydrogen was never
/// modeled), the flag falls back to the bonding atom instead.
pub fn parseBranchLinks(allocator: Allocator, mdl: *Model, block: *const cif.Block, lookup: *const AtomLookup) !void {
    const loop = block.findLoop("_pdbx_entity_branch_link.link_id") orelse return;

    const col_entity = loop.findTag("_pdbx_entity_branch_link.entity_id") orelse return error.MissingCoordinateField;
    const col_num1 = loop.findTag("_pdbx_entity_branch_link.entity_branch_list_num_1") orelse return error.MissingCoordinateField;
    const col_num2 = loop.findTag("_pdbx_entity_branch_link.entity_branch_list_num_2") orelse return error.MissingCoordinateField;
    const col_atom1 = loop.findTag("_pdbx_entity_branch_link.atom_id_1") orelse return error.MissingCoordinateField;
    const col_atom2 = loop.findTag("_pdbx_entity_branch_link.atom_id_2") orelse return error.MissingCoordinateField;
    const col_leaving1 = loop.findTag("_pdbx_entity_branch_link.leaving_atom_id_1") orelse return error.MissingCoordinateField;
    const col_leaving2 = loop.findTag("_pdbx_entity_branch_link.leaving_atom_id_2") orelse return error.MissingCoordinateField;

    // Build entity_id -> [asym_id] mapping from mdl.chains
    var entity_to_asyms = std.StringHashMap(std.ArrayListUnmanaged([]const u8)).init(allocator);
    defer {
        var it = entity_to_asyms.valueIterator();
        while (it.next()) |list| {
            list.deinit(allocator);
        }
        entity_to_asyms.deinit();
    }

    for (mdl.chains.items) |*chain| {
        const entity_id = chain.entityIdSlice();
        const asym_id = chain.labelSlice();
        if (entity_id.len == 0) continue;

        const gop = try entity_to_asyms.getOrPut(entity_id);
        if (!gop.found_existing) {
            gop.value_ptr.* = .empty;
        }
        try gop.value_ptr.append(allocator, asym_id);
    }

    const nrows = loop.length();
    for (0..nrows) |row| {
        const entity_id = cif.asString(loop.val(row, col_entity) orelse continue);
        const num1_str = cif.asString(loop.val(row, col_num1) orelse continue);
        const num2_str = cif.asString(loop.val(row, col_num2) orelse continue);
        const atom_name1 = cif.asString(loop.val(row, col_atom1) orelse continue);
        const atom_name2 = cif.asString(loop.val(row, col_atom2) orelse continue);
        const leaving_name1 = cif.asString(loop.val(row, col_leaving1) orelse continue);
        const leaving_name2 = cif.asString(loop.val(row, col_leaving2) orelse continue);

        const asyms = entity_to_asyms.get(entity_id) orelse continue;

        for (asyms.items) |asym_id| {
            // Resolve bonding atoms using num (= auth_seq_id for branched entities)
            const idx1 = lookup.get(.{
                .label_asym_id = asym_id,
                .seq_id = num1_str,
                .atom_name = atom_name1,
            }) orelse continue;
            const idx2 = lookup.get(.{
                .label_asym_id = asym_id,
                .seq_id = num2_str,
                .atom_name = atom_name2,
            }) orelse continue;

            try mdl.bonds.append(mdl.allocator, bond_mod.Bond{
                .atom_1 = idx1,
                .atom_2 = idx2,
                .order = .single,
                .source = .branch_link,
            });

            // Set bonded_inter_residue on leaving atoms.
            // If the leaving atom is absent (e.g. hydrogen not in model),
            // fall back to the bonding atom.
            const leaving_idx1 = lookup.get(.{
                .label_asym_id = asym_id,
                .seq_id = num1_str,
                .atom_name = leaving_name1,
            }) orelse idx1;
            mdl.atoms.items[leaving_idx1].flags.bonded_inter_residue = true;

            const leaving_idx2 = lookup.get(.{
                .label_asym_id = asym_id,
                .seq_id = num2_str,
                .atom_name = leaving_name2,
            }) orelse idx2;
            mdl.atoms.items[leaving_idx2].flags.bonded_inter_residue = true;
        }
    }
}

// ── Inline Component Parsing ──────────────────────────────────────────────────

/// Helper: find atom index in a slice by name.
fn findAtomIdx(atoms: []const ccd_mod.CompAtom, name: []const u8) ?u16 {
    for (atoms, 0..) |*a, i| {
        if (std.mem.eql(u8, a.nameSlice(), name)) return @intCast(i);
    }
    return null;
}

/// Temporary bond record storing atom names (strings) before index resolution.
const RawBond = struct {
    atom1: []const u8,
    atom2: []const u8,
    order: ccd_mod.BondOrder,
    aromatic: bool,
};

/// Parse inline `_chem_comp_atom` and `_chem_comp_bond` loops from a structure
/// CIF block into a `ccd_mod.ComponentDict`.
///
/// Returns `null` when neither loop is present in the block.
/// The returned ComponentDict owns all its memory.
pub fn parseInlineComponents(allocator: Allocator, block: *const cif.Block) !?ccd_mod.ComponentDict {
    // Check that at least the bond loop is present; if not, return null.
    const bond_loop = block.findLoop("_chem_comp_bond.comp_id") orelse return null;

    // ── Phase 1: group atoms by comp_id ──────────────────────────────────────

    // comp_id -> list of CompAtom (temporary, unowned strings)
    var atom_groups = std.StringHashMap(std.ArrayListUnmanaged(ccd_mod.CompAtom)).init(allocator);
    defer {
        var it = atom_groups.valueIterator();
        while (it.next()) |list| list.deinit(allocator);
        atom_groups.deinit();
    }

    if (block.findLoop("_chem_comp_atom.comp_id")) |atom_loop| {
        const col_comp = atom_loop.findTag("_chem_comp_atom.comp_id") orelse return null;
        const col_atom = atom_loop.findTag("_chem_comp_atom.atom_id");
        const col_sym = atom_loop.findTag("_chem_comp_atom.type_symbol");
        const col_charge = atom_loop.findTag("_chem_comp_atom.charge");
        const col_leaving = atom_loop.findTag("_chem_comp_atom.pdbx_leaving_atom_flag");
        const col_arom = atom_loop.findTag("_chem_comp_atom.pdbx_aromatic_flag");

        const nrows = atom_loop.length();
        for (0..nrows) |row| {
            const comp_id = cif.asString(atom_loop.val(row, col_comp) orelse continue);

            var atom = ccd_mod.CompAtom{};

            if (col_atom) |c| {
                if (atom_loop.val(row, c)) |v| {
                    const s = cif.asString(v);
                    const len = @min(s.len, 4);
                    atom.name_len = @intCast(len);
                    @memcpy(atom.name[0..len], s[0..len]);
                }
            }

            if (col_sym) |c| {
                if (atom_loop.val(row, c)) |v| {
                    const s = cif.asString(v);
                    const len = @min(s.len, 2);
                    @memcpy(atom.element_symbol[0..len], s[0..len]);
                }
            }

            if (col_charge) |c| {
                if (atom_loop.val(row, c)) |v| {
                    atom.charge = cif.value.asIntOr(i8, v, 0);
                }
            }

            if (col_leaving) |c| {
                if (atom_loop.val(row, c)) |v| {
                    atom.leaving = std.ascii.eqlIgnoreCase(cif.asString(v), "Y");
                }
            }

            if (col_arom) |c| {
                if (atom_loop.val(row, c)) |v| {
                    atom.aromatic = std.ascii.eqlIgnoreCase(cif.asString(v), "Y");
                }
            }

            const gop = try atom_groups.getOrPut(comp_id);
            if (!gop.found_existing) gop.value_ptr.* = .{};
            try gop.value_ptr.append(allocator, atom);
        }
    }

    // ── Phase 2: group bonds by comp_id ──────────────────────────────────────

    // comp_id -> list of RawBond (strings are slices into the CIF source)
    var bond_groups = std.StringHashMap(std.ArrayListUnmanaged(RawBond)).init(allocator);
    defer {
        var it = bond_groups.valueIterator();
        while (it.next()) |list| list.deinit(allocator);
        bond_groups.deinit();
    }

    {
        const col_comp = bond_loop.findTag("_chem_comp_bond.comp_id") orelse return null;
        const col_a1 = bond_loop.findTag("_chem_comp_bond.atom_id_1");
        const col_a2 = bond_loop.findTag("_chem_comp_bond.atom_id_2");
        const col_order = bond_loop.findTag("_chem_comp_bond.value_order");
        const col_arom = bond_loop.findTag("_chem_comp_bond.pdbx_aromatic_flag");

        const nrows = bond_loop.length();
        for (0..nrows) |row| {
            const comp_id = cif.asString(bond_loop.val(row, col_comp) orelse continue);

            const atom1 = if (col_a1) |c|
                cif.asString(bond_loop.val(row, c) orelse continue)
            else
                continue;

            const atom2 = if (col_a2) |c|
                cif.asString(bond_loop.val(row, c) orelse continue)
            else
                continue;

            const order: ccd_mod.BondOrder = if (col_order) |c|
                if (bond_loop.val(row, c)) |v| ccd_mod.BondOrder.fromString(cif.asString(v)) else .unknown
            else
                .unknown;

            const aromatic: bool = if (col_arom) |c|
                if (bond_loop.val(row, c)) |v| std.ascii.eqlIgnoreCase(cif.asString(v), "Y") else false
            else
                false;

            const gop = try bond_groups.getOrPut(comp_id);
            if (!gop.found_existing) gop.value_ptr.* = .{};
            try gop.value_ptr.append(allocator, RawBond{
                .atom1 = atom1,
                .atom2 = atom2,
                .order = order,
                .aromatic = aromatic,
            });
        }
    }

    // ── Phase 3: build ComponentDict ─────────────────────────────────────────

    var dict = ccd_mod.ComponentDict{
        .components = std.StringHashMap(ccd_mod.Component).init(allocator),
        .allocator = allocator,
    };
    errdefer dict.deinit();

    // Collect all unique comp_ids from both groups.
    var all_ids = std.StringHashMap(void).init(allocator);
    defer all_ids.deinit();
    {
        var it = atom_groups.keyIterator();
        while (it.next()) |k| try all_ids.put(k.*, {});
    }
    {
        var it = bond_groups.keyIterator();
        while (it.next()) |k| try all_ids.put(k.*, {});
    }

    var id_it = all_ids.keyIterator();
    while (id_it.next()) |id_ptr| {
        const comp_id = id_ptr.*;

        // Atoms slice (may be empty if no atom loop).
        const raw_atoms: []const ccd_mod.CompAtom = if (atom_groups.get(comp_id)) |list|
            list.items
        else
            &[_]ccd_mod.CompAtom{};

        const owned_atoms = try allocator.dupe(ccd_mod.CompAtom, raw_atoms);
        errdefer allocator.free(owned_atoms);

        // Resolve bonds to index-based CompBond.
        var resolved_bonds = std.ArrayListUnmanaged(ccd_mod.CompBond){};
        defer resolved_bonds.deinit(allocator);

        if (bond_groups.get(comp_id)) |raw_bond_list| {
            for (raw_bond_list.items) |rb| {
                const idx1 = findAtomIdx(owned_atoms, rb.atom1) orelse continue;
                const idx2 = findAtomIdx(owned_atoms, rb.atom2) orelse continue;
                try resolved_bonds.append(allocator, ccd_mod.CompBond{
                    .atom_idx_1 = idx1,
                    .atom_idx_2 = idx2,
                    .order = rb.order,
                    .aromatic = rb.aromatic,
                });
            }
        }

        const owned_bonds = try allocator.dupe(ccd_mod.CompBond, resolved_bonds.items);
        errdefer allocator.free(owned_bonds);

        const owned_id = try allocator.dupe(u8, comp_id);
        errdefer allocator.free(owned_id);

        const owned_type = try allocator.dupe(u8, "");
        errdefer allocator.free(owned_type);

        const key = try allocator.dupe(u8, comp_id);
        errdefer allocator.free(key);

        try dict.components.put(key, ccd_mod.Component{
            .comp_id = owned_id,
            .comp_type = owned_type,
            .atoms = owned_atoms,
            .bonds = owned_bonds,
        });
    }

    return dict;
}

// ── Tests ─────────────────────────────────────────────────────────────────────

const testing = std.testing;

test "parse tiny mmCIF" {
    const source = @embedFile("test_data/tiny.cif");
    var mdl = try parseModel(testing.allocator, source);
    defer mdl.deinit();

    try testing.expectEqual(@as(usize, 5), mdl.atoms.items.len);
    try testing.expectEqual(@as(usize, 1), mdl.residues.items.len);
    try testing.expectEqual(@as(usize, 1), mdl.chains.items.len);
    try testing.expectApproxEqAbs(mdl.atoms.items[0].pos.x, 1.0, 1e-3);
    try testing.expectApproxEqAbs(mdl.atoms.items[0].pos.y, 2.0, 1e-3);
    try testing.expectEqualStrings("ALA", mdl.residues.items[0].compIdSlice());
}

test "parse multi-chain multi-residue mmCIF" {
    const source = @embedFile("test_data/multi_chain.cif");
    var mdl = try parseModel(testing.allocator, source);
    defer mdl.deinit();

    // 11 atoms total
    try testing.expectEqual(@as(usize, 11), mdl.atoms.items.len);
    // 3 residues: ALA(A/1), GLY(A/2), VAL(B/1)
    try testing.expectEqual(@as(usize, 3), mdl.residues.items.len);
    // 2 chains: A and B
    try testing.expectEqual(@as(usize, 2), mdl.chains.items.len);

    // Chain A: residues 0-1, atoms 0-7
    try testing.expectEqualStrings("A", mdl.chains.items[0].labelSlice());
    try testing.expectEqual(@as(u32, 0), mdl.chains.items[0].residue_start);
    try testing.expectEqual(@as(u32, 2), mdl.chains.items[0].residue_end);

    // Chain B: residue 2, atoms 8-10
    try testing.expectEqualStrings("B", mdl.chains.items[1].labelSlice());
    try testing.expectEqual(@as(u32, 2), mdl.chains.items[1].residue_start);
    try testing.expectEqual(@as(u32, 3), mdl.chains.items[1].residue_end);

    // Residue 0: ALA, atoms 0-3
    try testing.expectEqualStrings("ALA", mdl.residues.items[0].compIdSlice());
    try testing.expectEqual(@as(u32, 0), mdl.residues.items[0].atom_start);
    try testing.expectEqual(@as(u32, 4), mdl.residues.items[0].atom_end);

    // Residue 1: GLY, atoms 4-7
    try testing.expectEqualStrings("GLY", mdl.residues.items[1].compIdSlice());
    try testing.expectEqual(@as(u32, 4), mdl.residues.items[1].atom_start);
    try testing.expectEqual(@as(u32, 8), mdl.residues.items[1].atom_end);

    // Residue 2: VAL, atoms 8-10
    try testing.expectEqualStrings("VAL", mdl.residues.items[2].compIdSlice());
    try testing.expectEqual(@as(u32, 8), mdl.residues.items[2].atom_start);
    try testing.expectEqual(@as(u32, 11), mdl.residues.items[2].atom_end);

    // Last atom coordinate check
    try testing.expectApproxEqAbs(mdl.atoms.items[10].pos.x, 13.0, 1e-3);
}

test "error: no atom_site loop" {
    const source = "data_EMPTY\n_entry.id EMPTY\n";
    const result = parseModel(testing.allocator, source);
    try testing.expectError(MmcifError.NoAtomSiteLoop, result);
}

test "error: missing coordinate field" {
    const source =
        \\data_BAD
        \\loop_
        \\_atom_site.label_atom_id
        \\CA
    ;
    const result = parseModel(testing.allocator, source);
    // No Cartn_x → NoAtomSiteLoop (findLoop won't find it)
    try testing.expectError(MmcifError.NoAtomSiteLoop, result);
}

test "error: invalid coordinate value" {
    const source =
        \\data_BAD
        \\loop_
        \\_atom_site.Cartn_x
        \\_atom_site.Cartn_y
        \\_atom_site.Cartn_z
        \\1.0 2.0 abc
    ;
    const result = parseModel(testing.allocator, source);
    try testing.expectError(MmcifError.InvalidCoordinateValue, result);
}

test "parse chain break from pdbx_poly_seq_scheme" {
    const source = @embedFile("test_data/gap_chain.cif");
    var mdl = try parseModel(testing.allocator, source);
    defer mdl.deinit();

    try testing.expectEqual(@as(usize, 2), mdl.residues.items.len);
    try testing.expectEqual(@as(usize, 1), mdl.chains.items.len);

    // First residue: no chain break before
    try testing.expect(!mdl.residues.items[0].is_chain_break_before);
    // Second residue: chain break before (seq_id 2 is unobserved)
    try testing.expect(mdl.residues.items[1].is_chain_break_before);
}

test "parse without pdbx_poly_seq_scheme is backward compatible" {
    const source = @embedFile("test_data/tiny.cif");
    var mdl = try parseModel(testing.allocator, source);
    defer mdl.deinit();

    try testing.expect(!mdl.residues.items[0].is_chain_break_before);
}

test "model reports zero unobs atoms for tiny.cif" {
    const source = @embedFile("test_data/tiny.cif");
    var mdl = try parseModel(testing.allocator, source);
    defer mdl.deinit();

    try testing.expectEqual(@as(u32, 0), mdl.n_unobs_atoms);
}

test "insertion code splits residues with same seq_id" {
    const source = @embedFile("test_data/ins_code.cif");
    var mdl = try parseModel(testing.allocator, source);
    defer mdl.deinit();

    // Same seq_id=1 but different ins_codes (blank vs 'A') → 2 residues
    try testing.expectEqual(@as(usize, 2), mdl.residues.items.len);
    try testing.expectEqual(@as(usize, 1), mdl.chains.items.len);

    // First residue: ALA with ins_code=' '
    try testing.expectEqualStrings("ALA", mdl.residues.items[0].compIdSlice());
    try testing.expectEqual(@as(u8, ' '), mdl.residues.items[0].ins_code);

    // Second residue: GLY with ins_code='A'
    try testing.expectEqualStrings("GLY", mdl.residues.items[1].compIdSlice());
    try testing.expectEqual(@as(u8, 'A'), mdl.residues.items[1].ins_code);
}

test "buildAtomLookup resolves atom indices" {
    const source = @embedFile("test_data/disulfide.cif");
    var doc = cif.readString(testing.allocator, source) catch unreachable;
    defer doc.deinit();
    const block = &doc.blocks.items[0];

    var mdl = try parseModel(testing.allocator, source);
    defer mdl.deinit();

    var lookup = try buildAtomLookup(testing.allocator, block);
    defer lookup.deinit();

    // SG of CYS residue 1 (seq_id=1) should be atom index 5
    const sg1 = lookup.get(.{ .label_asym_id = "A", .seq_id = "1", .atom_name = "SG" });
    try testing.expect(sg1 != null);
    try testing.expectEqual(@as(u32, 5), sg1.?);

    // SG of CYS residue 2 (seq_id=2) should be atom index 11
    const sg2 = lookup.get(.{ .label_asym_id = "A", .seq_id = "2", .atom_name = "SG" });
    try testing.expect(sg2 != null);
    try testing.expectEqual(@as(u32, 11), sg2.?);
}

test "parseStructConn disulfide bond" {
    const source = @embedFile("test_data/disulfide.cif");
    var doc = cif.readString(testing.allocator, source) catch unreachable;
    defer doc.deinit();
    const block = &doc.blocks.items[0];

    var mdl = try parseModel(testing.allocator, source);
    defer mdl.deinit();

    var lookup = try buildAtomLookup(testing.allocator, block);
    defer lookup.deinit();

    try parseStructConn(&mdl, block, &lookup);

    // Should have 1 bond (disulfide SG-SG)
    try testing.expectEqual(@as(usize, 1), mdl.bonds.items.len);
    const bond = mdl.bonds.items[0];
    try testing.expectEqual(bond_mod.BondSource.struct_conn, bond.source);

    // Both SG atoms should have bonded_inter_residue flag
    try testing.expect(mdl.atoms.items[bond.atom_1].flags.bonded_inter_residue);
    try testing.expect(mdl.atoms.items[bond.atom_2].flags.bonded_inter_residue);
}

test "parseBranchLinks glycan bond" {
    const source = @embedFile("test_data/branch_link.cif");
    var doc = cif.readString(testing.allocator, source) catch unreachable;
    defer doc.deinit();
    const block = &doc.blocks.items[0];

    var mdl = try parseModel(testing.allocator, source);
    defer mdl.deinit();

    var lookup = try buildAtomLookup(testing.allocator, block);
    defer lookup.deinit();

    try parseBranchLinks(testing.allocator, &mdl, block, &lookup);

    // Should have 1 bond (NAG O4 — GAL C1)
    try testing.expectEqual(@as(usize, 1), mdl.bonds.items.len);
    const bond = mdl.bonds.items[0];
    try testing.expectEqual(bond_mod.BondSource.branch_link, bond.source);

    // Leaving atoms should have bonded_inter_residue flag.
    // leaving_atom_id_1=HO4 is absent in model → falls back to bonding atom O4 (index 1).
    // leaving_atom_id_2=O1 is present at index 3.
    try testing.expect(mdl.atoms.items[1].flags.bonded_inter_residue); // NAG O4 (leaving fallback)
    try testing.expect(mdl.atoms.items[3].flags.bonded_inter_residue); // GAL O1 (leaving)
}

test "parseInlineComponents returns ComponentDict" {
    const source = @embedFile("test_data/inline_comp.cif");
    var doc = cif.readString(testing.allocator, source) catch unreachable;
    defer doc.deinit();
    const block = &doc.blocks.items[0];

    var dict = try parseInlineComponents(testing.allocator, block);
    defer if (dict) |*d| d.deinit();

    try testing.expect(dict != null);
    const ala = dict.?.get("ALA");
    try testing.expect(ala != null);
    try testing.expectEqual(@as(usize, 11), ala.?.atoms.len);
    try testing.expectEqual(@as(usize, 10), ala.?.bonds.len);
}

test "parseInlineComponents returns null when no inline data" {
    const source = @embedFile("test_data/tiny.cif");
    var doc = cif.readString(testing.allocator, source) catch unreachable;
    defer doc.deinit();
    const block = &doc.blocks.items[0];

    var dict = try parseInlineComponents(testing.allocator, block);
    defer if (dict) |*d| d.deinit();

    try testing.expect(dict == null);
}

test "isCovalentConnType accepts covalent types" {
    try testing.expect(isCovalentConnType("disulf"));
    try testing.expect(isCovalentConnType("DISULF"));
    try testing.expect(isCovalentConnType("Disulf"));
    try testing.expect(isCovalentConnType("covale"));
    try testing.expect(isCovalentConnType("COVALE"));
    try testing.expect(isCovalentConnType("covale_base"));
    try testing.expect(isCovalentConnType("covale_phosph"));
}

test "isCovalentConnType rejects non-covalent types" {
    try testing.expect(!isCovalentConnType("hydrog"));
    try testing.expect(!isCovalentConnType("metalc"));
    try testing.expect(!isCovalentConnType("mismat"));
    try testing.expect(!isCovalentConnType("saltbr"));
    try testing.expect(!isCovalentConnType(""));
}

test "parseStructConn skips non-covalent connections" {
    const source = @embedFile("test_data/disulfide_with_hbond.cif");
    var doc = cif.readString(testing.allocator, source) catch unreachable;
    defer doc.deinit();
    const block = &doc.blocks.items[0];

    var mdl = try parseModel(testing.allocator, source);
    defer mdl.deinit();

    var lookup = try buildAtomLookup(testing.allocator, block);
    defer lookup.deinit();

    try parseStructConn(&mdl, block, &lookup);

    // Should have exactly 1 bond (disulfide SG-SG), hydrogen bond must be skipped
    try testing.expectEqual(@as(usize, 1), mdl.bonds.items.len);

    // Verify the bond is the disulfide (SG-SG), not the hydrogen bond (N-O)
    const bond = mdl.bonds.items[0];
    const a1 = mdl.atoms.items[bond.atom_1];
    const a2 = mdl.atoms.items[bond.atom_2];
    try testing.expectEqualStrings("SG", a1.nameSlice());
    try testing.expectEqualStrings("SG", a2.nameSlice());
}

test "parseModel stores auth_seq_id and falls back" {
    const source = @embedFile("test_data/entity_type.cif");
    var mdl = try parseModel(testing.allocator, source);
    defer mdl.deinit();

    // Polymer ALA: label_seq_id=1, auth_seq_id=1
    try testing.expectEqual(@as(i32, 1), mdl.residues.items[0].seq_id);
    try testing.expectEqual(@as(i32, 1), mdl.residues.items[0].auth_seq_id);

    // Non-polymer EDO: label_seq_id="." -> fallback to auth_seq_id=1
    try testing.expectEqual(@as(i32, 1), mdl.residues.items[1].seq_id);
    try testing.expectEqual(@as(i32, 1), mdl.residues.items[1].auth_seq_id);

    // Water HOH: label_seq_id="." -> fallback to auth_seq_id=1
    try testing.expectEqual(@as(i32, 1), mdl.residues.items[2].seq_id);

    // Branched NAG residue 1: label_seq_id="." -> fallback to auth_seq_id=1
    try testing.expectEqual(@as(i32, 1), mdl.residues.items[3].seq_id);
    try testing.expectEqual(@as(i32, 1), mdl.residues.items[3].auth_seq_id);

    // Branched NAG residue 2: label_seq_id="." -> fallback to auth_seq_id=2
    try testing.expectEqual(@as(i32, 2), mdl.residues.items[4].seq_id);
    try testing.expectEqual(@as(i32, 2), mdl.residues.items[4].auth_seq_id);
}

test "parseModel sets entity_type from _entity loop" {
    const source = @embedFile("test_data/entity_type.cif");
    var mdl = try parseModel(testing.allocator, source);
    defer mdl.deinit();

    // 5 residues: ALA, EDO, HOH, NAG(1), NAG(2)
    try testing.expectEqual(@as(usize, 5), mdl.residues.items.len);

    // ALA (entity 1) -> polymer
    try testing.expectEqual(EntityType.polymer, mdl.residues.items[0].entity_type);

    // EDO (entity 2) -> non_polymer
    try testing.expectEqual(EntityType.non_polymer, mdl.residues.items[1].entity_type);

    // HOH (entity 3) -> water
    try testing.expectEqual(EntityType.water, mdl.residues.items[2].entity_type);

    // NAG (entity 4) -> branched
    try testing.expectEqual(EntityType.branched, mdl.residues.items[3].entity_type);
    try testing.expectEqual(EntityType.branched, mdl.residues.items[4].entity_type);
}

test "parseModel without _entity loop keeps entity_type unknown" {
    const source = @embedFile("test_data/tiny.cif");
    var mdl = try parseModel(testing.allocator, source);
    defer mdl.deinit();

    // tiny.cif has no _entity loop — all residues should stay .unknown
    for (mdl.residues.items) |res| {
        try testing.expectEqual(EntityType.unknown, res.entity_type);
    }
}

test "entityTypeFromString maps correctly" {
    try testing.expectEqual(EntityType.polymer, entityTypeFromString("polymer"));
    try testing.expectEqual(EntityType.polymer, entityTypeFromString("POLYMER"));
    try testing.expectEqual(EntityType.non_polymer, entityTypeFromString("non-polymer"));
    try testing.expectEqual(EntityType.branched, entityTypeFromString("branched"));
    try testing.expectEqual(EntityType.water, entityTypeFromString("water"));
    try testing.expectEqual(EntityType.unknown, entityTypeFromString("macrolide"));
    try testing.expectEqual(EntityType.unknown, entityTypeFromString(""));
}
