//! mmCIF extraction: parses _atom_site loop data into a Model.

const std = @import("std");
const Allocator = std.mem.Allocator;

const cif = @import("cif.zig");
const model_mod = @import("model.zig");
const element = @import("element.zig");

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
        }

        // Residue boundary detection (includes insertion code)
        const new_residue = !in_residue or
            !std.mem.eql(u8, seq_id, cur_seq_id) or
            !std.mem.eql(u8, comp_id, cur_comp_id) or
            !std.mem.eql(u8, ins_code_str, cur_ins_code);

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
            new_res.seq_id = cif.value.asIntOr(i32, seq_id, 0);
            new_res.ins_code = if (ins_code_str.len > 0) ins_code_str[0] else ' ';
            new_res.chain_idx = @intCast(mdl.chains.items.len - 1);
            new_res.atom_start = res_atom_start;
            new_res.atom_end = res_atom_start; // will be updated
            try mdl.residues.append(mdl.allocator, new_res);
            cur_seq_id = seq_id;
            cur_comp_id = comp_id;
            cur_ins_code = ins_code_str;
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

    return mdl;
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
