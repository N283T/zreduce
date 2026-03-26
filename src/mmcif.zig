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
};

/// Parse an mmCIF source string and extract all _atom_site records into a Model.
pub fn parseModel(allocator: Allocator, source: []const u8) MmcifError!Model {
    var doc = cif.readString(allocator, source) catch return MmcifError.OutOfMemory;
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
    var res_atom_start: u32 = 0;
    var in_chain = false;
    var in_residue = false;

    for (0..nrows) |row| {
        const x_str = loop.val(row, cols.cartn_x.?) orelse "0";
        const y_str = loop.val(row, cols.cartn_y.?) orelse "0";
        const z_str = loop.val(row, cols.cartn_z.?) orelse "0";

        const x = cif.asFloatOr(x_str, 0.0);
        const y = cif.asFloatOr(y_str, 0.0);
        const z = cif.asFloatOr(z_str, 0.0);

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
            new_chain.residue_start = @intCast(mdl.residues.items.len);
            try mdl.chains.append(mdl.allocator, new_chain);
            cur_label_asym_id = label_asym;
            in_chain = true;
            // Force new residue
            cur_seq_id = "";
            cur_comp_id = "";
        }

        // Residue boundary detection
        const new_residue = !in_residue or
            !std.mem.eql(u8, seq_id, cur_seq_id) or
            !std.mem.eql(u8, comp_id, cur_comp_id);

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
            new_res.chain_idx = @intCast(mdl.chains.items.len - 1);
            new_res.atom_start = res_atom_start;
            new_res.atom_end = res_atom_start; // will be updated
            try mdl.residues.append(mdl.allocator, new_res);
            cur_seq_id = seq_id;
            cur_comp_id = comp_id;
            in_residue = true;
        }

        // Build atom
        var atom = Atom{
            .pos = .{ .x = x, .y = y, .z = z },
        };
        atom.setName(atom_name);
        atom.element_type = element.elementFromSymbol(type_sym);
        atom.is_hydrogen = type_sym.len > 0 and (type_sym[0] == 'H' or type_sym[0] == 'h');
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
