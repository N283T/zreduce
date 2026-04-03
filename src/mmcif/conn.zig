//! Inter-residue connection parsing: _struct_conn and _pdbx_entity_branch_link.

const std = @import("std");
const Allocator = std.mem.Allocator;

const cif = @import("../cif.zig");
const model_mod = @import("../model.zig");
const bond_mod = @import("../model/bond.zig");

const Model = model_mod.Model;

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
    const loop = block.findLoop("_atom_site.Cartn_x") orelse return error.NoAtomSiteLoop;

    const col_asym = loop.findTag("_atom_site.label_asym_id") orelse return error.MissingRequiredField;
    const col_seq = loop.findTag("_atom_site.label_seq_id") orelse return error.MissingRequiredField;
    const col_atom = loop.findTag("_atom_site.label_atom_id") orelse return error.MissingRequiredField;
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

/// Build an AtomLookup for a specific row range within the _atom_site loop.
/// Values stored are local indices: (row - row_start), so they map to the
/// corresponding Model's atom array. Used for multi-model support.
pub fn buildAtomLookupForRange(allocator: Allocator, block: *const cif.Block, row_start: u32, row_end: u32) !AtomLookup {
    const loop = block.findLoop("_atom_site.Cartn_x") orelse return error.NoAtomSiteLoop;

    const col_asym = loop.findTag("_atom_site.label_asym_id") orelse return error.MissingRequiredField;
    const col_seq = loop.findTag("_atom_site.label_seq_id") orelse return error.MissingRequiredField;
    const col_atom = loop.findTag("_atom_site.label_atom_id") orelse return error.MissingRequiredField;
    const col_auth_seq = loop.findTag("_atom_site.auth_seq_id");

    var lookup = AtomLookup.initContext(allocator, AtomLookupContext{});
    errdefer lookup.deinit();

    const n = row_end - row_start;
    try lookup.ensureTotalCapacity(@intCast(n));

    for (row_start..row_end) |row| {
        const asym = cif.asString(loop.val(row, col_asym) orelse continue);
        const seq = cif.asString(loop.val(row, col_seq) orelse continue);
        const atom_name = cif.asString(loop.val(row, col_atom) orelse continue);

        const idx: u32 = @intCast(row - row_start);

        const key = AtomLookupKey{
            .label_asym_id = asym,
            .seq_id = seq,
            .atom_name = atom_name,
        };
        const gop = try lookup.getOrPut(key);
        if (!gop.found_existing) {
            gop.value_ptr.* = idx;
        }

        if (seq.len == 0) {
            if (col_auth_seq) |c| {
                const auth_seq = cif.asString(loop.val(row, c) orelse ".");
                if (auth_seq.len > 0) {
                    const auth_key = AtomLookupKey{
                        .label_asym_id = asym,
                        .seq_id = auth_seq,
                        .atom_name = atom_name,
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

/// Returns true if the connection type should be treated as a structural bond
/// for hydrogen placement purposes: covale*, disulf, or metalc.
/// metalc (metal coordination) bonds are included so that coordinating atoms
/// (e.g. CYS SG, HIS NE2) receive the bonded_inter_residue flag and are
/// excluded from hydrogen placement.
pub fn isCovalentConnType(conn_type: []const u8) bool {
    var buf: [32]u8 = undefined;
    const len = @min(conn_type.len, buf.len);
    for (0..len) |i| {
        buf[i] = std.ascii.toLower(conn_type[i]);
    }
    const lower = buf[0..len];
    if (std.mem.startsWith(u8, lower, "covale")) return true;
    if (std.mem.eql(u8, lower, "disulf")) return true;
    if (std.mem.eql(u8, lower, "metalc")) return true; // metal coordination bonds
    return false;
}

/// Parse _struct_conn loop and add inter-residue bonds to Model.bonds.
/// Sets bonded_inter_residue = true on both partner atoms.
pub fn parseStructConn(mdl: *Model, block: *const cif.Block, lookup: *const AtomLookup) !void {
    const sc = block.findLoop("_struct_conn.conn_type_id") orelse return;

    const col_type = sc.findTag("_struct_conn.conn_type_id") orelse return error.MissingRequiredField;
    const col_asym1 = sc.findTag("_struct_conn.ptnr1_label_asym_id") orelse return error.MissingRequiredField;
    const col_seq1 = sc.findTag("_struct_conn.ptnr1_label_seq_id") orelse return error.MissingRequiredField;
    const col_atom1 = sc.findTag("_struct_conn.ptnr1_label_atom_id") orelse return error.MissingRequiredField;
    const col_asym2 = sc.findTag("_struct_conn.ptnr2_label_asym_id") orelse return error.MissingRequiredField;
    const col_seq2 = sc.findTag("_struct_conn.ptnr2_label_seq_id") orelse return error.MissingRequiredField;
    const col_atom2 = sc.findTag("_struct_conn.ptnr2_label_atom_id") orelse return error.MissingRequiredField;
    const col_sym1 = sc.findTag("_struct_conn.ptnr1_symmetry");
    const col_sym2 = sc.findTag("_struct_conn.ptnr2_symmetry");
    const col_order = sc.findTag("_struct_conn.pdbx_value_order");
    // auth_seq_id columns: fallback for branched entities where label_seq_id is "."
    const col_auth_seq1 = sc.findTag("_struct_conn.ptnr1_auth_seq_id");
    const col_auth_seq2 = sc.findTag("_struct_conn.ptnr2_auth_seq_id");

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

        // For branched entities, label_seq_id is "." (parsed as "" by asString).
        // Fall back to auth_seq_id to disambiguate residues in the same chain.
        // Keep label_asym_id — buildAtomLookup indexes by (label_asym_id, auth_seq_id).
        const eff_seq1 = if (seq1.len == 0 and col_auth_seq1 != null)
            cif.asString(sc.val(row, col_auth_seq1.?) orelse continue)
        else
            seq1;
        const eff_seq2 = if (seq2.len == 0 and col_auth_seq2 != null)
            cif.asString(sc.val(row, col_auth_seq2.?) orelse continue)
        else
            seq2;

        const idx1 = lookup.get(.{ .label_asym_id = asym1, .seq_id = eff_seq1, .atom_name = atom1 }) orelse continue;
        const idx2 = lookup.get(.{ .label_asym_id = asym2, .seq_id = eff_seq2, .atom_name = atom2 }) orelse continue;

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

/// Parse _pdbx_entity_branch_link loop and add branch bonds to Model.bonds.
/// Sets bonded_inter_residue = true on the leaving atoms of each bond partner.
/// If the leaving atom is not present in the model (e.g. the hydrogen was never
/// modeled), the flag falls back to the bonding atom instead.
pub fn parseBranchLinks(allocator: Allocator, mdl: *Model, block: *const cif.Block, lookup: *const AtomLookup) !void {
    const loop = block.findLoop("_pdbx_entity_branch_link.link_id") orelse return;

    const col_entity = loop.findTag("_pdbx_entity_branch_link.entity_id") orelse return error.MissingRequiredField;
    const col_num1 = loop.findTag("_pdbx_entity_branch_link.entity_branch_list_num_1") orelse return error.MissingRequiredField;
    const col_num2 = loop.findTag("_pdbx_entity_branch_link.entity_branch_list_num_2") orelse return error.MissingRequiredField;
    const col_atom1 = loop.findTag("_pdbx_entity_branch_link.atom_id_1") orelse return error.MissingRequiredField;
    const col_atom2 = loop.findTag("_pdbx_entity_branch_link.atom_id_2") orelse return error.MissingRequiredField;
    const col_leaving1 = loop.findTag("_pdbx_entity_branch_link.leaving_atom_id_1") orelse return error.MissingRequiredField;
    const col_leaving2 = loop.findTag("_pdbx_entity_branch_link.leaving_atom_id_2") orelse return error.MissingRequiredField;

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

// ── Tests ─────────────────────────────────────────────────────────────────────

const testing = std.testing;

test "buildAtomLookup resolves atom indices" {
    const mmcif = @import("../mmcif.zig");
    const source = @embedFile("../test_data/disulfide.cif");
    var doc = cif.readString(testing.allocator, source) catch unreachable;
    defer doc.deinit();
    const block = &doc.blocks.items[0];

    var mdl = try mmcif.parseModel(testing.allocator, source);
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
    const mmcif = @import("../mmcif.zig");
    const source = @embedFile("../test_data/disulfide.cif");
    var doc = cif.readString(testing.allocator, source) catch unreachable;
    defer doc.deinit();
    const block = &doc.blocks.items[0];

    var mdl = try mmcif.parseModel(testing.allocator, source);
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
    const mmcif = @import("../mmcif.zig");
    const source = @embedFile("../test_data/branch_link.cif");
    var doc = cif.readString(testing.allocator, source) catch unreachable;
    defer doc.deinit();
    const block = &doc.blocks.items[0];

    var mdl = try mmcif.parseModel(testing.allocator, source);
    defer mdl.deinit();

    var lookup = try buildAtomLookup(testing.allocator, block);
    defer lookup.deinit();

    try parseBranchLinks(testing.allocator, &mdl, block, &lookup);

    // Should have 1 bond (NAG O4 — GAL C1)
    try testing.expectEqual(@as(usize, 1), mdl.bonds.items.len);
    const bond = mdl.bonds.items[0];
    try testing.expectEqual(bond_mod.BondSource.branch_link, bond.source);

    // Leaving atoms should have bonded_inter_residue flag.
    // leaving_atom_id_1=HO4 is absent in model -> falls back to bonding atom O4 (index 1).
    // leaving_atom_id_2=O1 is present at index 3.
    try testing.expect(mdl.atoms.items[1].flags.bonded_inter_residue); // NAG O4 (leaving fallback)
    try testing.expect(mdl.atoms.items[3].flags.bonded_inter_residue); // GAL O1 (leaving)
}

test "isCovalentConnType accepts covalent types" {
    try testing.expect(isCovalentConnType("disulf"));
    try testing.expect(isCovalentConnType("DISULF"));
    try testing.expect(isCovalentConnType("Disulf"));
    try testing.expect(isCovalentConnType("covale"));
    try testing.expect(isCovalentConnType("COVALE"));
    try testing.expect(isCovalentConnType("covale_base"));
    try testing.expect(isCovalentConnType("covale_phosph"));
    try testing.expect(isCovalentConnType("metalc"));
    try testing.expect(isCovalentConnType("METALC"));
    try testing.expect(isCovalentConnType("Metalc"));
}

test "isCovalentConnType rejects non-covalent types" {
    try testing.expect(!isCovalentConnType("hydrog"));
    try testing.expect(!isCovalentConnType("mismat"));
    try testing.expect(!isCovalentConnType("saltbr"));
    try testing.expect(!isCovalentConnType(""));
}

test "parseStructConn skips non-covalent connections" {
    const mmcif = @import("../mmcif.zig");
    const source = @embedFile("../test_data/disulfide_with_hbond.cif");
    var doc = cif.readString(testing.allocator, source) catch unreachable;
    defer doc.deinit();
    const block = &doc.blocks.items[0];

    var mdl = try mmcif.parseModel(testing.allocator, source);
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
