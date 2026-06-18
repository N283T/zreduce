//! CCD (Chemical Component Dictionary) streaming parser.
//! Builds a ComponentDict from CIF source WITHOUT constructing a full Document.

const std = @import("std");
const Allocator = std.mem.Allocator;

const cif = @import("cif.zig");
const Tokenizer = cif.tokenizer.Tokenizer;
const Token = cif.tokenizer.Token;
const TokenType = cif.tokenizer.TokenType;

// ---------------------------------------------------------------------------
// Public types
// ---------------------------------------------------------------------------

pub const BondOrder = enum(u3) {
    single,
    double,
    triple,
    aromatic,
    delocalized,
    unknown,

    pub fn fromString(s: []const u8) BondOrder {
        if (std.ascii.eqlIgnoreCase(s, "SING")) return .single;
        if (std.ascii.eqlIgnoreCase(s, "DOUB")) return .double;
        if (std.ascii.eqlIgnoreCase(s, "TRIP")) return .triple;
        if (std.ascii.eqlIgnoreCase(s, "AROM")) return .aromatic;
        if (std.ascii.eqlIgnoreCase(s, "DELO")) return .delocalized;
        return .unknown;
    }
};

pub const CompAtom = struct {
    name: [4]u8 = .{ ' ', ' ', ' ', ' ' },
    name_len: u4 = 0,
    element_symbol: [2]u8 = .{ ' ', ' ' },
    charge: i8 = 0,
    leaving: bool = false,
    aromatic: bool = false,
    ideal_x: f32 = 0.0,
    ideal_y: f32 = 0.0,
    ideal_z: f32 = 0.0,

    pub fn nameSlice(self: *const CompAtom) []const u8 {
        return self.name[0..@min(@as(usize, self.name_len), 4)];
    }
};

pub const CompBond = struct {
    atom_idx_1: u16,
    atom_idx_2: u16,
    order: BondOrder,
    aromatic: bool = false,
};

pub const Component = struct {
    comp_id: []const u8,
    comp_type: []const u8,
    atoms: []CompAtom,
    bonds: []CompBond,
};

pub const ComponentDict = struct {
    components: std.StringHashMap(Component),
    allocator: Allocator,

    pub fn deinit(self: *ComponentDict) void {
        var iter = self.components.iterator();
        while (iter.next()) |entry| {
            self.allocator.free(entry.key_ptr.*);
            // comp_id shares memory with key — already freed above.
            self.allocator.free(entry.value_ptr.comp_type);
            self.allocator.free(entry.value_ptr.atoms);
            self.allocator.free(entry.value_ptr.bonds);
        }
        self.components.deinit();
    }

    pub fn get(self: *const ComponentDict, comp_id: []const u8) ?Component {
        return self.components.get(comp_id);
    }
};

// ---------------------------------------------------------------------------
// Internal builder state
// ---------------------------------------------------------------------------

const LoopKind = enum { atom, bond, other };

const AtomColIdx = struct {
    atom_id: ?usize = null,
    type_symbol: ?usize = null,
    charge: ?usize = null,
    leaving: ?usize = null,
    x_ideal: ?usize = null,
    y_ideal: ?usize = null,
    z_ideal: ?usize = null,
    aromatic: ?usize = null,
};

const BondColIdx = struct {
    atom_id_1: ?usize = null,
    atom_id_2: ?usize = null,
    value_order: ?usize = null,
    aromatic: ?usize = null,
};

const Builder = struct {
    allocator: Allocator,
    comp_id: ?[]const u8 = null, // slice into source
    comp_type: []const u8 = "", // slice into source
    atoms: std.ArrayListUnmanaged(CompAtom) = .empty,
    bonds: std.ArrayListUnmanaged(CompBond) = .empty,

    fn deinit(self: *Builder) void {
        self.atoms.deinit(self.allocator);
        self.bonds.deinit(self.allocator);
    }

    fn reset(self: *Builder) void {
        self.comp_id = null;
        self.comp_type = "";
        self.atoms.clearRetainingCapacity();
        self.bonds.clearRetainingCapacity();
    }

    /// Find atom index by name. Returns null if not found.
    fn findAtom(self: *const Builder, name: []const u8) ?u16 {
        for (self.atoms.items, 0..) |*a, i| {
            if (std.mem.eql(u8, a.nameSlice(), name)) {
                return @intCast(i);
            }
        }
        return null;
    }

    /// Flush current component into dict. Only stores if comp_id is set and atoms exist.
    fn flush(self: *Builder, dict: *ComponentDict) !void {
        const id = self.comp_id orelse return;
        if (self.atoms.items.len == 0) return;

        // Single allocation shared as both HashMap key and Component.comp_id.
        const key = try dict.allocator.dupe(u8, id);
        errdefer dict.allocator.free(key);

        const owned_type = try dict.allocator.dupe(u8, self.comp_type);
        errdefer dict.allocator.free(owned_type);

        const owned_atoms = try dict.allocator.dupe(CompAtom, self.atoms.items);
        errdefer dict.allocator.free(owned_atoms);

        const owned_bonds = try dict.allocator.dupe(CompBond, self.bonds.items);
        errdefer dict.allocator.free(owned_bonds);

        const comp = Component{
            .comp_id = key,
            .comp_type = owned_type,
            .atoms = owned_atoms,
            .bonds = owned_bonds,
        };

        // If a duplicate exists, free the old entry first.
        if (dict.components.fetchRemove(key)) |old| {
            dict.allocator.free(old.key);
            // comp_id shares memory with key — already freed above.
            dict.allocator.free(old.value.comp_type);
            dict.allocator.free(old.value.atoms);
            dict.allocator.free(old.value.bonds);
        }

        try dict.components.put(key, comp);
    }
};

// ---------------------------------------------------------------------------
// Tag matching helpers
// ---------------------------------------------------------------------------

fn tagEq(tag: []const u8, expected: []const u8) bool {
    return std.ascii.eqlIgnoreCase(tag, expected);
}

fn setElement(atom: *CompAtom, s: []const u8) void {
    const len = @min(s.len, 2);
    atom.element_symbol = .{ ' ', ' ' };
    @memcpy(atom.element_symbol[0..len], s[0..len]);
}

// ---------------------------------------------------------------------------
// Main parser
// ---------------------------------------------------------------------------

pub fn parseComponentDict(allocator: Allocator, source: []const u8) !ComponentDict {
    var dict = ComponentDict{
        .components = std.StringHashMap(Component).init(allocator),
        .allocator = allocator,
    };
    errdefer dict.deinit();

    var builder = Builder{ .allocator = allocator };
    defer builder.deinit();

    var tok = try Tokenizer.init(source);
    var pending: ?Token = null;

    while (true) {
        const token = if (pending) |p| blk: {
            pending = null;
            break :blk p;
        } else tok.next();

        switch (token.type) {
            .eof => {
                try builder.flush(&dict);
                break;
            },

            .data => {
                try builder.flush(&dict);
                builder.reset();
                const text = token.text(source);
                builder.comp_id = if (text.len > 5) text[5..] else "";
            },

            .tag => {
                const tag = token.text(source);
                if (tagEq(tag, "_chem_comp.type")) {
                    const val = tok.next();
                    if (val.type == .value) {
                        builder.comp_type = val.text(source);
                    }
                } else {
                    // Consume the value following any other lone tag
                    _ = tok.next();
                }
            },

            .loop => {
                pending = try parseLoop(&tok, source, &builder);
            },

            .save_begin, .save_end => {},
            .invalid => return error.InvalidCifToken,
            .value => return error.UnexpectedValue,
        }
    }

    return dict;
}

/// Parse a loop_ block. Returns the first non-value token after the loop (may be eof).
fn parseLoop(tok: *Tokenizer, source: []const u8, builder: *Builder) !Token {
    // --- Phase 1: collect tags ---
    var tags = std.ArrayListUnmanaged([]const u8).empty;
    defer tags.deinit(builder.allocator);

    var first_non_tag: Token = undefined;
    while (true) {
        const t = tok.next();
        if (t.type == .tag) {
            try tags.append(builder.allocator, t.text(source));
        } else {
            first_non_tag = t;
            break;
        }
    }

    if (tags.items.len == 0) return first_non_tag;

    // --- Phase 2: identify loop kind ---
    const kind = detectLoopKind(tags.items);

    // Build column index maps
    var atom_cols = AtomColIdx{};
    var bond_cols = BondColIdx{};

    if (kind == .atom) {
        for (tags.items, 0..) |tag, i| {
            if (tagEq(tag, "_chem_comp_atom.atom_id")) atom_cols.atom_id = i;
            if (tagEq(tag, "_chem_comp_atom.type_symbol")) atom_cols.type_symbol = i;
            if (tagEq(tag, "_chem_comp_atom.charge")) atom_cols.charge = i;
            if (tagEq(tag, "_chem_comp_atom.pdbx_leaving_atom_flag")) atom_cols.leaving = i;
            if (tagEq(tag, "_chem_comp_atom.pdbx_aromatic_flag")) atom_cols.aromatic = i;
            if (tagEq(tag, "_chem_comp_atom.pdbx_model_cartn_x_ideal")) atom_cols.x_ideal = i;
            if (tagEq(tag, "_chem_comp_atom.pdbx_model_cartn_y_ideal")) atom_cols.y_ideal = i;
            if (tagEq(tag, "_chem_comp_atom.pdbx_model_cartn_z_ideal")) atom_cols.z_ideal = i;
        }
    } else if (kind == .bond) {
        for (tags.items, 0..) |tag, i| {
            if (tagEq(tag, "_chem_comp_bond.atom_id_1")) bond_cols.atom_id_1 = i;
            if (tagEq(tag, "_chem_comp_bond.atom_id_2")) bond_cols.atom_id_2 = i;
            if (tagEq(tag, "_chem_comp_bond.value_order")) bond_cols.value_order = i;
            if (tagEq(tag, "_chem_comp_bond.pdbx_aromatic_flag")) bond_cols.aromatic = i;
        }
    }

    const ncols = tags.items.len;

    // --- Phase 3: consume values row by row ---
    // Collect one row at a time (ncols values)
    var row = std.ArrayListUnmanaged([]const u8).empty;
    defer row.deinit(builder.allocator);

    var next_tok = first_non_tag;
    while (next_tok.type == .value) {
        // Build a row
        row.clearRetainingCapacity();
        try row.append(builder.allocator, next_tok.text(source));
        var col: usize = 1;
        while (col < ncols) : (col += 1) {
            const t = tok.next();
            if (t.type != .value) {
                // Incomplete row — treat as loop end
                return t;
            }
            try row.append(builder.allocator, t.text(source));
        }

        // Process row
        switch (kind) {
            .atom => try processAtomRow(builder, row.items, atom_cols),
            .bond => try processBondRow(builder, row.items, bond_cols),
            .other => {},
        }

        next_tok = tok.next();
    }

    return next_tok;
}

fn detectLoopKind(tags: []const []const u8) LoopKind {
    for (tags) |tag| {
        if (tag.len >= 16 and std.ascii.startsWithIgnoreCase(tag, "_chem_comp_atom.")) return .atom;
        if (tag.len >= 16 and std.ascii.startsWithIgnoreCase(tag, "_chem_comp_bond.")) return .bond;
    }
    return .other;
}

fn processAtomRow(builder: *Builder, row: []const []const u8, cols: AtomColIdx) !void {
    var atom = CompAtom{};

    if (cols.atom_id) |c| {
        if (c < row.len) {
            const s = row[c];
            const len = @min(s.len, 4);
            atom.name_len = @intCast(len);
            @memcpy(atom.name[0..len], s[0..len]);
        }
    }

    if (cols.type_symbol) |c| {
        if (c < row.len) {
            setElement(&atom, row[c]);
        }
    }

    if (cols.charge) |c| {
        if (c < row.len) {
            atom.charge = cif.value.asIntOr(i8, row[c], 0);
        }
    }

    if (cols.leaving) |c| {
        if (c < row.len) {
            atom.leaving = std.ascii.eqlIgnoreCase(row[c], "Y");
        }
    }

    if (cols.aromatic) |c| {
        if (c < row.len) {
            atom.aromatic = std.ascii.eqlIgnoreCase(row[c], "Y");
        }
    }

    if (cols.x_ideal) |c| {
        if (c < row.len) atom.ideal_x = cif.asFloatOr(row[c], 0.0);
    }
    if (cols.y_ideal) |c| {
        if (c < row.len) atom.ideal_y = cif.asFloatOr(row[c], 0.0);
    }
    if (cols.z_ideal) |c| {
        if (c < row.len) atom.ideal_z = cif.asFloatOr(row[c], 0.0);
    }

    try builder.atoms.append(builder.allocator, atom);
}

fn processBondRow(builder: *Builder, row: []const []const u8, cols: BondColIdx) !void {
    const idx1 = if (cols.atom_id_1) |c|
        if (c < row.len) builder.findAtom(row[c]) else null
    else
        null;

    const idx2 = if (cols.atom_id_2) |c|
        if (c < row.len) builder.findAtom(row[c]) else null
    else
        null;

    // Skip bonds where we can't resolve atom indices
    if (idx1 == null or idx2 == null) return;

    const order: BondOrder = if (cols.value_order) |c|
        if (c < row.len) BondOrder.fromString(row[c]) else .unknown
    else
        .unknown;

    const aromatic: bool = if (cols.aromatic) |c|
        if (c < row.len) std.ascii.eqlIgnoreCase(row[c], "Y") else false
    else
        false;

    const bond = CompBond{
        .atom_idx_1 = idx1.?,
        .atom_idx_2 = idx2.?,
        .order = order,
        .aromatic = aromatic,
    };

    try builder.bonds.append(builder.allocator, bond);
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

const testing = std.testing;

test "parse CCD fragment" {
    const source =
        \\data_ALA
        \\_chem_comp.type 'L-peptide linking'
        \\loop_
        \\_chem_comp_atom.atom_id
        \\_chem_comp_atom.type_symbol
        \\_chem_comp_atom.charge
        \\_chem_comp_atom.pdbx_leaving_atom_flag
        \\N   N 0 N
        \\CA  C 0 N
        \\C   C 0 N
        \\O   O 0 N
        \\CB  C 0 N
        \\H   H 0 N
        \\HA  H 0 N
        \\HB1 H 0 N
        \\HB2 H 0 N
        \\HB3 H 0 N
        \\OXT O 0 Y
    ;
    var dict = try parseComponentDict(std.testing.allocator, source);
    defer dict.deinit();

    const ala = dict.get("ALA");
    try testing.expect(ala != null);
    try testing.expectEqual(@as(usize, 11), ala.?.atoms.len);
}

test "BondOrder.fromString" {
    try testing.expectEqual(BondOrder.single, BondOrder.fromString("SING"));
    try testing.expectEqual(BondOrder.double, BondOrder.fromString("DOUB"));
    try testing.expectEqual(BondOrder.triple, BondOrder.fromString("TRIP"));
    try testing.expectEqual(BondOrder.aromatic, BondOrder.fromString("AROM"));
    try testing.expectEqual(BondOrder.delocalized, BondOrder.fromString("DELO"));
    try testing.expectEqual(BondOrder.unknown, BondOrder.fromString("???"));
}

test "parse CCD with bonds" {
    const source =
        \\data_GLY
        \\_chem_comp.type 'L-peptide linking'
        \\loop_
        \\_chem_comp_atom.atom_id
        \\_chem_comp_atom.type_symbol
        \\_chem_comp_atom.charge
        \\_chem_comp_atom.pdbx_leaving_atom_flag
        \\N N 0 N
        \\CA C 0 N
        \\C C 0 N
        \\loop_
        \\_chem_comp_bond.atom_id_1
        \\_chem_comp_bond.atom_id_2
        \\_chem_comp_bond.value_order
        \\_chem_comp_bond.pdbx_aromatic_flag
        \\N CA SING N
        \\CA C SING N
    ;
    var dict = try parseComponentDict(std.testing.allocator, source);
    defer dict.deinit();

    const gly = dict.get("GLY");
    try testing.expect(gly != null);
    try testing.expectEqual(@as(usize, 3), gly.?.atoms.len);
    try testing.expectEqual(@as(usize, 2), gly.?.bonds.len);
    try testing.expectEqual(BondOrder.single, gly.?.bonds[0].order);
}

test "comp_type stored correctly" {
    const source =
        \\data_ATP
        \\_chem_comp.type 'non-polymer'
        \\loop_
        \\_chem_comp_atom.atom_id
        \\_chem_comp_atom.type_symbol
        \\_chem_comp_atom.charge
        \\_chem_comp_atom.pdbx_leaving_atom_flag
        \\P P 0 N
    ;
    var dict = try parseComponentDict(std.testing.allocator, source);
    defer dict.deinit();

    const atp = dict.get("ATP");
    try testing.expect(atp != null);
    try testing.expectEqualStrings("non-polymer", atp.?.comp_type);
}
