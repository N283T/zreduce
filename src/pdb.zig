//! PDB format parser: parses ATOM/HETATM records into a Model.
//!
//! Only MODEL 1 is parsed when MODEL/ENDMDL records are present.
//! CONECT and MASTER records are dropped entirely.
//! All other non-coordinate records are kept as raw_line passthrough.

const std = @import("std");
const Allocator = std.mem.Allocator;

const model_mod = @import("model.zig");
const element = @import("element.zig");
const math = @import("math.zig");

const Model = model_mod.Model;
const Atom = model_mod.Atom;
const Residue = model_mod.Residue;
const Chain = model_mod.Chain;

pub const PdbError = error{
    OutOfMemory,
    InvalidCoordinateValue,
};

/// A parsed PDB record: either a coordinate atom or a raw passthrough line.
pub const PdbRecord = union(enum) {
    atom_site,
    raw_line: []const u8,
};

/// Full parse result including the model and passthrough records.
pub const PdbParseResult = struct {
    model: Model,
    records: std.ArrayListUnmanaged(PdbRecord),
    source: []const u8,

    pub fn deinit(self: *PdbParseResult, allocator: Allocator) void {
        self.model.deinit();
        self.records.deinit(allocator);
    }
};

// ---------------------------------------------------------------------------
// Standard polymer residue names (20 AA + special + nucleotides)
// ---------------------------------------------------------------------------

const polymer_comp_ids = std.StaticStringMap(void).initComptime(.{
    // 20 standard amino acids
    .{ "ALA", {} }, .{ "ARG", {} }, .{ "ASN", {} }, .{ "ASP", {} },
    .{ "CYS", {} }, .{ "GLN", {} }, .{ "GLU", {} }, .{ "GLY", {} },
    .{ "HIS", {} }, .{ "ILE", {} }, .{ "LEU", {} }, .{ "LYS", {} },
    .{ "MET", {} }, .{ "PHE", {} }, .{ "PRO", {} }, .{ "SER", {} },
    .{ "THR", {} }, .{ "TRP", {} }, .{ "TYR", {} }, .{ "VAL", {} },
    // Special amino acids
    .{ "MSE", {} }, .{ "SEC", {} }, .{ "PYL", {} },
    // Standard nucleotides (RNA)
    .{ "A",   {} }, .{ "C",   {} }, .{ "G",   {} }, .{ "U",   {} },
    // Standard nucleotides (DNA)
    .{ "DA",  {} }, .{ "DC",  {} }, .{ "DG",  {} }, .{ "DT",  {} },
    .{ "DI",  {} },
});

fn entityTypeFromCompId(comp_id: []const u8) model_mod.residue.EntityType {
    if (std.mem.eql(u8, comp_id, "HOH") or
        std.mem.eql(u8, comp_id, "WAT") or
        std.mem.eql(u8, comp_id, "H2O"))
    {
        return .water;
    }
    if (polymer_comp_ids.has(comp_id)) return .polymer;
    return .non_polymer;
}

// ---------------------------------------------------------------------------
// Line parsing helpers
// ---------------------------------------------------------------------------

/// Return line[start..end] if the line is long enough, else return "".
fn safeSlice(line: []const u8, start: usize, end: usize) []const u8 {
    if (start >= line.len) return "";
    const e = @min(end, line.len);
    return line[start..e];
}

/// Parse a float from a slice, trimming whitespace. Returns 0 on failure.
fn parseFloat(s: []const u8) f32 {
    const trimmed = std.mem.trim(u8, s, " ");
    if (trimmed.len == 0) return 0.0;
    return std.fmt.parseFloat(f32, trimmed) catch 0.0;
}

/// Parse an integer from a slice, trimming whitespace. Returns 0 on failure.
fn parseInt(comptime T: type, s: []const u8) T {
    const trimmed = std.mem.trim(u8, s, " ");
    if (trimmed.len == 0) return 0;
    return std.fmt.parseInt(T, trimmed, 10) catch 0;
}

/// Infer element from PDB atom name when columns 76-77 are absent.
/// PDB atom names are 4 chars (cols 12-15). The element is typically the
/// first non-digit, non-space character when the name starts in col 13
/// (1-char element), or the first two chars for 4-char names starting at
/// col 12 (e.g. " CA " -> C, "FE1 " -> Fe).
fn inferElementFromName(raw_name: []const u8) element.AtomType {
    // raw_name is the 4-char field from cols 12-15 (may have leading/trailing spaces)
    const trimmed = std.mem.trim(u8, raw_name, " ");
    if (trimmed.len == 0) return .unknown;

    // Strip leading digits (e.g. "1H", "2HB" in hydrogen names)
    var start: usize = 0;
    while (start < trimmed.len and std.ascii.isDigit(trimmed[start])) : (start += 1) {}
    if (start >= trimmed.len) return .unknown;

    // The element symbol is at most 2 chars starting at `start`
    const sym_end = @min(start + 2, trimmed.len);
    return element.elementFromSymbol(trimmed[start..sym_end]);
}

// ---------------------------------------------------------------------------
// Core parser state
// ---------------------------------------------------------------------------

const ParseState = struct {
    mdl: Model,
    allocator: Allocator,

    // Current chain tracking
    cur_chain_id: u8 = 0,
    cur_chain_idx: u32 = 0,
    in_chain: bool = false,
    next_entity_id: u32 = 1,

    // Current residue tracking
    cur_seq_id: i32 = 0,
    cur_ins_code: u8 = 0,
    cur_comp_id: [5]u8 = undefined,
    cur_comp_id_len: usize = 0,
    in_residue: bool = false,
    res_atom_start: u32 = 0,

    fn init(allocator: Allocator) ParseState {
        return .{
            .mdl = Model.init(allocator),
            .allocator = allocator,
        };
    }

    fn deinit(self: *ParseState) void {
        self.mdl.deinit();
    }

    fn curCompIdSlice(self: *const ParseState) []const u8 {
        return self.cur_comp_id[0..self.cur_comp_id_len];
    }

    fn setCompId(self: *ParseState, s: []const u8) void {
        const n = @min(s.len, 5);
        self.cur_comp_id_len = n;
        for (0..n) |i| self.cur_comp_id[i] = s[i];
    }

    /// Flush the in-progress residue (set atom_end) and append it.
    fn flushResidue(self: *ParseState) !void {
        if (!self.in_residue) return;
        const atom_end: u32 = @intCast(self.mdl.atoms.items.len);
        var res = Residue{
            .chain_idx = self.cur_chain_idx,
            .seq_id = self.cur_seq_id,
            .auth_seq_id = self.cur_seq_id,
            .ins_code = self.cur_ins_code,
            .atom_start = self.res_atom_start,
            .atom_end = atom_end,
            .entity_type = entityTypeFromCompId(self.curCompIdSlice()),
        };
        res.setCompId(self.curCompIdSlice());
        try self.mdl.residues.append(self.allocator, res);
        self.in_residue = false;
    }

    /// Flush the in-progress chain (set residue_end) and append it.
    fn flushChain(self: *ParseState) !void {
        if (!self.in_chain) return;
        try self.flushResidue();
        const res_end: u32 = @intCast(self.mdl.residues.items.len);
        self.mdl.chains.items[self.cur_chain_idx].residue_end = res_end;
        self.in_chain = false;
    }

    /// Open a new chain for `chain_id`. Does not flush the old chain; caller
    /// must call flushChain() first when transitioning between chains.
    fn openChain(self: *ParseState, chain_id: u8) !void {
        const chain_str = [_]u8{chain_id};
        var entity_buf: [4]u8 = undefined;
        const entity_str = std.fmt.bufPrint(&entity_buf, "{d}", .{self.next_entity_id}) catch "1";
        self.next_entity_id += 1;

        var ch = Chain{
            .residue_start = @intCast(self.mdl.residues.items.len),
            .residue_end = @intCast(self.mdl.residues.items.len),
        };
        ch.setLabelAsymId(&chain_str);
        ch.setAuthAsymId(&chain_str);
        ch.setEntityId(entity_str);

        self.cur_chain_idx = @intCast(self.mdl.chains.items.len);
        self.cur_chain_id = chain_id;
        self.in_chain = true;
        try self.mdl.chains.append(self.allocator, ch);
    }

    /// Process a single ATOM/HETATM line.
    fn processAtomLine(self: *ParseState, line: []const u8) !void {
        // --- Extract fixed-width fields ---
        const raw_name = safeSlice(line, 12, 16); // 4 chars, space-padded
        const altloc_ch = if (line.len > 16) line[16] else ' ';
        const comp_id_raw = std.mem.trim(u8, safeSlice(line, 17, 20), " ");
        const chain_id_ch = if (line.len > 21) line[21] else ' ';
        const seq_id_str = safeSlice(line, 22, 26);
        const ins_code_ch: u8 = if (line.len > 26 and line[26] != ' ') line[26] else ' ';
        const x_str = safeSlice(line, 30, 38);
        const y_str = safeSlice(line, 38, 46);
        const z_str = safeSlice(line, 46, 54);
        const occ_str = safeSlice(line, 54, 60);
        const bfac_str = safeSlice(line, 60, 66);
        const elem_raw = safeSlice(line, 76, 78);

        const seq_id = parseInt(i32, seq_id_str);
        const x = parseFloat(x_str);
        const y = parseFloat(y_str);
        const z = parseFloat(z_str);

        // Validate coordinates by checking the strings are non-empty
        if (std.mem.trim(u8, x_str, " ").len == 0) return PdbError.InvalidCoordinateValue;
        if (std.mem.trim(u8, y_str, " ").len == 0) return PdbError.InvalidCoordinateValue;
        if (std.mem.trim(u8, z_str, " ").len == 0) return PdbError.InvalidCoordinateValue;

        // --- Chain transition ---
        if (!self.in_chain or chain_id_ch != self.cur_chain_id) {
            try self.flushChain();
            try self.openChain(chain_id_ch);
            self.in_residue = false;
        }

        // --- Residue transition ---
        const same_residue = self.in_residue and
            self.cur_seq_id == seq_id and
            self.cur_ins_code == ins_code_ch and
            std.mem.eql(u8, self.curCompIdSlice(), comp_id_raw);

        if (!same_residue) {
            try self.flushResidue();
            self.cur_seq_id = seq_id;
            self.cur_ins_code = ins_code_ch;
            self.setCompId(comp_id_raw);
            self.res_atom_start = @intCast(self.mdl.atoms.items.len);
            self.in_residue = true;
        }

        // --- Determine element ---
        const elem_trimmed = std.mem.trim(u8, elem_raw, " ");
        const atom_type = if (elem_trimmed.len > 0)
            element.elementFromSymbol(elem_trimmed)
        else
            inferElementFromName(raw_name);

        const is_h = switch (atom_type) {
            .H, .Har, .Hpol, .Ha_p, .HOd => true,
            else => false,
        };
        const vdw = atom_type.info().explicit_radius;

        // --- Trim atom name ---
        const name_trimmed = std.mem.trim(u8, raw_name, " ");

        var atm = Atom{
            .pos = .{ .x = x, .y = y, .z = z },
            .element_type = atom_type,
            .residue_idx = @intCast(self.mdl.residues.items.len),
            .altloc = altloc_ch,
            .occupancy = parseFloat(occ_str),
            .b_factor = parseFloat(bfac_str),
            .is_hydrogen = is_h,
            .vdw_radius = vdw,
        };
        atm.setName(name_trimmed);

        try self.mdl.atoms.append(self.allocator, atm);
    }

    /// Finalise after all lines are processed.
    fn finalize(self: *ParseState) !void {
        try self.flushChain();
        self.mdl.original_atom_count = @intCast(self.mdl.atoms.items.len);
        // Apply chain-break detection based on seq_id gaps
        try detectChainBreaks(&self.mdl);
    }
};

/// Detect chain breaks: for consecutive polymer residues in the same chain,
/// mark is_chain_break_before when seq_id gap > 1.
fn detectChainBreaks(mdl: *Model) !void {
    if (mdl.residues.items.len < 2) return;
    for (mdl.chains.items) |chain| {
        const res_slice = mdl.residues.items[chain.residue_start..chain.residue_end];
        if (res_slice.len < 2) continue;
        var prev_seq: i32 = res_slice[0].seq_id;
        var prev_is_polymer = res_slice[0].entity_type == .polymer;
        for (res_slice[1..]) |*res| {
            const is_poly = res.entity_type == .polymer;
            if (is_poly and prev_is_polymer and (res.seq_id - prev_seq) > 1) {
                res.is_chain_break_before = true;
            }
            prev_seq = res.seq_id;
            prev_is_polymer = is_poly;
        }
    }
}

// ---------------------------------------------------------------------------
// Public API
// ---------------------------------------------------------------------------

/// Strip DOS line ending from a line slice (remove trailing \r).
fn stripCR(line: []const u8) []const u8 {
    if (line.len > 0 and line[line.len - 1] == '\r') return line[0 .. line.len - 1];
    return line;
}

/// Parse PDB source, keeping passthrough records. Returns a PdbParseResult.
pub fn parse(allocator: Allocator, source: []const u8) PdbError!PdbParseResult {
    var state = ParseState.init(allocator);
    errdefer state.deinit();

    var records = std.ArrayListUnmanaged(PdbRecord){};
    errdefer records.deinit(allocator);

    var in_model_block = false;
    var model_done = false;

    var iter = std.mem.splitScalar(u8, source, '\n');
    while (iter.next()) |raw_line| {
        const line = stripCR(raw_line);
        if (line.len == 0) continue;

        // Record type is cols 0-5 (up to 6 chars), left-justified
        const rec_type = safeSlice(line, 0, 6);

        // Handle MODEL/ENDMDL — only parse MODEL 1
        if (std.mem.startsWith(u8, rec_type, "MODEL ") or std.mem.eql(u8, std.mem.trim(u8, rec_type, " "), "MODEL")) {
            if (!model_done and !in_model_block) {
                in_model_block = true;
            } else {
                // Skip subsequent models; still passthrough the line
                try records.append(allocator, .{ .raw_line = line });
            }
            continue;
        }
        if (std.mem.startsWith(u8, rec_type, "ENDMDL")) {
            if (in_model_block) {
                in_model_block = false;
                model_done = true;
            }
            try records.append(allocator, .{ .raw_line = line });
            continue;
        }

        // If we've finished MODEL 1 and there are more atoms, skip them
        if (model_done and
            (std.mem.startsWith(u8, rec_type, "ATOM  ") or
             std.mem.startsWith(u8, rec_type, "HETATM")))
        {
            continue;
        }

        // Drop CONECT and MASTER records entirely
        if (std.mem.startsWith(u8, rec_type, "CONECT") or
            std.mem.startsWith(u8, rec_type, "MASTER"))
        {
            continue;
        }

        // TER record: flush current chain
        if (std.mem.startsWith(u8, rec_type, "TER")) {
            try state.flushChain();
            try records.append(allocator, .{ .raw_line = line });
            continue;
        }

        // Coordinate records
        if (std.mem.startsWith(u8, rec_type, "ATOM  ") or
            std.mem.startsWith(u8, rec_type, "HETATM"))
        {
            try state.processAtomLine(line);
            try records.append(allocator, .atom_site);
            continue;
        }

        // Everything else: passthrough
        try records.append(allocator, .{ .raw_line = line });
    }

    try state.finalize();

    return PdbParseResult{
        .model = state.mdl,
        .records = records,
        .source = source,
    };
}

/// Parse PDB source into a Model, discarding passthrough records.
pub fn parseModel(allocator: Allocator, source: []const u8) PdbError!Model {
    var result = try parse(allocator, source);
    result.records.deinit(allocator);
    return result.model;
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

const testing = std.testing;

test "parse tiny PDB" {
    const source = @embedFile("test_data/tiny.pdb");
    var mdl = try parseModel(testing.allocator, source);
    defer mdl.deinit();
    try testing.expectEqual(@as(usize, 5), mdl.atoms.items.len);
    try testing.expectEqual(@as(usize, 1), mdl.residues.items.len);
    try testing.expectEqual(@as(usize, 1), mdl.chains.items.len);
    try testing.expectApproxEqAbs(mdl.atoms.items[0].pos.x, 1.0, 1e-3);
    try testing.expectApproxEqAbs(mdl.atoms.items[0].pos.y, 2.0, 1e-3);
    try testing.expectEqualStrings("ALA", mdl.residues.items[0].compIdSlice());
    try testing.expectEqualStrings("N", mdl.atoms.items[0].nameSlice());
}

test "parse multi-chain PDB" {
    const source = @embedFile("test_data/multi_chain.pdb");
    var mdl = try parseModel(testing.allocator, source);
    defer mdl.deinit();
    try testing.expectEqual(@as(usize, 12), mdl.atoms.items.len);
    try testing.expectEqual(@as(usize, 3), mdl.residues.items.len);
    try testing.expectEqual(@as(usize, 2), mdl.chains.items.len);
    try testing.expectEqualStrings("A", mdl.chains.items[0].labelSlice());
    try testing.expectEqual(@as(u32, 0), mdl.chains.items[0].residue_start);
    try testing.expectEqual(@as(u32, 2), mdl.chains.items[0].residue_end);
    try testing.expectEqualStrings("B", mdl.chains.items[1].labelSlice());
    try testing.expectEqual(@as(u32, 2), mdl.chains.items[1].residue_start);
    try testing.expectEqual(@as(u32, 3), mdl.chains.items[1].residue_end);
    try testing.expectApproxEqAbs(mdl.atoms.items[11].pos.x, 13.0, 1e-3);
}

test "parse HETATM and entity types" {
    const source = @embedFile("test_data/hetatm.pdb");
    var mdl = try parseModel(testing.allocator, source);
    defer mdl.deinit();
    try testing.expectEqual(@as(usize, 9), mdl.atoms.items.len);
    try testing.expectEqual(@as(usize, 3), mdl.residues.items.len);
    // ALA -> polymer
    try testing.expectEqual(model_mod.residue.EntityType.polymer, mdl.residues.items[0].entity_type);
    // EDO -> non_polymer
    try testing.expectEqual(model_mod.residue.EntityType.non_polymer, mdl.residues.items[1].entity_type);
    // HOH -> water
    try testing.expectEqual(model_mod.residue.EntityType.water, mdl.residues.items[2].entity_type);
}
