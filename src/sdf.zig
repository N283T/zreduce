//! SDF/MOL V2000 format parser.
//! Produces ccd.ComponentDict structs compatible with the CCD dictionary system.
//!
//! Supports:
//!   - Single and multi-molecule SDF files (molecules separated by $$$$)
//!   - MOL V2000 atom and bond blocks
//!   - M  CHG property lines for formal charges
//!   - Malformed molecules are skipped with a warning log

const std = @import("std");
const Allocator = std.mem.Allocator;

const ccd = @import("ccd.zig");
const CompAtom = ccd.CompAtom;
const CompBond = ccd.CompBond;
const BondOrder = ccd.BondOrder;
const Component = ccd.Component;
const ComponentDict = ccd.ComponentDict;

// ---------------------------------------------------------------------------
// Public API
// ---------------------------------------------------------------------------

/// Parse SDF file content into a ComponentDict.
/// Multiple molecules separated by $$$$ are supported.
/// The molecule name (line 1) is used as comp_id key.
/// Molecules with empty/whitespace names are skipped.
/// Malformed molecules are skipped with a warning.
pub fn parseSdf(allocator: Allocator, source: []const u8) !ComponentDict {
    var dict = ComponentDict{
        .components = std.StringHashMap(Component).init(allocator),
        .allocator = allocator,
    };
    errdefer dict.deinit();

    var remaining = source;
    while (remaining.len > 0) {
        const mol_end = std.mem.indexOf(u8, remaining, "$$$$") orelse remaining.len;
        const mol_block = remaining[0..mol_end];

        // Advance past the $$$$ separator (and trailing newline if present)
        if (mol_end < remaining.len) {
            var skip = mol_end + 4;
            if (skip < remaining.len and remaining[skip] == '\r') skip += 1;
            if (skip < remaining.len and remaining[skip] == '\n') skip += 1;
            remaining = remaining[skip..];
        } else {
            remaining = remaining[remaining.len..];
        }

        parseMolBlock(allocator, mol_block, &dict) catch |err| {
            std.log.warn("sdf: skipping malformed molecule block: {s}", .{@errorName(err)});
        };
    }

    return dict;
}

// ---------------------------------------------------------------------------
// Internal: parse one MOL V2000 block
// ---------------------------------------------------------------------------

const ParseError = error{
    EmptyName,
    MissingCountsLine,
    InvalidCountsLine,
    TooManyAtoms,
    TooManyBonds,
    TruncatedAtomBlock,
    TruncatedBondBlock,
};

fn parseMolBlock(allocator: Allocator, block: []const u8, dict: *ComponentDict) !void {
    var lines = std.mem.splitScalar(u8, block, '\n');

    // --- Line 1: molecule name (comp_id) ---
    const name_raw = lines.next() orelse return error.EmptyName;
    const name = std.mem.trim(u8, stripCR(name_raw), " \t");
    if (name.len == 0) return error.EmptyName;

    // --- Line 2: program/timestamp (ignored) ---
    _ = lines.next();
    // --- Line 3: comment (ignored) ---
    _ = lines.next();

    // --- Counts line ---
    const counts_raw = lines.next() orelse return error.MissingCountsLine;
    const counts_line = stripCR(counts_raw);
    if (counts_line.len < 6) return error.InvalidCountsLine;

    const n_atoms = parseFixedInt(counts_line, 0, 3) catch return error.InvalidCountsLine;
    const n_bonds = parseFixedInt(counts_line, 3, 6) catch return error.InvalidCountsLine;

    if (n_atoms > 999) return error.TooManyAtoms;
    if (n_bonds > 999) return error.TooManyBonds;

    // --- Atom block ---
    var atoms = std.ArrayListUnmanaged(CompAtom){};
    defer atoms.deinit(allocator);

    // Track per-element counters for name generation.
    // Keys are short element strings stored as fixed [2]u8, not heap-allocated.
    var elem_counts = std.AutoHashMap([2]u8, u16).init(allocator);
    defer elem_counts.deinit();

    var atom_i: usize = 0;
    while (atom_i < n_atoms) : (atom_i += 1) {
        const raw = lines.next() orelse return error.TruncatedAtomBlock;
        const line = stripCR(raw);
        const atom = try parseAtomLine(line, &elem_counts);
        try atoms.append(allocator, atom);
    }

    // --- Bond block ---
    var bonds = std.ArrayListUnmanaged(CompBond){};
    defer bonds.deinit(allocator);

    var bond_i: usize = 0;
    while (bond_i < n_bonds) : (bond_i += 1) {
        const raw = lines.next() orelse return error.TruncatedBondBlock;
        const line = stripCR(raw);
        const bond = parseBondLine(line) catch continue; // skip malformed bonds
        if (bond.atom_idx_1 >= atoms.items.len or bond.atom_idx_2 >= atoms.items.len) continue;
        try bonds.append(allocator, bond);
    }

    // --- Properties block: M  CHG ---
    while (lines.next()) |raw| {
        const line = stripCR(raw);
        if (std.mem.startsWith(u8, line, "M  CHG")) {
            applyCharges(line, atoms.items);
        }
    }

    // --- Flush into dict ---
    const key = try dict.allocator.dupe(u8, name);
    errdefer dict.allocator.free(key);

    const owned_type = try dict.allocator.dupe(u8, "non-polymer");
    errdefer dict.allocator.free(owned_type);

    const owned_atoms = try dict.allocator.dupe(CompAtom, atoms.items);
    errdefer dict.allocator.free(owned_atoms);

    const owned_bonds = try dict.allocator.dupe(CompBond, bonds.items);
    errdefer dict.allocator.free(owned_bonds);

    // Replace duplicate if present
    if (dict.components.fetchRemove(key)) |old| {
        dict.allocator.free(old.key);
        dict.allocator.free(old.value.comp_type);
        dict.allocator.free(old.value.atoms);
        dict.allocator.free(old.value.bonds);
    }

    try dict.components.put(key, Component{
        .comp_id = key,
        .comp_type = owned_type,
        .atoms = owned_atoms,
        .bonds = owned_bonds,
    });
}

// ---------------------------------------------------------------------------
// Atom line parsing
// ---------------------------------------------------------------------------

/// Parse one atom line from the V2000 atom block.
/// Coordinates are in columns 0-9, 10-19, 20-29; symbol in 31-33 (0-based).
/// Generates a unique PDB-style 4-char name: e.g. "C1  " for first carbon.
fn parseAtomLine(
    line: []const u8,
    elem_counts: *std.AutoHashMap([2]u8, u16),
) !CompAtom {
    if (line.len < 34) return error.TruncatedAtomBlock;

    const x = parseF32Field(line[0..10]);
    const y = parseF32Field(line[10..20]);
    const z = parseF32Field(line[20..30]);

    // Symbol is at column 31-33 (1-char gap after z-field)
    const sym_raw = std.mem.trim(u8, line[31..@min(line.len, 34)], " \t");
    const sym_len = @min(sym_raw.len, 2);

    var atom = CompAtom{
        .ideal_x = x,
        .ideal_y = y,
        .ideal_z = z,
    };

    // Fill element_symbol (left-aligned, space-padded)
    atom.element_symbol = .{ ' ', ' ' };
    @memcpy(atom.element_symbol[0..sym_len], sym_raw[0..sym_len]);

    // Build a [2]u8 key using uppercase element letters
    var elem_key: [2]u8 = .{ ' ', ' ' };
    for (sym_raw[0..sym_len], 0..) |ch, i| {
        elem_key[i] = std.ascii.toUpper(ch);
    }

    const count_entry = try elem_counts.getOrPut(elem_key);
    if (!count_entry.found_existing) {
        count_entry.value_ptr.* = 0;
    }
    count_entry.value_ptr.* += 1;
    const idx = count_entry.value_ptr.*;

    // Build name string: e.g. "C1", "N12" — left-justified in 4 chars
    var name_buf: [8]u8 = undefined;
    // Render only the non-space part of elem_key as the element prefix
    const elem_str = std.mem.trimRight(u8, &elem_key, " ");
    const name_str = std.fmt.bufPrint(&name_buf, "{s}{d}", .{ elem_str, idx }) catch "X";

    const copy_len = @min(name_str.len, 4);
    atom.name = .{ ' ', ' ', ' ', ' ' };
    @memcpy(atom.name[0..copy_len], name_str[0..copy_len]);
    atom.name_len = @intCast(copy_len);

    return atom;
}

// ---------------------------------------------------------------------------
// Bond line parsing
// ---------------------------------------------------------------------------

fn parseBondLine(line: []const u8) !CompBond {
    if (line.len < 9) return error.TruncatedBondBlock;

    const a1_1based = parseFixedInt(line, 0, 3) catch return error.TruncatedBondBlock;
    const a2_1based = parseFixedInt(line, 3, 6) catch return error.TruncatedBondBlock;
    const bond_type = parseFixedInt(line, 6, 9) catch return error.TruncatedBondBlock;

    if (a1_1based == 0 or a2_1based == 0) return error.InvalidCountsLine;

    const order: BondOrder = switch (bond_type) {
        1 => .single,
        2 => .double,
        3 => .triple,
        4 => .aromatic,
        else => .unknown,
    };

    return CompBond{
        .atom_idx_1 = @intCast(a1_1based - 1),
        .atom_idx_2 = @intCast(a2_1based - 1),
        .order = order,
        .aromatic = order == .aromatic,
    };
}

// ---------------------------------------------------------------------------
// M  CHG property line
// ---------------------------------------------------------------------------

/// Apply formal charges from an `M  CHG  n  a1 c1  a2 c2 ...` line.
/// Atom indices in the line are 1-based.
fn applyCharges(line: []const u8, atoms: []CompAtom) void {
    // Format: "M  CHG  n  a1 c1  a2 c2 ..."
    // After "M  CHG" there is a count field then pairs.
    // We tokenize the remainder by whitespace for robustness.
    const rest = if (line.len > 6) line[6..] else return;

    var it = std.mem.tokenizeAny(u8, rest, " \t");

    // First token: count of pairs
    const count_str = it.next() orelse return;
    const count = std.fmt.parseInt(usize, count_str, 10) catch return;

    var pair: usize = 0;
    while (pair < count) : (pair += 1) {
        const atom_str = it.next() orelse return;
        const charge_str = it.next() orelse return;

        const atom_1based = std.fmt.parseInt(usize, atom_str, 10) catch continue;
        const charge = std.fmt.parseInt(i8, charge_str, 10) catch continue;

        if (atom_1based == 0 or atom_1based > atoms.len) continue;
        atoms[atom_1based - 1].charge = charge;
    }
}

// ---------------------------------------------------------------------------
// Low-level field helpers
// ---------------------------------------------------------------------------

/// Strip trailing carriage return.
fn stripCR(s: []const u8) []const u8 {
    if (s.len > 0 and s[s.len - 1] == '\r') return s[0 .. s.len - 1];
    return s;
}

/// Parse a fixed-width integer field from `line[start..end]` (right-justified, space-padded).
fn parseFixedInt(line: []const u8, start: usize, end: usize) !usize {
    if (end > line.len) return error.InvalidCountsLine;
    const field = std.mem.trim(u8, line[start..end], " \t");
    return std.fmt.parseInt(usize, field, 10) catch error.InvalidCountsLine;
}

/// Parse a floating-point field, returning 0.0 on failure.
fn parseF32Field(field: []const u8) f32 {
    const trimmed = std.mem.trim(u8, field, " \t");
    return std.fmt.parseFloat(f32, trimmed) catch 0.0;
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

const testing = std.testing;

// Methanol: C + O + H*4 (CH3OH)
//   Atoms: C, O, H, H, H, H
//   Bonds: C-O single, C-H x3, O-H x1
const methanol_sdf =
    \\methanol
    \\  test
    \\  comment
    \\  6  5  0  0  0  0  0  0  0  0999 V2000
    \\   -0.6622    0.5342    0.0000 C   0  0  0  0  0  0
    \\    0.7707    0.5342    0.0000 O   0  0  0  0  0  0
    \\   -1.0387   -0.4789    0.0000 H   0  0  0  0  0  0
    \\   -1.0387    1.0441    0.8776 H   0  0  0  0  0  0
    \\   -1.0387    1.0441   -0.8776 H   0  0  0  0  0  0
    \\    1.1321    0.0000    0.7550 H   0  0  0  0  0  0
    \\  1  2  1  0
    \\  1  3  1  0
    \\  1  4  1  0
    \\  1  5  1  0
    \\  2  6  1  0
    \\M  END
    \\$$$$
;

test "parse single molecule (methanol)" {
    var dict = try parseSdf(testing.allocator, methanol_sdf);
    defer dict.deinit();

    const mol = dict.get("methanol");
    try testing.expect(mol != null);

    const m = mol.?;
    try testing.expectEqual(@as(usize, 6), m.atoms.len);
    try testing.expectEqual(@as(usize, 5), m.bonds.len);
    try testing.expectEqualStrings("non-polymer", m.comp_type);
    try testing.expectEqualStrings("methanol", m.comp_id);
}

test "atom names are generated correctly" {
    var dict = try parseSdf(testing.allocator, methanol_sdf);
    defer dict.deinit();

    const m = dict.get("methanol").?;

    // atom 0: first carbon → "C1  "
    try testing.expectEqualSlices(u8, "C1", m.atoms[0].nameSlice());
    // atom 1: first oxygen → "O1  "
    try testing.expectEqualSlices(u8, "O1", m.atoms[1].nameSlice());
    // atom 2: first hydrogen → "H1  "
    try testing.expectEqualSlices(u8, "H1", m.atoms[2].nameSlice());
    // atom 5: fourth hydrogen → "H4  "
    try testing.expectEqualSlices(u8, "H4", m.atoms[5].nameSlice());
}

test "element symbols are stored correctly" {
    var dict = try parseSdf(testing.allocator, methanol_sdf);
    defer dict.deinit();

    const m = dict.get("methanol").?;
    try testing.expectEqual([2]u8{ 'C', ' ' }, m.atoms[0].element_symbol);
    try testing.expectEqual([2]u8{ 'O', ' ' }, m.atoms[1].element_symbol);
    try testing.expectEqual([2]u8{ 'H', ' ' }, m.atoms[2].element_symbol);
}

test "coordinates are stored correctly" {
    var dict = try parseSdf(testing.allocator, methanol_sdf);
    defer dict.deinit();

    const m = dict.get("methanol").?;
    try testing.expectApproxEqAbs(@as(f32, -0.6622), m.atoms[0].ideal_x, 1e-3);
    try testing.expectApproxEqAbs(@as(f32, 0.5342), m.atoms[0].ideal_y, 1e-3);
    try testing.expectApproxEqAbs(@as(f32, 0.0), m.atoms[0].ideal_z, 1e-3);
}

test "bond orders are mapped correctly" {
    const sdf =
        \\ethanediol
        \\  test
        \\
        \\  3  2  0  0  0  0  0  0  0  0999 V2000
        \\    0.0000    0.0000    0.0000 C   0  0
        \\    1.5000    0.0000    0.0000 C   0  0
        \\    3.0000    0.0000    0.0000 O   0  0
        \\  1  2  2  0
        \\  2  3  3  0
        \\M  END
        \\$$$$
    ;
    var dict = try parseSdf(testing.allocator, sdf);
    defer dict.deinit();

    const m = dict.get("ethanediol").?;
    try testing.expectEqual(BondOrder.double, m.bonds[0].order);
    try testing.expectEqual(BondOrder.triple, m.bonds[1].order);
}

test "aromatic bond type 4" {
    const sdf =
        \\benzene
        \\  test
        \\
        \\  2  1  0  0  0  0  0  0  0  0999 V2000
        \\    0.0000    0.0000    0.0000 C   0  0
        \\    1.4000    0.0000    0.0000 C   0  0
        \\  1  2  4  0
        \\M  END
        \\$$$$
    ;
    var dict = try parseSdf(testing.allocator, sdf);
    defer dict.deinit();

    const m = dict.get("benzene").?;
    try testing.expectEqual(BondOrder.aromatic, m.bonds[0].order);
    try testing.expect(m.bonds[0].aromatic);
}

test "M CHG formal charges" {
    const sdf =
        \\charged
        \\  test
        \\
        \\  2  0  0  0  0  0  0  0  0  0999 V2000
        \\    0.0000    0.0000    0.0000 N   0  0
        \\    1.4000    0.0000    0.0000 O   0  0
        \\M  CHG  2  1  1  2 -1
        \\M  END
        \\$$$$
    ;
    var dict = try parseSdf(testing.allocator, sdf);
    defer dict.deinit();

    const m = dict.get("charged").?;
    try testing.expectEqual(@as(i8, 1), m.atoms[0].charge);
    try testing.expectEqual(@as(i8, -1), m.atoms[1].charge);
}

test "parse multiple molecules" {
    const sdf =
        \\mol_a
        \\  test
        \\
        \\  1  0  0  0  0  0  0  0  0  0999 V2000
        \\    0.0000    0.0000    0.0000 C   0  0
        \\M  END
        \\$$$$
        \\mol_b
        \\  test
        \\
        \\  2  1  0  0  0  0  0  0  0  0999 V2000
        \\    0.0000    0.0000    0.0000 N   0  0
        \\    1.0000    0.0000    0.0000 O   0  0
        \\  1  2  1  0
        \\M  END
        \\$$$$
    ;
    var dict = try parseSdf(testing.allocator, sdf);
    defer dict.deinit();

    try testing.expect(dict.get("mol_a") != null);
    try testing.expect(dict.get("mol_b") != null);

    const a = dict.get("mol_a").?;
    try testing.expectEqual(@as(usize, 1), a.atoms.len);
    try testing.expectEqual(@as(usize, 0), a.bonds.len);

    const b = dict.get("mol_b").?;
    try testing.expectEqual(@as(usize, 2), b.atoms.len);
    try testing.expectEqual(@as(usize, 1), b.bonds.len);
}

test "empty input returns empty dict" {
    var dict = try parseSdf(testing.allocator, "");
    defer dict.deinit();
    try testing.expectEqual(@as(usize, 0), dict.components.count());
}

test "whitespace-only name is skipped" {
    const sdf =
        \\
        \\  test
        \\
        \\  1  0  0  0  0  0  0  0  0  0999 V2000
        \\    0.0000    0.0000    0.0000 C   0  0
        \\M  END
        \\$$$$
    ;
    var dict = try parseSdf(testing.allocator, sdf);
    defer dict.deinit();
    try testing.expectEqual(@as(usize, 0), dict.components.count());
}

test "malformed counts line is skipped" {
    const sdf =
        \\bad_mol
        \\  test
        \\
        \\  NOTANUMBER V2000
        \\M  END
        \\$$$$
        \\good_mol
        \\  test
        \\
        \\  1  0  0  0  0  0  0  0  0  0999 V2000
        \\    0.0000    0.0000    0.0000 C   0  0
        \\M  END
        \\$$$$
    ;
    var dict = try parseSdf(testing.allocator, sdf);
    defer dict.deinit();

    // bad_mol skipped, good_mol parsed
    try testing.expect(dict.get("bad_mol") == null);
    try testing.expect(dict.get("good_mol") != null);
}

test "duplicate molecule names: last wins" {
    const sdf =
        \\dup
        \\  test
        \\
        \\  1  0  0  0  0  0  0  0  0  0999 V2000
        \\    0.0000    0.0000    0.0000 C   0  0
        \\M  END
        \\$$$$
        \\dup
        \\  test
        \\
        \\  2  0  0  0  0  0  0  0  0  0999 V2000
        \\    0.0000    0.0000    0.0000 N   0  0
        \\    1.0000    0.0000    0.0000 O   0  0
        \\M  END
        \\$$$$
    ;
    var dict = try parseSdf(testing.allocator, sdf);
    defer dict.deinit();

    const m = dict.get("dup").?;
    // Second molecule (2 atoms) replaces the first (1 atom)
    try testing.expectEqual(@as(usize, 2), m.atoms.len);
}

test "bond referencing out-of-range atom is silently skipped" {
    const sdf =
        \\range_test
        \\  test
        \\
        \\  2  2  0  0  0  0  0  0  0  0999 V2000
        \\    0.0000    0.0000    0.0000 C   0  0
        \\    1.0000    0.0000    0.0000 C   0  0
        \\  1  2  1  0
        \\  1 99  1  0
        \\M  END
        \\$$$$
    ;
    var dict = try parseSdf(testing.allocator, sdf);
    defer dict.deinit();

    const m = dict.get("range_test").?;
    // Only the valid bond (1-2) should be kept
    try testing.expectEqual(@as(usize, 1), m.bonds.len);
}

test "SDF without trailing separator" {
    const sdf =
        \\nosep
        \\  test
        \\
        \\  1  0  0  0  0  0  0  0  0  0999 V2000
        \\    0.0000    0.0000    0.0000 O   0  0
        \\M  END
    ;
    var dict = try parseSdf(testing.allocator, sdf);
    defer dict.deinit();

    const m = dict.get("nosep").?;
    try testing.expectEqual(@as(usize, 1), m.atoms.len);
}
