# PDB Format Support Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Enable `zreduce run input.pdb -o output.pdb` with full hydrogen placement and optimization, matching the existing mmCIF pipeline.

**Architecture:** Add a PDB parser (`src/pdb.zig`) that produces the same `Model` struct as `mmcif.parseModel()`, and a PDB writer (`src/writer/pdb_writer.zig`) that outputs fixed-width PDB ATOM/HETATM lines with passthrough of non-atom records. The run pipeline detects format by file extension and branches at parse/write stages only — all placement and optimization code is shared.

**Tech Stack:** Zig 0.15, existing Model/Chain/Residue/Atom types, existing gzip module, existing format helpers (`writeFixedFloat3`, `elementSymbol`, etc.)

---

### File Structure

| Action | Path | Responsibility |
|--------|------|---------------|
| Create | `src/pdb.zig` | PDB parser: ATOM/HETATM fixed-width → Model, passthrough record storage |
| Create | `src/writer/pdb_writer.zig` | PDB writer: Model → fixed-width ATOM/HETATM, passthrough, serial renumbering |
| Create | `src/test_data/tiny.pdb` | Minimal PDB test fixture (ALA, 5 atoms) |
| Create | `src/test_data/multi_chain.pdb` | Multi-chain PDB test fixture (ALA+GLY+VAL) |
| Create | `src/test_data/hetatm.pdb` | HETATM + water test fixture |
| Modify | `src/root.zig` | Add `pub const pdb = @import("pdb.zig");` |
| Modify | `src/writer.zig` | Add `pub const pdb_writer = @import("writer/pdb_writer.zig");` |
| Modify | `src/run.zig` | Format detection, PDB parse/write branches |
| Modify | `src/main.zig` | Update help text |

---

### Task 1: PDB Test Fixtures

Create minimal PDB files as test data. These mirror the existing mmCIF test fixtures.

**Files:**
- Create: `src/test_data/tiny.pdb`
- Create: `src/test_data/multi_chain.pdb`
- Create: `src/test_data/hetatm.pdb`

- [ ] **Step 1: Create tiny.pdb (single ALA residue, 5 atoms)**

```
HEADER    TEST STRUCTURE
ATOM      1  N   ALA A   1       1.000   2.000   3.000  1.00 10.00           N
ATOM      2  CA  ALA A   1       2.000   3.000   4.000  1.00 10.00           C
ATOM      3  C   ALA A   1       3.000   4.000   5.000  1.00 10.00           C
ATOM      4  O   ALA A   1       4.000   5.000   6.000  1.00 10.00           O
ATOM      5  CB  ALA A   1       2.500   2.500   3.500  1.00 10.00           C
TER       6      ALA A   1
END
```

Write to `src/test_data/tiny.pdb`. Coordinates match `tiny.cif` exactly.

- [ ] **Step 2: Create multi_chain.pdb (2 chains, 3 residues, 11 atoms)**

```
HEADER    MULTI CHAIN TEST
ATOM      1  N   ALA A   1       1.000   2.000   3.000  1.00 10.00           N
ATOM      2  CA  ALA A   1       2.000   3.000   4.000  1.00 10.00           C
ATOM      3  C   ALA A   1       3.000   4.000   5.000  1.00 10.00           C
ATOM      4  CB  ALA A   1       2.500   2.500   3.500  1.00 10.00           C
ATOM      5  N   GLY A   2       4.000   5.000   6.000  1.00 10.00           N
ATOM      6  CA  GLY A   2       5.000   6.000   7.000  1.00 10.00           C
ATOM      7  C   GLY A   2       6.000   7.000   8.000  1.00 10.00           C
ATOM      8  O   GLY A   2       7.000   8.000   9.000  1.00 10.00           O
TER       9      GLY A   2
ATOM     10  N   VAL B   1      10.000  11.000  12.000  1.00 10.00           N
ATOM     11  CA  VAL B   1      11.000  12.000  13.000  1.00 10.00           C
ATOM     12  CB  VAL B   1      13.000  13.000  13.000  1.00 10.00           C
TER      13      VAL B   1
END
```

Write to `src/test_data/multi_chain.pdb`. Coordinates match `multi_chain.cif`.

- [ ] **Step 3: Create hetatm.pdb (polymer + HETATM ligand + water)**

```
HEADER    HETATM TEST
ATOM      1  N   ALA A   1       1.000   2.000   3.000  1.00 10.00           N
ATOM      2  CA  ALA A   1       2.000   3.000   4.000  1.00 10.00           C
ATOM      3  C   ALA A   1       3.000   4.000   5.000  1.00 10.00           C
ATOM      4  O   ALA A   1       4.000   5.000   6.000  1.00 10.00           O
ATOM      5  CB  ALA A   1       2.500   2.500   3.500  1.00 10.00           C
TER       6      ALA A   1
HETATM    7  C1  EDO B   1      10.000  10.000  10.000  1.00 20.00           C
HETATM    8  O1  EDO B   1      11.000  10.000  10.000  1.00 20.00           O
HETATM    9  O   HOH C   1      20.000  20.000  20.000  1.00 30.00           O
END
```

Write to `src/test_data/hetatm.pdb`.

- [ ] **Step 4: Commit test fixtures**

```bash
git add src/test_data/tiny.pdb src/test_data/multi_chain.pdb src/test_data/hetatm.pdb
git commit -m "test: add PDB format test fixtures"
```

---

### Task 2: PDB Parser — Core ATOM/HETATM Parsing

Parse PDB fixed-width ATOM/HETATM records into the existing `Model` struct.

**Files:**
- Create: `src/pdb.zig`
- Modify: `src/root.zig`

- [ ] **Step 1: Write failing test for single-residue PDB parsing**

In `src/pdb.zig`, add the test at the end of the file (the implementation will go above):

```zig
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
```

- [ ] **Step 2: Run test to verify it fails**

Run: `zig build test --summary all -Dtest-filter="parse tiny PDB" 2>&1 | tail -5`
Expected: compilation error (parseModel not defined)

- [ ] **Step 3: Implement PDB parser core**

In `src/pdb.zig`, implement:

```zig
//! PDB format parser: fixed-width ATOM/HETATM records into Model.

const std = @import("std");
const Allocator = std.mem.Allocator;
const model_mod = @import("model.zig");
const element = @import("element.zig");

const Model = model_mod.Model;
const Atom = model_mod.Atom;
const Residue = model_mod.Residue;
const Chain = model_mod.Chain;

pub const PdbError = error{
    InvalidCoordinateValue,
    OutOfMemory,
    EmptyFile,
};

/// A stored line from the PDB file for passthrough output.
pub const PdbRecord = union(enum) {
    /// Placeholder marking where ATOM/HETATM/TER records for a chain go.
    atom_site,
    /// Any non-ATOM line preserved verbatim (HEADER, REMARK, HELIX, etc.)
    raw_line: []const u8,
};

/// Result of parsing a PDB file: the model plus passthrough records.
pub const PdbParseResult = struct {
    model: Model,
    records: std.ArrayListUnmanaged(PdbRecord),
    source: []const u8, // kept alive for raw_line slices

    pub fn deinit(self: *PdbParseResult) void {
        self.model.deinit();
        self.records.deinit(self.model.allocator);
    }
};

/// Standard amino acid and nucleotide comp_ids for entity type heuristic.
fn isStandardPolymerComp(comp_id: []const u8) bool {
    const standard_aa = [_][]const u8{
        "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY",
        "HIS", "ILE", "LEU", "LYS", "MET", "PHE", "PRO", "SER",
        "THR", "TRP", "TYR", "VAL", "MSE", "SEC", "PYL",
    };
    const standard_na = [_][]const u8{
        "A", "C", "G", "U", "DA", "DC", "DG", "DT", "DI",
    };
    for (standard_aa) |aa| {
        if (std.ascii.eqlIgnoreCase(comp_id, aa)) return true;
    }
    for (standard_na) |na| {
        if (std.ascii.eqlIgnoreCase(comp_id, na)) return true;
    }
    return false;
}

fn trimRight(s: []const u8) []const u8 {
    var end = s.len;
    while (end > 0 and s[end - 1] == ' ') end -= 1;
    return s[0..end];
}

fn trimLeft(s: []const u8) []const u8 {
    var start: usize = 0;
    while (start < s.len and s[start] == ' ') start += 1;
    return s[start..];
}

fn trim(s: []const u8) []const u8 {
    return trimLeft(trimRight(s));
}

/// Parse a PDB-format integer field, ignoring leading/trailing spaces.
fn parseIntField(field: []const u8) ?i32 {
    const trimmed = trim(field);
    if (trimmed.len == 0) return null;
    return std.fmt.parseInt(i32, trimmed, 10) catch null;
}

/// Parse a PDB-format float field, ignoring leading/trailing spaces.
fn parseFloatField(field: []const u8) ?f32 {
    const trimmed = trim(field);
    if (trimmed.len == 0) return null;
    return std.fmt.parseFloat(f32, trimmed) catch null;
}

/// Safely extract a substring from a line, returning empty if out of bounds.
fn safeSlice(line: []const u8, start: usize, end: usize) []const u8 {
    const s = @min(start, line.len);
    const e = @min(end, line.len);
    if (s >= e) return "";
    return line[s..e];
}

/// Parse a PDB source string into a Model.
/// Only MODEL 1 is extracted from multi-model files.
pub fn parseModel(allocator: Allocator, source: []const u8) PdbError!Model {
    const result = try parse(allocator, source);
    // Caller only wants the model — free the records list.
    var records = result.records;
    records.deinit(allocator);
    return result.model;
}

/// Parse a PDB source string into a Model plus passthrough records.
pub fn parse(allocator: Allocator, source: []const u8) PdbError!PdbParseResult {
    var mdl = Model.init(allocator);
    errdefer mdl.deinit();

    var records = std.ArrayListUnmanaged(PdbRecord).empty;
    errdefer records.deinit(allocator);

    // State tracking
    var cur_chain_id: u8 = 0; // 0 = no chain yet
    var cur_seq_num: i32 = -999999;
    var cur_ins_code: u8 = 0;
    var cur_comp_id: [5]u8 = .{ 0, 0, 0, 0, 0 };
    var in_chain = false;
    var in_residue = false;
    var atom_site_marker_added = false;

    // Multi-model: track whether we are inside MODEL 1 or should skip.
    var seen_model = false;
    var in_active_model = true; // true until we see MODEL with number != 1

    var line_iter = std.mem.splitScalar(u8, source, '\n');
    while (line_iter.next()) |raw_line| {
        // Strip trailing CR for DOS line endings
        const line = if (raw_line.len > 0 and raw_line[raw_line.len - 1] == '\r')
            raw_line[0 .. raw_line.len - 1]
        else
            raw_line;

        if (line.len < 6) {
            try records.append(allocator, .{ .raw_line = line });
            continue;
        }

        const rec_type = safeSlice(line, 0, 6);

        // MODEL record handling
        if (std.mem.startsWith(u8, rec_type, "MODEL")) {
            seen_model = true;
            const model_num = parseIntField(safeSlice(line, 6, 14));
            in_active_model = (model_num == null or model_num.? == 1);
            if (in_active_model) {
                try records.append(allocator, .{ .raw_line = line });
            }
            continue;
        }
        if (std.mem.startsWith(u8, rec_type, "ENDMDL")) {
            if (in_active_model and seen_model) {
                try records.append(allocator, .{ .raw_line = line });
                // Done with MODEL 1 — skip remaining models
                in_active_model = false;
            }
            continue;
        }

        // Skip non-MODEL-1 atoms in multi-model files
        if (seen_model and !in_active_model) continue;

        const is_atom = std.mem.eql(u8, rec_type, "ATOM  ");
        const is_hetatm = std.mem.eql(u8, rec_type, "HETATM");

        if (is_atom or is_hetatm) {
            // Add atom_site marker on first ATOM/HETATM
            if (!atom_site_marker_added) {
                try records.append(allocator, .atom_site);
                atom_site_marker_added = true;
            }

            // Parse fixed-width columns (0-indexed):
            //  6-10  serial (ignored, will renumber)
            // 12-15  atom name
            // 16     altloc
            // 17-19  residue name
            // 21     chain ID
            // 22-25  residue seq number
            // 26     insertion code
            // 30-37  x
            // 38-45  y
            // 46-53  z
            // 54-59  occupancy
            // 60-65  b-factor
            // 76-77  element symbol

            const atom_name_raw = safeSlice(line, 12, 16);
            const altloc_raw = if (line.len > 16) line[16] else ' ';
            const comp_id_raw = trim(safeSlice(line, 17, 20));
            const chain_id: u8 = if (line.len > 21) line[21] else ' ';
            const seq_num = parseIntField(safeSlice(line, 22, 26)) orelse 0;
            const ins_code: u8 = if (line.len > 26) line[26] else ' ';

            const x = parseFloatField(safeSlice(line, 30, 38)) orelse return PdbError.InvalidCoordinateValue;
            const y = parseFloatField(safeSlice(line, 38, 46)) orelse return PdbError.InvalidCoordinateValue;
            const z = parseFloatField(safeSlice(line, 46, 54)) orelse return PdbError.InvalidCoordinateValue;

            const occ = parseFloatField(safeSlice(line, 54, 60)) orelse 1.0;
            const bfac = parseFloatField(safeSlice(line, 60, 66)) orelse 0.0;

            // Element symbol: columns 76-77 (right-justified)
            const elem_raw = trim(safeSlice(line, 76, 78));

            // Chain boundary: new chain ID
            if (!in_chain or chain_id != cur_chain_id) {
                // Close previous residue
                if (in_residue) {
                    const atom_end: u32 = @intCast(mdl.atoms.items.len);
                    mdl.residues.items[mdl.residues.items.len - 1].atom_end = atom_end;
                    in_residue = false;
                }
                // Close previous chain
                if (in_chain) {
                    mdl.chains.items[mdl.chains.items.len - 1].residue_end = @intCast(mdl.residues.items.len);
                }
                // Start new chain
                var new_chain = Chain{};
                const chain_str = &[_]u8{chain_id};
                new_chain.setLabelAsymId(chain_str);
                new_chain.setAuthAsymId(chain_str);
                // Assign sequential entity_id
                var entity_buf: [4]u8 = undefined;
                const entity_str = std.fmt.bufPrint(&entity_buf, "{d}", .{mdl.chains.items.len + 1}) catch "1";
                new_chain.setEntityId(entity_str);
                new_chain.residue_start = @intCast(mdl.residues.items.len);
                try mdl.chains.append(allocator, new_chain);
                cur_chain_id = chain_id;
                in_chain = true;
                // Force new residue
                cur_seq_num = -999999;
                cur_ins_code = 0;
                @memset(&cur_comp_id, 0);
            }

            // Residue boundary: seq_num + ins_code + comp_id change
            var comp_id_buf: [5]u8 = .{ ' ', ' ', ' ', ' ', ' ' };
            const comp_len = @min(comp_id_raw.len, 5);
            @memcpy(comp_id_buf[0..comp_len], comp_id_raw[0..comp_len]);

            const new_residue = !in_residue or
                seq_num != cur_seq_num or
                ins_code != cur_ins_code or
                !std.mem.eql(u8, &comp_id_buf, &cur_comp_id);

            if (new_residue) {
                if (in_residue) {
                    mdl.residues.items[mdl.residues.items.len - 1].atom_end = @intCast(mdl.atoms.items.len);
                }
                var new_res = Residue{};
                new_res.setCompId(comp_id_raw);
                new_res.seq_id = seq_num;
                new_res.auth_seq_id = seq_num;
                new_res.ins_code = if (ins_code != ' ') ins_code else ' ';
                new_res.chain_idx = @intCast(mdl.chains.items.len - 1);
                new_res.atom_start = @intCast(mdl.atoms.items.len);
                new_res.atom_end = @intCast(mdl.atoms.items.len);
                try mdl.residues.append(allocator, new_res);
                cur_seq_num = seq_num;
                cur_ins_code = ins_code;
                cur_comp_id = comp_id_buf;
                in_residue = true;
            }

            // Build atom
            var atom = Atom{ .pos = .{ .x = x, .y = y, .z = z } };

            // PDB atom name: columns 12-15, preserve as-is (space-padded).
            // Trim for setName (which stores trimmed).
            const atom_name_trimmed = trim(atom_name_raw);
            atom.setName(atom_name_trimmed);

            // Element: prefer column 76-77, fall back to first non-space in name.
            const elem_str = if (elem_raw.len > 0) elem_raw else blk: {
                // Fallback: first alpha chars of atom name
                const t = trimLeft(atom_name_raw);
                if (t.len >= 2 and std.ascii.isAlphabetic(t[0]) and std.ascii.isAlphabetic(t[1])) {
                    break :blk t[0..2];
                } else if (t.len >= 1 and std.ascii.isAlphabetic(t[0])) {
                    break :blk t[0..1];
                }
                break :blk "X";
            };
            atom.element_type = element.elementFromSymbol(elem_str);
            atom.is_hydrogen = switch (atom.element_type) {
                .H, .Har, .Hpol, .Ha_p, .HOd => true,
                else => false,
            };

            atom.occupancy = occ;
            atom.b_factor = bfac;
            atom.serial = @intCast(mdl.atoms.items.len);
            atom.residue_idx = @intCast(mdl.residues.items.len - 1);

            const altloc = if (altloc_raw == ' ' or altloc_raw == '.') @as(u8, ' ') else altloc_raw;
            atom.altloc = altloc;

            try mdl.atoms.append(allocator, atom);
            continue;
        }

        // TER record: just continue (chain boundary detected by chain_id change)
        if (std.mem.startsWith(u8, rec_type, "TER")) {
            // Don't store TER — we regenerate them on output
            continue;
        }

        // CONECT/MASTER: drop
        if (std.mem.startsWith(u8, rec_type, "CONECT") or
            std.mem.startsWith(u8, rec_type, "MASTER"))
        {
            continue;
        }

        // END record: stop parsing
        if (std.mem.eql(u8, rec_type, "END   ") or
            std.mem.eql(u8, trimRight(rec_type), "END"))
        {
            try records.append(allocator, .{ .raw_line = line });
            break;
        }

        // All other records: passthrough
        try records.append(allocator, .{ .raw_line = line });
    }

    // Close final residue
    if (in_residue) {
        mdl.residues.items[mdl.residues.items.len - 1].atom_end = @intCast(mdl.atoms.items.len);
    }

    // Close final chain
    if (in_chain) {
        mdl.chains.items[mdl.chains.items.len - 1].residue_end = @intCast(mdl.residues.items.len);
    }

    mdl.original_atom_count = @intCast(mdl.atoms.items.len);

    // Apply entity type heuristic
    for (mdl.residues.items) |*res| {
        const cid = res.compIdSlice();
        if (std.ascii.eqlIgnoreCase(cid, "HOH") or std.ascii.eqlIgnoreCase(cid, "DOD")) {
            res.entity_type = .water;
        } else if (isStandardPolymerComp(cid)) {
            res.entity_type = .polymer;
        } else {
            res.entity_type = .non_polymer;
        }
    }

    // Detect chain breaks from sequence number gaps within polymer chains.
    for (mdl.chains.items) |chain| {
        const res_slice = mdl.residues.items[chain.residue_start..chain.residue_end];
        if (res_slice.len < 2) continue;
        for (1..res_slice.len) |i| {
            if (res_slice[i].entity_type != .polymer) continue;
            if (res_slice[i - 1].entity_type != .polymer) continue;
            const gap = res_slice[i].seq_id - res_slice[i - 1].seq_id;
            if (gap > 1) {
                res_slice[i].is_chain_break_before = true;
            }
        }
    }

    return .{
        .model = mdl,
        .records = records,
        .source = source,
    };
}
```

- [ ] **Step 4: Run test to verify it passes**

Run: `zig build test --summary all -Dtest-filter="parse tiny PDB" 2>&1 | tail -5`
Expected: PASS

- [ ] **Step 5: Write and run multi-chain test**

Add to `src/pdb.zig`:

```zig
test "parse multi-chain PDB" {
    const source = @embedFile("test_data/multi_chain.pdb");
    var mdl = try parseModel(testing.allocator, source);
    defer mdl.deinit();

    try testing.expectEqual(@as(usize, 12), mdl.atoms.items.len);
    try testing.expectEqual(@as(usize, 3), mdl.residues.items.len);
    try testing.expectEqual(@as(usize, 2), mdl.chains.items.len);

    // Chain A: 2 residues
    try testing.expectEqualStrings("A", mdl.chains.items[0].labelSlice());
    try testing.expectEqual(@as(u32, 0), mdl.chains.items[0].residue_start);
    try testing.expectEqual(@as(u32, 2), mdl.chains.items[0].residue_end);

    // Chain B: 1 residue
    try testing.expectEqualStrings("B", mdl.chains.items[1].labelSlice());
    try testing.expectEqual(@as(u32, 2), mdl.chains.items[1].residue_start);
    try testing.expectEqual(@as(u32, 3), mdl.chains.items[1].residue_end);

    // Last atom coordinate
    try testing.expectApproxEqAbs(mdl.atoms.items[11].pos.x, 13.0, 1e-3);
}
```

Run: `zig build test --summary all -Dtest-filter="parse multi-chain PDB" 2>&1 | tail -5`
Expected: PASS

- [ ] **Step 6: Write and run HETATM/entity type test**

Add to `src/pdb.zig`:

```zig
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
```

Run: `zig build test --summary all -Dtest-filter="parse HETATM" 2>&1 | tail -5`
Expected: PASS

- [ ] **Step 7: Register in root.zig**

In `src/root.zig`, add after the `mmcif` import:

```zig
pub const pdb = @import("pdb.zig");
```

- [ ] **Step 8: Run full test suite to verify no regressions**

Run: `zig build test --summary all 2>&1 | tail -10`
Expected: All existing tests pass, plus new PDB tests.

- [ ] **Step 9: Commit**

```bash
git add src/pdb.zig src/root.zig
git commit -m "feat: add PDB format parser (ATOM/HETATM/TER/MODEL)"
```

---

### Task 3: PDB Writer

Write Model atoms back to PDB fixed-width format with passthrough of non-atom records.

**Files:**
- Create: `src/writer/pdb_writer.zig`
- Modify: `src/writer.zig`

- [ ] **Step 1: Write failing test for PDB writer**

Create `src/writer/pdb_writer.zig` with the test first:

```zig
test "write tiny PDB round-trip" {
    const pdb_mod = @import("../pdb.zig");
    const source = @embedFile("../test_data/tiny.pdb");
    var result = try pdb_mod.parse(testing.allocator, source);
    defer result.deinit();

    var buf = std.ArrayList(u8).empty;
    defer buf.deinit(testing.allocator);
    const w = buf.writer(testing.allocator);
    try writeModel(w, &result.model, result.records.items, .hydrogen);

    const output = buf.items;

    // Should contain HEADER passthrough
    try testing.expect(std.mem.indexOf(u8, output, "HEADER") != null);
    // Should contain 5 ATOM lines
    var atom_count: usize = 0;
    var line_it = std.mem.splitScalar(u8, output, '\n');
    while (line_it.next()) |line| {
        if (line.len >= 4 and std.mem.eql(u8, line[0..4], "ATOM")) atom_count += 1;
    }
    try testing.expectEqual(@as(usize, 5), atom_count);
    // Should end with END
    try testing.expect(std.mem.indexOf(u8, output, "END") != null);
}
```

- [ ] **Step 2: Implement PDB writer**

Add the implementation above the test in `src/writer/pdb_writer.zig`:

```zig
//! PDB format writer: outputs Model atoms as fixed-width ATOM/HETATM records
//! with passthrough of non-atom records from the original PDB file.

const std = @import("std");
const model_mod = @import("../model.zig");
const Model = model_mod.Model;
const Atom = model_mod.Atom;
const Residue = model_mod.Residue;
const Chain = model_mod.Chain;
const pdb_mod = @import("../pdb.zig");
const PdbRecord = pdb_mod.PdbRecord;
const element = @import("../element.zig");
const format = @import("format.zig");
const mover_mod = @import("../optimize/mover.zig");
const place = @import("../place.zig");

const elementSymbol = format.elementSymbol;

/// Write the model to PDB format with passthrough records.
pub fn writeModel(
    writer: anytype,
    model: *const Model,
    records: []const PdbRecord,
    output_isotope: place.OutputIsotope,
) !void {
    for (records) |record| {
        switch (record) {
            .raw_line => |line| {
                try writer.writeAll(line);
                try writer.writeByte('\n');
            },
            .atom_site => {
                try writeAtomRecords(writer, model, output_isotope);
            },
        }
    }
}

/// Write all ATOM/HETATM records grouped by chain, with TER records between chains.
fn writeAtomRecords(writer: anytype, model: *const Model, output_isotope: place.OutputIsotope) !void {
    // Pre-index added H atoms by residue_idx (same approach as mmcif_writer).
    const n_residues = model.residues.items.len;
    const allocator = model.allocator;

    const added_counts = try allocator.alloc(u32, n_residues + 1);
    defer allocator.free(added_counts);
    @memset(added_counts, 0);
    for (model.atoms.items) |atom| {
        if (!atom.is_added) continue;
        added_counts[atom.residue_idx] += 1;
    }

    const added_offsets = try allocator.alloc(u32, n_residues + 1);
    defer allocator.free(added_offsets);
    added_offsets[0] = 0;
    for (0..n_residues) |r| {
        added_offsets[r + 1] = added_offsets[r] + added_counts[r];
    }
    const total_added = added_offsets[n_residues];

    const added_indices = try allocator.alloc(u32, total_added);
    defer allocator.free(added_indices);
    @memset(added_counts, 0);
    for (model.atoms.items, 0..) |atom, atom_idx| {
        if (!atom.is_added) continue;
        const r = atom.residue_idx;
        added_indices[added_offsets[r] + added_counts[r]] = @intCast(atom_idx);
        added_counts[r] += 1;
    }

    var serial: u32 = 1;

    for (model.chains.items) |chain| {
        for (model.residues.items[chain.residue_start..chain.residue_end], chain.residue_start..) |res, res_idx| {
            const is_hetatm = (res.entity_type == .non_polymer or
                res.entity_type == .water or
                res.entity_type == .unknown);

            // Original heavy atoms
            for (model.atoms.items[res.atom_start..res.atom_end]) |atom| {
                try writeAtomLine(writer, atom, res, chain, serial, is_hetatm, output_isotope);
                serial += 1;
            }

            // Added H atoms
            const h_start = added_offsets[res_idx];
            const h_end = added_offsets[res_idx + 1];
            for (added_indices[h_start..h_end]) |atom_idx| {
                const atom = model.atoms.items[atom_idx];
                if (mover_mod.isAbsentH(atom)) continue;
                try writeAtomLine(writer, atom, res, chain, serial, is_hetatm, output_isotope);
                serial += 1;
            }
        }

        // TER after each chain
        if (chain.residue_end > chain.residue_start) {
            const last_res = model.residues.items[chain.residue_end - 1];
            try writeTerLine(writer, last_res, chain, serial);
            serial += 1;
        }
    }
}

/// Write a single ATOM or HETATM line in PDB fixed-width format.
fn writeAtomLine(
    writer: anytype,
    atom: Atom,
    res: Residue,
    chain: Chain,
    serial: u32,
    is_hetatm: bool,
    output_isotope: place.OutputIsotope,
) !void {
    // Record type (columns 1-6)
    if (is_hetatm and !atom.is_added) {
        try writer.writeAll("HETATM");
    } else {
        try writer.writeAll("ATOM  ");
    }

    // Serial number (columns 7-11, right-justified)
    try writer.print("{d: >5}", .{serial});

    // Space (column 12)
    try writer.writeByte(' ');

    // Atom name (columns 13-16, left-justified in 4-char field)
    // PDB convention: 1-char elements start at col 14, 2-char at col 13
    const name = atom.nameSlice();
    const elem_sym = atomTypeSymbol(atom, output_isotope);
    if (name.len < 4 and elem_sym.len == 1) {
        try writer.writeByte(' ');
        try writer.writeAll(name);
        const pad = 3 - name.len;
        for (0..pad) |_| try writer.writeByte(' ');
    } else {
        try writer.writeAll(name);
        if (name.len < 4) {
            const pad = 4 - name.len;
            for (0..pad) |_| try writer.writeByte(' ');
        }
    }

    // Altloc (column 17)
    if (atom.altloc == ' ') {
        try writer.writeByte(' ');
    } else {
        try writer.writeByte(atom.altloc);
    }

    // Residue name (columns 18-20, right-justified in 3-char field)
    const comp = res.compIdSlice();
    if (comp.len < 3) {
        const pad = 3 - comp.len;
        for (0..pad) |_| try writer.writeByte(' ');
    }
    try writer.writeAll(comp[0..@min(comp.len, 3)]);

    // Space (column 21)
    try writer.writeByte(' ');

    // Chain ID (column 22)
    const chain_id = chain.authSlice();
    if (chain_id.len > 0) {
        try writer.writeByte(chain_id[0]);
    } else {
        try writer.writeByte(' ');
    }

    // Residue sequence number (columns 23-26, right-justified)
    try writer.print("{d: >4}", .{res.auth_seq_id});

    // Insertion code (column 27)
    if (res.ins_code != ' ' and res.ins_code != 0) {
        try writer.writeByte(res.ins_code);
    } else {
        try writer.writeByte(' ');
    }

    // Spaces (columns 28-30)
    try writer.writeAll("   ");

    // Coordinates (columns 31-54: x, y, z each 8.3f)
    try writePdbFloat83(writer, atom.pos.x);
    try writePdbFloat83(writer, atom.pos.y);
    try writePdbFloat83(writer, atom.pos.z);

    // Occupancy (columns 55-60: 6.2f)
    try writePdbFloat62(writer, atom.occupancy);

    // B-factor (columns 61-66: 6.2f)
    try writePdbFloat62(writer, atom.b_factor);

    // Spaces (columns 67-76)
    try writer.writeAll("          ");

    // Element symbol (columns 77-78, right-justified)
    if (elem_sym.len == 1) {
        try writer.writeByte(' ');
        try writer.writeAll(elem_sym);
    } else {
        try writer.writeAll(elem_sym[0..2]);
    }

    // Newline
    try writer.writeByte('\n');
}

/// Write a TER record.
fn writeTerLine(writer: anytype, last_res: Residue, chain: Chain, serial: u32) !void {
    try writer.writeAll("TER   ");
    try writer.print("{d: >5}", .{serial});
    try writer.writeAll("      "); // cols 12-16 + col 17
    const comp = last_res.compIdSlice();
    if (comp.len < 3) {
        const pad = 3 - comp.len;
        for (0..pad) |_| try writer.writeByte(' ');
    }
    try writer.writeAll(comp[0..@min(comp.len, 3)]);
    try writer.writeByte(' ');
    const chain_id = chain.authSlice();
    if (chain_id.len > 0) {
        try writer.writeByte(chain_id[0]);
    } else {
        try writer.writeByte(' ');
    }
    try writer.print("{d: >4}", .{last_res.auth_seq_id});
    if (last_res.ins_code != ' ' and last_res.ins_code != 0) {
        try writer.writeByte(last_res.ins_code);
    }
    try writer.writeByte('\n');
}

/// Get the element symbol string, respecting isotope settings.
fn atomTypeSymbol(atom: Atom, output_isotope: place.OutputIsotope) []const u8 {
    if (atom.is_hydrogen and atom.is_added and output_isotope == .deuterium) {
        return "D";
    }
    return elementSymbol(atom.element_type);
}

/// Write a float in PDB 8.3f format (right-justified, 8 chars total).
fn writePdbFloat83(writer: anytype, val: f32) !void {
    var buf: [16]u8 = undefined;
    var fbs = std.io.fixedBufferStream(&buf);
    const w = fbs.writer();
    format.writeFixedFloat3(w, val) catch {
        try writer.writeAll("   0.000");
        return;
    };
    const s = fbs.getWritten();
    if (s.len < 8) {
        const pad = 8 - s.len;
        for (0..pad) |_| try writer.writeByte(' ');
    }
    try writer.writeAll(s);
}

/// Write a float in PDB 6.2f format (right-justified, 6 chars total).
fn writePdbFloat62(writer: anytype, val: f32) !void {
    var buf: [16]u8 = undefined;
    var fbs = std.io.fixedBufferStream(&buf);
    const w = fbs.writer();
    format.writeFixedFloat2(w, val) catch {
        try writer.writeAll("  0.00");
        return;
    };
    const s = fbs.getWritten();
    if (s.len < 6) {
        const pad = 6 - s.len;
        for (0..pad) |_| try writer.writeByte(' ');
    }
    try writer.writeAll(s);
}
```

- [ ] **Step 3: Register in writer.zig**

In `src/writer.zig`, add:

```zig
pub const pdb_writer = @import("writer/pdb_writer.zig");
```

- [ ] **Step 4: Run test**

Run: `zig build test --summary all -Dtest-filter="write tiny PDB" 2>&1 | tail -10`
Expected: PASS

- [ ] **Step 5: Write round-trip coordinate accuracy test**

Add to `src/writer/pdb_writer.zig`:

```zig
test "PDB writer preserves coordinates" {
    const pdb_mod = @import("../pdb.zig");
    const source = @embedFile("../test_data/tiny.pdb");
    var result = try pdb_mod.parse(testing.allocator, source);
    defer result.deinit();

    var buf = std.ArrayList(u8).empty;
    defer buf.deinit(testing.allocator);
    const w = buf.writer(testing.allocator);
    try writeModel(w, &result.model, result.records.items, .hydrogen);

    // Re-parse the output
    var mdl2 = try pdb_mod.parseModel(testing.allocator, buf.items);
    defer mdl2.deinit();

    try testing.expectEqual(result.model.atoms.items.len, mdl2.atoms.items.len);
    for (result.model.atoms.items, mdl2.atoms.items) |a1, a2| {
        try testing.expectApproxEqAbs(a1.pos.x, a2.pos.x, 1e-3);
        try testing.expectApproxEqAbs(a1.pos.y, a2.pos.y, 1e-3);
        try testing.expectApproxEqAbs(a1.pos.z, a2.pos.z, 1e-3);
    }
}
```

Run: `zig build test --summary all -Dtest-filter="PDB writer preserves" 2>&1 | tail -5`
Expected: PASS

- [ ] **Step 6: Run full test suite**

Run: `zig build test --summary all 2>&1 | tail -10`
Expected: All tests pass.

- [ ] **Step 7: Commit**

```bash
git add src/writer/pdb_writer.zig src/writer.zig
git commit -m "feat: add PDB format writer with passthrough support"
```

---

### Task 4: Format Detection and Pipeline Integration

Wire PDB parser and writer into `run.zig` to enable end-to-end `zreduce run input.pdb -o output.pdb`.

**Files:**
- Modify: `src/run.zig`
- Modify: `src/main.zig`

- [ ] **Step 1: Add format detection helper to run.zig**

Add a format enum and detection function near the top of `src/run.zig` (after the imports):

```zig
pub const InputFormat = enum {
    mmcif,
    pdb,
};

/// Detect input format from file extension.
pub fn detectFormat(path: []const u8) InputFormat {
    const lower_path = path; // extensions are typically already lowercase
    if (std.mem.endsWith(u8, lower_path, ".pdb") or
        std.mem.endsWith(u8, lower_path, ".pdb.gz") or
        std.mem.endsWith(u8, lower_path, ".ent") or
        std.mem.endsWith(u8, lower_path, ".ent.gz"))
    {
        return .pdb;
    }
    return .mmcif;
}
```

- [ ] **Step 2: Add PDB processing branch in processFile**

Refactor `processFile` in `src/run.zig` to branch based on detected format. The key changes:

1. After reading the source, detect format.
2. For PDB: use `pdb.parse()` instead of `cif.readString()` + `mmcif.parseModel()`.
3. Skip mmCIF-specific steps for PDB (struct_conn, branch_links, inline_comp).
4. For PDB output: use `pdb_writer.writeModel()` instead of `mmcif_writer`.

The processFile function should be modified to add an `else` branch after format detection. The PDB path skips: CIF document parsing, struct_conn, branch_links, inline components, atom_lookup. It does share: applyChemistry, addHydrogens, generateMovers, optimize, validate.

For the write step, detect output format from the output path extension (or inherit from input format if no output path).

- [ ] **Step 3: Update ProcessConfig to store format info**

Add to `ProcessConfig`:

```zig
format: zreduce.run.InputFormat = .mmcif,
```

This is set by the caller based on `detectFormat(input_path)`.

- [ ] **Step 4: Update main.zig runSubcommand to set format**

In `runSubcommand` in `src/main.zig`, add format detection when building `proc_config`:

```zig
.format = zreduce.run.detectFormat(config.input_path),
```

- [ ] **Step 5: Update help text in main.zig**

Change the usage strings:
- `printUsage`: change "Hydrogen placement for mmCIF structures" to "Hydrogen placement for mmCIF/PDB structures"
- `printRunUsage`: change `<input.cif>` to `<input.cif|input.pdb>`
- `printBatchUsage`: similar update

- [ ] **Step 6: Run full test suite**

Run: `zig build test --summary all 2>&1 | tail -10`
Expected: All tests pass.

- [ ] **Step 7: Build and smoke test with a real PDB file**

Run: `zig build -Doptimize=ReleaseFast`

If a real PDB file is available, test:
```bash
./zig-out/bin/zreduce run test.pdb -o test_out.pdb
```

Verify output has ATOM lines with H atoms added, passthrough records preserved.

- [ ] **Step 8: Commit**

```bash
git add src/run.zig src/main.zig
git commit -m "feat: integrate PDB format into run pipeline with format auto-detection"
```

---

### Task 5: Gzip Support for PDB

Ensure `.pdb.gz` works for both input and output.

**Files:**
- Modify: `src/run.zig` (already handles .gz in readFile; need to handle .pdb.gz output)

- [ ] **Step 1: Verify input .pdb.gz already works**

The existing `readFile()` in `run.zig` already decompresses `.gz` files. The format detection `detectFormat()` already matches `.pdb.gz`. This should work out of the box.

Write a test to confirm:

Create a gzipped PDB in a test, parse it, verify atoms parsed correctly.

- [ ] **Step 2: Verify output .pdb.gz works**

The output path branching in step 4.2 should already use `GzipWriter` when the output path ends in `.gz`. Verify the PDB writer path handles this correctly — the writer accepts an `anytype` writer, so GzipWriter's AnyWriter should work.

- [ ] **Step 3: Run full test suite**

Run: `zig build test --summary all 2>&1 | tail -10`
Expected: All tests pass.

- [ ] **Step 4: Commit (if any changes were needed)**

```bash
git add -A
git commit -m "feat: verify gzip support for PDB format I/O"
```

---

### Task 6: mmCIF Regression Verification

Verify the mmCIF pipeline is completely unaffected by the PDB changes.

**Files:** None (testing only)

- [ ] **Step 1: Run full test suite**

Run: `zig build test --summary all 2>&1 | tail -10`
Expected: All 253+ existing tests pass. No regressions.

- [ ] **Step 2: Build optimized binary**

Run: `zig build -Doptimize=ReleaseFast`

- [ ] **Step 3: Benchmark small structure (mmCIF)**

Run zreduce on a small AF model (309-residue) and verify timing is in the ~0.03s range:

```bash
time ./zig-out/bin/zreduce run examples/data/<small_af_model>.cif -o /tmp/test_out.cif
```

Compare output with a pre-change baseline if available.

- [ ] **Step 4: Benchmark large structure (mmCIF)**

Run zreduce on a larger structure (2339-residue) and verify timing is in the ~1.0s range:

```bash
time ./zig-out/bin/zreduce run examples/data/<large_model>.cif -o /tmp/test_large_out.cif
```

- [ ] **Step 5: Verify mmCIF output is identical**

Diff the mmCIF output against a known-good baseline to confirm byte-identical results:

```bash
diff /tmp/test_out.cif /tmp/baseline_out.cif
```

If no baseline exists, at minimum verify the output parses correctly and has the expected number of atoms.

- [ ] **Step 6: Document results**

Note any timing differences. If all is well, no commit needed — this is verification only.

---

### Task 7: Integration Test with Real PDB

End-to-end test with a real-world PDB file to validate the full pipeline.

**Files:** None (testing only)

- [ ] **Step 1: Download or locate a real PDB file**

Use a well-known small protein structure (e.g., 1CRN, crambin). If not available locally, download from RCSB:

```bash
curl -o /tmp/1crn.pdb "https://files.rcsb.org/download/1CRN.pdb"
```

- [ ] **Step 2: Run zreduce on the PDB file**

```bash
./zig-out/bin/zreduce run /tmp/1crn.pdb -o /tmp/1crn_h.pdb
```

Verify:
- No errors/crashes
- Output file is valid PDB format
- H atoms were added (count ATOM lines with element H)
- Original non-ATOM records are preserved (HEADER, REMARK, etc.)
- TER records separate chains
- Serial numbers are sequential

- [ ] **Step 3: Compare H count with mmCIF path**

If an mmCIF version of the same structure is available, compare the number of placed H atoms between PDB and mmCIF paths. They should be identical (same placement engine).

- [ ] **Step 4: Test with CCD dictionary**

```bash
./zig-out/bin/zreduce run /tmp/1crn.pdb -d /path/to/components.cif -o /tmp/1crn_h_ccd.pdb
```

Verify it works with external CCD.

- [ ] **Step 5: Create PR**

After all tests pass:

```bash
git push -u origin feature/pdb-format-support
gh pr create --title "feat: add PDB format input/output support" --body "$(cat <<'EOF'
## Summary
- Add PDB format parser (ATOM/HETATM fixed-width, MODEL 1 only, TER/chain detection)
- Add PDB writer (passthrough non-atom records, serial renumbering, TER generation)
- Auto-detect format by file extension (.pdb/.pdb.gz/.ent/.ent.gz)
- Entity type heuristic for PDB (polymer/non-polymer/water from comp_id)
- Chain break detection from sequence number gaps
- gzip support for .pdb.gz input/output
- All existing mmCIF tests pass unchanged

## Test plan
- [ ] `zig build test --summary all` passes
- [ ] Small AF model benchmark unchanged (~0.03s)
- [ ] Large model benchmark unchanged (~1.0s)
- [ ] Real PDB file (1CRN) processes correctly
- [ ] PDB output preserves passthrough records
- [ ] H atom count matches mmCIF path for same structure
- [ ] .pdb.gz round-trip works
EOF
)"
```
