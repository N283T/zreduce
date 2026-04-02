# compile-dict Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add a `compile-dict` subcommand that pre-compiles CCD `components.cif` into a fast-loading binary format, and auto-detect binary vs CIF in existing `-d` flag.

**Architecture:** New `src/ccd_binary.zig` module handles binary serialization/deserialization of `ComponentDict`. `main.zig` gains a `compile-dict` subcommand and auto-detection logic. `batch.zig` gains the same auto-detection. No changes to `ccd.zig` or placement logic.

**Tech Stack:** Zig, `std.io.FixedBufferStream` for reader/writer abstraction.

**Spec:** `docs/superpowers/specs/2026-04-02-compile-dict-design.md`

---

## File Map

| File | Action | Responsibility |
|------|--------|---------------|
| `src/ccd_binary.zig` | Create | Binary format write/read/detect |
| `src/root.zig` | Modify (line 9) | Re-export `ccd_binary` |
| `src/main.zig` | Modify (lines 181-264) | `compile-dict` subcommand + auto-detect in `runSubcommand` |
| `src/batch.zig` | Modify (lines 506-512) | Auto-detect in `run()` |

---

### Task 1: `isBinaryDict` — magic byte detection

**Files:**
- Create: `src/ccd_binary.zig`

- [ ] **Step 1: Write the failing tests**

Create `src/ccd_binary.zig` with tests only:

```zig
const std = @import("std");
const ccd = @import("ccd.zig");

pub const MAGIC = [4]u8{ 'Z', 'R', 'D', 'C' };
pub const FORMAT_VERSION: u8 = 1;
pub const HEADER_SIZE: usize = 12; // 4 magic + 1 version + 3 reserved + 4 count

pub fn isBinaryDict(data: []const u8) bool {
    _ = data;
    return false; // stub
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

const testing = std.testing;

test "isBinaryDict: valid magic" {
    const header = [_]u8{ 'Z', 'R', 'D', 'C', 1, 0, 0, 0, 0, 0, 0, 0 };
    try testing.expect(isBinaryDict(&header));
}

test "isBinaryDict: wrong magic" {
    const header = [_]u8{ 'X', 'Y', 'Z', 'W', 1, 0, 0, 0, 0, 0, 0, 0 };
    try testing.expect(!isBinaryDict(&header));
}

test "isBinaryDict: empty data" {
    const data = [_]u8{};
    try testing.expect(!isBinaryDict(&data));
}

test "isBinaryDict: too short" {
    const data = [_]u8{ 'Z', 'R' };
    try testing.expect(!isBinaryDict(&data));
}
```

- [ ] **Step 2: Register module in root.zig**

In `src/root.zig`, add after the `ccd` line (line 9):

```zig
pub const ccd_binary = @import("ccd_binary.zig");
```

- [ ] **Step 3: Run tests to verify they fail**

Run: `zig build test --summary all 2>&1 | grep -E 'ccd_binary|FAIL|PASS'`
Expected: `isBinaryDict: valid magic` fails (stub returns false).

- [ ] **Step 4: Implement `isBinaryDict`**

Replace the stub in `src/ccd_binary.zig`:

```zig
pub fn isBinaryDict(data: []const u8) bool {
    if (data.len < HEADER_SIZE) return false;
    return std.mem.eql(u8, data[0..4], &MAGIC);
}
```

- [ ] **Step 5: Run tests to verify they pass**

Run: `zig build test --summary all 2>&1 | grep -E 'ccd_binary|FAIL|PASS'`
Expected: All 4 `isBinaryDict` tests PASS.

- [ ] **Step 6: Commit**

```bash
git add src/ccd_binary.zig src/root.zig
git commit -m "feat: add ccd_binary module with isBinaryDict magic byte detection"
```

---

### Task 2: `writeDict` — binary serialization

**Files:**
- Modify: `src/ccd_binary.zig`

- [ ] **Step 1: Write the failing round-trip test**

Add to `src/ccd_binary.zig`:

```zig
test "writeDict + readDict: round-trip" {
    const allocator = testing.allocator;

    // Build a small ComponentDict by hand
    var dict = ccd.ComponentDict{
        .components = std.StringHashMap(ccd.Component).init(allocator),
        .allocator = allocator,
    };
    defer dict.deinit();

    // Component "TST" with 2 atoms and 1 bond
    const comp_id = try allocator.dupe(u8, "TST");
    const comp_type = try allocator.dupe(u8, "non-polymer");
    var atoms = try allocator.alloc(ccd.CompAtom, 2);
    atoms[0] = ccd.CompAtom{
        .name = [4]u8{ ' ', 'C', ' ', ' ' },
        .name_len = 1,
        .element_symbol = [2]u8{ 'C', ' ' },
        .charge = 0,
        .leaving = false,
        .aromatic = true,
        .ideal_x = 1.5,
        .ideal_y = -2.3,
        .ideal_z = 0.0,
    };
    atoms[1] = ccd.CompAtom{
        .name = [4]u8{ ' ', 'H', '1', ' ' },
        .name_len = 2,
        .element_symbol = [2]u8{ 'H', ' ' },
        .charge = -1,
        .leaving = true,
        .aromatic = false,
        .ideal_x = 3.0,
        .ideal_y = 0.1,
        .ideal_z = -1.0,
    };
    var bonds = try allocator.alloc(ccd.CompBond, 1);
    bonds[0] = ccd.CompBond{
        .atom_idx_1 = 0,
        .atom_idx_2 = 1,
        .order = .single,
        .aromatic = false,
    };

    const key = try allocator.dupe(u8, "TST");
    try dict.components.put(key, ccd.Component{
        .comp_id = comp_id,
        .comp_type = comp_type,
        .atoms = atoms,
        .bonds = bonds,
    });

    // Serialize
    var buf: [4096]u8 = undefined;
    var fbs = std.io.fixedBufferStream(&buf);
    try writeDict(fbs.writer(), &dict);

    // Verify header magic
    try testing.expectEqualSlices(u8, &MAGIC, buf[0..4]);
    try testing.expectEqual(@as(u8, FORMAT_VERSION), buf[4]);
}

test "writeDict: empty dictionary" {
    const allocator = testing.allocator;
    var dict = ccd.ComponentDict{
        .components = std.StringHashMap(ccd.Component).init(allocator),
        .allocator = allocator,
    };
    defer dict.deinit();

    var buf: [64]u8 = undefined;
    var fbs = std.io.fixedBufferStream(&buf);
    try writeDict(fbs.writer(), &dict);

    // Header: magic(4) + version(1) + reserved(3) + count(4) = 12 bytes
    try testing.expectEqual(@as(usize, HEADER_SIZE), fbs.pos);
    // Count should be 0
    try testing.expectEqual(@as(u32, 0), std.mem.readInt(u32, buf[8..12], .little));
}
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `zig build test --summary all 2>&1 | grep -E 'round-trip|empty dict|FAIL'`
Expected: FAIL — `writeDict` not defined yet.

- [ ] **Step 3: Implement `writeDict`**

Add to `src/ccd_binary.zig` before the tests section:

```zig
const Allocator = std.mem.Allocator;
const ComponentDict = ccd.ComponentDict;
const Component = ccd.Component;
const CompAtom = ccd.CompAtom;
const CompBond = ccd.CompBond;
const BondOrder = ccd.BondOrder;

// ---------------------------------------------------------------------------
// Binary serialization
// ---------------------------------------------------------------------------

/// Packed CompAtom for binary format (24 bytes).
const PackedAtom = extern struct {
    name: [4]u8,
    name_len: u8,
    element_symbol: [2]u8,
    charge: i8,
    flags: u8, // bit0=leaving, bit1=aromatic
    _pad: [3]u8 = .{ 0, 0, 0 },
    ideal_x: f32,
    ideal_y: f32,
    ideal_z: f32,

    comptime {
        std.debug.assert(@sizeOf(PackedAtom) == 24);
    }

    fn fromCompAtom(a: CompAtom) PackedAtom {
        return .{
            .name = a.name,
            .name_len = @as(u8, a.name_len),
            .element_symbol = a.element_symbol,
            .charge = a.charge,
            .flags = (@as(u8, @intFromBool(a.leaving))) | (@as(u8, @intFromBool(a.aromatic)) << 1),
            .ideal_x = a.ideal_x,
            .ideal_y = a.ideal_y,
            .ideal_z = a.ideal_z,
        };
    }

    fn toCompAtom(self: PackedAtom) CompAtom {
        return .{
            .name = self.name,
            .name_len = @intCast(self.name_len & 0x0F),
            .element_symbol = self.element_symbol,
            .charge = self.charge,
            .leaving = (self.flags & 1) != 0,
            .aromatic = (self.flags & 2) != 0,
            .ideal_x = self.ideal_x,
            .ideal_y = self.ideal_y,
            .ideal_z = self.ideal_z,
        };
    }
};

/// Packed CompBond for binary format (6 bytes).
const PackedBond = extern struct {
    atom_idx_1: u16 align(1),
    atom_idx_2: u16 align(1),
    order: u8,
    flags: u8, // bit0=aromatic

    comptime {
        std.debug.assert(@sizeOf(PackedBond) == 6);
    }

    fn fromCompBond(b: CompBond) PackedBond {
        return .{
            .atom_idx_1 = b.atom_idx_1,
            .atom_idx_2 = b.atom_idx_2,
            .order = @intFromEnum(b.order),
            .flags = @intFromBool(b.aromatic),
        };
    }

    fn toCompBond(self: PackedBond) CompBond {
        return .{
            .atom_idx_1 = self.atom_idx_1,
            .atom_idx_2 = self.atom_idx_2,
            .order = @enumFromInt(self.order),
            .aromatic = (self.flags & 1) != 0,
        };
    }
};

/// Write ComponentDict in binary format.
pub fn writeDict(writer: anytype, dict: *const ComponentDict) !void {
    // Header
    try writer.writeAll(&MAGIC);
    try writer.writeByte(FORMAT_VERSION);
    try writer.writeAll(&[3]u8{ 0, 0, 0 }); // reserved
    try writer.writeInt(u32, @intCast(dict.components.count()), .little);

    // Components (iteration order is non-deterministic but that's fine)
    var iter = dict.components.iterator();
    while (iter.next()) |entry| {
        const comp = entry.value_ptr.*;
        // comp_id
        try writer.writeByte(@intCast(comp.comp_id.len));
        try writer.writeAll(comp.comp_id);
        // comp_type
        try writer.writeByte(@intCast(comp.comp_type.len));
        try writer.writeAll(comp.comp_type);
        // atoms
        try writer.writeInt(u16, @intCast(comp.atoms.len), .little);
        for (comp.atoms) |atom| {
            const packed = PackedAtom.fromCompAtom(atom);
            try writer.writeAll(std.mem.asBytes(&packed));
        }
        // bonds
        try writer.writeInt(u16, @intCast(comp.bonds.len), .little);
        for (comp.bonds) |bond| {
            const packed = PackedBond.fromCompBond(bond);
            try writer.writeAll(std.mem.asBytes(&packed));
        }
    }
}
```

- [ ] **Step 4: Run tests to verify they pass**

Run: `zig build test --summary all 2>&1 | grep -E 'round-trip|empty dict|FAIL|PASS'`
Expected: `writeDict + readDict: round-trip` still fails (readDict not yet), but `writeDict: empty dictionary` PASSES.

- [ ] **Step 5: Commit**

```bash
git add src/ccd_binary.zig
git commit -m "feat: implement writeDict for binary CCD serialization"
```

---

### Task 3: `readDict` — binary deserialization

**Files:**
- Modify: `src/ccd_binary.zig`

- [ ] **Step 1: Write additional failing tests**

Add to `src/ccd_binary.zig`:

```zig
test "readDict: version mismatch" {
    const allocator = testing.allocator;
    var buf = [_]u8{ 'Z', 'R', 'D', 'C', 99, 0, 0, 0, 0, 0, 0, 0 };
    var fbs = std.io.fixedBufferStream(@as([]const u8, &buf));
    const result = readDict(allocator, fbs.reader());
    try testing.expectError(error.UnsupportedVersion, result);
}

test "readDict: invalid magic" {
    const allocator = testing.allocator;
    var buf = [_]u8{ 'X', 'X', 'X', 'X', 1, 0, 0, 0, 0, 0, 0, 0 };
    var fbs = std.io.fixedBufferStream(@as([]const u8, &buf));
    const result = readDict(allocator, fbs.reader());
    try testing.expectError(error.InvalidMagic, result);
}

test "readDict: truncated data" {
    const allocator = testing.allocator;
    // Valid header claiming 1 component, but no component data follows
    var buf: [HEADER_SIZE]u8 = undefined;
    @memcpy(buf[0..4], &MAGIC);
    buf[4] = FORMAT_VERSION;
    buf[5] = 0;
    buf[6] = 0;
    buf[7] = 0;
    std.mem.writeInt(u32, buf[8..12], 1, .little);
    var fbs = std.io.fixedBufferStream(@as([]const u8, &buf));
    const result = readDict(allocator, fbs.reader());
    try testing.expectError(error.UnexpectedEof, result);
}
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `zig build test --summary all 2>&1 | grep -E 'version mismatch|invalid magic|truncated|FAIL'`
Expected: FAIL — `readDict` not defined yet.

- [ ] **Step 3: Implement `readDict`**

Add to `src/ccd_binary.zig` after `writeDict`:

```zig
const ReadError = error{
    InvalidMagic,
    UnsupportedVersion,
    UnexpectedEof,
    OutOfMemory,
};

/// Read binary format into a ComponentDict.
pub fn readDict(allocator: Allocator, reader: anytype) ReadError!ComponentDict {
    // Header
    var header: [HEADER_SIZE]u8 = undefined;
    reader.readNoEof(&header) catch return error.UnexpectedEof;

    if (!std.mem.eql(u8, header[0..4], &MAGIC)) return error.InvalidMagic;
    if (header[4] != FORMAT_VERSION) return error.UnsupportedVersion;

    const count = std.mem.readInt(u32, header[8..12], .little);

    var dict = ComponentDict{
        .components = std.StringHashMap(Component).init(allocator),
        .allocator = allocator,
    };
    errdefer dict.deinit();

    try dict.components.ensureTotalCapacity(count);

    var i: u32 = 0;
    while (i < count) : (i += 1) {
        // comp_id
        const id_len = reader.readByte() catch return error.UnexpectedEof;
        const comp_id = allocator.alloc(u8, id_len) catch return error.OutOfMemory;
        errdefer allocator.free(comp_id);
        reader.readNoEof(comp_id) catch return error.UnexpectedEof;

        // comp_type
        const type_len = reader.readByte() catch return error.UnexpectedEof;
        const comp_type = allocator.alloc(u8, type_len) catch return error.OutOfMemory;
        errdefer allocator.free(comp_type);
        reader.readNoEof(comp_type) catch return error.UnexpectedEof;

        // atoms
        var atom_count_buf: [2]u8 = undefined;
        reader.readNoEof(&atom_count_buf) catch return error.UnexpectedEof;
        const atom_count = std.mem.readInt(u16, &atom_count_buf, .little);
        const atoms = allocator.alloc(CompAtom, atom_count) catch return error.OutOfMemory;
        errdefer allocator.free(atoms);
        for (atoms) |*atom| {
            var packed: PackedAtom = undefined;
            reader.readNoEof(std.mem.asBytes(&packed)) catch return error.UnexpectedEof;
            atom.* = packed.toCompAtom();
        }

        // bonds
        var bond_count_buf: [2]u8 = undefined;
        reader.readNoEof(&bond_count_buf) catch return error.UnexpectedEof;
        const bond_count = std.mem.readInt(u16, &bond_count_buf, .little);
        const bonds = allocator.alloc(CompBond, bond_count) catch return error.OutOfMemory;
        errdefer allocator.free(bonds);
        for (bonds) |*bond| {
            var packed: PackedBond = undefined;
            reader.readNoEof(std.mem.asBytes(&packed)) catch return error.UnexpectedEof;
            bond.* = packed.toCompBond();
        }

        // HashMap key
        const key = allocator.dupe(u8, comp_id) catch return error.OutOfMemory;
        errdefer allocator.free(key);

        dict.components.putAssumeCapacity(key, Component{
            .comp_id = comp_id,
            .comp_type = comp_type,
            .atoms = atoms,
            .bonds = bonds,
        });
    }

    return dict;
}
```

- [ ] **Step 4: Run tests to verify all pass**

Run: `zig build test --summary all 2>&1 | grep -E 'ccd_binary|round-trip|FAIL|PASS'`
Expected: ALL tests pass, including the round-trip test from Task 2.

- [ ] **Step 5: Enhance round-trip test to verify field values**

The round-trip test from Task 2 only checks the header. Add field verification to the end of the `"writeDict + readDict: round-trip"` test:

```zig
    // Deserialize
    var read_fbs = std.io.fixedBufferStream(buf[0..fbs.pos]);
    var result = try readDict(allocator, read_fbs.reader());
    defer result.deinit();

    try testing.expectEqual(@as(u32, 1), result.components.count());
    const comp = result.get("TST").?;
    try testing.expectEqualStrings("TST", comp.comp_id);
    try testing.expectEqualStrings("non-polymer", comp.comp_type);
    try testing.expectEqual(@as(usize, 2), comp.atoms.len);
    try testing.expectEqual(@as(usize, 1), comp.bonds.len);

    // Verify atom fields
    const a0 = comp.atoms[0];
    try testing.expectEqualSlices(u8, "C", a0.nameSlice());
    try testing.expectEqual([2]u8{ 'C', ' ' }, a0.element_symbol);
    try testing.expectEqual(@as(i8, 0), a0.charge);
    try testing.expect(!a0.leaving);
    try testing.expect(a0.aromatic);
    try testing.expectApproxEqAbs(@as(f32, 1.5), a0.ideal_x, 1e-6);
    try testing.expectApproxEqAbs(@as(f32, -2.3), a0.ideal_y, 1e-6);

    const a1 = comp.atoms[1];
    try testing.expectEqual(@as(i8, -1), a1.charge);
    try testing.expect(a1.leaving);
    try testing.expect(!a1.aromatic);

    // Verify bond fields
    const b0 = comp.bonds[0];
    try testing.expectEqual(@as(u16, 0), b0.atom_idx_1);
    try testing.expectEqual(@as(u16, 1), b0.atom_idx_2);
    try testing.expectEqual(ccd.BondOrder.single, b0.order);
    try testing.expect(!b0.aromatic);
```

- [ ] **Step 6: Run tests to verify all pass**

Run: `zig build test --summary all 2>&1 | grep -E 'round-trip|FAIL|PASS'`
Expected: PASS

- [ ] **Step 7: Commit**

```bash
git add src/ccd_binary.zig
git commit -m "feat: implement readDict for binary CCD deserialization"
```

---

### Task 4: Auto-detect in `main.zig` and `batch.zig`

**Files:**
- Modify: `src/main.zig` (lines 457-468 in `runSubcommand`)
- Modify: `src/batch.zig` (lines 506-512 in `run()`)

- [ ] **Step 1: Extract dict loading helper**

The dict loading logic is duplicated in `main.zig:runSubcommand` and `batch.zig:run`. Instead of duplicating the auto-detect logic, add a helper in `ccd_binary.zig`:

```zig
/// Load a ComponentDict from either binary (.zdict) or CIF text format.
/// Auto-detects based on magic bytes.
pub fn loadDict(allocator: Allocator, data: []const u8) !ComponentDict {
    if (isBinaryDict(data)) {
        var fbs = std.io.fixedBufferStream(data);
        return readDict(allocator, fbs.reader());
    }
    return ccd.parseComponentDict(allocator, data);
}
```

- [ ] **Step 2: Add test for `loadDict` with CIF input**

```zig
test "loadDict: CIF text input" {
    const allocator = testing.allocator;
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
    ;
    var result = try loadDict(allocator, source);
    defer result.deinit();
    try testing.expect(result.get("ALA") != null);
}

test "loadDict: binary input round-trip" {
    const allocator = testing.allocator;

    // First create a binary dict
    var dict = ccd.ComponentDict{
        .components = std.StringHashMap(ccd.Component).init(allocator),
        .allocator = allocator,
    };
    defer dict.deinit();

    const comp_id = try allocator.dupe(u8, "XX");
    const comp_type = try allocator.dupe(u8, "other");
    const atoms = try allocator.alloc(ccd.CompAtom, 0);
    const bonds = try allocator.alloc(ccd.CompBond, 0);
    const key = try allocator.dupe(u8, "XX");
    try dict.components.put(key, ccd.Component{
        .comp_id = comp_id,
        .comp_type = comp_type,
        .atoms = atoms,
        .bonds = bonds,
    });

    var buf: [4096]u8 = undefined;
    var fbs = std.io.fixedBufferStream(&buf);
    try writeDict(fbs.writer(), &dict);

    // Load via auto-detect
    var result = try loadDict(allocator, buf[0..fbs.pos]);
    defer result.deinit();
    try testing.expect(result.get("XX") != null);
}
```

- [ ] **Step 3: Run tests**

Run: `zig build test --summary all 2>&1 | grep -E 'loadDict|FAIL|PASS'`
Expected: PASS

- [ ] **Step 4: Update `main.zig` `runSubcommand`**

Replace lines 458-468 in `src/main.zig`:

```zig
    // Load CCD dictionary (once, before processFile)
    var ccd_dict: ?zreduce.ccd.ComponentDict = null;
    if (config.dict_path) |dict_path| {
        const dict_source = zreduce.run.readFile(allocator, dict_path) catch |err| {
            std.debug.print("Error: cannot read dictionary '{s}': {s}\n", .{ dict_path, @errorName(err) });
            std.process.exit(1);
        };
        defer allocator.free(dict_source);
        ccd_dict = zreduce.ccd_binary.loadDict(allocator, dict_source) catch |err| {
            std.debug.print("Error: failed to load CCD dictionary: {s}\n", .{@errorName(err)});
            std.process.exit(1);
        };
    }
```

- [ ] **Step 5: Update `batch.zig` `run()`**

Replace lines 508-511 in `src/batch.zig`:

```zig
    if (config.dict_path) |dict_path| {
        const dict_source = try run_mod.readFile(allocator, dict_path);
        defer allocator.free(dict_source);
        ccd_dict = try zreduce.ccd_binary.loadDict(allocator, dict_source);
    }
```

- [ ] **Step 6: Build and run full test suite**

Run: `zig build test --summary all`
Expected: All tests pass (existing + new).

- [ ] **Step 7: Commit**

```bash
git add src/ccd_binary.zig src/main.zig src/batch.zig
git commit -m "feat: auto-detect binary vs CIF dictionary format in run and batch"
```

---

### Task 5: `compile-dict` subcommand

**Files:**
- Modify: `src/main.zig`

- [ ] **Step 1: Add `compile-dict` subcommand dispatch**

In `src/main.zig`, in the `main()` function, add the new subcommand branch after the `"batch"` check (around line 258):

```zig
    } else if (std.mem.eql(u8, subcmd, "compile-dict")) {
        compileDictSubcommand(allocator, args[2..]);
    } else {
```

- [ ] **Step 2: Add argument parser and handler**

Add these functions to `src/main.zig`:

```zig
const CompileDictConfig = struct {
    input_path: []const u8,
    output_path: ?[]const u8 = null,
};

fn parseCompileDictArgs(args: []const []const u8) ?CompileDictConfig {
    var config = CompileDictConfig{ .input_path = undefined };
    var input_set = false;
    var i: usize = 0;

    while (i < args.len) : (i += 1) {
        const arg = args[i];
        if (std.mem.eql(u8, arg, "-h") or std.mem.eql(u8, arg, "--help")) {
            printCompileDictUsage();
            return null;
        } else if (std.mem.eql(u8, arg, "-o") or std.mem.eql(u8, arg, "--output")) {
            i += 1;
            if (i >= args.len) {
                std.debug.print("Error: {s} requires a path argument\n", .{arg});
                std.process.exit(1);
            }
            config.output_path = args[i];
        } else if (arg.len > 0 and arg[0] == '-') {
            std.debug.print("Error: unknown option '{s}'\n", .{arg});
            std.process.exit(1);
        } else {
            if (input_set) {
                std.debug.print("Error: unexpected argument '{s}'\n", .{arg});
                std.process.exit(1);
            }
            config.input_path = arg;
            input_set = true;
        }
    }

    if (!input_set) {
        std.debug.print("Error: missing input CIF path\n", .{});
        printCompileDictUsage();
        std.process.exit(1);
    }

    if (config.output_path == null) {
        std.debug.print("Error: -o/--output is required\n", .{});
        std.process.exit(1);
    }

    return config;
}

fn printCompileDictUsage() void {
    std.debug.print(
        \\USAGE:
        \\    zreduce compile-dict [OPTIONS] <input.cif>
        \\
        \\OPTIONS:
        \\    -h, --help         Show this help message
        \\    -o, --output PATH  Output binary dictionary file (required)
        \\
    , .{});
}

fn compileDictSubcommand(allocator: Allocator, args: []const []const u8) void {
    const config = parseCompileDictArgs(args) orelse return;

    // Read input
    const source = zreduce.run.readFile(allocator, config.input_path) catch |err| {
        std.debug.print("Error: cannot read '{s}': {s}\n", .{ config.input_path, @errorName(err) });
        std.process.exit(1);
    };
    defer allocator.free(source);

    // Reject if already binary
    if (zreduce.ccd_binary.isBinaryDict(source)) {
        std.debug.print("Error: input is already a compiled dictionary\n", .{});
        std.process.exit(1);
    }

    // Parse CIF
    var dict = zreduce.ccd.parseComponentDict(allocator, source) catch |err| {
        std.debug.print("Error: failed to parse CCD dictionary: {s}\n", .{@errorName(err)});
        std.process.exit(1);
    };
    defer dict.deinit();

    // Write binary
    const output_path = config.output_path.?;
    const file = std.fs.cwd().createFile(output_path, .{}) catch |err| {
        std.debug.print("Error: cannot create '{s}': {s}\n", .{ output_path, @errorName(err) });
        std.process.exit(1);
    };
    defer file.close();

    var buf_writer = std.io.bufferedWriter(file.writer());
    zreduce.ccd_binary.writeDict(buf_writer.writer(), &dict) catch |err| {
        std.debug.print("Error: failed to write binary dictionary: {s}\n", .{@errorName(err)});
        std.process.exit(1);
    };
    buf_writer.flush() catch |err| {
        std.debug.print("Error: failed to flush output: {s}\n", .{@errorName(err)});
        std.process.exit(1);
    };

    std.debug.print("Compiled {d} components to '{s}'\n", .{ dict.components.count(), output_path });
}
```

- [ ] **Step 3: Update usage text**

In `printUsage()`, add the new subcommand to the COMMANDS list:

```zig
        \\COMMANDS:
        \\    run           Process a single mmCIF file
        \\    batch         Process all mmCIF files in a directory
        \\    compile-dict  Pre-compile CCD dictionary to binary format
```

- [ ] **Step 4: Build and verify**

Run: `zig build -Doptimize=ReleaseFast`
Expected: Compiles without errors.

Run: `./zig-out/bin/zreduce compile-dict --help`
Expected: Prints compile-dict usage.

Run: `./zig-out/bin/zreduce --help`
Expected: Shows `compile-dict` in the commands list.

- [ ] **Step 5: Commit**

```bash
git add src/main.zig
git commit -m "feat: add compile-dict subcommand for pre-compiling CCD to binary"
```

---

### Task 6: Integration test — full round-trip with real CCD

**Files:** None (manual verification)

- [ ] **Step 1: Run full test suite**

Run: `zig build test --summary all`
Expected: All tests pass.

- [ ] **Step 2: Compile the real CCD**

Run: `./zig-out/bin/zreduce compile-dict /Users/nagaet/data/ccd/components.cif -o /tmp/components.zdict`
Expected: Prints component count and output path. Note the file sizes:

```bash
ls -lh /Users/nagaet/data/ccd/components.cif /tmp/components.zdict
```

- [ ] **Step 3: Verify binary dict produces identical results**

Pick a test structure that uses CCD (non-standard residues):

```bash
# With CIF dictionary
./zig-out/bin/zreduce run examples/data/fold_test2.cif -d /Users/nagaet/data/ccd/components.cif -o /tmp/out_cif.cif

# With binary dictionary
./zig-out/bin/zreduce run examples/data/fold_test2.cif -d /tmp/components.zdict -o /tmp/out_zdict.cif

# Compare
diff /tmp/out_cif.cif /tmp/out_zdict.cif
```

Expected: No differences.

- [ ] **Step 4: Verify binary-already-compiled rejection**

```bash
./zig-out/bin/zreduce compile-dict /tmp/components.zdict -o /tmp/test.zdict
```

Expected: `Error: input is already a compiled dictionary`

- [ ] **Step 5: Commit any final fixes if needed**

```bash
git add -A
git commit -m "test: verify compile-dict integration with real CCD data"
```

(Skip commit if no changes were needed.)
