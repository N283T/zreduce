# Bond Graph Parsing Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Parse `_struct_conn`, `_chem_comp_bond`/`_chem_comp_atom` (inline), and `_pdbx_entity_branch_link` from mmCIF files to build a bond graph and skip H placement on leaving atoms.

**Architecture:** New parsing functions in `mmcif.zig` extract bond data from three mmCIF categories. An atom lookup map resolves atom identities to Model indices. A `bonded_inter_residue` flag on `AtomFlags` tells the placer to skip H placement on bonded atoms. Inline `_chem_comp_bond`/`_chem_comp_atom` produces a `ComponentDict` (same type as CCD) that takes priority over external CCD dictionary.

**Tech Stack:** Zig, existing CIF tokenizer/parser

---

## File Map

| File | Action | Responsibility |
|------|--------|----------------|
| `src/element.zig` | Modify | Add `bonded_inter_residue` bit to `AtomFlags` |
| `src/mmcif.zig` | Modify | Add `parseStructConn`, `parseBranchLinks`, `parseInlineComponents`, atom lookup builder |
| `src/model/chain.zig` | Modify | Store `entity_id` as string (currently `u16`, needs change for branch_link entity mapping) |
| `src/place/placer.zig` | Modify | Check `bonded_inter_residue` flag, accept inline dict priority |
| `src/run.zig` | Modify | Wire new parsing steps into pipeline |
| `src/batch.zig` | Modify | Wire new parsing steps (same as run.zig) |
| `src/test_data/disulfide.cif` | Create | Test fixture: two CYS residues with `_struct_conn` disulfide |
| `src/test_data/branch_link.cif` | Create | Test fixture: sugar with `_pdbx_entity_branch_link` |
| `src/test_data/inline_comp.cif` | Create | Test fixture: structure with inline `_chem_comp_atom`/`_chem_comp_bond` |

---

### Task 1: Add `bonded_inter_residue` flag to AtomFlags

**Files:**
- Modify: `src/element.zig:5-14`

- [ ] **Step 1: Write the failing test**

Add to `src/element.zig` test section:

```zig
test "bonded_inter_residue flag" {
    var flags = AtomFlags{};
    try std.testing.expect(!flags.bonded_inter_residue);
    flags.bonded_inter_residue = true;
    try std.testing.expect(flags.bonded_inter_residue);

    // mergeFlags preserves bonded_inter_residue
    const a = AtomFlags{ .donor = true };
    const b = AtomFlags{ .bonded_inter_residue = true };
    const merged = mergeFlags(a, b);
    try std.testing.expect(merged.donor);
    try std.testing.expect(merged.bonded_inter_residue);
}
```

- [ ] **Step 2: Run test to verify it fails**

Run: `zig build test --summary all 2>&1 | grep "bonded_inter_residue"`
Expected: compilation error — `bonded_inter_residue` not a field of `AtomFlags`

- [ ] **Step 3: Add `bonded_inter_residue` field to `AtomFlags`**

In `src/element.zig`, modify `AtomFlags`:

```zig
pub const AtomFlags = packed struct {
    donor: bool = false,
    acceptor: bool = false,
    aromatic: bool = false,
    positive: bool = false,
    negative: bool = false,
    metallic: bool = false,
    hb_only_dummy: bool = false,
    bonded_inter_residue: bool = false,
};
```

Remove the `_padding: u1 = 0` field and replace with `bonded_inter_residue: bool = false`. The packed struct remains 8 bits (1 byte) since we had a padding bit.

- [ ] **Step 4: Run test to verify it passes**

Run: `zig build test --summary all 2>&1 | tail -5`
Expected: all tests pass

- [ ] **Step 5: Commit**

```bash
git add src/element.zig
git commit -m "feat: add bonded_inter_residue flag to AtomFlags"
```

---

### Task 2: Store entity_id as string in Chain

**Files:**
- Modify: `src/model/chain.zig`
- Modify: `src/mmcif.zig` (entity_id parsing from `_atom_site.label_entity_id`)

The current `Chain.entity_id` is a `u16`. For `_pdbx_entity_branch_link`, we need to match entity_id strings from the loop to chain entity_ids. Change to a fixed-size string (like asym_id).

- [ ] **Step 1: Write the failing test**

Add to `src/model/chain.zig`:

```zig
test "Chain entity_id string" {
    var c = Chain{};
    c.setEntityId("2");
    try std.testing.expectEqualStrings("2", c.entityIdSlice());
}
```

- [ ] **Step 2: Run test to verify it fails**

Run: `zig build test --summary all 2>&1 | grep "entity_id"`
Expected: compilation error — `setEntityId` not found

- [ ] **Step 3: Implement entity_id as fixed-size string**

In `src/model/chain.zig`, replace `entity_id: u16 = 0` with:

```zig
entity_id: [4]u8 = .{ ' ', ' ', ' ', ' ' },
entity_id_len: u4 = 0,

pub fn entityIdSlice(self: *const Chain) []const u8 {
    return self.entity_id[0..@min(@as(usize, self.entity_id_len), 4)];
}

pub fn setEntityId(self: *Chain, id: []const u8) void {
    const len: u4 = @intCast(@min(id.len, 4));
    self.entity_id = .{ ' ', ' ', ' ', ' ' };
    for (0..len) |i| self.entity_id[i] = id[i];
    self.entity_id_len = len;
}
```

- [ ] **Step 4: Update mmcif.zig to parse `_atom_site.label_entity_id`**

In `src/mmcif.zig`, add to `AtomSiteColumns`:

```zig
label_entity_id: ?usize = null,
```

In the column mapping section:

```zig
cols.label_entity_id = loop.findTag("_atom_site.label_entity_id");
```

In the chain construction section, after `new_chain.setAuthAsymId(auth_asym)`:

```zig
const entity_id = if (cols.label_entity_id) |c| cif.asString(loop.val(row, c) orelse ".") else "";
new_chain.setEntityId(entity_id);
```

- [ ] **Step 5: Fix any compilation errors from old entity_id usage**

Search the codebase for `entity_id` usage and update any references from `u16` to the new string API. Currently only `chain.zig` defines `entity_id` and it's not referenced elsewhere (the field was unused).

- [ ] **Step 6: Run all tests**

Run: `zig build test --summary all`
Expected: all tests pass

- [ ] **Step 7: Commit**

```bash
git add src/model/chain.zig src/mmcif.zig
git commit -m "feat: store entity_id as string in Chain, parse from _atom_site"
```

---

### Task 3: Create test fixtures

**Files:**
- Create: `src/test_data/disulfide.cif`
- Create: `src/test_data/branch_link.cif`
- Create: `src/test_data/inline_comp.cif`

- [ ] **Step 1: Create disulfide test fixture**

Create `src/test_data/disulfide.cif` — two CYS residues in chain A with a `_struct_conn` disulfide bond between SG atoms:

```cif
data_DISULFIDE
#
loop_
_atom_site.group_PDB
_atom_site.id
_atom_site.type_symbol
_atom_site.label_atom_id
_atom_site.label_comp_id
_atom_site.label_asym_id
_atom_site.label_entity_id
_atom_site.label_seq_id
_atom_site.Cartn_x
_atom_site.Cartn_y
_atom_site.Cartn_z
_atom_site.occupancy
_atom_site.B_iso_or_equiv
_atom_site.label_alt_id
_atom_site.pdbx_PDB_ins_code
ATOM 1  N N   CYS A 1 1  1.0  0.0  0.0 1.0 0.0 . .
ATOM 2  C CA  CYS A 1 1  2.0  0.0  0.0 1.0 0.0 . .
ATOM 3  C C   CYS A 1 1  3.0  0.0  0.0 1.0 0.0 . .
ATOM 4  O O   CYS A 1 1  3.5  1.0  0.0 1.0 0.0 . .
ATOM 5  C CB  CYS A 1 1  2.0  1.0  1.0 1.0 0.0 . .
ATOM 6  S SG  CYS A 1 1  2.0  2.0  2.0 1.0 0.0 . .
ATOM 7  N N   CYS A 1 2  5.0  0.0  0.0 1.0 0.0 . .
ATOM 8  C CA  CYS A 1 2  6.0  0.0  0.0 1.0 0.0 . .
ATOM 9  C C   CYS A 1 2  7.0  0.0  0.0 1.0 0.0 . .
ATOM 10 O O   CYS A 1 2  7.5  1.0  0.0 1.0 0.0 . .
ATOM 11 C CB  CYS A 1 2  6.0  1.0  1.0 1.0 0.0 . .
ATOM 12 S SG  CYS A 1 2  4.0  2.0  2.0 1.0 0.0 . .
#
loop_
_struct_conn.id
_struct_conn.conn_type_id
_struct_conn.ptnr1_label_asym_id
_struct_conn.ptnr1_label_comp_id
_struct_conn.ptnr1_label_seq_id
_struct_conn.ptnr1_label_atom_id
_struct_conn.ptnr1_symmetry
_struct_conn.ptnr2_label_asym_id
_struct_conn.ptnr2_label_comp_id
_struct_conn.ptnr2_label_seq_id
_struct_conn.ptnr2_label_atom_id
_struct_conn.ptnr2_symmetry
_struct_conn.pdbx_value_order
disulf1 disulf A CYS 1 SG 1_555 A CYS 2 SG 1_555 SING
#
```

- [ ] **Step 2: Create branch_link test fixture**

Create `src/test_data/branch_link.cif` — two sugar residues with a branch link:

```cif
data_BRANCH
#
loop_
_atom_site.group_PDB
_atom_site.id
_atom_site.type_symbol
_atom_site.label_atom_id
_atom_site.label_comp_id
_atom_site.label_asym_id
_atom_site.label_entity_id
_atom_site.label_seq_id
_atom_site.Cartn_x
_atom_site.Cartn_y
_atom_site.Cartn_z
_atom_site.occupancy
_atom_site.B_iso_or_equiv
_atom_site.auth_seq_id
_atom_site.label_alt_id
_atom_site.pdbx_PDB_ins_code
HETATM 1  C C1  NAG B 2 . 1.0 0.0 0.0 1.0 0.0 1 . .
HETATM 2  O O4  NAG B 2 . 2.0 0.0 0.0 1.0 0.0 1 . .
HETATM 3  C C1  GAL B 2 . 3.0 0.0 0.0 1.0 0.0 2 . .
HETATM 4  O O1  GAL B 2 . 4.0 0.0 0.0 1.0 0.0 2 . .
#
loop_
_pdbx_entity_branch_link.link_id
_pdbx_entity_branch_link.entity_id
_pdbx_entity_branch_link.entity_branch_list_num_1
_pdbx_entity_branch_link.comp_id_1
_pdbx_entity_branch_link.atom_id_1
_pdbx_entity_branch_link.leaving_atom_id_1
_pdbx_entity_branch_link.entity_branch_list_num_2
_pdbx_entity_branch_link.comp_id_2
_pdbx_entity_branch_link.atom_id_2
_pdbx_entity_branch_link.leaving_atom_id_2
1 2 1 NAG O4 HO4 2 GAL C1 O1
#
```

- [ ] **Step 3: Create inline component test fixture**

Create `src/test_data/inline_comp.cif` — structure with inline `_chem_comp_atom` and `_chem_comp_bond`:

```cif
data_INLINE
#
loop_
_atom_site.group_PDB
_atom_site.id
_atom_site.type_symbol
_atom_site.label_atom_id
_atom_site.label_comp_id
_atom_site.label_asym_id
_atom_site.label_entity_id
_atom_site.label_seq_id
_atom_site.Cartn_x
_atom_site.Cartn_y
_atom_site.Cartn_z
_atom_site.occupancy
_atom_site.B_iso_or_equiv
_atom_site.label_alt_id
_atom_site.pdbx_PDB_ins_code
ATOM 1  N N   ALA A 1 1  1.458 0.000  0.000 1.0 0.0 . .
ATOM 2  C CA  ALA A 1 1  2.009 1.420  0.000 1.0 0.0 . .
ATOM 3  C C   ALA A 1 1  3.534 1.424  0.000 1.0 0.0 . .
ATOM 4  O O   ALA A 1 1  4.139 0.366 -0.002 1.0 0.0 . .
ATOM 5  C CB  ALA A 1 1  1.498 2.138  1.248 1.0 0.0 . .
#
loop_
_chem_comp_atom.comp_id
_chem_comp_atom.atom_id
_chem_comp_atom.type_symbol
_chem_comp_atom.charge
_chem_comp_atom.pdbx_leaving_atom_flag
ALA N   N 0 N
ALA CA  C 0 N
ALA C   C 0 N
ALA O   O 0 N
ALA CB  C 0 N
ALA H   H 0 N
ALA HA  H 0 N
ALA HB1 H 0 N
ALA HB2 H 0 N
ALA HB3 H 0 N
ALA OXT O 0 Y
#
loop_
_chem_comp_bond.comp_id
_chem_comp_bond.atom_id_1
_chem_comp_bond.atom_id_2
_chem_comp_bond.value_order
_chem_comp_bond.pdbx_aromatic_flag
ALA N   CA  SING N
ALA CA  C   SING N
ALA CA  CB  SING N
ALA C   O   DOUB N
ALA N   H   SING N
ALA CA  HA  SING N
ALA CB  HB1 SING N
ALA CB  HB2 SING N
ALA CB  HB3 SING N
ALA C   OXT SING N
#
```

- [ ] **Step 4: Commit fixtures**

```bash
git add src/test_data/disulfide.cif src/test_data/branch_link.cif src/test_data/inline_comp.cif
git commit -m "test: add CIF fixtures for bond graph parsing"
```

---

### Task 4: Implement atom lookup builder in mmcif.zig

**Files:**
- Modify: `src/mmcif.zig`

Build a HashMap from `(label_asym_id, label_seq_id, atom_name)` to Model atom index. This is needed by `parseStructConn` and `parseBranchLinks`. Functions accept `*const cif.Block` (not raw source) to avoid re-parsing the CIF document multiple times — `run.zig` already has the parsed doc.

- [ ] **Step 1: Write the failing test**

Add to `src/mmcif.zig` test section:

```zig
test "buildAtomLookup resolves atom indices" {
    const source = @embedFile("test_data/disulfide.cif");

    // Parse the CIF document (needed for block access)
    var doc = cif.readString(testing.allocator, source) catch unreachable;
    defer doc.deinit();
    const block = &doc.blocks.items[0];

    var mdl = try parseModel(testing.allocator, source);
    defer mdl.deinit();

    var lookup = try buildAtomLookup(testing.allocator, block);
    defer lookup.deinit();

    // SG of CYS residue 1 (seq_id=1) should be atom index 5
    const sg1 = lookup.get(.{ .label_asym_id = "A", .label_seq_id = "1", .atom_name = "SG" });
    try testing.expect(sg1 != null);
    try testing.expectEqual(@as(u32, 5), sg1.?);

    // SG of CYS residue 2 (seq_id=2) should be atom index 11
    const sg2 = lookup.get(.{ .label_asym_id = "A", .label_seq_id = "2", .atom_name = "SG" });
    try testing.expect(sg2 != null);
    try testing.expectEqual(@as(u32, 11), sg2.?);
}
```

- [ ] **Step 2: Run test to verify it fails**

Run: `zig build test --summary all 2>&1 | grep "buildAtomLookup"`
Expected: compilation error — `buildAtomLookup` not found

- [ ] **Step 3: Implement AtomLookupKey and buildAtomLookup**

Add to `src/mmcif.zig`:

```zig
/// Key for looking up atoms by identity (chain, seq_id, atom_name).
pub const AtomLookupKey = struct {
    label_asym_id: []const u8,
    label_seq_id: []const u8,
    atom_name: []const u8,
};

const AtomLookupContext = struct {
    pub fn hash(_: AtomLookupContext, key: AtomLookupKey) u64 {
        var h = std.hash.Wyhash.init(0);
        h.update(key.label_asym_id);
        h.update("\x00");
        h.update(key.label_seq_id);
        h.update("\x00");
        h.update(key.atom_name);
        return h.final();
    }

    pub fn eql(_: AtomLookupContext, a: AtomLookupKey, b: AtomLookupKey) bool {
        return std.mem.eql(u8, a.label_asym_id, b.label_asym_id) and
            std.mem.eql(u8, a.label_seq_id, b.label_seq_id) and
            std.mem.eql(u8, a.atom_name, b.atom_name);
    }
};

pub const AtomLookup = std.HashMap(AtomLookupKey, u32, AtomLookupContext, 80);

/// Build an atom lookup from a CIF block's _atom_site loop.
/// Keys use the CIF source strings (label_asym_id, label_seq_id, atom_name).
/// For altloc atoms, only the first occurrence is indexed (blank altloc preferred).
/// Also registers auth_seq_id entries for branched entities where label_seq_id is ".".
pub fn buildAtomLookup(allocator: Allocator, block: *const cif.Block) !AtomLookup {
    const loop = block.findLoop("_atom_site.Cartn_x") orelse return error.NoAtomSiteLoop;

    const col_asym = loop.findTag("_atom_site.label_asym_id") orelse return error.MissingCoordinateField;
    const col_seq = loop.findTag("_atom_site.label_seq_id") orelse return error.MissingCoordinateField;
    const col_atom = loop.findTag("_atom_site.label_atom_id") orelse return error.MissingCoordinateField;
    const col_alt = loop.findTag("_atom_site.label_alt_id");
    const col_auth_seq = loop.findTag("_atom_site.auth_seq_id");

    var lookup = AtomLookup.init(allocator);
    errdefer lookup.deinit();

    const nrows = loop.length();
    for (0..nrows) |row| {
        const asym = cif.asString(loop.val(row, col_asym) orelse continue);
        const seq = cif.asString(loop.val(row, col_seq) orelse continue);
        const atom_name = cif.asString(loop.val(row, col_atom) orelse continue);

        // For altloc: only index the first occurrence of each (asym, seq, name) tuple.
        // This means blank altloc or the first conformer seen.
        const key = AtomLookupKey{
            .label_asym_id = asym,
            .label_seq_id = seq,
            .atom_name = atom_name,
        };

        if (col_alt) |c| {
            const alt = cif.asString(loop.val(row, c) orelse ".");
            if (alt.len > 0) {
                // Non-blank altloc: only index if not already present
                if (lookup.contains(key)) continue;
            }
        }

        const atom_idx: u32 = @intCast(row);
        // Row index == Model atom index because parseModel processes all rows in order
        try lookup.put(key, atom_idx);

        // Also register auth_seq_id for branched entities (label_seq_id is "." there)
        if (col_auth_seq) |c| {
            const auth_seq = cif.asString(loop.val(row, c) orelse ".");
            if (auth_seq.len > 0) {
                const auth_key = AtomLookupKey{
                    .label_asym_id = asym,
                    .label_seq_id = auth_seq,
                    .atom_name = atom_name,
                };
                const r = try lookup.getOrPut(auth_key);
                if (!r.found_existing) {
                    r.value_ptr.* = atom_idx;
                }
            }
        }
    }

    return lookup;
}
```

**Important note:** Row index == Model atom index because `parseModel` appends atoms in CIF row order without skipping. This invariant is verified in the test.

- [ ] **Step 4: Run test to verify it passes**

Run: `zig build test --summary all 2>&1 | tail -5`
Expected: all tests pass

- [ ] **Step 5: Commit**

```bash
git add src/mmcif.zig
git commit -m "feat: add buildAtomLookup for atom identity resolution"
```

---

### Task 5: Implement `parseStructConn`

**Files:**
- Modify: `src/mmcif.zig`

- [ ] **Step 1: Write the failing test**

Add to `src/mmcif.zig`:

```zig
test "parseStructConn disulfide bond" {
    const source = @embedFile("test_data/disulfide.cif");
    var doc = cif.readString(testing.allocator, source) catch unreachable;
    defer doc.deinit();
    const block = &doc.blocks.items[0];

    var mdl = try parseModel(testing.allocator, source);
    defer mdl.deinit();

    try parseStructConn(testing.allocator, &mdl, block);

    // Should have 1 bond (disulfide SG-SG)
    try testing.expectEqual(@as(usize, 1), mdl.bonds.items.len);
    const bond = mdl.bonds.items[0];
    try testing.expectEqual(@import("model/bond.zig").BondSource.struct_conn, bond.source);

    // Both SG atoms should have bonded_inter_residue flag
    try testing.expect(mdl.atoms.items[bond.atom_1].flags.bonded_inter_residue);
    try testing.expect(mdl.atoms.items[bond.atom_2].flags.bonded_inter_residue);
}
```

- [ ] **Step 2: Run test to verify it fails**

Run: `zig build test --summary all 2>&1 | grep "parseStructConn"`
Expected: compilation error — `parseStructConn` not found

- [ ] **Step 3: Implement parseStructConn**

Add to `src/mmcif.zig`:

```zig
const bond_mod = @import("model/bond.zig");

/// Parse _struct_conn loop and add inter-residue bonds to the Model.
/// Sets `bonded_inter_residue` flag on both partner atoms.
/// Accepts a pre-parsed CIF block to avoid re-parsing the document.
pub fn parseStructConn(allocator: Allocator, mdl: *Model, block: *const cif.Block) !void {
    const sc = block.findLoop("_struct_conn.conn_type_id") orelse return;

    const col_type = sc.findTag("_struct_conn.conn_type_id") orelse return;
    const col_asym1 = sc.findTag("_struct_conn.ptnr1_label_asym_id") orelse return;
    const col_seq1 = sc.findTag("_struct_conn.ptnr1_label_seq_id") orelse return;
    const col_atom1 = sc.findTag("_struct_conn.ptnr1_label_atom_id") orelse return;
    const col_asym2 = sc.findTag("_struct_conn.ptnr2_label_asym_id") orelse return;
    const col_seq2 = sc.findTag("_struct_conn.ptnr2_label_seq_id") orelse return;
    const col_atom2 = sc.findTag("_struct_conn.ptnr2_label_atom_id") orelse return;
    const col_sym1 = sc.findTag("_struct_conn.ptnr1_symmetry");
    const col_sym2 = sc.findTag("_struct_conn.ptnr2_symmetry");
    const col_order = sc.findTag("_struct_conn.pdbx_value_order");

    // Build atom lookup from the same block
    var lookup = try buildAtomLookup(allocator, block);
    defer lookup.deinit();

    for (0..sc.length()) |row| {
        const conn_type = cif.asString(sc.val(row, col_type) orelse continue);

        // Only process covalent bond types (covale*, disulf)
        if (!isCovalentConnType(conn_type)) continue;

        // Skip inter-symmetry bonds
        if (col_sym1 != null and col_sym2 != null) {
            const sym1 = cif.asString(sc.val(row, col_sym1.?) orelse "1_555");
            const sym2 = cif.asString(sc.val(row, col_sym2.?) orelse "1_555");
            if (!std.mem.eql(u8, sym1, sym2)) continue;
        }

        const asym1 = cif.asString(sc.val(row, col_asym1) orelse continue);
        const seq1 = cif.asString(sc.val(row, col_seq1) orelse continue);
        const name1 = cif.asString(sc.val(row, col_atom1) orelse continue);
        const asym2 = cif.asString(sc.val(row, col_asym2) orelse continue);
        const seq2 = cif.asString(sc.val(row, col_seq2) orelse continue);
        const name2 = cif.asString(sc.val(row, col_atom2) orelse continue);

        const a1 = lookup.get(.{
            .label_asym_id = asym1,
            .label_seq_id = seq1,
            .atom_name = name1,
        }) orelse continue;

        const a2 = lookup.get(.{
            .label_asym_id = asym2,
            .label_seq_id = seq2,
            .atom_name = name2,
        }) orelse continue;

        const order = if (col_order) |co|
            bond_mod.BondOrder.fromString(cif.asString(sc.val(row, co) orelse "SING"))
        else
            bond_mod.BondOrder.single;

        try mdl.bonds.append(mdl.allocator, bond_mod.Bond{
            .atom_1 = a1,
            .atom_2 = a2,
            .order = order,
            .source = .struct_conn,
        });

        // Set bonded_inter_residue flag on both atoms
        mdl.atoms.items[a1].flags.bonded_inter_residue = true;
        mdl.atoms.items[a2].flags.bonded_inter_residue = true;
    }
}

fn isCovalentConnType(conn_type: []const u8) bool {
    var buf: [16]u8 = undefined;
    const len = @min(conn_type.len, buf.len);
    for (0..len) |i| {
        buf[i] = std.ascii.toLower(conn_type[i]);
    }
    const lower = buf[0..len];
    if (std.mem.startsWith(u8, lower, "covale")) return true;
    if (std.mem.eql(u8, lower, "disulf")) return true;
    return false;
}
```

- [ ] **Step 4: Run tests to verify they pass**

Run: `zig build test --summary all`
Expected: all tests pass

- [ ] **Step 5: Commit**

```bash
git add src/mmcif.zig
git commit -m "feat: implement parseStructConn for inter-residue bonds"
```

---

### Task 6: Implement `parseBranchLinks`

**Files:**
- Modify: `src/mmcif.zig`

- [ ] **Step 1: Write the failing test**

Add to `src/mmcif.zig`:

```zig
test "parseBranchLinks glycan bond" {
    const source = @embedFile("test_data/branch_link.cif");
    var doc = cif.readString(testing.allocator, source) catch unreachable;
    defer doc.deinit();
    const block = &doc.blocks.items[0];

    var mdl = try parseModel(testing.allocator, source);
    defer mdl.deinit();

    try parseBranchLinks(testing.allocator, &mdl, block);

    // Should have 1 bond (NAG O4 — GAL C1)
    try testing.expectEqual(@as(usize, 1), mdl.bonds.items.len);
    const bond = mdl.bonds.items[0];
    try testing.expectEqual(bond_mod.BondSource.branch_link, bond.source);

    // Leaving atoms (O4 of NAG and O1 of GAL) should have bonded_inter_residue flag
    // O4 is atom index 1, O1 is atom index 3
    try testing.expect(mdl.atoms.items[1].flags.bonded_inter_residue);  // NAG O4
    try testing.expect(mdl.atoms.items[3].flags.bonded_inter_residue);  // GAL O1
}
```

- [ ] **Step 2: Run test to verify it fails**

Run: `zig build test --summary all 2>&1 | grep "parseBranchLinks"`
Expected: compilation error

- [ ] **Step 3: Implement parseBranchLinks**

Add to `src/mmcif.zig`:

```zig
/// Parse _pdbx_entity_branch_link and add branch bonds to the Model.
/// Sets `bonded_inter_residue` flag on the leaving atoms (not the bonding atoms).
/// Accepts a pre-parsed CIF block to avoid re-parsing.
pub fn parseBranchLinks(allocator: Allocator, mdl: *Model, block: *const cif.Block) !void {
    const bl = block.findLoop("_pdbx_entity_branch_link.link_id") orelse return;

    const col_eid = bl.findTag("_pdbx_entity_branch_link.entity_id") orelse return;
    const col_num1 = bl.findTag("_pdbx_entity_branch_link.entity_branch_list_num_1") orelse return;
    const col_comp1 = bl.findTag("_pdbx_entity_branch_link.comp_id_1") orelse return;
    const col_atom1 = bl.findTag("_pdbx_entity_branch_link.atom_id_1") orelse return;
    const col_leaving1 = bl.findTag("_pdbx_entity_branch_link.leaving_atom_id_1") orelse return;
    const col_num2 = bl.findTag("_pdbx_entity_branch_link.entity_branch_list_num_2") orelse return;
    const col_comp2 = bl.findTag("_pdbx_entity_branch_link.comp_id_2") orelse return;
    const col_atom2 = bl.findTag("_pdbx_entity_branch_link.atom_id_2") orelse return;
    const col_leaving2 = bl.findTag("_pdbx_entity_branch_link.leaving_atom_id_2") orelse return;

    // Build atom lookup (includes auth_seq_id entries for branched entities)
    var lookup = try buildAtomLookup(allocator, block);
    defer lookup.deinit();

    // Build entity_id -> [asym_id] mapping
    var entity_to_asyms = std.StringHashMap(std.ArrayListUnmanaged([]const u8)).init(allocator);
    defer {
        var it = entity_to_asyms.iterator();
        while (it.next()) |entry| entry.value_ptr.deinit(allocator);
        entity_to_asyms.deinit();
    }
    for (mdl.chains.items) |chain| {
        const eid = chain.entityIdSlice();
        if (eid.len == 0) continue;
        const entry = try entity_to_asyms.getOrPut(eid);
        if (!entry.found_existing) entry.value_ptr.* = .{};
        try entry.value_ptr.append(allocator, chain.labelSlice());
    }

    for (0..bl.length()) |row| {
        const entity_id = cif.asString(bl.val(row, col_eid) orelse continue);
        const num1 = cif.asString(bl.val(row, col_num1) orelse continue);
        const comp1 = cif.asString(bl.val(row, col_comp1) orelse continue);
        const atom1_name = cif.asString(bl.val(row, col_atom1) orelse continue);
        const leaving1 = cif.asString(bl.val(row, col_leaving1) orelse continue);
        const num2 = cif.asString(bl.val(row, col_num2) orelse continue);
        const comp2 = cif.asString(bl.val(row, col_comp2) orelse continue);
        const atom2_name = cif.asString(bl.val(row, col_atom2) orelse continue);
        const leaving2 = cif.asString(bl.val(row, col_leaving2) orelse continue);

        const asyms = entity_to_asyms.get(entity_id) orelse continue;

        for (asyms.items) |asym_id| {
            // Resolve bonding atoms
            const a1 = lookup.get(.{
                .label_asym_id = asym_id,
                .label_seq_id = num1,
                .atom_name = atom1_name,
            }) orelse continue;

            const a2 = lookup.get(.{
                .label_asym_id = asym_id,
                .label_seq_id = num2,
                .atom_name = atom2_name,
            }) orelse continue;

            try mdl.bonds.append(mdl.allocator, Bond{
                .atom_1 = a1,
                .atom_2 = a2,
                .order = .single,
                .source = .branch_link,
            });

            // Set bonded_inter_residue on LEAVING atoms (not the bonding atoms)
            if (leaving1.len > 0) {
                if (lookup.get(.{
                    .label_asym_id = asym_id,
                    .label_seq_id = num1,
                    .atom_name = leaving1,
                })) |leaving_idx| {
                    mdl.atoms.items[leaving_idx].flags.bonded_inter_residue = true;
                }
            }
            if (leaving2.len > 0) {
                if (lookup.get(.{
                    .label_asym_id = asym_id,
                    .label_seq_id = num2,
                    .atom_name = leaving2,
                })) |leaving_idx| {
                    mdl.atoms.items[leaving_idx].flags.bonded_inter_residue = true;
                }
            }
        }
    }
}
```

- [ ] **Step 4: Run tests to verify they pass**

Run: `zig build test --summary all`
Expected: all tests pass

- [ ] **Step 5: Commit**

```bash
git add src/mmcif.zig
git commit -m "feat: implement parseBranchLinks for glycan inter-residue bonds"
```

---

### Task 7: Implement `parseInlineComponents`

**Files:**
- Modify: `src/mmcif.zig`

Parse inline `_chem_comp_atom` and `_chem_comp_bond` from the structure file into a `ccd.ComponentDict`. This reuses the same types as the external CCD parser.

- [ ] **Step 1: Write the failing test**

Add to `src/mmcif.zig`:

```zig
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
    try testing.expectEqual(@as(usize, 11), ala.?.atoms.len);  // 5 heavy + 5 H + OXT
    try testing.expectEqual(@as(usize, 10), ala.?.bonds.len);   // 10 bonds
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
```

- [ ] **Step 2: Run test to verify it fails**

Expected: compilation error — `parseInlineComponents` not found

- [ ] **Step 3: Implement parseInlineComponents**

Add to `src/mmcif.zig`:

```zig
const ccd_mod = @import("ccd.zig");

/// Parse inline _chem_comp_atom and _chem_comp_bond from a CIF block.
/// Returns a ComponentDict compatible with CCD dictionary, or null if not present.
/// Accepts a pre-parsed CIF block to avoid re-parsing.
pub fn parseInlineComponents(allocator: Allocator, block: *const cif.Block) !?ccd_mod.ComponentDict {
    // Check if _chem_comp_bond loop exists
    const bond_loop = block.findLoop("_chem_comp_bond.comp_id") orelse return null;

    const bcol_comp = bond_loop.findTag("_chem_comp_bond.comp_id") orelse return null;
    const bcol_a1 = bond_loop.findTag("_chem_comp_bond.atom_id_1") orelse return null;
    const bcol_a2 = bond_loop.findTag("_chem_comp_bond.atom_id_2") orelse return null;
    const bcol_order = bond_loop.findTag("_chem_comp_bond.value_order");
    const bcol_arom = bond_loop.findTag("_chem_comp_bond.pdbx_aromatic_flag");

    // Optional _chem_comp_atom loop (from same block)
    const atom_loop = block.findLoop("_chem_comp_atom.comp_id");

    var dict = ccd_mod.ComponentDict{
        .components = std.StringHashMap(ccd_mod.Component).init(allocator),
        .allocator = allocator,
    };
    errdefer dict.deinit();

    // Phase 1: Group atoms by comp_id
    const AtomList = std.ArrayListUnmanaged(ccd_mod.CompAtom);
    var atom_groups = std.StringHashMap(AtomList).init(allocator);
    defer {
        var ait = atom_groups.iterator();
        while (ait.next()) |entry| entry.value_ptr.deinit(allocator);
        atom_groups.deinit();
    }

    if (atom_loop) |al| {
        const acol_comp = al.findTag("_chem_comp_atom.comp_id") orelse return null;
        const acol_id = al.findTag("_chem_comp_atom.atom_id") orelse return null;
        const acol_type = al.findTag("_chem_comp_atom.type_symbol");
        const acol_charge = al.findTag("_chem_comp_atom.charge");
        const acol_leaving = al.findTag("_chem_comp_atom.pdbx_leaving_atom_flag");

        for (0..al.length()) |row| {
            const comp = cif.asString(al.val(row, acol_comp) orelse continue);
            const atom_id = cif.asString(al.val(row, acol_id) orelse continue);

            var atom = ccd_mod.CompAtom{};
            const len = @min(atom_id.len, 4);
            atom.name_len = @intCast(len);
            @memcpy(atom.name[0..len], atom_id[0..len]);

            if (acol_type) |c| {
                const sym = cif.asString(al.val(row, c) orelse "");
                const slen = @min(sym.len, 2);
                atom.element_symbol = .{ ' ', ' ' };
                @memcpy(atom.element_symbol[0..slen], sym[0..slen]);
            }

            if (acol_charge) |c| {
                atom.charge = cif.value.asIntOr(i8, al.val(row, c) orelse "0", 0);
            }

            if (acol_leaving) |c| {
                const lv = cif.asString(al.val(row, c) orelse "N");
                atom.leaving = lv.len > 0 and (lv[0] == 'Y' or lv[0] == 'y');
            }

            const entry = try atom_groups.getOrPut(comp);
            if (!entry.found_existing) entry.value_ptr.* = .{};
            try entry.value_ptr.append(allocator, atom);
        }
    }

    // Phase 2: Group bonds by comp_id, resolve to atom indices
    const BondEntry = struct { atom1: []const u8, atom2: []const u8, order: ccd_mod.BondOrder, aromatic: bool };
    const BondList = std.ArrayListUnmanaged(BondEntry);
    var bond_groups = std.StringHashMap(BondList).init(allocator);
    defer {
        var bit = bond_groups.iterator();
        while (bit.next()) |entry| entry.value_ptr.deinit(allocator);
        bond_groups.deinit();
    }

    for (0..bond_loop.length()) |row| {
        const comp = cif.asString(bond_loop.val(row, bcol_comp) orelse continue);
        const a1 = cif.asString(bond_loop.val(row, bcol_a1) orelse continue);
        const a2 = cif.asString(bond_loop.val(row, bcol_a2) orelse continue);
        const order = if (bcol_order) |c|
            ccd_mod.BondOrder.fromString(cif.asString(bond_loop.val(row, c) orelse "SING"))
        else
            ccd_mod.BondOrder.unknown;
        const aromatic = if (bcol_arom) |c| blk: {
            const s = cif.asString(bond_loop.val(row, c) orelse "N");
            break :blk s.len > 0 and (s[0] == 'Y' or s[0] == 'y');
        } else false;

        const entry = try bond_groups.getOrPut(comp);
        if (!entry.found_existing) entry.value_ptr.* = .{};
        try entry.value_ptr.append(allocator, .{
            .atom1 = a1, .atom2 = a2, .order = order, .aromatic = aromatic,
        });
    }

    // Phase 3: Build components
    // Collect all comp_ids from both groups
    var all_comp_ids = std.StringHashMap(void).init(allocator);
    defer all_comp_ids.deinit();
    {
        var ait = atom_groups.iterator();
        while (ait.next()) |entry| try all_comp_ids.put(entry.key_ptr.*, {});
        var bit = bond_groups.iterator();
        while (bit.next()) |entry| try all_comp_ids.put(entry.key_ptr.*, {});
    }

    var cit = all_comp_ids.iterator();
    while (cit.next()) |entry| {
        const comp_id = entry.key_ptr.*;

        // Get or create atom list
        const atoms = if (atom_groups.getPtr(comp_id)) |al|
            try allocator.dupe(ccd_mod.CompAtom, al.items)
        else
            try allocator.alloc(ccd_mod.CompAtom, 0);
        errdefer allocator.free(atoms);

        // Resolve bonds to atom indices
        var comp_bonds = std.ArrayListUnmanaged(ccd_mod.CompBond){};
        defer comp_bonds.deinit(allocator);

        if (bond_groups.getPtr(comp_id)) |bl| {
            for (bl.items) |be| {
                const idx1 = findAtomIdx(atoms, be.atom1) orelse continue;
                const idx2 = findAtomIdx(atoms, be.atom2) orelse continue;
                try comp_bonds.append(allocator, .{
                    .atom_idx_1 = idx1,
                    .atom_idx_2 = idx2,
                    .order = be.order,
                    .aromatic = be.aromatic,
                });
            }
        }

        const owned_bonds = try allocator.dupe(ccd_mod.CompBond, comp_bonds.items);
        errdefer allocator.free(owned_bonds);

        const owned_id = try allocator.dupe(u8, comp_id);
        errdefer allocator.free(owned_id);

        const key = try allocator.dupe(u8, comp_id);
        errdefer allocator.free(key);

        try dict.components.put(key, ccd_mod.Component{
            .comp_id = owned_id,
            .comp_type = try allocator.dupe(u8, ""),
            .atoms = atoms,
            .bonds = owned_bonds,
        });
    }

    return dict;
}

fn findAtomIdx(atoms: []const ccd_mod.CompAtom, name: []const u8) ?u16 {
    for (atoms, 0..) |*a, i| {
        if (std.mem.eql(u8, a.nameSlice(), name)) return @intCast(i);
    }
    return null;
}
```

- [ ] **Step 4: Run tests**

Run: `zig build test --summary all`
Expected: all tests pass

- [ ] **Step 5: Commit**

```bash
git add src/mmcif.zig
git commit -m "feat: implement parseInlineComponents for inline CCD data"
```

---

### Task 8: Integrate into placer.zig — skip H on bonded atoms

**Files:**
- Modify: `src/place/placer.zig`

- [ ] **Step 1: Write the failing test**

Add to `src/place/placer.zig` (or in mmcif.zig integration test section):

```zig
test "addHydrogens skips H on bonded_inter_residue atom" {
    // Parse disulfide fixture
    const source = @embedFile("../test_data/disulfide.cif");
    var mdl = try @import("../mmcif.zig").parseModel(testing.allocator, source);
    defer mdl.deinit();

    // Manually set bonded_inter_residue on SG atoms (index 5 and 11)
    mdl.atoms.items[5].flags.bonded_inter_residue = true;
    mdl.atoms.items[11].flags.bonded_inter_residue = true;

    // Apply chemistry and place hydrogens
    applyChemistry(&mdl);
    const result = try addHydrogens(&mdl, null);

    // Verify no HG was placed on either CYS SG
    for (mdl.atoms.items) |atom| {
        if (atom.is_added) {
            try testing.expect(!std.mem.eql(u8, atom.nameSlice(), "HG"));
        }
    }
    _ = result;
}
```

- [ ] **Step 2: Run test to verify it fails**

The test should fail because currently `addHydrogens` does not check the flag. HG may or may not be placed depending on whether CYS plans include it and whether distance-based neighbor finding succeeds. We need to verify the test can differentiate.

- [ ] **Step 3: Add bonded_inter_residue check to executePlan**

In `src/place/placer.zig`, at the start of `executePlan`, after resolving the base atom, add:

```zig
fn executePlan(mdl: *Model, res: Residue, res_idx: u32, plan: *const standard.PlacementPlan, bonds: ?[]const topology.BondEntry, target_altloc: u8) !bool {
    // Resolve parent heavy atom (connected[0]) for metadata and position
    const base_atom = findAtom(mdl, res, plan.connected[0], target_altloc) orelse return false;

    // Skip H placement if parent atom is involved in inter-residue bond
    if (base_atom.flags.bonded_inter_residue) return false;

    // ... rest of function unchanged
```

Also add the same check in the CCD-derived placement path (the `else if (ccd_dict)` branch in `addHydrogens`). Before calling `executePlan` for CCD-derived plans:

```zig
// In the ccd_derive loop, skip plans whose parent is bonded
for (plans) |plan| {
    // Check if parent atom has bonded_inter_residue flag
    const parent_atom = findAtom(mdl, res, plan.connected[0], ' ');
    if (parent_atom != null and parent_atom.?.flags.bonded_inter_residue) {
        result.n_skipped += 1;
        continue;
    }
    if (try executePlan(mdl, res, @intCast(res_idx), &plan, null, ' ')) {
        result.n_placed += 1;
    } else {
        result.n_skipped += 1;
    }
}
```

- [ ] **Step 4: Run all tests**

Run: `zig build test --summary all`
Expected: all tests pass (including existing 253+)

- [ ] **Step 5: Commit**

```bash
git add src/place/placer.zig
git commit -m "feat: skip H placement on bonded_inter_residue atoms"
```

---

### Task 9: Wire into run.zig pipeline

**Files:**
- Modify: `src/run.zig`
- Modify: `src/batch.zig`

- [ ] **Step 1: Update processFile in run.zig**

In `src/run.zig`, modify `processFile` to add inline components and bond parsing between steps 3 and 4:

```zig
pub fn processFile(allocator: Allocator, config: ProcessConfig) !ProcessResult {
    // 1. Read input mmCIF
    const source = try readFile(allocator, config.input_path);
    defer allocator.free(source);

    // 2. Parse CIF document
    var doc = try zreduce.cif.readString(allocator, source);
    defer doc.deinit();

    // 3. Extract model from CIF
    var mdl = try zreduce.mmcif.parseModel(allocator, source);
    defer mdl.deinit();

    // 3a. Parse inline component dictionary (takes priority over external CCD)
    const block = &doc.blocks.items[0];
    var inline_dict = try zreduce.mmcif.parseInlineComponents(allocator, block);
    defer if (inline_dict) |*d| d.deinit();

    // 3b. Parse inter-residue bonds (_struct_conn)
    try zreduce.mmcif.parseStructConn(allocator, &mdl, block);

    // 3c. Parse branch links (_pdbx_entity_branch_link)
    try zreduce.mmcif.parseBranchLinks(allocator, &mdl, block);

    // Effective dictionary: inline > external CCD
    const effective_dict: ?*const zreduce.ccd.ComponentDict = if (inline_dict) |*d|
        d
    else
        config.dict;

    // 4. Apply chemistry annotations
    zreduce.place.applyChemistry(&mdl);

    // 5. Place hydrogens (use effective dict)
    const place_result = try zreduce.place.addHydrogens(
        &mdl,
        effective_dict,
    );

    // ... rest unchanged, but also pass effective_dict to generateMovers
```

Also update the `generateMovers` call to use `effective_dict`:

```zig
    const gen_result = try zreduce.optimize.generateMovers(
        allocator,
        &mdl,
        config.no_flip,
        effective_dict,
    );
```

- [ ] **Step 2: Update batch.zig similarly**

In `src/batch.zig`, apply the same changes in the batch worker function. The batch worker should also call `parseInlineComponents`, `parseStructConn`, and `parseBranchLinks`, and use the effective dictionary.

Check if batch.zig calls run.processFile or has its own pipeline:

```bash
grep -n "processFile\|parseModel" src/batch.zig
```

If batch.zig calls `run.processFile`, it's already covered. If it has its own pipeline, duplicate the changes.

- [ ] **Step 3: Run all tests**

Run: `zig build test --summary all`
Expected: all tests pass

- [ ] **Step 4: Build and run a quick smoke test**

```bash
zig build -Doptimize=ReleaseFast
# Test with a real structure that has disulfide bonds
./zig-out/bin/zreduce run examples/data/AF-P0A7N4-F1-model_v4.cif -o /tmp/test_output.cif
```

- [ ] **Step 5: Commit**

```bash
git add src/run.zig src/batch.zig
git commit -m "feat: wire bond graph parsing into processing pipeline"
```

---

### Task 10: Integration test with real data

**Files:**
- Modify: `src/integration_test.zig` (or wherever integration tests live)

- [ ] **Step 1: Add integration test for struct_conn skipping**

Verify end-to-end that a structure with disulfide bonds does not place HG on bonded SG:

```zig
test "integration: disulfide SG-HG not placed" {
    const source = @embedFile("test_data/disulfide.cif");
    var doc = cif.readString(testing.allocator, source) catch unreachable;
    defer doc.deinit();
    const block = &doc.blocks.items[0];

    var mdl = try mmcif.parseModel(testing.allocator, source);
    defer mdl.deinit();

    try mmcif.parseStructConn(testing.allocator, &mdl, block);
    place.applyChemistry(&mdl);
    _ = try place.addHydrogens(&mdl, null);

    // Check no HG atom was placed
    for (mdl.atoms.items) |atom| {
        if (atom.is_added and std.mem.eql(u8, atom.nameSlice(), "HG")) {
            return error.UnexpectedHGPlacement;
        }
    }
}
```

- [ ] **Step 2: Add integration test for inline component dict**

```zig
test "integration: inline dict used for placement" {
    const source = @embedFile("test_data/inline_comp.cif");
    var doc = cif.readString(testing.allocator, source) catch unreachable;
    defer doc.deinit();
    const block = &doc.blocks.items[0];

    var mdl = try mmcif.parseModel(testing.allocator, source);
    defer mdl.deinit();

    var inline_dict = try mmcif.parseInlineComponents(testing.allocator, block);
    defer if (inline_dict) |*d| d.deinit();

    place.applyChemistry(&mdl);
    const result = try place.addHydrogens(&mdl, if (inline_dict) |*d| d else null);

    // ALA should have hydrogens placed (H, HA, HB1, HB2, HB3)
    try std.testing.expect(result.n_placed > 0);
}
```

- [ ] **Step 3: Run all tests**

Run: `zig build test --summary all`
Expected: all tests pass

- [ ] **Step 4: Commit**

```bash
git add src/integration_test.zig
git commit -m "test: add integration tests for bond graph and inline dict"
```

---

### Task 11: Verify all existing tests still pass and clean up

- [ ] **Step 1: Run full test suite**

```bash
zig build test --summary all
```

Expected: all 253+ tests plus new tests pass.

- [ ] **Step 2: Build release and test with real structures**

```bash
zig build -Doptimize=ReleaseFast
./zig-out/bin/zreduce run examples/data/AF-P0A7N4-F1-model_v4.cif -o /tmp/test1.cif
```

Compare output with previous version to ensure no regressions in standard cases.

- [ ] **Step 3: Test with glycoprotein (if available)**

```bash
# Download 5fyl if not already available
./zig-out/bin/zreduce run examples/data/5fyl.cif -o /tmp/test_glyco.cif
```

- [ ] **Step 4: Final commit if any cleanup needed**

```bash
git add -A
git commit -m "chore: final cleanup for bond graph feature"
```
