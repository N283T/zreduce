# Chain-Break Detection + Unobserved Atom Handling — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Detect sequence gaps from `_pdbx_poly_seq_scheme` so residues after chain breaks are treated as N-terminal (NH3+), and parse `_pdbx_unobs_or_zero_occ_atoms` for diagnostic awareness of missing atoms.

**Architecture:** Add `is_chain_break_before: bool` to Residue. After `_atom_site` parsing in `parseModel`, optionally parse `_pdbx_poly_seq_scheme` to detect gaps and mark residues. Extend `is_nterm`/`is_cterm` detection in `placer.zig` and `applyChemistry` to include chain-break awareness. Parse `_pdbx_unobs_or_zero_occ_atoms` into a lightweight set on Model for future diagnostic use.

**Tech Stack:** Zig, zig test

---

### Task 1: Add `is_chain_break_before` to Residue

**Files:**
- Modify: `src/model/residue.zig`

- [ ] **Step 1: Add field**

Add to the Residue struct (after `entity_type`):

```zig
    is_chain_break_before: bool = false,
```

- [ ] **Step 2: Run tests**

Run: `cd /Users/nagaet/zreduce && zig build test --summary all 2>&1 | tail -20`
Expected: All tests pass (default value false preserves behavior).

- [ ] **Step 3: Commit**

```bash
git add src/model/residue.zig
git commit -m "feat: add is_chain_break_before field to Residue"
```

---

### Task 2: Parse `_pdbx_poly_seq_scheme` for chain-break detection

**Files:**
- Modify: `src/mmcif.zig`

- [ ] **Step 1: Create test fixture with sequence gap**

Create `src/test_data/gap_chain.cif` — two ALA residues (seq_id 1 and 3) with a gap at seq_id 2:

```
data_GAP
#
_entry.id GAP
#
loop_
_atom_site.group_PDB
_atom_site.id
_atom_site.type_symbol
_atom_site.label_atom_id
_atom_site.label_comp_id
_atom_site.label_asym_id
_atom_site.label_seq_id
_atom_site.Cartn_x
_atom_site.Cartn_y
_atom_site.Cartn_z
_atom_site.occupancy
_atom_site.B_iso_or_equiv
_atom_site.label_alt_id
_atom_site.auth_asym_id
ATOM 1  N N   ALA A 1 1.000 2.000 3.000 1.00 10.0 . A
ATOM 2  C CA  ALA A 1 2.000 3.000 4.000 1.00 10.0 . A
ATOM 3  C C   ALA A 1 3.000 4.000 5.000 1.00 10.0 . A
ATOM 4  O O   ALA A 1 4.000 5.000 6.000 1.00 10.0 . A
ATOM 5  C CB  ALA A 1 2.500 2.500 3.500 1.00 10.0 . A
ATOM 6  N N   ALA A 3 11.000 12.000 13.000 1.00 10.0 . A
ATOM 7  C CA  ALA A 3 12.000 13.000 14.000 1.00 10.0 . A
ATOM 8  C C   ALA A 3 13.000 14.000 15.000 1.00 10.0 . A
ATOM 9  O O   ALA A 3 14.000 15.000 16.000 1.00 10.0 . A
ATOM 10 C CB  ALA A 3 12.500 12.500 13.500 1.00 10.0 . A
#
loop_
_pdbx_poly_seq_scheme.asym_id
_pdbx_poly_seq_scheme.entity_id
_pdbx_poly_seq_scheme.seq_id
_pdbx_poly_seq_scheme.mon_id
_pdbx_poly_seq_scheme.pdb_seq_num
_pdbx_poly_seq_scheme.auth_seq_num
_pdbx_poly_seq_scheme.pdb_strand_id
A 1 1 ALA 1 1 A
A 1 2 ALA 2 ? A
A 1 3 ALA 3 3 A
#
```

Row 2 has `auth_seq_num = ?` → residue at seq_id 2 is unobserved → chain break between seq_id 1 and 3.

- [ ] **Step 2: Write the failing test**

Add to `src/mmcif.zig` test section:

```zig
test "parse chain break from pdbx_poly_seq_scheme" {
    const source = @embedFile("test_data/gap_chain.cif");
    var mdl = try parseModel(testing.allocator, source);
    defer mdl.deinit();

    // Two residues in one chain
    try testing.expectEqual(@as(usize, 2), mdl.residues.items.len);
    try testing.expectEqual(@as(usize, 1), mdl.chains.items.len);

    // First residue (seq_id 1): no chain break before
    try testing.expect(!mdl.residues.items[0].is_chain_break_before);
    // Second residue (seq_id 3): chain break before (seq_id 2 is unobserved)
    try testing.expect(mdl.residues.items[1].is_chain_break_before);
}
```

- [ ] **Step 3: Run test to verify it fails**

Expected: FAIL — `is_chain_break_before` is always false.

- [ ] **Step 4: Implement `_pdbx_poly_seq_scheme` parsing**

In `parseModel`, after closing the final chain (line ~198) and before `return mdl;`, add:

```zig
    // Parse _pdbx_poly_seq_scheme for chain-break detection (optional)
    if (block.findLoop("_pdbx_poly_seq_scheme.seq_id")) |pss| {
        const col_asym = pss.findTag("_pdbx_poly_seq_scheme.asym_id");
        const col_seq = pss.findTag("_pdbx_poly_seq_scheme.seq_id");
        const col_auth_seq = pss.findTag("_pdbx_poly_seq_scheme.auth_seq_num");

        if (col_asym != null and col_seq != null and col_auth_seq != null) {
            detectChainBreaks(&mdl, pss, col_asym.?, col_seq.?, col_auth_seq.?);
        }
    }
```

Add the helper function before `parseModel`:

```zig
/// Scan _pdbx_poly_seq_scheme to detect sequence gaps and mark residues
/// that follow an unobserved residue (auth_seq_num = '?' or '.').
fn detectChainBreaks(
    mdl: *Model,
    pss: *const cif.types.Loop,
    col_asym: usize,
    col_seq: usize,
    col_auth_seq: usize,
) void {
    // For each chain, track whether we've seen a gap since the last observed residue.
    // When we encounter an observed residue after a gap, find the matching model
    // residue and set is_chain_break_before = true.

    const nrows = pss.length();
    var prev_asym: []const u8 = "";
    var gap_pending = false;

    for (0..nrows) |row| {
        const asym = cif.asString(pss.val(row, col_asym) orelse continue);
        const seq_str = pss.val(row, col_seq) orelse continue;
        const auth_seq = cif.asString(pss.val(row, col_auth_seq) orelse "?");

        // Reset gap tracking on chain change
        if (!std.mem.eql(u8, asym, prev_asym)) {
            gap_pending = false;
            prev_asym = asym;
        }

        // Check if this scheme row is unobserved
        if (auth_seq.len == 0 or std.mem.eql(u8, auth_seq, "?")) {
            gap_pending = true;
            continue;
        }

        // This is an observed row. If there was a gap before it, find
        // the corresponding model residue and mark it.
        if (gap_pending) {
            const seq_id = cif.value.asIntOr(i32, seq_str, 0);
            // Find the model residue matching (asym_id, seq_id)
            for (mdl.residues.items) |*res| {
                const chain = mdl.chains.items[res.chain_idx];
                const chain_asym = chain.labelAsymIdSlice();
                if (std.mem.eql(u8, chain_asym, asym) and res.seq_id == seq_id) {
                    res.is_chain_break_before = true;
                    break;
                }
            }
            gap_pending = false;
        }
    }
}
```

Note: Check what method gives the chain label_asym_id. The Chain struct stores it as `label_asym_id: [4]u8`. There should be a `labelAsymIdSlice()` or similar. If not, you may need to add one or compare the raw bytes. Read `src/model/chain.zig` to find the accessor.

- [ ] **Step 5: Run tests**

Run: `cd /Users/nagaet/zreduce && zig build test --summary all 2>&1 | tail -20`
Expected: All tests pass.

- [ ] **Step 6: Commit**

```bash
git add src/mmcif.zig src/test_data/gap_chain.cif
git commit -m "feat: parse _pdbx_poly_seq_scheme for chain-break detection"
```

---

### Task 3: Extend is_nterm/is_cterm with chain-break awareness

**Files:**
- Modify: `src/place/placer.zig`

- [ ] **Step 1: Write the failing test**

```zig
test "chain break residue gets NH3+ placement" {
    const source = @embedFile("../test_data/gap_chain.cif");
    var mdl = try mmcif.parseModel(testing.allocator, source);
    defer mdl.deinit();

    applyChemistry(&mdl);
    _ = try addHydrogens(&mdl, null);

    // Second residue (seq_id 3) should be treated as N-terminal after chain break
    // → should have H1, H2, H3 (NH3+) instead of single backbone H
    var h1_count: u32 = 0;
    for (mdl.atoms.items) |atom| {
        if (atom.is_added and atom.residue_idx == 1 and std.mem.eql(u8, atom.nameSlice(), "H1")) {
            h1_count += 1;
        }
    }
    try testing.expectEqual(@as(u32, 1), h1_count);
}
```

- [ ] **Step 2: Run test to verify it fails**

Expected: FAIL — H1 not placed on second residue (treated as internal, gets backbone H instead).

- [ ] **Step 3: Update `addHydrogens` N-terminal detection**

In `addHydrogens`, change:
```zig
const is_nterm = (res_idx == chain.residue_start);
```
to:
```zig
const is_nterm = (res_idx == chain.residue_start) or res.is_chain_break_before;
```

- [ ] **Step 4: Update `applyChemistry` terminal detection**

In `applyChemistry`, update:
```zig
const is_nterm = (res_idx == chain.residue_start);
const is_cterm = (res_idx == chain.residue_end - 1);
```
to:
```zig
const is_nterm = (res_idx == chain.residue_start) or res.is_chain_break_before;
const is_cterm = (res_idx == chain.residue_end - 1) or
    (res_idx + 1 < n_residues and mdl.residues.items[res_idx + 1].is_chain_break_before);
```

Note: `n_residues` is already available as the loop bound in `applyChemistry`.

- [ ] **Step 5: Run tests**

Run: `cd /Users/nagaet/zreduce && zig build test --summary all 2>&1 | tail -20`
Expected: All tests pass.

- [ ] **Step 6: Add test for C-terminal before chain break**

```zig
test "residue before chain break gets C-terminal charge" {
    const source = @embedFile("../test_data/gap_chain.cif");
    var mdl = try mmcif.parseModel(testing.allocator, source);
    defer mdl.deinit();

    applyChemistry(&mdl);

    // First residue (seq_id 1): should get negative flag on O (C-terminal before gap)
    const res0 = mdl.residues.items[0];
    const atoms0 = mdl.atoms.items[res0.atom_start..res0.atom_end];
    for (atoms0) |atom| {
        if (std.mem.eql(u8, atom.nameSlice(), "O")) {
            try testing.expect(atom.flags.negative);
        }
    }
}
```

- [ ] **Step 7: Run tests, commit**

```bash
git add src/place/placer.zig
git commit -m "feat: extend N/C-terminal detection with chain-break awareness"
```

---

### Task 4: Parse `_pdbx_unobs_or_zero_occ_atoms` (diagnostic)

**Files:**
- Modify: `src/mmcif.zig`
- Modify: `src/model/model.zig` (add unobs count field)

- [ ] **Step 1: Add unobs atom count to Model**

In `src/model/model.zig`, add a simple counter to Model:

```zig
    n_unobs_atoms: u32 = 0,
```

This is a lightweight approach — just track the count for now. A full hash set can be added later if diagnostic per-atom queries are needed.

- [ ] **Step 2: Parse the loop in `parseModel`**

After the `_pdbx_poly_seq_scheme` parsing block, add:

```zig
    // Parse _pdbx_unobs_or_zero_occ_atoms count (optional, diagnostic)
    if (block.findLoop("_pdbx_unobs_or_zero_occ_atoms.label_atom_id")) |unobs| {
        mdl.n_unobs_atoms = @intCast(unobs.length());
    }
```

- [ ] **Step 3: Add test**

```zig
test "parse without pdbx_poly_seq_scheme is backward compatible" {
    const source = @embedFile("test_data/tiny.cif");
    var mdl = try parseModel(testing.allocator, source);
    defer mdl.deinit();

    // No chain breaks when scheme is absent
    try testing.expect(!mdl.residues.items[0].is_chain_break_before);
    try testing.expectEqual(@as(u32, 0), mdl.n_unobs_atoms);
}
```

- [ ] **Step 4: Run tests, commit**

```bash
git add src/mmcif.zig src/model/model.zig
git commit -m "feat: parse _pdbx_unobs_or_zero_occ_atoms count for diagnostics"
```

---

### Task 5: Regression test and full verification

- [ ] **Step 1: Run full test suite**

Run: `cd /Users/nagaet/zreduce && zig build test --summary all 2>&1`
Expected: All tests pass.

- [ ] **Step 2: Run full build**

Run: `cd /Users/nagaet/zreduce && zig build 2>&1`
Expected: Clean build.

- [ ] **Step 3: Verify existing tests unchanged**

All existing tests use CIF files without `_pdbx_poly_seq_scheme`, so `is_chain_break_before` defaults to false and behavior is preserved.

- [ ] **Step 4: Commit if adjustments needed**

```bash
git add -A
git commit -m "test: verify full regression for chain-break detection"
```
