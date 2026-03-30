# Entity Type + seq_id Fallback Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Parse `_entity.type` to set `Residue.entity_type`, add `auth_seq_id` fallback when `label_seq_id` is `.` for non-polymer/branched entities.

**Architecture:** Add `branched` to `EntityType` enum and `auth_seq_id` to `Residue`. In `parseModel`, parse `auth_seq_id` column, apply label→auth fallback for seq_id, then post-process entity types from `_entity` loop.

**Tech Stack:** Zig, existing CIF parser

---

## File Map

| File | Action | Responsibility |
|------|--------|----------------|
| `src/model/residue.zig` | Modify | Add `branched` to EntityType, add `auth_seq_id` field |
| `src/mmcif.zig` | Modify | Parse `auth_seq_id`, implement fallback, parse `_entity` loop |
| `src/test_data/entity_type.cif` | Create | Test fixture with polymer + non-polymer + water + branched |

---

### Task 1: Add `branched` to EntityType and `auth_seq_id` to Residue

**Files:**
- Modify: `src/model/residue.zig`

- [ ] **Step 1: Write the failing test**

Add to `src/model/residue.zig`:

```zig
test "EntityType has branched variant" {
    const e: EntityType = .branched;
    try std.testing.expect(e == .branched);
}

test "Residue auth_seq_id default" {
    const r = Residue{};
    try std.testing.expectEqual(@as(i32, 0), r.auth_seq_id);
}
```

- [ ] **Step 2: Run test to verify it fails**

Run: `zig build test --summary all 2>&1 | grep -E "branched|auth_seq"`
Expected: compilation error

- [ ] **Step 3: Implement changes**

In `src/model/residue.zig`:

Change EntityType:
```zig
pub const EntityType = enum { polymer, non_polymer, branched, water, unknown };
```

Add field to Residue (after `seq_id`):
```zig
seq_id: i32 = 0,
auth_seq_id: i32 = 0,
```

- [ ] **Step 4: Run all tests**

Run: `zig build test --summary all`
Expected: all tests pass

- [ ] **Step 5: Commit**

```bash
git add src/model/residue.zig
git commit -m "feat: add branched EntityType and auth_seq_id to Residue"
```

---

### Task 2: Create test fixture with entity types

**Files:**
- Create: `src/test_data/entity_type.cif`

- [ ] **Step 1: Create fixture**

Create `src/test_data/entity_type.cif` — structure with polymer (ALA), non-polymer (EDO), and water (HOH), each in different chains with different entity types. Branched entity (NAG) with `label_seq_id = .` and `auth_seq_id` numbering:

```cif
data_ENTITY
#
loop_
_entity.id
_entity.type
1 polymer
2 non-polymer
3 water
4 branched
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
_atom_site.auth_seq_id
_atom_site.Cartn_x
_atom_site.Cartn_y
_atom_site.Cartn_z
_atom_site.occupancy
_atom_site.B_iso_or_equiv
_atom_site.label_alt_id
_atom_site.pdbx_PDB_ins_code
ATOM   1  N N   ALA A 1 1  1  1.0 0.0 0.0 1.0 0.0 . .
ATOM   2  C CA  ALA A 1 1  1  2.0 0.0 0.0 1.0 0.0 . .
ATOM   3  C C   ALA A 1 1  1  3.0 0.0 0.0 1.0 0.0 . .
ATOM   4  O O   ALA A 1 1  1  4.0 0.0 0.0 1.0 0.0 . .
ATOM   5  C CB  ALA A 1 1  1  2.5 1.0 0.0 1.0 0.0 . .
HETATM 6  C C1  EDO B 2 .  1  5.0 0.0 0.0 1.0 0.0 . .
HETATM 7  O O1  EDO B 2 .  1  6.0 0.0 0.0 1.0 0.0 . .
HETATM 8  O O   HOH C 3 .  1  7.0 0.0 0.0 1.0 0.0 . .
HETATM 9  C C1  NAG D 4 .  1  8.0 0.0 0.0 1.0 0.0 . .
HETATM 10 O O4  NAG D 4 .  1  9.0 0.0 0.0 1.0 0.0 . .
HETATM 11 C C1  NAG D 4 .  2 10.0 0.0 0.0 1.0 0.0 . .
HETATM 12 O O1  NAG D 4 .  2 11.0 0.0 0.0 1.0 0.0 . .
#
```

Key points:
- Entity 1 (chain A): polymer, `label_seq_id=1`
- Entity 2 (chain B): non-polymer EDO, `label_seq_id=.`, `auth_seq_id=1`
- Entity 3 (chain C): water, `label_seq_id=.`, `auth_seq_id=1`
- Entity 4 (chain D): branched NAG, `label_seq_id=.`, `auth_seq_id=1,2` (two residues)

- [ ] **Step 2: Commit**

```bash
git add src/test_data/entity_type.cif
git commit -m "test: add entity type fixture with polymer/non-polymer/water/branched"
```

---

### Task 3: Parse `auth_seq_id` and implement seq_id fallback

**Files:**
- Modify: `src/mmcif.zig`

- [ ] **Step 1: Write the failing test**

Add to `src/mmcif.zig`:

```zig
test "parseModel stores auth_seq_id" {
    const source = @embedFile("test_data/entity_type.cif");
    var mdl = try parseModel(testing.allocator, source);
    defer mdl.deinit();

    // Polymer ALA: label_seq_id=1, auth_seq_id=1
    try testing.expectEqual(@as(i32, 1), mdl.residues.items[0].seq_id);
    try testing.expectEqual(@as(i32, 1), mdl.residues.items[0].auth_seq_id);

    // Non-polymer EDO: label_seq_id="." -> fallback to auth_seq_id=1
    try testing.expectEqual(@as(i32, 1), mdl.residues.items[1].seq_id);
    try testing.expectEqual(@as(i32, 1), mdl.residues.items[1].auth_seq_id);

    // Water HOH: label_seq_id="." -> fallback to auth_seq_id=1
    try testing.expectEqual(@as(i32, 1), mdl.residues.items[2].seq_id);

    // Branched NAG residue 1: label_seq_id="." -> fallback to auth_seq_id=1
    try testing.expectEqual(@as(i32, 1), mdl.residues.items[3].seq_id);
    try testing.expectEqual(@as(i32, 1), mdl.residues.items[3].auth_seq_id);

    // Branched NAG residue 2: label_seq_id="." -> fallback to auth_seq_id=2
    try testing.expectEqual(@as(i32, 2), mdl.residues.items[4].seq_id);
    try testing.expectEqual(@as(i32, 2), mdl.residues.items[4].auth_seq_id);
}
```

- [ ] **Step 2: Run test to verify it fails**

Expected: `auth_seq_id` not a field (compile error), or seq_id values wrong

- [ ] **Step 3: Implement auth_seq_id parsing and fallback**

In `src/mmcif.zig`, add to `AtomSiteColumns`:
```zig
auth_seq_id: ?usize = null,
```

In column mapping (after existing mappings):
```zig
cols.auth_seq_id = loop.findTag("_atom_site.auth_seq_id");
```

In the field extraction section (around line 160), add after `seq_id`:
```zig
const auth_seq = if (cols.auth_seq_id) |c| cif.asString(loop.val(row, c) orelse ".") else "";
```

In residue creation (around line 216), change:
```zig
// Old:
new_res.seq_id = cif.value.asIntOr(i32, seq_id, 0);

// New: fallback to auth_seq_id when label_seq_id is empty (non-polymer/branched)
const label_seq = cif.value.asIntOr(i32, seq_id, 0);
const auth_seq_int = cif.value.asIntOr(i32, auth_seq, 0);
new_res.seq_id = if (seq_id.len == 0) auth_seq_int else label_seq;
new_res.auth_seq_id = auth_seq_int;
```

Note: `cif.asString` converts `.` and `?` to empty string, so `seq_id.len == 0` catches both.

Also update residue boundary detection to include `auth_seq` in the comparison:
```zig
// Around line 200-203, add auth_seq to tracking
const new_residue = !in_residue or
    !std.mem.eql(u8, seq_id, cur_seq_id) or
    !std.mem.eql(u8, comp_id, cur_comp_id) or
    !std.mem.eql(u8, ins_code_str, cur_ins_code) or
    !std.mem.eql(u8, auth_seq, cur_auth_seq);
```

Add `cur_auth_seq` state variable alongside existing state vars:
```zig
var cur_auth_seq: []const u8 = "";
```

And update it in the new residue block:
```zig
cur_auth_seq = auth_seq;
```

Reset it on chain change:
```zig
cur_auth_seq = "";
```

- [ ] **Step 4: Run all tests**

Run: `zig build test --summary all`
Expected: all tests pass

- [ ] **Step 5: Commit**

```bash
git add src/mmcif.zig
git commit -m "feat: parse auth_seq_id with label->auth fallback for seq_id"
```

---

### Task 4: Parse `_entity` loop and set entity_type

**Files:**
- Modify: `src/mmcif.zig`

- [ ] **Step 1: Write the failing test**

Add to `src/mmcif.zig`:

```zig
test "parseModel sets entity_type from _entity loop" {
    const source = @embedFile("test_data/entity_type.cif");
    var mdl = try parseModel(testing.allocator, source);
    defer mdl.deinit();

    // 5 residues: ALA, EDO, HOH, NAG(1), NAG(2)
    try testing.expectEqual(@as(usize, 5), mdl.residues.items.len);

    // ALA (entity 1) -> polymer
    try testing.expectEqual(Residue.EntityType.polymer, mdl.residues.items[0].entity_type);

    // EDO (entity 2) -> non_polymer
    try testing.expectEqual(Residue.EntityType.non_polymer, mdl.residues.items[1].entity_type);

    // HOH (entity 3) -> water
    try testing.expectEqual(Residue.EntityType.water, mdl.residues.items[2].entity_type);

    // NAG (entity 4) -> branched
    try testing.expectEqual(Residue.EntityType.branched, mdl.residues.items[3].entity_type);
    try testing.expectEqual(Residue.EntityType.branched, mdl.residues.items[4].entity_type);
}
```

- [ ] **Step 2: Run test to verify it fails**

Expected: entity_type is `.unknown` for all residues

- [ ] **Step 3: Implement entity type parsing**

In `src/mmcif.zig`, add a helper function before `parseModel`:

```zig
/// Parse _entity loop and set entity_type on residues based on chain entity_id.
fn applyEntityTypes(mdl: *Model, block: *const cif.Block) void {
    const ent = block.findLoop("_entity.id") orelse return;
    const col_id = ent.findTag("_entity.id") orelse return;
    const col_type = ent.findTag("_entity.type") orelse return;

    for (mdl.residues.items) |*res| {
        const chain = mdl.chains.items[res.chain_idx];
        const entity_id = chain.entityIdSlice();
        if (entity_id.len == 0) continue;

        // Look up entity type
        for (0..ent.length()) |row| {
            const eid = cif.asString(ent.val(row, col_id) orelse continue);
            if (std.mem.eql(u8, eid, entity_id)) {
                const etype = cif.asString(ent.val(row, col_type) orelse continue);
                res.entity_type = entityTypeFromString(etype);
                break;
            }
        }
    }
}

fn entityTypeFromString(s: []const u8) Residue.EntityType {
    if (std.ascii.eqlIgnoreCase(s, "polymer")) return .polymer;
    if (std.ascii.eqlIgnoreCase(s, "non-polymer")) return .non_polymer;
    if (std.ascii.eqlIgnoreCase(s, "branched")) return .branched;
    if (std.ascii.eqlIgnoreCase(s, "water")) return .water;
    return .unknown;
}
```

Then in `parseModel`, change from taking `source: []const u8` to also doing entity parsing. Currently `parseModel` parses the CIF internally. Add the entity parsing after the existing `_pdbx_poly_seq_scheme` and `_pdbx_unobs_or_zero_occ_atoms` blocks (around line 278):

```zig
// Parse _entity loop for entity types (optional)
applyEntityTypes(&mdl, block);
```

IMPORTANT: `parseModel` currently creates its own `doc` internally and defers `doc.deinit()`. The `block` reference is available inside `parseModel` as `&doc.blocks.items[0]`. Add the `applyEntityTypes` call before `return mdl;`, after the existing optional loop parsing.

- [ ] **Step 4: Run all tests**

Run: `zig build test --summary all`
Expected: all tests pass

- [ ] **Step 5: Commit**

```bash
git add src/mmcif.zig
git commit -m "feat: parse _entity.type and set Residue.entity_type (#65)"
```

---

### Task 5: Verify backward compatibility and existing tests

- [ ] **Step 1: Run full test suite**

```bash
zig build test --summary all
```

Expected: all 292+ tests pass.

- [ ] **Step 2: Build release and smoke test**

```bash
zig build -Doptimize=ReleaseFast
./zig-out/bin/zreduce run examples/data/AF-P0A9J6-F1-model_v6.cif -o /tmp/test_entity.cif
```

Expected: works correctly. AlphaFold models may not have `_entity` loop, so entity_type stays `.unknown` — this is fine.

- [ ] **Step 3: Commit if any cleanup needed**

```bash
git add -A
git commit -m "chore: verify backward compatibility for entity type + seq_id fallback"
```
