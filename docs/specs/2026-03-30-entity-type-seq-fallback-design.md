# Entity Type Detection + seq_id Fallback Design

**Date:** 2026-03-30
**Issues:** #65, #62
**Status:** Approved

## Overview

Parse `_entity.type` to set `Residue.entity_type`, and implement `label_seq_id` → `auth_seq_id` fallback for non-polymer/branched entities where `label_seq_id` is `.`.

## Changes

### 1. EntityType Enum

Add `branched` variant:

```
pub const EntityType = enum { polymer, non_polymer, branched, water, unknown };
```

### 2. Residue Struct

Add `auth_seq_id: i32` field. Both `seq_id` (from label) and `auth_seq_id` (from auth) are always stored.

For non-polymer/branched entities where `label_seq_id` is `.`, `seq_id` receives the `auth_seq_id` value as fallback.

### 3. Entity Type Parsing (mmcif.zig)

After `parseModel` constructs the Model, parse `_entity` loop:
- Map `_entity.id` → `_entity.type` string
- For each Residue, look up its chain's `entity_id` in the map
- Set `entity_type` based on the type string: `"polymer"` → `.polymer`, `"non-polymer"` → `.non_polymer`, `"branched"` → `.branched`, `"water"` → `.water`

### 4. seq_id Fallback (mmcif.zig)

In `parseModel`, during residue creation:
- Always parse and store `auth_seq_id`
- When `label_seq_id` is `.` or `?` (empty after `cif.asString`): use `auth_seq_id` as `seq_id`
- When `label_seq_id` is a valid integer: use it as `seq_id`

### 5. AtomSiteColumns

Add `auth_seq_id: ?usize = null` column mapping.

### 6. Impact on Existing Code

- `buildAtomLookup`: Already registers auth_seq_id entries. No change needed.
- `mmcif_writer`: Already uses `entity_type` for ATOM/HETATM. Now gets correct values.
- `detectChainBreaks`: Uses `res.seq_id` — now correct for non-polymers.
- `parseStructConn`/`parseBranchLinks`: Unaffected (use AtomLookup, not Residue.seq_id).

## Reference

- atomworks: `src/atomworks/io/utils/bonds.py` lines 424-439 (polymer→label, non-polymer→auth)
- atomworks: `src/atomworks/io/transforms/categories.py` lines 57-136 (entity type parsing)
