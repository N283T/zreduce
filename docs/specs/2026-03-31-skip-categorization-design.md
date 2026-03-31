# Issue #93: Categorize executePlan Skip Reasons

## Problem

`executePlan` returns `bool` — all skip reasons collapse into a single `n_skipped` counter. This made the nucleotide plan bug (#91, 455 silently skipped H) invisible.

## Design: Typed Skip Reasons (Approach A)

### New Types

```zig
pub const PlaceResult = enum {
    placed,
    existing_h,        // H already present in residue
    inter_residue,     // parent atom bonded_inter_residue (disulfide, glycosidic)
    missing_parent,    // parent heavy atom (connected[0]) not found
    missing_ref,       // reference neighbor or geometric lookup failed
};

pub const PlacementResult = struct {
    n_placed: u32 = 0,
    n_skipped_existing: u32 = 0,
    n_skipped_inter_residue: u32 = 0,
    n_skipped_missing_ref: u32 = 0,
    n_residues: u32 = 0,
};
```

### Changes

1. **`executePlan`**: return `!PlaceResult` instead of `!bool`
   - L195 parent not found → `.missing_parent`
   - L199 inter-residue → `.inter_residue`
   - L208 H exists → `.existing_h`
   - All other `return false` (missing neighbor/geometric) → `.missing_ref`
   - Success → `.placed`

2. **`addHydrogens`**: switch on PlaceResult to increment correct counter

3. **`PlacementResult`**: replace `n_skipped` with 3 categorized counters.
   Add `totalSkipped()` helper for backward compat in output.

4. **Consumers** (`main.zig`, `run.zig`, `batch.zig`): update to use categorized counters.
   - `main.zig`: show breakdown when `--validate` or verbose
   - `run.zig`: adapt ProcessResult
   - `batch.zig`: use `n_placed` (already does, `n_skipped` unused in JSONL)

5. **Warning**: when `n_skipped_missing_ref > 0`, emit stderr warning with count.

6. **Tests**: update existing count tests to use new fields.

### Non-goals

- Per-residue/per-atom warning messages (too verbose; save for `--validate` diagnostics later)
- Changing nterm placement skip tracking (separate functions, keep as-is for now)
