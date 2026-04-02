# compile-dict: Pre-compiled CCD Binary Dictionary

**Issue**: #52
**Date**: 2026-04-02
**Status**: Approved

## Problem

Loading `components.cif` (~1GB) at runtime requires full text tokenization and parsing,
which is the startup bottleneck for both `run` and `batch` modes. A pre-compiled binary
format eliminates this parsing overhead.

## Usage

```bash
# One-time compilation
zreduce compile-dict components.cif -o components.zdict

# Use compiled dictionary (auto-detected)
zreduce run input.cif -d components.zdict -o output.cif
zreduce batch input_dir/ -d components.zdict -o output_dir/
```

## Design Decisions

| Decision | Choice | Rationale |
|----------|--------|-----------|
| Data granularity | Component level (CompAtom + CompBond) | CCD parse is the bottleneck; derivePlans is cheap and depends on runtime residue state |
| Format detection | Magic bytes + recommended `.zdict` extension | Robust detection without extension dependency; extension for UX clarity |
| Versioning | Header version byte, reject on mismatch | Component structure is stable; re-run compile-dict on mismatch |
| Input scope | Full CCD only, no filtering | YAGNI; binary is already compact enough |
| Subcommand name | `compile-dict` | Unambiguous; avoids collision with future compile targets |
| Serialization | Custom flat binary | Zero dependencies; aligns with full-load-into-HashMap usage pattern |

## Binary Format

```
Offset  Size     Field
────────────────────────────────────────────
0       4        Magic: "ZRDC" (0x5A 0x52 0x44 0x43)
4       1        Format version: 1
5       3        Reserved (zero-padded)
8       4        Component count (u32 little-endian)
12      ...      Component records (sequential)
```

### Component Record

```
[1B comp_id length] [N bytes comp_id]
[1B comp_type length] [N bytes comp_type]
[2B atom count (u16 LE)] [CompAtom x N]
[2B bond count (u16 LE)] [CompBond x N]
```

### CompAtom (24 bytes, fixed)

```
Offset  Size  Field
0       4     name ([4]u8, space-padded)
4       1     name_len (u4 stored as u8)
5       2     element_symbol ([2]u8)
7       1     charge (i8)
8       1     flags (bit0=leaving, bit1=aromatic)
9       3     padding (for f32 alignment)
12      4     ideal_x (f32 LE)
16      4     ideal_y (f32 LE)
20      4     ideal_z (f32 LE)
```

Total: 24 bytes.

### CompBond (6 bytes, fixed)

```
Offset  Size  Field
0       2     atom_idx_1 (u16 LE)
2       2     atom_idx_2 (u16 LE)
4       1     order (u8, BondOrder enum value)
5       1     flags (bit0=aromatic)
```

## Module Structure

### New file: `src/ccd_binary.zig`

```zig
pub const MAGIC = "ZRDC";
pub const FORMAT_VERSION: u8 = 1;

/// Write a ComponentDict to binary format.
pub fn writeDict(writer: anytype, dict: *const ComponentDict) !void

/// Read binary format into a ComponentDict.
pub fn readDict(allocator: Allocator, reader: anytype) !ComponentDict

/// Check if data starts with ZRDC magic bytes.
pub fn isBinaryDict(header: []const u8) bool
```

### Changes to existing files

**`src/main.zig`**:
- Add `compile-dict` subcommand parsing and dispatch
- In `runSubcommand`: after `readFile`, check `isBinaryDict` to choose read path
- compile-dict subcommand: read CIF, parse, write binary, print stats

**`src/batch.zig`**:
- Same `isBinaryDict` auto-detection for `-d` flag

**`src/root.zig`**:
- Re-export `ccd_binary` module

## Error Handling

| Scenario | Behavior |
|----------|----------|
| Magic bytes mismatch | `error.InvalidMagic` — "not a valid zreduce dictionary file" |
| Version mismatch | `error.UnsupportedVersion` — "dictionary format version N not supported, please re-run compile-dict" |
| Truncated data | `error.UnexpectedEof` |
| compile-dict: `-o` missing | Error: output path required |
| compile-dict: input is already binary | Error: "input is already a compiled dictionary" |
| compile-dict: CIF parse failure | Existing error path |

No fallback or retry logic. Errors terminate with a clear message.

## Testing

### Unit tests (`src/ccd_binary.zig`)

1. **Round-trip**: build ComponentDict by hand → writeDict → readDict → verify all fields match
2. **`isBinaryDict`**: valid magic, invalid magic, empty, too-short data
3. **Version mismatch**: tampered version byte → `error.UnsupportedVersion`
4. **Empty dictionary**: 0 components round-trip
5. **Truncated data**: partial binary → `error.UnexpectedEof`

### Integration test (manual)

- `components.cif` → `compile-dict` → `.zdict` → `run`/`batch` output matches CIF-direct results
