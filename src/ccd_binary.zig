//! Binary format for CCD (Chemical Component Dictionary).
//! Provides fast load/save of ComponentDict without CIF text parsing.
//!
//! Format:
//!   Header: MAGIC(4) + version(1) + reserved(3) + component_count(u32 LE)
//!   Per component:
//!     comp_id_len(u8) + comp_id(bytes)
//!     comp_type_len(u8) + comp_type(bytes)
//!     atom_count(u16 LE) + PackedAtom * atom_count
//!     bond_count(u16 LE) + PackedBond * bond_count

const std = @import("std");
const Allocator = std.mem.Allocator;
const ccd = @import("ccd.zig");

pub const ComponentDict = ccd.ComponentDict;
pub const Component = ccd.Component;
pub const CompAtom = ccd.CompAtom;
pub const CompBond = ccd.CompBond;

// ---------------------------------------------------------------------------
// Constants
// ---------------------------------------------------------------------------

pub const MAGIC: [4]u8 = .{ 'Z', 'R', 'D', 'C' };
pub const FORMAT_VERSION: u8 = 1;
pub const HEADER_SIZE: usize = 12; // 4 (magic) + 1 (version) + 3 (reserved) + 4 (count)

// ---------------------------------------------------------------------------
// isBinaryDict
// ---------------------------------------------------------------------------

/// Returns true if data begins with the binary dict magic bytes and is at
/// least HEADER_SIZE bytes long.
pub fn isBinaryDict(data: []const u8) bool {
    if (data.len < HEADER_SIZE) return false;
    return std.mem.eql(u8, data[0..4], &MAGIC);
}

// ---------------------------------------------------------------------------
// Tests: isBinaryDict
// ---------------------------------------------------------------------------

const testing = std.testing;

test "isBinaryDict: valid magic" {
    var header: [HEADER_SIZE]u8 = undefined;
    @memcpy(header[0..4], &MAGIC);
    header[4] = FORMAT_VERSION;
    @memset(header[5..], 0);
    try testing.expect(isBinaryDict(&header));
}

test "isBinaryDict: wrong magic" {
    var header: [HEADER_SIZE]u8 = undefined;
    @memcpy(header[0..4], "XXXX");
    @memset(header[4..], 0);
    try testing.expect(!isBinaryDict(&header));
}

test "isBinaryDict: empty data" {
    try testing.expect(!isBinaryDict(&[_]u8{}));
}

test "isBinaryDict: too short" {
    try testing.expect(!isBinaryDict(&[_]u8{ 'Z', 'R' }));
}
