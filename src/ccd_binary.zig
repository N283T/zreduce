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
// PackedAtom  (extern struct, 24 bytes)
// ---------------------------------------------------------------------------

pub const PackedAtom = extern struct {
    name: [4]u8,
    name_len: u8, // stores u4 value; mask with 0x0F on read
    element_symbol: [2]u8,
    charge: i8,
    /// bit0 = leaving, bit1 = aromatic
    flags: u8,
    _pad: [3]u8,
    ideal_x: f32,
    ideal_y: f32,
    ideal_z: f32,

    comptime {
        std.debug.assert(@sizeOf(PackedAtom) == 24);
    }

    pub fn fromCompAtom(a: CompAtom) PackedAtom {
        var flags: u8 = 0;
        if (a.leaving) flags |= 0x01;
        if (a.aromatic) flags |= 0x02;
        return .{
            .name = a.name,
            .name_len = @as(u8, a.name_len),
            .element_symbol = a.element_symbol,
            .charge = a.charge,
            .flags = flags,
            ._pad = .{ 0, 0, 0 },
            .ideal_x = a.ideal_x,
            .ideal_y = a.ideal_y,
            .ideal_z = a.ideal_z,
        };
    }

    pub fn toCompAtom(self: PackedAtom) CompAtom {
        return .{
            .name = self.name,
            .name_len = @truncate(self.name_len & 0x0F),
            .element_symbol = self.element_symbol,
            .charge = self.charge,
            .leaving = (self.flags & 0x01) != 0,
            .aromatic = (self.flags & 0x02) != 0,
            .ideal_x = self.ideal_x,
            .ideal_y = self.ideal_y,
            .ideal_z = self.ideal_z,
        };
    }
};

// ---------------------------------------------------------------------------
// PackedBond  (extern struct, 6 bytes)
// ---------------------------------------------------------------------------

pub const PackedBond = extern struct {
    atom_idx_1: u16,
    atom_idx_2: u16,
    order: u8, // BondOrder tag value
    /// bit0 = aromatic
    flags: u8,

    comptime {
        std.debug.assert(@sizeOf(PackedBond) == 6);
    }

    pub fn fromCompBond(b: CompBond) PackedBond {
        var flags: u8 = 0;
        if (b.aromatic) flags |= 0x01;
        return .{
            .atom_idx_1 = b.atom_idx_1,
            .atom_idx_2 = b.atom_idx_2,
            .order = @intFromEnum(b.order),
            .flags = flags,
        };
    }

    pub fn toCompBond(self: PackedBond) CompBond {
        return .{
            .atom_idx_1 = self.atom_idx_1,
            .atom_idx_2 = self.atom_idx_2,
            .order = @enumFromInt(self.order),
            .aromatic = (self.flags & 0x01) != 0,
        };
    }
};

// ---------------------------------------------------------------------------
// writeDict
// ---------------------------------------------------------------------------

/// Serialize dict to writer in binary format.
pub fn writeDict(writer: anytype, dict: *const ComponentDict) !void {
    // Header
    try writer.writeAll(&MAGIC);
    try writer.writeByte(FORMAT_VERSION);
    try writer.writeAll(&[_]u8{ 0, 0, 0 }); // reserved
    const count: u32 = @intCast(dict.components.count());
    try writer.writeInt(u32, count, .little);

    // Components
    var iter = dict.components.iterator();
    while (iter.next()) |entry| {
        const comp = entry.value_ptr.*;

        // comp_id
        const id_len: u8 = @intCast(comp.comp_id.len);
        try writer.writeByte(id_len);
        try writer.writeAll(comp.comp_id);

        // comp_type
        const type_len: u8 = @intCast(comp.comp_type.len);
        try writer.writeByte(type_len);
        try writer.writeAll(comp.comp_type);

        // atoms
        const atom_count: u16 = @intCast(comp.atoms.len);
        try writer.writeInt(u16, atom_count, .little);
        for (comp.atoms) |atom| {
            const pa = PackedAtom.fromCompAtom(atom);
            try writer.writeAll(std.mem.asBytes(&pa));
        }

        // bonds
        const bond_count: u16 = @intCast(comp.bonds.len);
        try writer.writeInt(u16, bond_count, .little);
        for (comp.bonds) |bond| {
            const pb = PackedBond.fromCompBond(bond);
            try writer.writeAll(std.mem.asBytes(&pb));
        }
    }
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

test "writeDict: empty dictionary" {
    const allocator = testing.allocator;

    var dict = ComponentDict{
        .components = std.StringHashMap(Component).init(allocator),
        .allocator = allocator,
    };
    defer dict.deinit();

    var buf: std.ArrayListUnmanaged(u8) = .empty;
    defer buf.deinit(allocator);
    try writeDict(buf.writer(allocator), &dict);

    try testing.expectEqual(HEADER_SIZE, buf.items.len);
    // Verify component count field is 0
    const count = std.mem.readInt(u32, buf.items[8..12], .little);
    try testing.expectEqual(@as(u32, 0), count);
}
