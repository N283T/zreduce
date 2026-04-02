//! CIF module: re-exports all CIF sub-modules and convenience aliases.

pub const types = @import("cif/types.zig");
pub const char_table = @import("cif/char_table.zig");
pub const tokenizer = @import("cif/tokenizer.zig");
pub const parser = @import("cif/parser.zig");
pub const value = @import("cif/value.zig");

// Convenience re-exports
pub const Document = types.Document;
pub const Block = types.Block;
pub const Loop = types.Loop;
pub const readString = parser.readString;
pub const isNull = value.isNull;
pub const asFloat = value.asFloat;
pub const asFloatOr = value.asFloatOr;
pub const asInt = value.asInt;
pub const asIntOr = value.asIntOr;
pub const asString = value.asString;

test {
    @import("std").testing.refAllDecls(@This());
}
