//! zreduce: Hydrogen placement and optimization for mmCIF structures.

pub const math = @import("math.zig");
pub const element = @import("element.zig");

test {
    @import("std").testing.refAllDecls(@This());
    // CIF sub-modules (full pub export happens in Task 5)
    _ = @import("cif/char_table.zig");
    _ = @import("cif/tokenizer.zig");
}
