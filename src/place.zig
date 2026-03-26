//! Place module: hydrogen placement geometry and utilities.

pub const geometry = @import("place/geometry.zig");

test {
    @import("std").testing.refAllDecls(@This());
}
