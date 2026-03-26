//! Place module: hydrogen placement geometry and utilities.

pub const geometry = @import("place/geometry.zig");
pub const standard = @import("place/standard.zig");

test {
    @import("std").testing.refAllDecls(@This());
}
