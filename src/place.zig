//! Place module: hydrogen placement geometry and utilities.

pub const geometry = @import("place/geometry.zig");
pub const standard = @import("place/standard.zig");
pub const topology = @import("place/topology.zig");
pub const het = @import("place/het.zig");
pub const placer = @import("place/placer.zig");

// Convenience re-exports
pub const addHydrogens = placer.addHydrogens;
pub const PlacementResult = placer.PlacementResult;
pub const derivePlans = het.derivePlans;

test {
    @import("std").testing.refAllDecls(@This());
}
