//! Place module: hydrogen placement geometry and utilities.

pub const chemistry = @import("place/chemistry.zig");
pub const geometry = @import("place/geometry.zig");
pub const standard = @import("place/standard.zig");
pub const nucleotide = @import("place/nucleotide.zig");
pub const modified = @import("place/modified.zig");
pub const topology = @import("place/topology.zig");
pub const ccd_derive = @import("place/ccd_derive.zig");
pub const protonation = @import("place/protonation.zig");
pub const placer = @import("place/placer.zig");

// Convenience re-exports
pub const addHydrogens = placer.addHydrogens;
pub const addHydrogensWithConfig = placer.addHydrogensWithConfig;
pub const applyChemistry = placer.applyChemistry;
pub const applyChemistryWithConfig = placer.applyChemistryWithConfig;
pub const PlacementResult = placer.PlacementResult;
pub const PlacementConfig = placer.PlacementConfig;
pub const WaterConfig = placer.WaterConfig;
pub const ProtonationOverrides = protonation.ProtonationOverrides;
pub const derivePlans = ccd_derive.derivePlans;

test {
    @import("std").testing.refAllDecls(@This());
}
