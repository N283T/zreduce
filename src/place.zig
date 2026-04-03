//! Place module: hydrogen placement geometry and utilities.

pub const chemistry = @import("place/chemistry.zig");
pub const bond_policy = @import("place/bond_policy.zig");
pub const geometry = @import("place/geometry.zig");
pub const standard = @import("place/standard.zig");
pub const nucleotide = @import("place/nucleotide.zig");
pub const modified = @import("place/modified.zig");
pub const topology = @import("place/topology.zig");
pub const ccd_derive = @import("place/ccd_derive.zig");
pub const protonation = @import("place/protonation.zig");
pub const placer = @import("place/placer.zig");
pub const lookup = @import("place/lookup.zig");
pub const terminal_placement = @import("place/terminal.zig");
pub const execute = @import("place/execute.zig");
pub const water = @import("place/water.zig");

// Convenience re-exports
pub const addHydrogens = placer.addHydrogens;
pub const addHydrogensWithConfig = placer.addHydrogensWithConfig;
pub const applyChemistry = placer.applyChemistry;
pub const applyChemistryWithConfig = placer.applyChemistryWithConfig;
pub const PlacementResult = placer.PlacementResult;
pub const PlacementConfig = placer.PlacementConfig;
pub const WaterConfig = placer.WaterConfig;
pub const BondLengthMode = bond_policy.BondLengthMode;
pub const OutputIsotope = bond_policy.OutputIsotope;
pub const BondPolicy = bond_policy.BondPolicy;
pub const ProtonationOverrides = protonation.ProtonationOverrides;
pub const derivePlans = ccd_derive.derivePlans;

test {
    @import("std").testing.refAllDecls(@This());
    _ = @import("place/placer_test.zig");
}
