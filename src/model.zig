//! Molecular model data structures.
//! Re-exports all model submodules and provides convenience type aliases.

pub const atom = @import("model/atom.zig");
pub const residue = @import("model/residue.zig");
pub const chain = @import("model/chain.zig");
pub const bond = @import("model/bond.zig");
pub const model = @import("model/model.zig");
pub const neighbor = @import("model/neighbor.zig");

// Convenience re-exports
pub const Atom = model.Atom;
pub const Residue = model.Residue;
pub const Chain = model.Chain;
pub const Bond = model.Bond;
pub const Model = model.Model;
pub const CellList = neighbor.CellList;

test {
    @import("std").testing.refAllDecls(@This());
}
