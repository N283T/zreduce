const std = @import("std");
const element = @import("../element.zig");

pub const BondLengthMode = enum {
    xray,
    neutron,
};

pub const OutputIsotope = enum {
    hydrogen,
    deuterium,
};

pub const BondPolicy = struct {
    mode: BondLengthMode = .neutron,
    output_isotope: OutputIsotope = .hydrogen,
};

pub fn adjustedBondLength(mode: BondLengthMode, base_len: f32, parent_type: element.AtomType, h_type: element.AtomType) f32 {
    return switch (mode) {
        .neutron => base_len,
        .xray => switch (parent_type) {
            .O => 0.84,
            .S, .Se => 1.20,
            .N, .Nacc => 0.86,
            else => switch (h_type) {
                .Har, .Ha_p => 0.93,
                else => 0.98,
            },
        },
    };
}

const testing = std.testing;

test "neutron mode preserves base bond length" {
    const base: f32 = 1.01;
    try testing.expectApproxEqAbs(base, adjustedBondLength(.neutron, base, .N, .Hpol), 1e-6);
    try testing.expectApproxEqAbs(base, adjustedBondLength(.neutron, base, .C, .Har), 1e-6);
}

test "xray mode uses parent/type specific target lengths" {
    try testing.expectApproxEqAbs(@as(f32, 0.84), adjustedBondLength(.xray, 1.00, .O, .Hpol), 1e-6);
    try testing.expectApproxEqAbs(@as(f32, 0.86), adjustedBondLength(.xray, 1.00, .N, .Hpol), 1e-6);
    try testing.expectApproxEqAbs(@as(f32, 1.20), adjustedBondLength(.xray, 1.00, .S, .Hpol), 1e-6);
    try testing.expectApproxEqAbs(@as(f32, 0.93), adjustedBondLength(.xray, 1.00, .C, .Har), 1e-6);
    try testing.expectApproxEqAbs(@as(f32, 0.98), adjustedBondLength(.xray, 1.00, .C, .Hpol), 1e-6);
}
