//! JSON log writer: outputs optimization results as JSON.

const std = @import("std");
const testing = std.testing;
const model_mod = @import("../model.zig");
const mover_mod = @import("../optimize/mover.zig");
const bond_policy_mod = @import("../place/bond_policy.zig");

/// Write JSON log of optimization results.
pub fn writeLog(
    out_writer: anytype,
    version: []const u8,
    input_file: []const u8,
    n_hydrogens: u32,
    bond_policy: bond_policy_mod.BondPolicy,
    movers: []const mover_mod.Mover,
    residues: []const model_mod.Residue,
    chains: []const model_mod.Chain,
) !void {
    try out_writer.writeAll("{\n");
    try out_writer.writeAll("  \"version\": \"");
    try writeJsonString(out_writer, version);
    try out_writer.writeAll("\",\n");
    try out_writer.writeAll("  \"input\": \"");
    try writeJsonString(out_writer, input_file);
    try out_writer.writeAll("\",\n");
    try out_writer.print("  \"hydrogens_added\": {d},\n", .{n_hydrogens});
    try out_writer.print("  \"bond_mode\": \"{s}\",\n", .{bondModeStr(bond_policy.mode)});
    try out_writer.print("  \"output_isotope\": \"{s}\",\n", .{outputIsotopeStr(bond_policy.output_isotope)});
    try out_writer.writeAll("  \"movers\": [\n");

    for (movers, 0..) |m, i| {
        const res = residues[m.residue_idx];
        const chain = chains[res.chain_idx];
        const kind_str = moverKindStr(m.kind);

        try out_writer.print("    {{\"residue\": \"{s}.{s}.{d}\", \"type\": \"{s}\", \"orientation\": {d}}}", .{
            chain.labelSlice(),
            res.compIdSlice(),
            res.seq_id,
            kind_str,
            m.best_orientation,
        });

        if (i < movers.len - 1) try out_writer.writeAll(",");
        try out_writer.writeAll("\n");
    }

    try out_writer.writeAll("  ]\n}\n");
}

/// Write a JSON-escaped string (without surrounding quotes).
fn writeJsonString(writer: anytype, s: []const u8) !void {
    for (s) |c| {
        switch (c) {
            '"' => try writer.writeAll("\\\""),
            '\\' => try writer.writeAll("\\\\"),
            '\n' => try writer.writeAll("\\n"),
            '\r' => try writer.writeAll("\\r"),
            '\t' => try writer.writeAll("\\t"),
            else => {
                if (c < 0x20) {
                    try writer.print("\\u{x:0>4}", .{c});
                } else {
                    try writer.writeByte(c);
                }
            },
        }
    }
}

fn moverKindStr(kind: mover_mod.MoverKind) []const u8 {
    return switch (kind) {
        .single_h_rotator => "single_h_rotator",
        .nh3_rotator => "nh3_rotator",
        .methyl_rotator => "methyl_rotator",
        .aromatic_methyl => "aromatic_methyl",
        .amide_flip => "amide_flip",
        .his_flip => "his_flip",
    };
}

fn bondModeStr(mode: bond_policy_mod.BondLengthMode) []const u8 {
    return switch (mode) {
        .xray => "xray",
        .neutron => "neutron",
    };
}

fn outputIsotopeStr(isotope: bond_policy_mod.OutputIsotope) []const u8 {
    return switch (isotope) {
        .hydrogen => "hydrogen",
        .deuterium => "deuterium",
    };
}

test "write JSON log" {
    var buf: std.ArrayListUnmanaged(u8) = .empty;
    defer buf.deinit(testing.allocator);

    try writeLog(buf.writer(testing.allocator), "0.1.0", "test.cif", 42, .{}, &.{}, &.{}, &.{});

    const output = buf.items;
    try testing.expect(std.mem.indexOf(u8, output, "\"version\": \"0.1.0\"") != null);
    try testing.expect(std.mem.indexOf(u8, output, "\"hydrogens_added\": 42") != null);
    try testing.expect(std.mem.indexOf(u8, output, "\"bond_mode\": \"neutron\"") != null);
    try testing.expect(std.mem.indexOf(u8, output, "\"output_isotope\": \"hydrogen\"") != null);
}
