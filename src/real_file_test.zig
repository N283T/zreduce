//! End-to-end tests with real PDB structures.
//!
//! Tests chain gap handling, insertion code handling, and full pipeline
//! correctness on deposited PDB structures (examples/data/*.cif.gz).
//!
//! Structures used:
//!   1rqf — 4 chains, mid-chain gaps in B/C (CK2 beta subunit)
//!   3rk2 — 4 chains, internal gaps in B/F (SNARE complex)
//!   6fys — 1 chain, 8 insertion code positions (nanobody)
//!   2cf8 — 3 chains, chymotrypsinogen insertion codes (thrombin)
//!   2hnt — 4 chains, both gaps AND insertion codes (gamma-thrombin)

const std = @import("std");
const testing = std.testing;
const zreduce = @import("root.zig");
const mmcif = zreduce.mmcif;
const place = zreduce.place;
const optimize = zreduce.optimize;
const cif = zreduce.cif;

/// Helper: run the full pipeline on a .cif.gz file, return model and placement stats.
/// Caller owns model and must deinit.
const PipelineResult = struct {
    model: zreduce.model.Model,
    n_placed: u32,
    n_residues: u32,
    n_skipped_existing: u32,
    n_skipped_inter_residue: u32,
    n_skipped_missing_ref: u32,
    n_movers: usize,

    pub fn deinit(self: *PipelineResult) void {
        self.model.deinit();
    }
};

/// Run the full pipeline without an external CCD dictionary.
/// Structures that need CCD-derived placement (non-standard ligands)
/// will have those residues skipped, which is acceptable for these tests.
fn runPipeline(allocator: std.mem.Allocator, path: []const u8) !PipelineResult {
    // 1. Read and decompress
    const source = try zreduce.run.readFile(allocator, path);
    defer allocator.free(source);

    // 2. Parse CIF document
    var doc = try cif.readString(allocator, source);
    defer doc.deinit();
    const block = &doc.blocks.items[0];

    // 3. Parse model
    var mdl = try mmcif.parseModel(allocator, source);
    errdefer mdl.deinit();

    // 4. Inline components
    var inline_dict = try mmcif.parseInlineComponents(allocator, block);
    defer if (inline_dict) |*d| d.deinit();

    // 5. Chemistry + bonds
    place.applyChemistry(&mdl);
    var atom_lookup = try mmcif.buildAtomLookup(allocator, block);
    defer atom_lookup.deinit();
    try mmcif.parseStructConn(&mdl, block, &atom_lookup);
    try mmcif.parseBranchLinks(allocator, &mdl, block, &atom_lookup);
    mmcif.flagLeavingAtoms(&mdl, if (inline_dict) |*d| d else null, null);

    // 6. Place hydrogens
    const place_result = try place.addHydrogens(&mdl, null, if (inline_dict) |*d| d else null);

    // 7. Optimize
    const gen_result = try optimize.generateMovers(
        allocator,
        &mdl,
        false,
        null,
        if (inline_dict) |*d| d else null,
        null,
    );
    var movers = gen_result.movers;
    defer {
        for (0..movers.len) |i| movers[i].deinit();
        allocator.free(movers);
    }

    if (movers.len > 0) {
        _ = try optimize.optimizer.optimize(allocator, movers, &mdl, .{ .n_threads = 1 });
    }

    // 8. Mark absent H
    for (mdl.atoms.items) |*atom| {
        if (optimize.mover.isAbsentH(atom.*)) atom.is_added = false;
    }

    // 9. Validate (should not crash on any real structure)
    var validation = try zreduce.validate.validateModel(allocator, &mdl);
    validation.deinit();

    return .{
        .model = mdl,
        .n_placed = place_result.n_placed,
        .n_residues = place_result.n_residues,
        .n_skipped_existing = place_result.n_skipped_existing,
        .n_skipped_inter_residue = place_result.n_skipped_inter_residue,
        .n_skipped_missing_ref = place_result.n_skipped_missing_ref,
        .n_movers = movers.len,
    };
}

/// Count added hydrogen atoms in the final model.
fn countFinalH(mdl: *const zreduce.model.Model) u32 {
    var count: u32 = 0;
    for (mdl.atoms.items) |atom| {
        if (atom.is_added and atom.is_hydrogen) count += 1;
    }
    return count;
}

/// Count residues with is_chain_break_before = true.
fn countChainBreaks(mdl: *const zreduce.model.Model) u32 {
    var count: u32 = 0;
    for (mdl.residues.items) |res| {
        if (res.is_chain_break_before) count += 1;
    }
    return count;
}

/// Check that the output is valid mmCIF (write and re-parse).
fn verifyOutputRoundTrip(allocator: std.mem.Allocator, mdl: *const zreduce.model.Model) !void {
    var buf: std.ArrayListUnmanaged(u8) = .empty;
    defer buf.deinit(allocator);
    try zreduce.writer.mmcif_writer.write(buf.writer(allocator), mdl, "TEST");

    // Re-parse to verify valid CIF
    var doc2 = try cif.readString(allocator, buf.items);
    defer doc2.deinit();
    try testing.expect(doc2.blocks.items.len > 0);
    const loop = doc2.blocks.items[0].findLoop("_atom_site.Cartn_x");
    try testing.expect(loop != null);
}

/// Verify added H bond lengths. Allows up to 1% outliers for real structures
/// (strained geometry, low-occupancy altlocs, etc.).
fn verifyBondLengths(mdl: *const zreduce.model.Model) !void {
    var n_checked: u32 = 0;
    var n_bad: u32 = 0;

    for (mdl.atoms.items) |atom| {
        if (!atom.is_added or !atom.is_hydrogen) continue;

        // Find nearest heavy atom in same residue as proxy for bond length.
        // Not exact (nearest != bonded parent) but sufficient for smoke testing.
        var min_dist: f32 = std.math.inf(f32);
        for (mdl.atoms.items) |other| {
            if (other.is_hydrogen) continue;
            if (other.residue_idx != atom.residue_idx) continue;
            const dist = atom.pos.distance(other.pos);
            if (dist < min_dist) min_dist = dist;
        }

        n_checked += 1;
        if (min_dist < 0.5 or min_dist > 1.6) {
            n_bad += 1;
        }
    }

    // Allow up to 1% outliers (strained real-world geometries)
    const max_bad = @max(n_checked / 100, 3);
    if (n_bad > max_bad) {
        std.debug.print(
            "\nToo many bad H bond lengths: {d}/{d} (max allowed: {d})\n",
            .{ n_bad, n_checked, max_bad },
        );
        return error.TestUnexpectedResult;
    }
}

// ============================================================
// Chain gap tests
// ============================================================

test "real-file: 1rqf — chain gaps, 4-chain CK2 beta" {
    // Paths are relative to project root (zig build test CWD).
    var result = try runPipeline(testing.allocator, "examples/data/1rqf.cif.gz");
    defer result.deinit();

    // 1rqf: ~1418 residues, ~10467 H placed (measured)
    try testing.expect(result.n_placed > 5000);
    try testing.expect(result.n_placed < 20000);
    try testing.expect(result.n_residues > 500);
    try testing.expect(result.n_movers > 0);

    // Mid-chain gaps in B/C
    const n_breaks = countChainBreaks(&result.model);
    try testing.expect(n_breaks > 0);

    const final_h = countFinalH(&result.model);
    try testing.expect(final_h > 5000);
    try testing.expect(final_h < 20000);

    try verifyBondLengths(&result.model);
    try verifyOutputRoundTrip(testing.allocator, &result.model);
}

test "real-file: 3rk2 — chain gaps, SNARE complex" {
    var result = try runPipeline(testing.allocator, "examples/data/3rk2.cif.gz");
    defer result.deinit();

    // 3rk2: ~505 residues, ~3734 H placed (measured)
    try testing.expect(result.n_placed > 1500);
    try testing.expect(result.n_placed < 8000);
    try testing.expect(result.n_residues > 200);
    try testing.expect(result.n_movers > 0);

    const n_breaks = countChainBreaks(&result.model);
    try testing.expect(n_breaks > 0);

    const final_h = countFinalH(&result.model);
    try testing.expect(final_h > 1500);
    try testing.expect(final_h < 8000);

    try verifyBondLengths(&result.model);
    try verifyOutputRoundTrip(testing.allocator, &result.model);
}

// ============================================================
// Insertion code tests
// ============================================================

test "real-file: 6fys — insertion codes, nanobody" {
    var result = try runPipeline(testing.allocator, "examples/data/6fys.cif.gz");
    defer result.deinit();

    // 6fys: ~1218 residues (with HOH), ~7320 H placed (measured)
    try testing.expect(result.n_placed > 3000);
    try testing.expect(result.n_placed < 15000);
    try testing.expect(result.n_residues > 100);
    try testing.expect(result.n_movers > 0);

    const final_h = countFinalH(&result.model);
    try testing.expect(final_h > 3000);
    try testing.expect(final_h < 15000);

    try verifyBondLengths(&result.model);
    try verifyOutputRoundTrip(testing.allocator, &result.model);
}

test "real-file: 2cf8 — insertion codes, thrombin 3-chain" {
    var result = try runPipeline(testing.allocator, "examples/data/2cf8.cif.gz");
    defer result.deinit();

    // 2cf8: ~684 residues, ~151 H placed (many HOH skipped without CCD)
    try testing.expect(result.n_placed > 50);
    try testing.expect(result.n_placed < 1000);
    try testing.expect(result.n_residues > 100);

    // Multi-chain: should have at least 3 chains
    try testing.expect(result.model.chains.items.len >= 3);

    const final_h = countFinalH(&result.model);
    try testing.expect(final_h > 50);
    try testing.expect(final_h < 1000);

    try verifyBondLengths(&result.model);
    try verifyOutputRoundTrip(testing.allocator, &result.model);
}

// ============================================================
// Both gaps and insertion codes
// ============================================================

test "real-file: 2hnt — gaps + insertion codes, gamma-thrombin" {
    var result = try runPipeline(testing.allocator, "examples/data/2hnt.cif.gz");
    defer result.deinit();

    // 2hnt: ~448 residues, ~2100 H placed (measured)
    try testing.expect(result.n_placed > 1000);
    try testing.expect(result.n_placed < 5000);
    try testing.expect(result.n_residues > 200);
    try testing.expect(result.n_movers > 0);

    // Should have chain breaks
    const n_breaks = countChainBreaks(&result.model);
    try testing.expect(n_breaks > 0);

    // Multi-chain
    try testing.expect(result.model.chains.items.len >= 2);

    const final_h = countFinalH(&result.model);
    try testing.expect(final_h > 1000);
    try testing.expect(final_h < 5000);

    try verifyBondLengths(&result.model);
    try verifyOutputRoundTrip(testing.allocator, &result.model);
}
