const std = @import("std");
const zreduce = @import("zreduce");
const build_options = @import("build_options");

const Allocator = std.mem.Allocator;

const Config = struct {
    input_path: []const u8,
    output_path: ?[]const u8 = null,
    dict_path: ?[]const u8 = null,
    json_path: ?[]const u8 = null,
    no_opt: bool = false,
    no_flip: bool = false,
    validate: bool = false,
};

fn parseArgs() ?Config {
    const args = std.process.argsAlloc(std.heap.page_allocator) catch {
        std.debug.print("Fatal: cannot read process arguments\n", .{});
        std.process.exit(2);
    };

    var config = Config{ .input_path = undefined };
    var input_set = false;
    var i: usize = 1;

    while (i < args.len) : (i += 1) {
        const arg = args[i];

        if (std.mem.eql(u8, arg, "-h") or std.mem.eql(u8, arg, "--help")) {
            printUsage(args[0]);
            return null; // exit 0
        } else if (std.mem.eql(u8, arg, "-V") or std.mem.eql(u8, arg, "--version")) {
            std.debug.print("zreduce {s}\n", .{build_options.version});
            return null; // exit 0
        } else if (std.mem.eql(u8, arg, "-o") or std.mem.eql(u8, arg, "--output")) {
            i += 1;
            if (i >= args.len) {
                std.debug.print("Error: {s} requires a path argument\n", .{arg});
                std.process.exit(1);
            }
            config.output_path = args[i];
        } else if (std.mem.eql(u8, arg, "-d") or std.mem.eql(u8, arg, "--dict")) {
            i += 1;
            if (i >= args.len) {
                std.debug.print("Error: {s} requires a path argument\n", .{arg});
                std.process.exit(1);
            }
            config.dict_path = args[i];
        } else if (std.mem.eql(u8, arg, "--json")) {
            i += 1;
            if (i >= args.len) {
                std.debug.print("Error: --json requires a path argument\n", .{});
                std.process.exit(1);
            }
            config.json_path = args[i];
        } else if (std.mem.eql(u8, arg, "--no-opt")) {
            config.no_opt = true;
        } else if (std.mem.eql(u8, arg, "--no-flip")) {
            config.no_flip = true;
        } else if (std.mem.eql(u8, arg, "--validate")) {
            config.validate = true;
        } else if (arg.len > 0 and arg[0] == '-') {
            std.debug.print("Error: unknown option '{s}'\n", .{arg});
            std.process.exit(1);
        } else {
            if (input_set) {
                std.debug.print("Error: unexpected positional argument '{s}'\n", .{arg});
                std.process.exit(1);
            }
            config.input_path = arg;
            input_set = true;
        }
    }

    if (!input_set) {
        std.debug.print("Error: missing input file\n", .{});
        printUsage(args[0]);
        std.process.exit(1);
    }

    return config;
}

fn readFile(allocator: Allocator, path: []const u8) ![]const u8 {
    const file = try std.fs.cwd().openFile(path, .{});
    defer file.close();
    return try file.readToEndAlloc(allocator, 1024 * 1024 * 1024); // 1GB max
}

fn printUsage(program_name: []const u8) void {
    std.debug.print(
        \\zreduce {s} - Hydrogen placement for mmCIF structures
        \\
        \\USAGE:
        \\    {s} [OPTIONS] <input.cif> [-o output.cif]
        \\
        \\OPTIONS:
        \\    -h, --help         Show this help message
        \\    -V, --version      Show version
        \\    -d, --dict PATH    Path to components.cif for CCD HET groups
        \\    -o, --output PATH  Output file (default: stdout)
        \\    --json PATH        Write JSON log to file
        \\    --no-opt           Skip optimization (placement only)
        \\    --no-flip          Disable Asn/Gln/His flips
        \\    --validate         Print detailed validation diagnostics
        \\
    , .{ build_options.version, program_name });
}

pub fn main() !void {
    var gpa = std.heap.GeneralPurposeAllocator(.{}){};
    defer _ = gpa.deinit();
    const allocator = gpa.allocator();

    // 1. Parse args
    const config = parseArgs() orelse return;

    // 2. Read input mmCIF file
    const source = readFile(allocator, config.input_path) catch |err| {
        std.debug.print("Error: cannot read '{s}': {s}\n", .{ config.input_path, @errorName(err) });
        std.process.exit(1);
    };
    defer allocator.free(source);

    // 3. Parse CIF document (for preserving non-atom_site categories in output)
    var doc = zreduce.cif.readString(allocator, source) catch |err| {
        std.debug.print("Error: failed to parse CIF: {s}\n", .{@errorName(err)});
        std.process.exit(1);
    };
    defer doc.deinit();

    // 4. Extract model from CIF
    var mdl = zreduce.mmcif.parseModel(allocator, source) catch |err| {
        std.debug.print("Error: failed to parse mmCIF: {s}\n", .{@errorName(err)});
        std.process.exit(1);
    };
    defer mdl.deinit();

    // 5. Load CCD dictionary (optional)
    var ccd_dict: ?zreduce.ccd.ComponentDict = null;
    if (config.dict_path) |dict_path| {
        const dict_source = readFile(allocator, dict_path) catch |err| {
            std.debug.print("Error: cannot read dictionary '{s}': {s}\n", .{ dict_path, @errorName(err) });
            std.process.exit(1);
        };
        defer allocator.free(dict_source);
        ccd_dict = zreduce.ccd.parseComponentDict(allocator, dict_source) catch |err| {
            std.debug.print("Error: failed to parse CCD dictionary: {s}\n", .{@errorName(err)});
            std.process.exit(1);
        };
    }
    defer if (ccd_dict) |*d| d.deinit();

    // 5.5 Apply chemistry annotations to standard-residue heavy atoms
    zreduce.place.applyChemistry(&mdl);

    // 6. Place hydrogens
    const initial_count = mdl.atoms.items.len;
    const place_result = zreduce.place.addHydrogens(&mdl, if (ccd_dict) |*d| d else null) catch |err| {
        std.debug.print("Error: hydrogen placement failed: {s}\n", .{@errorName(err)});
        std.process.exit(1);
    };
    const n_added: u32 = @intCast(mdl.atoms.items.len - initial_count);

    // 7. Optimize (unless --no-opt)
    var movers: []zreduce.optimize.Mover = &.{};
    var movers_owned = false;
    defer {
        for (0..movers.len) |i| movers[i].deinit();
        if (movers_owned) allocator.free(movers);
    }

    if (!config.no_opt) {
        const gen_result = zreduce.optimize.generateMovers(allocator, &mdl, config.no_flip, if (ccd_dict) |*d| d else null) catch |err| {
            std.debug.print("Error: mover generation failed: {s}\n", .{@errorName(err)});
            std.process.exit(1);
        };
        movers = gen_result.movers;
        movers_owned = true;

        if (gen_result.n_skipped > 0) {
            std.debug.print("  Mover generation: {d} skipped (missing atoms or incomplete groups)\n", .{gen_result.n_skipped});
        }

        if (movers.len > 0) {
            const opt_config = zreduce.optimize.OptConfig{};
            const opt_result = zreduce.optimize.optimizer.optimize(allocator, movers, &mdl, opt_config) catch |err| {
                std.debug.print("Error: optimization failed: {s}\n", .{@errorName(err)});
                std.process.exit(1);
            };
            std.debug.print("  Movers: {d} ({d} singletons, {d} brute-force, {d} greedy)\n", .{
                movers.len,
                opt_result.n_singletons,
                opt_result.n_brute_force,
                opt_result.n_vertex_cut,
            });
        }
    }

    // 7.5 Validate model (always run, report issues — before sentinel removal)
    {
        var validation = zreduce.validate.validateModel(allocator, &mdl) catch |err| {
            std.debug.print("Error: validation failed: {s}\n", .{@errorName(err)});
            std.process.exit(1);
        };
        defer validation.deinit();

        if (!validation.ok()) {
            std.debug.print("  Validation: {d} issue(s) found\n", .{validation.issues.len});
            if (config.validate) {
                zreduce.validate.reportIssues(validation.issues, &mdl);
            }
        }
    }

    // 7.6 Remove absent H atoms (flipper sentinels) from the model
    for (mdl.atoms.items) |*atom| {
        if (zreduce.optimize.mover.isAbsentH(atom.*)) {
            atom.is_added = false; // mark as not placed — writer will skip
        }
    }

    // 8. Write output (preserving original CIF categories)
    var out_buf: [4096]u8 = undefined;
    if (config.output_path) |out_path| {
        const file = std.fs.cwd().createFile(out_path, .{}) catch |err| {
            std.debug.print("Error: cannot create output file '{s}': {s}\n", .{ out_path, @errorName(err) });
            std.process.exit(1);
        };
        defer file.close();
        var fw = file.writer(&out_buf);
        try zreduce.writer.mmcif_writer.writeWithDocument(&fw.interface, &mdl, &doc);
        try fw.interface.flush();
    } else {
        const stdout = std.fs.File.stdout();
        var sw = stdout.writer(&out_buf);
        try zreduce.writer.mmcif_writer.writeWithDocument(&sw.interface, &mdl, &doc);
        try sw.interface.flush();
    }

    // 9. Write JSON log (optional)
    if (config.json_path) |json_path| {
        var json_buf: [4096]u8 = undefined;
        const file = std.fs.cwd().createFile(json_path, .{}) catch |err| {
            std.debug.print("Error: cannot create JSON log '{s}': {s}\n", .{ json_path, @errorName(err) });
            std.process.exit(1);
        };
        defer file.close();
        var jw = file.writer(&json_buf);
        try zreduce.writer.json_writer.writeLog(
            &jw.interface,
            build_options.version,
            config.input_path,
            n_added,
            movers,
            mdl.residues.items,
            mdl.chains.items,
        );
        try jw.interface.flush();
    }

    // 10. Report summary to stderr
    std.debug.print("zreduce: placed {d} H atoms on {d} residues ({d} skipped)\n", .{
        place_result.n_placed,
        place_result.n_residues,
        place_result.n_skipped,
    });
}
