const std = @import("std");
const zreduce = @import("zreduce");
const build_options = @import("build_options");

const Allocator = std.mem.Allocator;

const RunConfig = struct {
    input_path: []const u8,
    output_path: ?[]const u8 = null,
    dict_path: ?[]const u8 = null,
    json_path: ?[]const u8 = null,
    no_opt: bool = false,
    no_flip: bool = false,
    validate: bool = false,
    water: zreduce.place.WaterConfig = .{},
    bond_policy: zreduce.place.BondPolicy = .{},
    protonation_path: ?[]const u8 = null,
    fix_path: ?[]const u8 = null,
    dump_movers_path: ?[]const u8 = null,
};

fn parseBondModeValue(s: []const u8) ?zreduce.place.BondLengthMode {
    if (std.ascii.eqlIgnoreCase(s, "neutron")) return .neutron;
    if (std.ascii.eqlIgnoreCase(s, "xray")) return .xray;
    return null;
}

fn parseOutputIsotopeValue(s: []const u8) ?zreduce.place.OutputIsotope {
    if (std.ascii.eqlIgnoreCase(s, "hydrogen") or std.ascii.eqlIgnoreCase(s, "h")) return .hydrogen;
    if (std.ascii.eqlIgnoreCase(s, "deuterium") or std.ascii.eqlIgnoreCase(s, "d")) return .deuterium;
    return null;
}

fn parseRunArgs(args: []const []const u8) ?RunConfig {
    var config = RunConfig{ .input_path = undefined };
    var input_set = false;
    var i: usize = 0;

    while (i < args.len) : (i += 1) {
        const arg = args[i];

        if (std.mem.eql(u8, arg, "-h") or std.mem.eql(u8, arg, "--help")) {
            printRunUsage();
            return null;
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
        } else if (std.mem.eql(u8, arg, "--protonation")) {
            i += 1;
            if (i >= args.len) {
                std.debug.print("Error: --protonation requires a path argument\n", .{});
                std.process.exit(1);
            }
            config.protonation_path = args[i];
        } else if (std.mem.eql(u8, arg, "--fix")) {
            i += 1;
            if (i >= args.len) {
                std.debug.print("Error: --fix requires a path argument\n", .{});
                std.process.exit(1);
            }
            config.fix_path = args[i];
        } else if (std.mem.eql(u8, arg, "--dump-movers")) {
            i += 1;
            if (i >= args.len) {
                std.debug.print("Error: --dump-movers requires a path argument\n", .{});
                std.process.exit(1);
            }
            config.dump_movers_path = args[i];
        } else if (std.mem.eql(u8, arg, "--no-opt")) {
            config.no_opt = true;
        } else if (std.mem.eql(u8, arg, "--no-flip")) {
            config.no_flip = true;
        } else if (std.mem.eql(u8, arg, "--validate")) {
            config.validate = true;
        } else if (std.mem.eql(u8, arg, "--water")) {
            config.water.enabled = true;
        } else if (std.mem.eql(u8, arg, "--no-water")) {
            config.water.enabled = false;
        } else if (std.mem.eql(u8, arg, "--water-phantom")) {
            config.water.enabled = true;
            config.water.phantom = true;
        } else if (std.mem.eql(u8, arg, "--water-occ-cutoff")) {
            i += 1;
            if (i >= args.len) {
                std.debug.print("Error: --water-occ-cutoff requires a numeric argument\n", .{});
                std.process.exit(1);
            }
            const val = std.fmt.parseFloat(f32, args[i]) catch {
                std.debug.print("Error: invalid water occupancy cutoff '{s}'\n", .{args[i]});
                std.process.exit(1);
            };
            if (!std.math.isFinite(val) or val < 0.0 or val > 1.0) {
                std.debug.print("Error: --water-occ-cutoff must be between 0.0 and 1.0\n", .{});
                std.process.exit(1);
            }
            config.water.occupancy_cutoff = val;
            config.water.enabled = true;
        } else if (std.mem.eql(u8, arg, "--water-b-cutoff")) {
            i += 1;
            if (i >= args.len) {
                std.debug.print("Error: --water-b-cutoff requires a numeric argument\n", .{});
                std.process.exit(1);
            }
            const val = std.fmt.parseFloat(f32, args[i]) catch {
                std.debug.print("Error: invalid water B-factor cutoff '{s}'\n", .{args[i]});
                std.process.exit(1);
            };
            if (!std.math.isFinite(val) or val < 0.0) {
                std.debug.print("Error: --water-b-cutoff must be a non-negative number\n", .{});
                std.process.exit(1);
            }
            config.water.b_factor_cutoff = val;
            config.water.enabled = true;
        } else if (std.mem.eql(u8, arg, "--bond-mode")) {
            i += 1;
            if (i >= args.len) {
                std.debug.print("Error: --bond-mode requires 'neutron' or 'xray'\n", .{});
                std.process.exit(1);
            }
            const parsed = parseBondModeValue(args[i]) orelse {
                std.debug.print("Error: invalid --bond-mode '{s}' (expected 'neutron' or 'xray')\n", .{args[i]});
                std.process.exit(1);
            };
            config.bond_policy.mode = parsed;
        } else if (std.mem.eql(u8, arg, "--isotope")) {
            i += 1;
            if (i >= args.len) {
                std.debug.print("Error: --isotope requires 'hydrogen'/'h' or 'deuterium'/'d'\n", .{});
                std.process.exit(1);
            }
            const parsed = parseOutputIsotopeValue(args[i]) orelse {
                std.debug.print("Error: invalid --isotope '{s}' (expected hydrogen|h|deuterium|d)\n", .{args[i]});
                std.process.exit(1);
            };
            config.bond_policy.output_isotope = parsed;
        } else if (std.mem.eql(u8, arg, "--deuterium")) {
            config.bond_policy.output_isotope = .deuterium;
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
        printRunUsage();
        std.process.exit(1);
    }

    return config;
}

fn printUsage(program_name: []const u8) void {
    std.debug.print(
        \\zreduce {s} - Hydrogen placement for mmCIF structures
        \\
        \\USAGE:
        \\    {s} <command> [OPTIONS] <args>
        \\
        \\COMMANDS:
        \\    run      Process a single mmCIF file
        \\    batch    Process all mmCIF files in a directory
        \\
        \\GLOBAL OPTIONS:
        \\    -h, --help       Show this help message
        \\    -V, --version    Show version
        \\
        \\Use '{s} <command> --help' for more information.
        \\
    , .{ build_options.version, program_name, program_name });
}

fn printRunUsage() void {
    std.debug.print(
        \\USAGE:
        \\    zreduce run [OPTIONS] <input.cif>
        \\
        \\OPTIONS:
        \\    -h, --help         Show this help message
        \\    -d, --dict PATH    CCD dictionary
        \\    -o, --output PATH  Output file (default: stdout)
        \\    --json PATH        Write JSON log to file
        \\    --protonation PATH  Residue protonation override file
        \\    --fix PATH         Force mover states from control file
        \\    --dump-movers PATH  Write available mover IDs/states to file
        \\    --no-opt           Skip optimization
        \\    --no-flip          Disable Asn/Gln/His flips
        \\    --validate         Print validation diagnostics
        \\    --water            Add water hydrogens (default: off)
        \\    --water-phantom    Allow zero-occupancy phantom water H when orientation is underdetermined
        \\    --water-occ-cutoff N  Skip waters with occupancy below N (default: 0.66)
        \\    --water-b-cutoff N    Skip waters with B-factor above N (default: 40.0)
        \\    --bond-mode MODE   Bond-length mode: neutron|xray (default: neutron)
        \\    --isotope NAME     Output isotope for added H: hydrogen|h|deuterium|d (default: hydrogen)
        \\    --deuterium        Shortcut for --isotope deuterium
        \\
    , .{});
}

pub fn main() !void {
    var gpa = std.heap.GeneralPurposeAllocator(.{}){};
    defer _ = gpa.deinit();
    const allocator = gpa.allocator();

    const args = std.process.argsAlloc(std.heap.page_allocator) catch {
        std.debug.print("Fatal: cannot read process arguments\n", .{});
        std.process.exit(2);
    };

    if (args.len < 2) {
        printUsage(args[0]);
        std.process.exit(1);
    }

    const subcmd = args[1];

    if (std.mem.eql(u8, subcmd, "-h") or std.mem.eql(u8, subcmd, "--help")) {
        printUsage(args[0]);
        return;
    }
    if (std.mem.eql(u8, subcmd, "-V") or std.mem.eql(u8, subcmd, "--version")) {
        std.debug.print("zreduce {s}\n", .{build_options.version});
        return;
    }

    if (std.mem.eql(u8, subcmd, "run")) {
        runSubcommand(allocator, args[2..]);
    } else if (std.mem.eql(u8, subcmd, "batch")) {
        batchSubcommand(allocator, args[2..]);
    } else {
        std.debug.print("Error: unknown subcommand '{s}'\n", .{subcmd});
        printUsage(args[0]);
        std.process.exit(1);
    }
}

fn parseBatchArgs(args: []const []const u8) ?zreduce.batch.BatchConfig {
    var config = zreduce.batch.BatchConfig{ .input_dir = undefined };
    var input_set = false;
    var i: usize = 0;

    while (i < args.len) : (i += 1) {
        const arg = args[i];
        if (std.mem.eql(u8, arg, "-h") or std.mem.eql(u8, arg, "--help")) {
            printBatchUsage();
            return null;
        } else if (std.mem.eql(u8, arg, "-o") or std.mem.eql(u8, arg, "--output")) {
            i += 1;
            if (i >= args.len) {
                std.debug.print("Error: {s} requires a path\n", .{arg});
                std.process.exit(1);
            }
            config.output_dir = args[i];
        } else if (std.mem.eql(u8, arg, "-d") or std.mem.eql(u8, arg, "--dict")) {
            i += 1;
            if (i >= args.len) {
                std.debug.print("Error: {s} requires a path\n", .{arg});
                std.process.exit(1);
            }
            config.dict_path = args[i];
        } else if (std.mem.eql(u8, arg, "--jsonl")) {
            i += 1;
            if (i >= args.len) {
                std.debug.print("Error: --jsonl requires a path\n", .{});
                std.process.exit(1);
            }
            config.jsonl_path = args[i];
        } else if (std.mem.eql(u8, arg, "--protonation")) {
            i += 1;
            if (i >= args.len) {
                std.debug.print("Error: --protonation requires a path\n", .{});
                std.process.exit(1);
            }
            config.protonation_path = args[i];
        } else if (std.mem.eql(u8, arg, "--fix")) {
            i += 1;
            if (i >= args.len) {
                std.debug.print("Error: --fix requires a path\n", .{});
                std.process.exit(1);
            }
            config.fix_path = args[i];
        } else if (std.mem.eql(u8, arg, "-j") or std.mem.eql(u8, arg, "--threads")) {
            i += 1;
            if (i >= args.len) {
                std.debug.print("Error: {s} requires a number\n", .{arg});
                std.process.exit(1);
            }
            config.n_threads = std.fmt.parseInt(u32, args[i], 10) catch {
                std.debug.print("Error: invalid thread count '{s}'\n", .{args[i]});
                std.process.exit(1);
            };
        } else if (std.mem.eql(u8, arg, "--no-opt")) {
            config.no_opt = true;
        } else if (std.mem.eql(u8, arg, "--no-flip")) {
            config.no_flip = true;
        } else if (std.mem.eql(u8, arg, "--quiet")) {
            config.quiet = true;
        } else if (std.mem.eql(u8, arg, "--water")) {
            config.water.enabled = true;
        } else if (std.mem.eql(u8, arg, "--no-water")) {
            config.water.enabled = false;
        } else if (std.mem.eql(u8, arg, "--water-phantom")) {
            config.water.enabled = true;
            config.water.phantom = true;
        } else if (std.mem.eql(u8, arg, "--water-occ-cutoff")) {
            i += 1;
            if (i >= args.len) {
                std.debug.print("Error: --water-occ-cutoff requires a numeric argument\n", .{});
                std.process.exit(1);
            }
            const val = std.fmt.parseFloat(f32, args[i]) catch {
                std.debug.print("Error: invalid water occupancy cutoff '{s}'\n", .{args[i]});
                std.process.exit(1);
            };
            if (!std.math.isFinite(val) or val < 0.0 or val > 1.0) {
                std.debug.print("Error: --water-occ-cutoff must be between 0.0 and 1.0\n", .{});
                std.process.exit(1);
            }
            config.water.occupancy_cutoff = val;
            config.water.enabled = true;
        } else if (std.mem.eql(u8, arg, "--water-b-cutoff")) {
            i += 1;
            if (i >= args.len) {
                std.debug.print("Error: --water-b-cutoff requires a numeric argument\n", .{});
                std.process.exit(1);
            }
            const val = std.fmt.parseFloat(f32, args[i]) catch {
                std.debug.print("Error: invalid water B-factor cutoff '{s}'\n", .{args[i]});
                std.process.exit(1);
            };
            if (!std.math.isFinite(val) or val < 0.0) {
                std.debug.print("Error: --water-b-cutoff must be a non-negative number\n", .{});
                std.process.exit(1);
            }
            config.water.b_factor_cutoff = val;
            config.water.enabled = true;
        } else if (std.mem.eql(u8, arg, "--bond-mode")) {
            i += 1;
            if (i >= args.len) {
                std.debug.print("Error: --bond-mode requires 'neutron' or 'xray'\n", .{});
                std.process.exit(1);
            }
            const parsed = parseBondModeValue(args[i]) orelse {
                std.debug.print("Error: invalid --bond-mode '{s}' (expected 'neutron' or 'xray')\n", .{args[i]});
                std.process.exit(1);
            };
            config.bond_policy.mode = parsed;
        } else if (std.mem.eql(u8, arg, "--isotope")) {
            i += 1;
            if (i >= args.len) {
                std.debug.print("Error: --isotope requires 'hydrogen'/'h' or 'deuterium'/'d'\n", .{});
                std.process.exit(1);
            }
            const parsed = parseOutputIsotopeValue(args[i]) orelse {
                std.debug.print("Error: invalid --isotope '{s}' (expected hydrogen|h|deuterium|d)\n", .{args[i]});
                std.process.exit(1);
            };
            config.bond_policy.output_isotope = parsed;
        } else if (std.mem.eql(u8, arg, "--deuterium")) {
            config.bond_policy.output_isotope = .deuterium;
        } else if (arg.len > 0 and arg[0] == '-') {
            std.debug.print("Error: unknown option '{s}'\n", .{arg});
            std.process.exit(1);
        } else {
            if (input_set) {
                std.debug.print("Error: unexpected argument '{s}'\n", .{arg});
                std.process.exit(1);
            }
            config.input_dir = arg;
            input_set = true;
        }
    }

    if (!input_set) {
        std.debug.print("Error: missing input directory\n", .{});
        printBatchUsage();
        std.process.exit(1);
    }

    config.json_version = build_options.version;
    return config;
}

fn printBatchUsage() void {
    std.debug.print(
        \\USAGE:
        \\    zreduce batch [OPTIONS] <input_dir>
        \\
        \\OPTIONS:
        \\    -h, --help         Show this help message
        \\    -d, --dict PATH    CCD dictionary (loaded once)
        \\    -o, --output PATH  Output directory (default: <input>_reduced/)
        \\    -j, --threads N    Thread count (default: auto-detect)
        \\    --jsonl PATH       Aggregated JSONL log file
        \\    --protonation PATH  Residue protonation override file
        \\    --fix PATH         Force mover states from control file
        \\    --no-opt           Skip optimization
        \\    --no-flip          Disable flips
        \\    --quiet            Suppress progress output
        \\    --water            Add water hydrogens (default: off)
        \\    --water-phantom    Allow zero-occupancy phantom water H when orientation is underdetermined
        \\    --water-occ-cutoff N  Skip waters with occupancy below N (default: 0.66)
        \\    --water-b-cutoff N    Skip waters with B-factor above N (default: 40.0)
        \\    --bond-mode MODE   Bond-length mode: neutron|xray (default: neutron)
        \\    --isotope NAME     Output isotope for added H: hydrogen|h|deuterium|d (default: hydrogen)
        \\    --deuterium        Shortcut for --isotope deuterium
        \\
    , .{});
}

fn batchSubcommand(allocator: Allocator, args: []const []const u8) void {
    const config = parseBatchArgs(args) orelse return;
    zreduce.batch.run(allocator, config) catch |err| {
        switch (err) {
            error.SomeFilesFailed, error.NoFilesFound => {},
            else => std.debug.print("Error: batch processing failed: {s}\n", .{@errorName(err)}),
        }
        std.process.exit(1);
    };
}

fn runSubcommand(allocator: Allocator, args: []const []const u8) void {
    const config = parseRunArgs(args) orelse return;

    // Load CCD dictionary (once, before processFile)
    var ccd_dict: ?zreduce.ccd.ComponentDict = null;
    if (config.dict_path) |dict_path| {
        const dict_source = zreduce.run.readFile(allocator, dict_path) catch |err| {
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

    const proc_config = zreduce.run.ProcessConfig{
        .input_path = config.input_path,
        .output_path = config.output_path,
        .dict = if (ccd_dict) |*d| d else null,
        .json_path = config.json_path,
        .json_version = build_options.version,
        .no_opt = config.no_opt,
        .no_flip = config.no_flip,
        .validate_flag = config.validate,
        .water = config.water,
        .bond_policy = config.bond_policy,
        .protonation_path = config.protonation_path,
        .fix_path = config.fix_path,
        .dump_movers_path = config.dump_movers_path,
    };

    const result = zreduce.run.processFile(allocator, proc_config) catch |err| {
        std.debug.print("Error: processing failed: {s}\n", .{@errorName(err)});
        std.process.exit(1);
    };

    if (result.n_movers > 0) {
        std.debug.print("  Movers: {d} ({d} singletons, {d} brute-force, {d} greedy)\n", .{
            result.n_movers,
            result.n_singletons,
            result.n_brute_force,
            result.n_vertex_cut,
        });
    }

    std.debug.print("zreduce: placed {d} H atoms on {d} residues ({d} skipped)\n", .{
        result.n_placed,
        result.n_residues,
        result.totalSkipped(),
    });
    if (result.n_skipped_missing_ref > 0) {
        std.debug.print("  warning: {d} H skipped due to missing reference atoms (potential plan bug)\n", .{result.n_skipped_missing_ref});
    }
}
