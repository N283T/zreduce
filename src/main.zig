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
};

fn readFile(allocator: Allocator, path: []const u8) ![]const u8 {
    const file = try std.fs.cwd().openFile(path, .{});
    defer file.close();
    return try file.readToEndAlloc(allocator, 1024 * 1024 * 1024);
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
        \\    --no-opt           Skip optimization
        \\    --no-flip          Disable Asn/Gln/His flips
        \\    --validate         Print validation diagnostics
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
        std.debug.print("Error: batch subcommand not yet implemented\n", .{});
        std.process.exit(1);
    } else {
        std.debug.print("Error: unknown subcommand '{s}'\n", .{subcmd});
        printUsage(args[0]);
        std.process.exit(1);
    }
}

fn runSubcommand(allocator: Allocator, args: []const []const u8) void {
    const config = parseRunArgs(args) orelse return;

    // Load CCD dictionary (once, before processFile)
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

    const proc_config = zreduce.run.ProcessConfig{
        .input_path = config.input_path,
        .output_path = config.output_path,
        .dict = if (ccd_dict) |*d| d else null,
        .json_path = config.json_path,
        .json_version = build_options.version,
        .no_opt = config.no_opt,
        .no_flip = config.no_flip,
        .validate_flag = config.validate,
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
        result.n_skipped,
    });
}
