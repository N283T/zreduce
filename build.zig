const std = @import("std");

const version = "0.1.0";

pub fn build(b: *std.Build) void {
    const target = b.standardTargetOptions(.{});
    const optimize = b.standardOptimizeOption(.{});

    // Library module
    const mod = b.addModule("zreduce", .{
        .root_source_file = b.path("src/root.zig"),
        .target = target,
    });

    // CLI executable
    const options = b.addOptions();
    options.addOption([]const u8, "version", version);

    const exe_module = b.createModule(.{
        .root_source_file = b.path("src/main.zig"),
        .target = target,
        .optimize = optimize,
        .imports = &.{
            .{ .name = "zreduce", .module = mod },
            .{ .name = "build_options", .module = options.createModule() },
        },
    });
    exe_module.link_libc = true;
    exe_module.linkSystemLibrary("z", .{});

    const exe = b.addExecutable(.{
        .name = "zreduce",
        .root_module = exe_module,
    });
    b.installArtifact(exe);

    // Run step
    const run_step = b.step("run", "Run zreduce");
    const run_cmd = b.addRunArtifact(exe);
    run_step.dependOn(&run_cmd.step);
    run_cmd.step.dependOn(b.getInstallStep());
    if (b.args) |args| {
        run_cmd.addArgs(args);
    }

    // Test step
    const mod_tests = b.addTest(.{ .root_module = mod });
    const exe_tests = b.addTest(.{ .root_module = exe.root_module });
    const test_step = b.step("test", "Run tests");
    test_step.dependOn(&b.addRunArtifact(mod_tests).step);
    test_step.dependOn(&b.addRunArtifact(exe_tests).step);
}
