const std = @import("std");

pub fn build(b: *std.Build) void {
    const target = b.standardTargetOptions(.{});
    const optimize = b.standardOptimizeOption(.{});

    const mod = b.addModule("eigencube", .{
        .root_source_file = b.path("src/eigencube.zig"),
        .target = target,
        .optimize = optimize,
    });

    const exe = b.addExecutable(.{
        .name = "eigencube",
        .root_source_file = b.path("src/main.zig"),
        .target = target,
        .optimize = optimize,
    });
    exe.root_module.addImport("eigencube", mod);
    b.installArtifact(exe);

    const run_cmd = b.addRunArtifact(exe);
    run_cmd.step.dependOn(b.getInstallStep());
    if (b.args) |args| {
        run_cmd.addArgs(args);
    }
    const run_step = b.step("run", "Run eigencube CLI");
    run_step.dependOn(&run_cmd.step);

    const tests = b.addTest(.{
        .root_source_file = b.path("tests/test_eigencube.zig"),
        .target = target,
        .optimize = optimize,
    });
    tests.root_module.addImport("eigencube", mod);
    const run_tests = b.addRunArtifact(tests);
    if (b.args) |args| {
        run_tests.addArgs(args);
    }
    const test_step = b.step("test", "Run all unit tests");
    test_step.dependOn(&run_tests.step);
}
