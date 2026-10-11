const std = @import("std");
const eigencube = @import("eigencube");

pub fn main() !void {
    var gpa = std.heap.GeneralPurposeAllocator(.{}){};
    defer _ = gpa.deinit();
    const allocator = gpa.allocator();

    var args = try std.process.argsWithAllocator(allocator);
    defer args.deinit();

    _ = args.next(); // skip exe name
    var seed: u64 = 42;
    if (args.next()) |seed_str| {
        if (std.mem.eql(u8, seed_str, "-h") or std.mem.eql(u8, seed_str, "--help")) {
            std.debug.print("Usage: zig build run -- [seed]\n", .{});
            std.debug.print("  [seed]   Solves a Rubik's cube scrambled from random seed (default: 42)\n", .{});
            return;
        }
        seed = std.fmt.parseInt(u64, seed_str, 10) catch {
            std.debug.print("Error: unrecognized option or invalid integer seed '{s}'.\n\n", .{seed_str});
            std.debug.print("Usage: zig build run -- [seed]\n", .{});
            std.debug.print("  [seed]   Solves a Rubik's cube scrambled from random seed (default: 42)\n", .{});
            std.process.exit(1);
        };
    }

    eigencube.log("Solving scrambled cube (seed={d})...", .{seed});
    const scrambled = eigencube.shuffle(100_000, seed, eigencube.SOLVED_CUBE);
    const moves = try eigencube.solve(allocator, scrambled);
    defer allocator.free(moves);

    std.debug.print("Solution length: {d}\n", .{moves.len});
}
