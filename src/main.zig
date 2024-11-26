const std = @import("std");
const root = @import("root.zig");

pub fn main() !void {
    // Create an allocator
    var gpa = std.heap.GeneralPurposeAllocator(.{}){};
    defer _ = gpa.deinit();
    const allocator = gpa.allocator();

    // Set up simulation parameters
    const numero_particulas: u32 = 100;
    const numero_dimensiones: u32 = 3;
    const lado: u32 = 1e6;
    const epsilon: f32 = 1.0;
    const sigma: f32 = 1.0;
    const dt: f32 = 1e-6;
    const iteraciones_max: u32 = 1e4;

    // Initialize the simulation
    var sim = try root.Simulacion.init(
        numero_particulas,
        numero_dimensiones,
        lado,
        .{ .posicion = .esquina, .velocidad = .al_azar },
        iteraciones_max,
        allocator,
        epsilon,
        sigma,
        dt
    );
    defer sim.deinit();

    // Run the simulation
    try sim.correr();

    // Export states to file
    try exportarEstados(&sim, "estados.txt");

    // Print some final statistics or results
    std.debug.print("Simulation completed with {} particles for {} iterations.\n", .{
        numero_particulas, iteraciones_max
    });
}

fn exportarEstados(sim: *root.Simulacion, filename: []const u8) !void {
    const file = try std.fs.cwd().createFile(filename, .{});
    defer file.close();

    var writer = file.writer();

    for (sim.estados) |estado| {
        // Write x coordinates
        for (estado.x) |x| {
            try writer.print("{d} ", .{x});
        }
        try writer.writeByte('\n');

        // Write y coordinates
        for (estado.y) |y| {
            try writer.print("{d} ", .{y});
        }
        try writer.writeByte('\n');

        // Write z coordinates (if 3D)
        if (estado.z) |z| {
            for (z) |z_val| {
                try writer.print("{d} ", .{z_val});
            }
        } else {
            // If 2D, write zeros for z
            for (0..estado.x.len) |_| {
                try writer.writeAll("0 ");
            }
        }
        try writer.writeByte('\n');
    }
}