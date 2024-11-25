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
    const lado: u32 = 10;
    const epsilon: f32 = 1.0;
    const sigma: f32 = 1.0;
    const dt: f32 = 0.01;
    const iteraciones_max: u32 = 1000;

    // Initialize the simulation
    var sim = try root.Simulacion.init(
        numero_particulas,
        numero_dimensiones,
        lado,
        .{ .posicion = .al_azar, .velocidad = .al_azar },
        iteraciones_max,
        allocator,
        epsilon,
        sigma,
        dt
    );
    defer sim.deinit();

    // Run the simulation
    try sim.correr();

    // Print some final statistics or results
    std.debug.print("Simulation completed with {} particles for {} iterations.\n", .{
        numero_particulas, iteraciones_max
    });
}