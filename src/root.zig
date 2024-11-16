const std = @import("std");
const testing = std.testing;

pub const EstrategiaInicializacion = enum {
    al_azar,
    esquina,
};

pub const Simulacion = struct {
    numero_particulas: u32,
    numero_dimensiones: u32,
    estrategia_inicializacion: EstrategiaInicializacion,
    iteraciones_max: u32,

    const Self = @This();

    pub fn init(
        numero_particulas: u32,
        numero_dimensiones: u32,
        estrategia_inicializacion: EstrategiaInicializacion,
        iteraciones_max: u32
    ) Self {
        return Self{
            .numero_particulas = numero_particulas,
            .numero_dimensiones = numero_dimensiones,
            .estrategia_inicializacion = estrategia_inicializacion,
            .iteraciones_max = iteraciones_max,
        };
    }

    pub fn correr(self: *Self) void {
        // Implementación de la simulación
        std.debug.print("Ejecutando simulación con {} partículas en {} dimensiones.\n", .{
            self.numero_particulas,
            self.numero_dimensiones,
        });
        // Aquí iría la lógica de la simulación
    }
};

test "Simulacion initialization" {
    const sim = Simulacion.init(100, 3, .al_azar, 1000);
    try testing.expectEqual(@as(u32, 100), sim.numero_particulas);
    try testing.expectEqual(@as(u32, 3), sim.numero_dimensiones);
    try testing.expectEqual(EstrategiaInicializacion.al_azar, sim.estrategia_inicializacion);
    try testing.expectEqual(@as(u32, 1000), sim.iteraciones_max);
}
