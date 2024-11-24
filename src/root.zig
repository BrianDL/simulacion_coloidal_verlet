const std = @import("std");
const testing = std.testing;

pub const EstratIniPosicion = enum { al_azar, esquina };

pub const EstratIniVelocidad = enum { cero, al_azar };

const EstrategiaInicializacion = struct {
    posicion: EstratIniPosicion,
    velocidad: EstratIniVelocidad,
};

const Estado = struct {
    x: []const f32,
    y: []const f32,
    z: ?[]const f32,
    vx: []const f32,
    vy: []const f32,
    vz: ?[]const f32,

    pub fn init(x: []const f32, y: []const f32, z: ?[]const f32, vx: []const f32, vy: []const f32, vz: ?[]const f32) Estado {
        return .{
            .x = x,
            .y = y,
            .z = z,
            .vx = vx,
            .vy = vy,
            .vz = vz,
        };
    }
};

pub const Simulacion = struct {
    estados: []Estado,
    lado: u32,
    estrategia_inicializacion: EstrategiaInicializacion,
    iteraciones_max: u32,
    allocator: std.mem.Allocator,

    const Self = @This();

    pub fn init(allocator: std.mem.Allocator, numero_particulas: u32, numero_dimensiones: u32, lado: u32, estrategia_inicializacion: EstrategiaInicializacion, iteraciones_max: u32) !Self {
        const x = try allocator.alloc(f32, numero_particulas);
        const y = try allocator.alloc(f32, numero_particulas);
        const z = if (numero_dimensiones == 3) try allocator.alloc(f32, numero_particulas) else null;
        const vx = try allocator.alloc(f32, numero_particulas);
        const vy = try allocator.alloc(f32, numero_particulas);
        const vz = if (numero_dimensiones == 3) try allocator.alloc(f32, numero_particulas) else null;

        // Inicializar posiciones
        switch (estrategia_inicializacion.posicion) {
            .al_azar => inicializarPosicionAlAzar(x, y, z, lado),
            .esquina => inicializarPosicionEsquina(x, y, z),
        }

        // Inicializar velocidades
        switch (estrategia_inicializacion.velocidad) {
            .cero => inicializarVelocidadCero(vx, vy, vz),
            .al_azar => inicializarVelocidadAlAzar(vx, vy, vz, 1.0),
        }

        const estado_inicial = Estado.init(x, y, z, vx, vy, vz);

        // Create the estados list with iteraciones_max length
        var estados = try allocator.alloc(Estado, iteraciones_max);
        estados[0] = estado_inicial;

        // The rest of the estados will be filled during the simulation

        return Self{
            .estados = estados,
            .lado = lado,
            .estrategia_inicializacion = estrategia_inicializacion,
            .iteraciones_max = iteraciones_max,
            .allocator = allocator,
        };
    }

    pub fn deinit(self: *Self) void {
        // Free the memory for the first state (which we initialized)
        self.allocator.free(self.estados[0].x);
        self.allocator.free(self.estados[0].y);
        if (self.estados[0].z) |z| self.allocator.free(z);

        self.allocator.free(self.estados[0].vx);
        self.allocator.free(self.estados[0].vy);
        if (self.estados[0].vz) |vz| self.allocator.free(vz);

        // Free the estados array itself
        self.allocator.free(self.estados);
    }

    fn inicializarPosicionAlAzar(x: []f32, y: []f32, z: ?[]f32, lado: u32) void {
        var rng = std.rand.DefaultPrng.init(blk: {
            var seed: u64 = undefined;
            std.crypto.random.bytes(std.mem.asBytes(&seed));
            break :blk seed;
        });
        const random = rng.random();
        const lado_flt = @as(f32, @floatFromInt(lado));

        for (x, y, 0..) |*xi, *yi, i| {
            xi.* = random.float(f32) * lado_flt;
            yi.* = random.float(f32) * lado_flt;
            if (z) |z_slice| {
                z_slice[i] = random.float(f32) * lado_flt;
            }
        }
    }

    fn inicializarPosicionEsquina(x: []f32, y: []f32, z: ?[]f32) void {
        @memset(x, 0);
        @memset(y, 0);
        if (z) |z_slice| {
            @memset(z_slice, 0);
        }
    }

    fn inicializarVelocidadCero(vx: []f32, vy: []f32, vz: ?[]f32) void {
        @memset(vx, 0);
        @memset(vy, 0);
        if (vz) |z_slice| {
            @memset(z_slice, 0);
        }
    }

    fn inicializarVelocidadAlAzar(vx: []f32, vy: []f32, vz: ?[]f32, velocidad_max: f32) void {
        var rng = std.rand.DefaultPrng.init(blk: {
            var seed: u64 = undefined;
            std.crypto.random.bytes(std.mem.asBytes(&seed));
            break :blk seed;
        });
        const random = rng.random();

        for (vx, vy, 0..) |*vxi, *vyi, i| {
            vxi.* = (random.float(f32) * 2 - 1) * velocidad_max;
            vyi.* = (random.float(f32) * 2 - 1) * velocidad_max;
            if (vz) |z_slice| {
                z_slice[i] = (random.float(f32) * 2 - 1) * velocidad_max;
            }
        }
    }

    pub fn correr(self: *Self) void {
        std.debug.print("Ejecutando simulación con {} partículas.\n", .{self.estados.len});
        // Aquí iría la lógica de la simulación
    }
};

// ///// INICIAN PRUEBAS UNITARIAS////////////

test "inicializarPosicionAlAzar" {
    const allocator = std.testing.allocator;
    const numero_particulas: u32 = 100;
    const numero_dimensiones: u32 = 3;
    const lado: u32 = 10;
    const iteraciones_max: u32 = 1;

    var sim = try Simulacion.init(allocator, numero_particulas, numero_dimensiones, lado, EstrategiaInicializacion{ .posicion = .al_azar, .velocidad = .cero }, iteraciones_max);
    defer sim.deinit();

    const estado = sim.estados[0];
    for (estado.x, estado.y, estado.z.?) |x, y, z| {
        try testing.expect(x >= 0 and x <= @as(f32, lado));
        try testing.expect(y >= 0 and y <= @as(f32, lado));
        try testing.expect(z >= 0 and z <= @as(f32, lado));
    }
}

test "inicializarPosicionEsquina" {
    const allocator = std.testing.allocator;
    const numero_particulas: u32 = 100;
    const numero_dimensiones: u32 = 3;
    const lado: u32 = 10;
    const iteraciones_max: u32 = 1;

    var sim = try Simulacion.init(allocator, numero_particulas, numero_dimensiones, lado, EstrategiaInicializacion{ .posicion = .esquina, .velocidad = .cero }, iteraciones_max);
    defer sim.deinit();

    const estado = sim.estados[0];
    for (estado.x, estado.y, estado.z.?) |x, y, z| {
        try testing.expectEqual(@as(f32, 0), x);
        try testing.expectEqual(@as(f32, 0), y);
        try testing.expectEqual(@as(f32, 0), z);
    }
}

test "inicializarVelocidadCero" {
    const allocator = std.testing.allocator;
    const numero_particulas: u32 = 100;
    const numero_dimensiones: u32 = 3;
    const lado: u32 = 10;
    const iteraciones_max: u32 = 1;

    var sim = try Simulacion.init(allocator, numero_particulas, numero_dimensiones, lado, EstrategiaInicializacion{ .posicion = .al_azar, .velocidad = .cero }, iteraciones_max);
    defer sim.deinit();

    const estado = sim.estados[0];
    for (estado.vx, estado.vy, estado.vz.?) |vx, vy, vz| {
        try testing.expectEqual(@as(f32, 0), vx);
        try testing.expectEqual(@as(f32, 0), vy);
        try testing.expectEqual(@as(f32, 0), vz);
    }
}

test "inicializarVelocidadAlAzar" {
    const allocator = std.testing.allocator;
    const numero_particulas: u32 = 100;
    const numero_dimensiones: u32 = 3;
    const lado: u32 = 10;
    const iteraciones_max: u32 = 1;
    const velocidad_max: f32 = 1.0;

    var sim = try Simulacion.init(allocator, numero_particulas, numero_dimensiones, lado, EstrategiaInicializacion{ .posicion = .al_azar, .velocidad = .al_azar }, iteraciones_max);
    defer sim.deinit();

    const estado = sim.estados[0];
    for (estado.vx, estado.vy, estado.vz.?) |vx, vy, vz| {
        try testing.expect(vx >= -velocidad_max and vx <= velocidad_max);
        try testing.expect(vy >= -velocidad_max and vy <= velocidad_max);
        try testing.expect(vz >= -velocidad_max and vz <= velocidad_max);
    }
}
