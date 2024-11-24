const std = @import("std");
const testing = std.testing;

const Vec3 = struct { x: f32, y: f32, z: f32 };

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

pub fn fuerzaLJ(sigma: f32, epsilon: f32, r: f32) f32 {
    const sigma_r = sigma / r;
    const sigma_r_6 = std.math.pow(f32, sigma_r, 6);
    const sigma_r_12 = sigma_r_6 * sigma_r_6;
    return 24.0 * epsilon * (2.0 * sigma_r_12 - sigma_r_6) / r;
}

fn calcularFuerzaLJEn(i: u32, estado: Estado, sigma: f32, epsilon: f32) Vec3 {
    var fx: f32 = 0;
    var fy: f32 = 0;
    var fz: f32 = 0;

    const xi = estado.x[i];
    const yi = estado.y[i];
    const zi = if (estado.z) |z| z[i] else 0;

    for (estado.x, estado.y, 0..) |xj, yj, j| {
        if (i == j) continue; // Skip self-interaction

        const dx = xj - xi;
        const dy = yj - yi;
        const dz = if (estado.z) |z| z[j] - zi else 0;

        const r_squared = dx * dx + dy * dy + dz * dz;
        const r = @sqrt(r_squared);

        const force_magnitude = fuerzaLJ(sigma, epsilon, r);

        // Calculate force components
        const force_factor = force_magnitude / r;
        fx += force_factor * dx;
        fy += force_factor * dy;
        if (estado.z != null) {
            fz += force_factor * dz;
        }
    }

    return Vec3{ .x = fx, .y = fy, .z = fz };
}

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

test "fuerzaLJ - Repulsive force at short distance" {
    const sigma: f32 = 1.0;
    const epsilon: f32 = 1.0;
    const r: f32 = 0.9 * sigma;
    const force = fuerzaLJ(sigma, epsilon, r);
    try testing.expect(force > 0);
}

test "fuerzaLJ - Attractive force at medium distance" {
    const sigma: f32 = 1.0;
    const epsilon: f32 = 1.0;
    const r: f32 = 1.5 * sigma;
    const force = fuerzaLJ(sigma, epsilon, r);
    try testing.expect(force < 0);
}

test "fuerzaLJ - Force approaches zero at large distance" {
    const sigma: f32 = 1.0;
    const epsilon: f32 = 1.0;
    const r: f32 = 1000.0 * sigma;
    const force = fuerzaLJ(sigma, epsilon, r);
    try testing.expectApproxEqAbs(@as(f32, 0.0), force, 1e-6);
}

test "fuerzaLJ - Force at equilibrium distance" {
    const sigma: f32 = 1.0;
    const epsilon: f32 = 1.0;
    const r: f32 = std.math.pow(f32, 2.0, 1.0/6.0) * sigma; // equilibrium distance
    const force = fuerzaLJ(sigma, epsilon, r);
    try testing.expectApproxEqAbs(@as(f32, 0.0), force, 1e-6);
}

test "calcularFuerzaLJEn - Diagonal configuration" {
    const allocator = std.testing.allocator;
    const numero_particulas: u32 = 2;

    // Allocate memory for positions and velocities
    var x = try allocator.alloc(f32, numero_particulas);
    defer allocator.free(x);
    var y = try allocator.alloc(f32, numero_particulas);
    defer allocator.free(y);
    var z = try allocator.alloc(f32, numero_particulas);
    defer allocator.free(z);
    var vx = try allocator.alloc(f32, numero_particulas);
    defer allocator.free(vx);
    var vy = try allocator.alloc(f32, numero_particulas);
    defer allocator.free(vy);
    var vz = try allocator.alloc(f32, numero_particulas);
    defer allocator.free(vz);

    // Initialize positions and velocities
    x[0] = 0.0; x[1] = 1.0;
    y[0] = 0.0; y[1] = 1.0;
    z[0] = 0.0; z[1] = 1.0;
    vx[0] = 0.0; vx[1] = 0.0;
    vy[0] = 0.0; vy[1] = 0.0;
    vz[0] = 0.0; vz[1] = 0.0;

    const estado = Estado.init(x, y, z, vx, vy, vz);
    
    // Test calcularFuerzaLJEn function
    const epsilon: f32 = 1.0;
    const sigma: f32 = 1.0;
    
    // Calculate forces for both particles
    const fuerza0 = calcularFuerzaLJEn(0, estado, sigma, epsilon);
    const fuerza1 = calcularFuerzaLJEn(1, estado, sigma, epsilon);

    // Expected force magnitude for particles at distance sqrt(3) with epsilon=1.0 and sigma=1.0
    const r = std.math.sqrt(3.0);
    const fuerza_esperada_magnitud = 24.0 * epsilon * (2.0 * std.math.pow(f32, sigma / r, 12.0) - std.math.pow(f32, sigma / r, 6.0)) / r;
    const fuerza_esperada_componente = fuerza_esperada_magnitud / std.math.sqrt(3.0);

    // Test force on particle 0
    try testing.expectApproxEqAbs(fuerza0.x, fuerza_esperada_componente, 1e-6);
    try testing.expectApproxEqAbs(fuerza0.y, fuerza_esperada_componente, 1e-6);
    try testing.expectApproxEqAbs(fuerza0.z, fuerza_esperada_componente, 1e-6);

    // Test force on particle 1 (should be equal and opposite)
    try testing.expectApproxEqAbs(fuerza1.x, -fuerza_esperada_componente, 1e-6);
    try testing.expectApproxEqAbs(fuerza1.y, -fuerza_esperada_componente, 1e-6);
    try testing.expectApproxEqAbs(fuerza1.z, -fuerza_esperada_componente, 1e-6);
}