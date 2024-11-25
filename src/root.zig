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
    
    fn deinit(self: *Estado, allocator: std.mem.Allocator) void {
        
        allocator.free(self.x);
        allocator.free(self.y);
        if (self.z) |z_slice| {
            allocator.free(z_slice);
        }

        allocator.free(self.vx);
        allocator.free(self.vy);
        if (self.vz) |vz_slice| {
            allocator.free(vz_slice);
        }
    }
};

pub fn fuerzaLJ(sigma: f32, epsilon: f32, r: f32) f32 {
    const sigma_r = sigma / r;
    const sigma_r_6 = std.math.pow(f32, sigma_r, 6);
    const sigma_r_12 = sigma_r_6 * sigma_r_6;
    return 24.0 * epsilon * (2.0 * sigma_r_12 - sigma_r_6) / r;
}

fn calcularFuerzaLJEn(i: usize, estado: Estado, sigma: f32, epsilon: f32) Vec3 {
    var fx: f32 = 0;
    var fy: f32 = 0;
    var fz: f32 = 0;

    const xi = estado.x[i];
    const yi = estado.y[i];
    const zi = if (estado.z) |z| z[i] else 0;

    for (estado.x, estado.y, 0..) |xj, yj, j| {
        if (i == j) continue; // Skip self-interaction

        const dx = xi - xj;
        const dy = yi - yj;
        const dz = if (estado.z) |z| zi - z[j] else 0;

        const r_squared = dx * dx + dy * dy + dz * dz;
        const r = @sqrt(r_squared);
        // std.debug.print("R: {any}\n", .{r});

        const force_magnitude = fuerzaLJ(sigma, epsilon, r);
        // std.debug.print("SIGMA_EPSILON: {any}, {any}\n", .{sigma, epsilon});
        // std.debug.print("FORCE_MAGNITUD: {any}\n", .{force_magnitude});

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

fn pasoVerlet(
    estado: Estado, 
    lado: f32,
    dt_param: ?f32,
    sigma_param: ?f32,
    epsilon_param: ?f32, 
    allocator: std.mem.Allocator,
    current_forces: ?[]Vec3
) !Estado {
    const dt = dt_param orelse 0.01;
    const sigma = sigma_param orelse 1.0;
    const epsilon = epsilon_param orelse 1.0;

    const n = estado.x.len;

    var nuevo_x = try allocator.alloc(f32, n);
    errdefer allocator.free(nuevo_x);
    var nuevo_y = try allocator.alloc(f32, n);
    errdefer allocator.free(nuevo_y);
    var nuevo_z = try allocator.alloc(f32, n);
    errdefer allocator.free(nuevo_z);
    var nuevo_vx = try allocator.alloc(f32, n);
    errdefer allocator.free(nuevo_vx);
    var nuevo_vy = try allocator.alloc(f32, n);
    errdefer allocator.free(nuevo_vy);
    var nuevo_vz = try allocator.alloc(f32, n);
    errdefer allocator.free(nuevo_vz);

    var fuerzas = try allocator.alloc(Vec3, n);
    defer allocator.free(fuerzas);

    // Use memoized forces if available, otherwise calculate
    if (current_forces) |forces| {
        std.mem.copyForwards(Vec3, fuerzas, forces);
    } else {
        for (0..n) |i| {
            fuerzas[i] = calcularFuerzaLJEn(i, estado, sigma, epsilon);
        }
        
        // std.debug.print("\nFUERZAS: {any}\n", .{fuerzas});
    }
    

    // Update positions and calculate half-step velocities
    for (0..n) |i| {
        const fuerza = fuerzas[i];
        
        // Calculate new positions
        var nuevo_x_temp = estado.x[i] + estado.vx[i] * dt + 0.5 * fuerza.x * dt * dt;
        var nuevo_y_temp = estado.y[i] + estado.vy[i] * dt + 0.5 * fuerza.y * dt * dt;
        var nuevo_z_temp = estado.z.?[i] + estado.vz.?[i] * dt + 0.5 * fuerza.z * dt * dt;

        // Check for boundary collisions and apply elastic reflection
        if (nuevo_x_temp >= lado or nuevo_x_temp <= 0) {
            nuevo_x_temp = estado.x[i];
            nuevo_vx[i] = -(estado.vx[i] + 0.5 * fuerza.x * dt);
        } else {
            nuevo_vx[i] = estado.vx[i] + 0.5 * fuerza.x * dt;
        }

        if (nuevo_y_temp >= lado or nuevo_y_temp <= 0) {
            nuevo_y_temp = estado.y[i];
            nuevo_vy[i] = -(estado.vy[i] + 0.5 * fuerza.y * dt);
        } else {
            nuevo_vy[i] = estado.vy[i] + 0.5 * fuerza.y * dt;
        }

        if (nuevo_z_temp >= lado or nuevo_z_temp <= 0) {
            nuevo_z_temp = estado.z.?[i];
            nuevo_vz[i] = -(estado.vz.?[i] + 0.5 * fuerza.z * dt);
        } else {
            nuevo_vz[i] = estado.vz.?[i] + 0.5 * fuerza.z * dt;
        }

        nuevo_x[i] = nuevo_x_temp;
        nuevo_y[i] = nuevo_y_temp;
        nuevo_z[i] = nuevo_z_temp;
    }

    // Calculate new forces and update velocities
    for (0..n) |i| {
        const nueva_fuerza = calcularFuerzaLJEn(i, Estado.init(nuevo_x, nuevo_y, nuevo_z, nuevo_vx, nuevo_vy, nuevo_vz), sigma, epsilon);
        
        nuevo_vx[i] += 0.5 * nueva_fuerza.x * dt;
        nuevo_vy[i] += 0.5 * nueva_fuerza.y * dt;
        nuevo_vz[i] += 0.5 * nueva_fuerza.z * dt;

        // Update forces for the next iteration if memoization is being used
        if (current_forces) |forces| {
            forces[i] = nueva_fuerza;
        }
    }

    return Estado.init(nuevo_x, nuevo_y, nuevo_z, nuevo_vx, nuevo_vy, nuevo_vz);
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
        for (0..self.estados.len) |i| {
            self.estados[i].deinit(self.allocator);
        }
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
    
    pub fn correr(self: *Self) !void {
        std.debug.print("Ejecutando simulación con {} partículas por {} iteraciones.\n", .{self.estados[0].x.len, self.iteraciones_max});

        const sigma: f32 = 1.0;  // You might want to make these parameters of Simulacion
        const epsilon: f32 = 1.0;
        const dt: f32 = 0.01;  // Time step, you might want to make this a parameter too

        // Allocate memory for memoized forces
        var current_forces = try self.allocator.alloc(Vec3, self.estados[0].x.len);
        defer self.allocator.free(current_forces);

        // Initialize forces for the first iteration
        for (0..self.estados[0].x.len) |i| {
            current_forces[i] = calcularFuerzaLJEn(i, self.estados[0], sigma, epsilon);
        }
        
        for (1..self.iteraciones_max) |i| {
            const estado_actual = self.estados[i - 1];
            const nuevo_estado = try pasoVerlet(
                estado_actual, self.lado,
                dt, sigma, epsilon,
                self.allocator,
                current_forces
            );
            self.estados[i] = nuevo_estado;

            // You could add some periodic output here, e.g.:
            if (i % 100 == 0) {
                std.debug.print("Completed iteration {}/{}\n", .{i, self.iteraciones_max});
            }
        }

        std.debug.print("Simulación completada.\n", .{});
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
    const fuerza_esperada_magnitud = -24.0 * epsilon * (2.0 * std.math.pow(f32, sigma / r, 12.0) - std.math.pow(f32, sigma / r, 6.0)) / r;
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

test "calcularFuerzaLJEn - Repulsive force at close distance" {
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
    x[0] = 0.0; x[1] = 0.8; // Particles are closer than sigma
    y[0] = 0.0; y[1] = 0.0;
    z[0] = 0.0; z[1] = 0.0;
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
    
    
    // Expected force magnitude for particles at distance 0.8 with epsilon=1.0 and sigma=1.0
    const r = 0.8;
    const fuerza_esperada_magnitud = 24.0 * epsilon * (2.0 * std.math.pow(f32, sigma / r, 12.0) - std.math.pow(f32, sigma / r, 6.0)) / r;

    // Test force on particle 0 (should be positive, pushing it away)
    try testing.expect(fuerza0.x < 0);
    try testing.expectApproxEqAbs(fuerza0.x, -fuerza_esperada_magnitud, 1e-6);
    try testing.expectApproxEqAbs(fuerza0.y, 0, 1e-6);
    try testing.expectApproxEqAbs(fuerza0.z, 0, 1e-6);

    // Test force on particle 1 (should be equal and opposite)
    try testing.expect(fuerza1.x > 0);
    try testing.expectApproxEqAbs(fuerza1.x, fuerza_esperada_magnitud, 1e-6);
    try testing.expectApproxEqAbs(fuerza1.y, 0, 1e-6);
    try testing.expectApproxEqAbs(fuerza1.z, 0, 1e-6);
}

test "pasoVerlet - Single particle, no forces" {
    const allocator = std.testing.allocator;
    const dt: f32 = 0.01;
    const sigma: f32 = 1.0;
    const epsilon: f32 = 1.0;
    const lado: f32 = 10.0;

    var x = try allocator.alloc(f32, 1);
    defer allocator.free(x);
    var y = try allocator.alloc(f32, 1);
    defer allocator.free(y);
    var z = try allocator.alloc(f32, 1);
    defer allocator.free(z);
    var vx = try allocator.alloc(f32, 1);
    defer allocator.free(vx);
    var vy = try allocator.alloc(f32, 1);
    defer allocator.free(vy);
    var vz = try allocator.alloc(f32, 1);
    defer allocator.free(vz);

    x[0] = 1.0; y[0] = 1.0; z[0] = 1.0;
    vx[0] = 1.0; vy[0] = 1.0; vz[0] = 1.0;

    const estado_inicial = Estado.init(x, y, z, vx, vy, vz);
    
    var estado_final = try pasoVerlet(estado_inicial, lado, dt, sigma, epsilon, allocator, null);
    defer estado_final.deinit(allocator);

    // Expected positions after dt
    const expected_pos = 1.0 + 1.0 * dt + 0.5 * 0.0 * dt * dt;
    // Expected velocities after dt (no forces, so they remain the same)
    const expected_vel = 1.0;

    try testing.expectApproxEqAbs(expected_pos, estado_final.x[0], 1e-6);
    try testing.expectApproxEqAbs(expected_pos, estado_final.y[0], 1e-6);
    
    const z_val = estado_final.z orelse unreachable;
    try testing.expectApproxEqAbs(expected_pos, z_val[0], 1e-6);
    
    try testing.expectApproxEqAbs(expected_vel, estado_final.vx[0], 1e-6);
    try testing.expectApproxEqAbs(expected_vel, estado_final.vy[0], 1e-6);
    
    const vz_val = estado_final.vz orelse unreachable;
    try testing.expectApproxEqAbs(expected_vel, vz_val[0], 1e-6);
}

test "pasoVerlet - Two particles, with forces" {
    const allocator = std.testing.allocator;
    const dt: f32 = 0.01;
    const sigma: f32 = 1.0;
    const epsilon: f32 = 1.0;
    const lado: f32 = 10.0;

    var x = try allocator.alloc(f32, 2);
    defer allocator.free(x);
    var y = try allocator.alloc(f32, 2);
    defer allocator.free(y);
    var z = try allocator.alloc(f32, 2);
    defer allocator.free(z);
    var vx = try allocator.alloc(f32, 2);
    defer allocator.free(vx);
    var vy = try allocator.alloc(f32, 2);
    defer allocator.free(vy);
    var vz = try allocator.alloc(f32, 2);
    defer allocator.free(vz);

    x[0] = 1.0; x[1] = 1.0 + 0.9 * sigma;
    y[0] = 0.0; y[1] = 0.0;
    z[0] = 0.0; z[1] = 0.0;
    vx[0] = 0.0; vx[1] = 0.0;
    vy[0] = 0.0; vy[1] = 0.0;
    vz[0] = 0.0; vz[1] = 0.0;

    const estado_inicial = Estado.init(x, y, z, vx, vy, vz);
    
    var estado_final = try pasoVerlet(estado_inicial, lado, dt, sigma, epsilon, allocator, null);
    defer estado_final.deinit(allocator);

    // Particles should move away from each other due to repulsive force
    try testing.expect(estado_final.x[0] < 1.0);
    try testing.expect(estado_final.x[1] > 1.0 + 0.9 * sigma);
    try testing.expectApproxEqAbs(0.0, estado_final.y[0], 1e-6);
    try testing.expectApproxEqAbs(0.0, estado_final.y[1], 1e-6);

    const z_val = estado_final.z orelse unreachable;
    try testing.expectApproxEqAbs(0.0, z_val[0], 1e-6);
    try testing.expectApproxEqAbs(0.0, z_val[1], 1e-6);
}

test "pasoVerlet - Memoization" {
    const allocator = std.testing.allocator;
    const dt: f32 = 0.01;
    const sigma: f32 = 1.0;
    const epsilon: f32 = 1.0;
    const lado: f32 = 10.0;

    var x = try allocator.alloc(f32, 2);
    defer allocator.free(x);
    var y = try allocator.alloc(f32, 2);
    defer allocator.free(y);
    var z = try allocator.alloc(f32, 2);
    defer allocator.free(z);
    var vx = try allocator.alloc(f32, 2);
    defer allocator.free(vx);
    var vy = try allocator.alloc(f32, 2);
    defer allocator.free(vy);
    var vz = try allocator.alloc(f32, 2);
    defer allocator.free(vz);

    x[0] = 0.0; x[1] = 1.0;
    y[0] = 0.0; y[1] = 0.0;
    z[0] = 0.0; z[1] = 0.0;
    vx[0] = 0.0; vx[1] = 0.0;
    vy[0] = 0.0; vy[1] = 0.0;
    vz[0] = 0.0; vz[1] = 0.0;

    const estado_inicial = Estado.init(x, y, z, vx, vy, vz);
    
    var current_forces = try allocator.alloc(Vec3, 2);
    defer allocator.free(current_forces);
    current_forces[0] = Vec3{ .x = 1.0, .y = 0.0, .z = 0.0 };
    current_forces[1] = Vec3{ .x = -1.0, .y = 0.0, .z = 0.0 };

    var estado_final = try pasoVerlet(estado_inicial, lado, dt, sigma, epsilon, allocator, current_forces);
    defer estado_final.deinit(allocator);

    // Check if the memoized forces were used
    const expected_pos_0 = 0.0 + 0.0 * dt + 0.5 * (1.0 / 1.0) * dt * dt;
    const expected_pos_1 = 1.0 + 0.0 * dt + 0.5 * (-1.0 / 1.0) * dt * dt;

    try testing.expectApproxEqAbs(expected_pos_0, estado_final.x[0], 1e-6);
    try testing.expectApproxEqAbs(expected_pos_1, estado_final.x[1], 1e-6);
}