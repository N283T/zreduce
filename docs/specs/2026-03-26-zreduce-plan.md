# zreduce Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Build a complete hydrogen placement and optimization tool for mmCIF structures in Zig, with full feature parity with original C++ reduce.

**Architecture:** Pipeline of specialized modules: CIF tokenizer/parser → mmCIF extraction → CCD loading → molecular model → geometric H placement → dot-sphere scoring → clique-based optimization → mmCIF/JSON output. Each module is independently testable with inline Zig tests.

**Tech Stack:** Zig (latest stable), zlib (for gzip CCD), no external Zig dependencies.

---

## File Structure

```
~/zreduce/
├── build.zig                     # Build configuration (exe + lib + tests)
├── build.zig.zon                 # Package metadata
├── src/
│   ├── main.zig                  # CLI entry point
│   ├── root.zig                  # Library re-exports
│   ├── math.zig                  # Vec3(T), rotation, dihedral, angle
│   ├── element.zig               # AtomType, VDW radii, flags (comptime tables)
│   ├── cif/
│   │   ├── char_table.zig        # Character classification LUT
│   │   ├── tokenizer.zig         # CIF token stream
│   │   ├── parser.zig            # Tokens → Document/Block/Loop
│   │   ├── types.zig             # Document, Block, Loop, Item, Pair
│   │   └── value.zig             # Value extraction (null check, as_float, as_int, as_string)
│   ├── cif.zig                   # CIF module re-exports
│   ├── mmcif.zig                 # atom_site + struct_conn → Model
│   ├── ccd.zig                   # Streaming CCD parser → ComponentDict
│   ├── model/
│   │   ├── atom.zig              # Atom struct
│   │   ├── residue.zig           # Residue struct
│   │   ├── chain.zig             # Chain struct
│   │   ├── bond.zig              # Bond struct, BondOrder
│   │   ├── model.zig             # Model aggregate + queries
│   │   └── neighbor.zig          # Spatial cell list
│   ├── model.zig                 # Model module re-exports
│   ├── place/
│   │   ├── geometry.zig          # Type 1-6 placement functions
│   │   ├── standard.zig          # Hardcoded plans for 20 AA + nucleic acids
│   │   ├── het.zig               # CCD-derived placement for HET groups
│   │   └── placer.zig            # Unified placement entry point
│   ├── place.zig                 # Place module re-exports
│   ├── optimize/
│   │   ├── dot_sphere.zig        # Dot sphere generation
│   │   ├── scorer.zig            # Bump/H-bond scoring
│   │   ├── mover.zig             # Mover interface + types
│   │   ├── rotator.zig           # OH/SH/NH3+/methyl movers
│   │   ├── flipper.zig           # Asn/Gln/His flip movers
│   │   ├── clique.zig            # Interaction graph + clique detection
│   │   └── optimizer.zig         # Brute-force + vertex-cut optimizer
│   ├── optimize.zig              # Optimize module re-exports
│   ├── writer/
│   │   ├── mmcif_writer.zig      # mmCIF output with H + flip corrections
│   │   └── json_writer.zig       # JSON log output
│   └── writer.zig                # Writer module re-exports
├── test_data/
│   ├── tiny.cif                  # Minimal mmCIF (1 residue, for unit tests)
│   ├── ala_dipeptide.cif         # Ala dipeptide (backbone H test)
│   ├── his_test.cif              # Structure with His (flip test)
│   └── het_test.cif              # Structure with HET group
└── docs/
    └── specs/
        ├── 2026-03-26-zreduce-design.md
        └── 2026-03-26-zreduce-plan.md
```

---

## Phase 1: Foundation

### Task 1: Project Scaffold + Build System

**Files:**
- Create: `~/zreduce/build.zig`
- Create: `~/zreduce/build.zig.zon`
- Create: `~/zreduce/src/main.zig`
- Create: `~/zreduce/src/root.zig`

- [ ] **Step 1: Initialize git repo**

```bash
cd ~/zreduce
git init
```

- [ ] **Step 2: Create build.zig.zon**

```zig
.{
    .name = .@"zreduce",
    .version = "0.1.0",
    .fingerprint = 0xDEADBEEF,
    .minimum_zig_version = "0.14.0",
    .paths = .{ "build.zig", "build.zig.zon", "src" },
}
```

- [ ] **Step 3: Create build.zig**

```zig
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
```

- [ ] **Step 4: Create minimal main.zig**

```zig
const std = @import("std");
const build_options = @import("build_options");

pub fn main() !void {
    const args = try std.process.argsAlloc(std.heap.page_allocator);
    defer std.process.argsFree(std.heap.page_allocator, args);

    if (args.len < 2) {
        printUsage(args[0]);
        std.process.exit(1);
    }

    const subcmd = args[1];
    if (std.mem.eql(u8, subcmd, "--help") or std.mem.eql(u8, subcmd, "-h")) {
        printUsage(args[0]);
        return;
    }
    if (std.mem.eql(u8, subcmd, "--version") or std.mem.eql(u8, subcmd, "-V")) {
        std.debug.print("zreduce {s}\n", .{build_options.version});
        return;
    }

    std.debug.print("Unknown command: {s}\n", .{subcmd});
    std.process.exit(1);
}

fn printUsage(program_name: []const u8) void {
    std.debug.print(
        \\zreduce {s} - Hydrogen placement for mmCIF structures
        \\
        \\USAGE:
        \\    {s} [OPTIONS] <input.cif> [-o output.cif]
        \\
        \\OPTIONS:
        \\    -h, --help       Show this help message
        \\    -V, --version    Show version
        \\    -d, --dict PATH  Path to components.cif[.gz]
        \\    -o, --output PATH  Output file (default: stdout)
        \\    --no-opt         Skip optimization (placement only)
        \\    --no-flip        Disable Asn/Gln/His flips
        \\
    , .{ build_options.version, program_name });
}
```

- [ ] **Step 5: Create root.zig**

```zig
//! zreduce: Hydrogen placement and optimization for mmCIF structures.

// Modules will be added as they are implemented.

test {
    @import("std").testing.refAllDecls(@This());
}
```

- [ ] **Step 6: Build and verify**

Run: `cd ~/zreduce && zig build`
Expected: Successful build, no errors.

Run: `zig build run -- --version`
Expected: `zreduce 0.1.0`

- [ ] **Step 7: Create .gitignore and commit**

```
.zig-cache/
zig-out/
zig-cache/
*.o
*.swp
```

```bash
cd ~/zreduce
git add .
git commit -m "feat: initial project scaffold with build system"
```

---

### Task 2: Math Library (Vec3, rotation, dihedral)

**Files:**
- Create: `~/zreduce/src/math.zig`

- [ ] **Step 1: Write failing tests for Vec3 basic operations**

```zig
// src/math.zig
const std = @import("std");
const testing = std.testing;

pub fn Vec3(comptime T: type) type {
    return struct { x: T, y: T, z: T };
}

test "Vec3 add" {
    const v1 = Vec3(f32){ .x = 1.0, .y = 2.0, .z = 3.0 };
    const v2 = Vec3(f32){ .x = 4.0, .y = 5.0, .z = 6.0 };
    const result = v1.add(v2);
    try testing.expectApproxEqAbs(result.x, 5.0, 1e-6);
    try testing.expectApproxEqAbs(result.y, 7.0, 1e-6);
    try testing.expectApproxEqAbs(result.z, 9.0, 1e-6);
}

test "Vec3 cross product" {
    const v1 = Vec3(f32){ .x = 1.0, .y = 0.0, .z = 0.0 };
    const v2 = Vec3(f32){ .x = 0.0, .y = 1.0, .z = 0.0 };
    const result = v1.cross(v2);
    try testing.expectApproxEqAbs(result.x, 0.0, 1e-6);
    try testing.expectApproxEqAbs(result.y, 0.0, 1e-6);
    try testing.expectApproxEqAbs(result.z, 1.0, 1e-6);
}

test "Vec3 normalize" {
    const v = Vec3(f32){ .x = 3.0, .y = 4.0, .z = 0.0 };
    const n = v.normalize();
    try testing.expectApproxEqAbs(n.length(), 1.0, 1e-6);
    try testing.expectApproxEqAbs(n.x, 0.6, 1e-6);
    try testing.expectApproxEqAbs(n.y, 0.8, 1e-6);
}
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `cd ~/zreduce && zig build test 2>&1 | head -20`
Expected: Compilation error — `add`, `cross`, `normalize` not defined.

- [ ] **Step 3: Implement Vec3 with all operations**

```zig
// src/math.zig
const std = @import("std");
const testing = std.testing;

pub fn Vec3(comptime T: type) type {
    return struct {
        const Self = @This();

        x: T,
        y: T,
        z: T,

        pub const zero = Self{ .x = 0, .y = 0, .z = 0 };

        pub fn add(self: Self, other: Self) Self {
            return .{ .x = self.x + other.x, .y = self.y + other.y, .z = self.z + other.z };
        }

        pub fn sub(self: Self, other: Self) Self {
            return .{ .x = self.x - other.x, .y = self.y - other.y, .z = self.z - other.z };
        }

        pub fn scale(self: Self, s: T) Self {
            return .{ .x = self.x * s, .y = self.y * s, .z = self.z * s };
        }

        /// Scale vector to given length. Returns zero if length is zero.
        pub fn scaleTo(self: Self, len: T) Self {
            const l = self.length();
            if (l < 1e-10) return zero;
            return self.scale(len / l);
        }

        pub fn dot(self: Self, other: Self) T {
            return self.x * other.x + self.y * other.y + self.z * other.z;
        }

        pub fn cross(self: Self, other: Self) Self {
            return .{
                .x = self.y * other.z - self.z * other.y,
                .y = self.z * other.x - self.x * other.z,
                .z = self.x * other.y - self.y * other.x,
            };
        }

        pub fn length(self: Self) T {
            return @sqrt(self.dot(self));
        }

        pub fn distance(self: Self, other: Self) T {
            return self.sub(other).length();
        }

        pub fn normalize(self: Self) Self {
            return self.scaleTo(1.0);
        }

        pub fn negate(self: Self) Self {
            return .{ .x = -self.x, .y = -self.y, .z = -self.z };
        }

        /// Convert from f64 to f32 or vice versa
        pub fn cast(self: Self, comptime U: type) Vec3(U) {
            return .{
                .x = @floatCast(self.x),
                .y = @floatCast(self.y),
                .z = @floatCast(self.z),
            };
        }
    };
}

/// Rotate point around axis (origin + direction) by angle in degrees.
pub fn rotateAroundAxis(
    comptime T: type,
    point: Vec3(T),
    origin: Vec3(T),
    axis_dir: Vec3(T),
    degrees: T,
) Vec3(T) {
    const rad = degrees * std.math.pi / 180.0;
    const cos_a = @cos(rad);
    const sin_a = @sin(rad);
    const u = axis_dir.normalize();
    const p = point.sub(origin);

    // Rodrigues' rotation formula: p' = p*cos + (u x p)*sin + u*(u.p)*(1-cos)
    const term1 = p.scale(cos_a);
    const term2 = u.cross(p).scale(sin_a);
    const term3 = u.scale(u.dot(p) * (1.0 - cos_a));

    return term1.add(term2).add(term3).add(origin);
}

/// Compute dihedral angle (in degrees) for 4 points a-b-c-d.
pub fn dihedralAngle(comptime T: type, a: Vec3(T), b: Vec3(T), c: Vec3(T), d: Vec3(T)) T {
    const b1 = b.sub(a);
    const b2 = c.sub(b);
    const b3 = d.sub(c);

    const n1 = b1.cross(b2);
    const n2 = b2.cross(b3);

    const m1 = n1.cross(b2.normalize());

    const x = n1.dot(n2);
    const y = m1.dot(n2);

    return std.math.atan2(y, x) * 180.0 / std.math.pi;
}

/// Compute angle (in degrees) at vertex b for points a-b-c.
pub fn angle(comptime T: type, a: Vec3(T), b: Vec3(T), c: Vec3(T)) T {
    const ba = a.sub(b).normalize();
    const bc = c.sub(b).normalize();
    const cos_val = @max(@as(T, -1.0), @min(@as(T, 1.0), ba.dot(bc)));
    return std.math.acos(cos_val) * 180.0 / std.math.pi;
}
```

- [ ] **Step 4: Add tests for rotation and dihedral**

```zig
test "rotateAroundAxis 90 degrees around Z" {
    const p = Vec3(f64){ .x = 1.0, .y = 0.0, .z = 0.0 };
    const origin = Vec3(f64).zero;
    const axis = Vec3(f64){ .x = 0.0, .y = 0.0, .z = 1.0 };
    const result = rotateAroundAxis(f64, p, origin, axis, 90.0);
    try testing.expectApproxEqAbs(result.x, 0.0, 1e-10);
    try testing.expectApproxEqAbs(result.y, 1.0, 1e-10);
    try testing.expectApproxEqAbs(result.z, 0.0, 1e-10);
}

test "rotateAroundAxis 180 degrees" {
    const p = Vec3(f64){ .x = 1.0, .y = 0.0, .z = 0.0 };
    const origin = Vec3(f64).zero;
    const axis = Vec3(f64){ .x = 0.0, .y = 0.0, .z = 1.0 };
    const result = rotateAroundAxis(f64, p, origin, axis, 180.0);
    try testing.expectApproxEqAbs(result.x, -1.0, 1e-10);
    try testing.expectApproxEqAbs(result.y, 0.0, 1e-10);
}

test "dihedral angle trans" {
    // Trans conformation: dihedral ~ 180 degrees
    const a = Vec3(f64){ .x = 1.0, .y = 1.0, .z = 0.0 };
    const b = Vec3(f64){ .x = 0.0, .y = 1.0, .z = 0.0 };
    const c = Vec3(f64){ .x = 0.0, .y = 0.0, .z = 0.0 };
    const d = Vec3(f64){ .x = 1.0, .y = 0.0, .z = 0.0 };
    const dh = dihedralAngle(f64, a, b, c, d);
    try testing.expectApproxEqAbs(@abs(dh), 180.0, 1e-6);
}

test "angle 90 degrees" {
    const a = Vec3(f64){ .x = 1.0, .y = 0.0, .z = 0.0 };
    const b = Vec3(f64).zero;
    const c = Vec3(f64){ .x = 0.0, .y = 1.0, .z = 0.0 };
    const ang = angle(f64, a, b, c);
    try testing.expectApproxEqAbs(ang, 90.0, 1e-6);
}
```

- [ ] **Step 5: Run tests**

Run: `cd ~/zreduce && zig build test`
Expected: All tests pass.

- [ ] **Step 6: Commit**

```bash
cd ~/zreduce && git add src/math.zig && git commit -m "feat: add Vec3 math library with rotation and dihedral"
```

---

### Task 3: Element Table

**Files:**
- Create: `~/zreduce/src/element.zig`

- [ ] **Step 1: Write failing test for element lookup**

```zig
const std = @import("std");
const testing = std.testing;

test "hydrogen VDW radius" {
    const h = AtomType.H;
    const info = h.info();
    try testing.expectApproxEqAbs(info.explicit_radius, 1.22, 1e-3);
}

test "polar hydrogen is donor" {
    const hpol = AtomType.Hpol;
    const info = hpol.info();
    try testing.expect(info.flags.donor);
    try testing.expectApproxEqAbs(info.explicit_radius, 1.05, 1e-3);
}

test "oxygen is acceptor" {
    const o = AtomType.O;
    const info = o.info();
    try testing.expect(info.flags.acceptor);
    try testing.expectApproxEqAbs(info.explicit_radius, 1.40, 1e-3);
}
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `cd ~/zreduce && zig build test 2>&1 | head -5`
Expected: Compilation error — `AtomType` not defined.

- [ ] **Step 3: Implement element table**

```zig
// src/element.zig
const std = @import("std");
const testing = std.testing;

pub const AtomFlags = packed struct {
    donor: bool = false,
    acceptor: bool = false,
    aromatic: bool = false,
    positive: bool = false,
    negative: bool = false,
    metallic: bool = false,
    hb_only_dummy: bool = false,
    _padding: u1 = 0,
};

pub const AtomTypeInfo = struct {
    explicit_radius: f32,
    implicit_radius: f32,
    covalent_radius: f32,
    flags: AtomFlags,
};

/// Atom types matching original reduce ElementInfo table.
/// VDW radii from Gavezzotti (1983) and Bondi (1964).
pub const AtomType = enum(u8) {
    // Hydrogen variants
    H,         // non-polar H (1.22 Å)
    Har,       // aromatic H (1.05 Å)
    Hpol,      // polar H, donor (1.05 Å)
    Ha_p,      // aromatic + polar H (1.05 Å)
    HOd,       // H-bond only dummy (1.05 Å)

    // Carbon variants
    C,         // sp3 carbon (1.70 Å)
    Car,       // aromatic carbon (1.75 Å, acceptor)
    C_eq_O,    // carbonyl carbon (1.65 Å)

    // Nitrogen variants
    N,         // nitrogen (1.55 Å)
    Nacc,      // nitrogen acceptor (1.55 Å)

    // Others
    O,         // oxygen (1.40 Å, acceptor)
    P,         // phosphorus (1.80 Å)
    S,         // sulfur (1.80 Å, acceptor)
    Se,        // selenium (1.90 Å)
    F,         // fluorine (1.30 Å, acceptor)
    Cl,        // chlorine (1.77 Å, acceptor)
    Br,        // bromine (1.95 Å, acceptor)
    I,         // iodine (2.10 Å, acceptor)

    // Metals
    Li, Na, Mg, K, Ca, Mn, Fe, Co, Ni, Cu, Zn, As, Rb, Sr, Mo, Ag, Cd, Sn, Cs, Ba, W, Pt, Au, Hg, Pb, U,

    unknown,

    const count = @typeInfo(AtomType).@"enum".fields.len;

    pub fn info(self: AtomType) AtomTypeInfo {
        return atom_type_table[@intFromEnum(self)];
    }
};

const D = AtomFlags{ .donor = true };
const A = AtomFlags{ .acceptor = true };
const DA = AtomFlags{ .donor = true, .acceptor = true };
const AR = AtomFlags{ .aromatic = true };
const ARA = AtomFlags{ .aromatic = true, .acceptor = true };
const M = AtomFlags{ .metallic = true };
const HBD = AtomFlags{ .donor = true, .hb_only_dummy = true };
const NONE = AtomFlags{};

/// Comptime-generated lookup table from original reduce constants.
const atom_type_table: [AtomType.count]AtomTypeInfo = blk: {
    var table: [AtomType.count]AtomTypeInfo = undefined;
    // Hydrogen variants
    table[@intFromEnum(AtomType.H)]     = .{ .explicit_radius = 1.22, .implicit_radius = 0.00, .covalent_radius = 0.30, .flags = NONE };
    table[@intFromEnum(AtomType.Har)]   = .{ .explicit_radius = 1.05, .implicit_radius = 0.00, .covalent_radius = 0.30, .flags = AR };
    table[@intFromEnum(AtomType.Hpol)]  = .{ .explicit_radius = 1.05, .implicit_radius = 0.00, .covalent_radius = 0.30, .flags = D };
    table[@intFromEnum(AtomType.Ha_p)]  = .{ .explicit_radius = 1.05, .implicit_radius = 0.00, .covalent_radius = 0.30, .flags = D };
    table[@intFromEnum(AtomType.HOd)]   = .{ .explicit_radius = 1.05, .implicit_radius = 0.00, .covalent_radius = 0.30, .flags = HBD };
    // Carbon variants
    table[@intFromEnum(AtomType.C)]      = .{ .explicit_radius = 1.70, .implicit_radius = 1.90, .covalent_radius = 0.77, .flags = NONE };
    table[@intFromEnum(AtomType.Car)]    = .{ .explicit_radius = 1.75, .implicit_radius = 1.75, .covalent_radius = 0.77, .flags = ARA };
    table[@intFromEnum(AtomType.C_eq_O)] = .{ .explicit_radius = 1.65, .implicit_radius = 1.65, .covalent_radius = 0.80, .flags = NONE };
    // Nitrogen variants
    table[@intFromEnum(AtomType.N)]     = .{ .explicit_radius = 1.55, .implicit_radius = 1.70, .covalent_radius = 0.70, .flags = NONE };
    table[@intFromEnum(AtomType.Nacc)]  = .{ .explicit_radius = 1.55, .implicit_radius = 1.70, .covalent_radius = 0.70, .flags = A };
    // Non-metals
    table[@intFromEnum(AtomType.O)]     = .{ .explicit_radius = 1.40, .implicit_radius = 1.50, .covalent_radius = 0.66, .flags = A };
    table[@intFromEnum(AtomType.P)]     = .{ .explicit_radius = 1.80, .implicit_radius = 1.80, .covalent_radius = 1.80, .flags = NONE };
    table[@intFromEnum(AtomType.S)]     = .{ .explicit_radius = 1.80, .implicit_radius = 1.90, .covalent_radius = 1.04, .flags = A };
    table[@intFromEnum(AtomType.Se)]    = .{ .explicit_radius = 1.90, .implicit_radius = 1.90, .covalent_radius = 1.17, .flags = A };
    // Halogens
    table[@intFromEnum(AtomType.F)]     = .{ .explicit_radius = 1.30, .implicit_radius = 1.30, .covalent_radius = 0.58, .flags = A };
    table[@intFromEnum(AtomType.Cl)]    = .{ .explicit_radius = 1.77, .implicit_radius = 1.77, .covalent_radius = 0.99, .flags = A };
    table[@intFromEnum(AtomType.Br)]    = .{ .explicit_radius = 1.95, .implicit_radius = 1.95, .covalent_radius = 1.14, .flags = A };
    table[@intFromEnum(AtomType.I)]     = .{ .explicit_radius = 2.10, .implicit_radius = 2.10, .covalent_radius = 1.33, .flags = A };
    // Metals (using ionic radii)
    table[@intFromEnum(AtomType.Li)]    = .{ .explicit_radius = 1.82, .implicit_radius = 1.82, .covalent_radius = 1.23, .flags = M };
    table[@intFromEnum(AtomType.Na)]    = .{ .explicit_radius = 2.27, .implicit_radius = 2.27, .covalent_radius = 1.54, .flags = M };
    table[@intFromEnum(AtomType.Mg)]    = .{ .explicit_radius = 1.73, .implicit_radius = 1.73, .covalent_radius = 1.36, .flags = M };
    table[@intFromEnum(AtomType.K)]     = .{ .explicit_radius = 2.75, .implicit_radius = 2.75, .covalent_radius = 1.96, .flags = M };
    table[@intFromEnum(AtomType.Ca)]    = .{ .explicit_radius = 2.31, .implicit_radius = 2.31, .covalent_radius = 1.74, .flags = M };
    table[@intFromEnum(AtomType.Mn)]    = .{ .explicit_radius = 1.73, .implicit_radius = 1.73, .covalent_radius = 1.17, .flags = M };
    table[@intFromEnum(AtomType.Fe)]    = .{ .explicit_radius = 1.65, .implicit_radius = 1.65, .covalent_radius = 1.17, .flags = M };
    table[@intFromEnum(AtomType.Co)]    = .{ .explicit_radius = 1.67, .implicit_radius = 1.67, .covalent_radius = 1.16, .flags = M };
    table[@intFromEnum(AtomType.Ni)]    = .{ .explicit_radius = 1.50, .implicit_radius = 1.50, .covalent_radius = 1.15, .flags = M };
    table[@intFromEnum(AtomType.Cu)]    = .{ .explicit_radius = 1.52, .implicit_radius = 1.52, .covalent_radius = 1.17, .flags = M };
    table[@intFromEnum(AtomType.Zn)]    = .{ .explicit_radius = 1.65, .implicit_radius = 1.65, .covalent_radius = 1.25, .flags = M };
    table[@intFromEnum(AtomType.As)]    = .{ .explicit_radius = 1.85, .implicit_radius = 1.85, .covalent_radius = 1.21, .flags = NONE };
    table[@intFromEnum(AtomType.Rb)]    = .{ .explicit_radius = 2.75, .implicit_radius = 2.75, .covalent_radius = 2.11, .flags = M };
    table[@intFromEnum(AtomType.Sr)]    = .{ .explicit_radius = 2.49, .implicit_radius = 2.49, .covalent_radius = 1.92, .flags = M };
    table[@intFromEnum(AtomType.Mo)]    = .{ .explicit_radius = 1.90, .implicit_radius = 1.90, .covalent_radius = 1.30, .flags = M };
    table[@intFromEnum(AtomType.Ag)]    = .{ .explicit_radius = 1.72, .implicit_radius = 1.72, .covalent_radius = 1.34, .flags = M };
    table[@intFromEnum(AtomType.Cd)]    = .{ .explicit_radius = 1.58, .implicit_radius = 1.58, .covalent_radius = 1.48, .flags = M };
    table[@intFromEnum(AtomType.Sn)]    = .{ .explicit_radius = 2.17, .implicit_radius = 2.17, .covalent_radius = 1.40, .flags = M };
    table[@intFromEnum(AtomType.Cs)]    = .{ .explicit_radius = 3.01, .implicit_radius = 3.01, .covalent_radius = 2.25, .flags = M };
    table[@intFromEnum(AtomType.Ba)]    = .{ .explicit_radius = 2.68, .implicit_radius = 2.68, .covalent_radius = 1.98, .flags = M };
    table[@intFromEnum(AtomType.W)]     = .{ .explicit_radius = 1.80, .implicit_radius = 1.80, .covalent_radius = 1.30, .flags = M };
    table[@intFromEnum(AtomType.Pt)]    = .{ .explicit_radius = 1.75, .implicit_radius = 1.75, .covalent_radius = 1.28, .flags = M };
    table[@intFromEnum(AtomType.Au)]    = .{ .explicit_radius = 1.66, .implicit_radius = 1.66, .covalent_radius = 1.34, .flags = M };
    table[@intFromEnum(AtomType.Hg)]    = .{ .explicit_radius = 1.55, .implicit_radius = 1.55, .covalent_radius = 1.49, .flags = M };
    table[@intFromEnum(AtomType.Pb)]    = .{ .explicit_radius = 2.02, .implicit_radius = 2.02, .covalent_radius = 1.47, .flags = M };
    table[@intFromEnum(AtomType.U)]     = .{ .explicit_radius = 1.86, .implicit_radius = 1.86, .covalent_radius = 1.42, .flags = M };
    table[@intFromEnum(AtomType.unknown)] = .{ .explicit_radius = 1.80, .implicit_radius = 1.80, .covalent_radius = 1.00, .flags = NONE };
    break :blk table;
};

/// Look up element by symbol string (e.g., "C", "N", "FE").
pub fn elementFromSymbol(symbol: []const u8) AtomType {
    if (symbol.len == 0) return .unknown;
    var buf: [2]u8 = .{ ' ', ' ' };
    buf[0] = std.ascii.toUpper(symbol[0]);
    if (symbol.len > 1) buf[1] = std.ascii.toLower(symbol[1]);

    const map = std.StaticStringMap(AtomType).initComptime(.{
        .{ "H ", .H },  .{ "C ", .C },  .{ "N ", .N },   .{ "O ", .O },
        .{ "P ", .P },  .{ "S ", .S },  .{ "Se", .Se },  .{ "F ", .F },
        .{ "Cl", .Cl }, .{ "Br", .Br }, .{ "I ", .I },
        .{ "Li", .Li }, .{ "Na", .Na }, .{ "Mg", .Mg },  .{ "K ", .K },
        .{ "Ca", .Ca }, .{ "Mn", .Mn }, .{ "Fe", .Fe },  .{ "Co", .Co },
        .{ "Ni", .Ni }, .{ "Cu", .Cu }, .{ "Zn", .Zn },  .{ "As", .As },
        .{ "Rb", .Rb }, .{ "Sr", .Sr }, .{ "Mo", .Mo },  .{ "Ag", .Ag },
        .{ "Cd", .Cd }, .{ "Sn", .Sn }, .{ "Cs", .Cs },  .{ "Ba", .Ba },
        .{ "W ", .W },  .{ "Pt", .Pt }, .{ "Au", .Au },  .{ "Hg", .Hg },
        .{ "Pb", .Pb }, .{ "U ", .U },
    });

    return map.get(buf[0..2]) orelse .unknown;
}
```

- [ ] **Step 4: Add element lookup tests and run**

```zig
test "elementFromSymbol" {
    try testing.expectEqual(AtomType.C, elementFromSymbol("C"));
    try testing.expectEqual(AtomType.Fe, elementFromSymbol("FE"));
    try testing.expectEqual(AtomType.Fe, elementFromSymbol("Fe"));
    try testing.expectEqual(AtomType.unknown, elementFromSymbol("Xx"));
}
```

Run: `cd ~/zreduce && zig build test`
Expected: All tests pass.

- [ ] **Step 5: Register in root.zig and commit**

Add to root.zig:
```zig
pub const math = @import("math.zig");
pub const element = @import("element.zig");
```

```bash
cd ~/zreduce && git add src/math.zig src/element.zig src/root.zig && git commit -m "feat: add math library and element table"
```

---

### Task 4: CIF Tokenizer

**Files:**
- Create: `~/zreduce/src/cif/char_table.zig`
- Create: `~/zreduce/src/cif/tokenizer.zig`

- [ ] **Step 1: Create char_table.zig**

Character classification for CIF format (whitespace, ordinary, non-blank).

```zig
// src/cif/char_table.zig

/// Character categories for CIF tokenization.
pub const CharType = enum(u2) {
    /// Whitespace: space, tab, newline, carriage return
    whitespace = 0,
    /// Ordinary printable character (valid in unquoted values)
    ordinary = 1,
    /// Special: ; # $ ' " _ (context-dependent meaning)
    special = 2,
};

/// Lookup table for ASCII characters 0-127.
pub const char_table: [128]CharType = blk: {
    var table: [128]CharType = .{.ordinary} ** 128;
    // Control characters → whitespace
    for (0..33) |i| table[i] = .whitespace;
    table[127] = .whitespace; // DEL
    // Explicit whitespace
    table[' '] = .whitespace;
    table['\t'] = .whitespace;
    table['\n'] = .whitespace;
    table['\r'] = .whitespace;
    // Special characters
    table['#'] = .special;
    table['$'] = .special;
    table['\''] = .special;
    table['"'] = .special;
    table['_'] = .special;
    table[';'] = .special;
    break :blk table;
};

/// Check if character is CIF whitespace
pub fn isWhitespace(c: u8) bool {
    if (c > 127) return false;
    return char_table[c] == .whitespace;
}

/// Check if character is ordinary (valid in unquoted values)
pub fn isOrdinary(c: u8) bool {
    if (c > 127) return true; // Non-ASCII treated as ordinary
    return char_table[c] == .ordinary;
}
```

- [ ] **Step 2: Write failing tokenizer test**

```zig
// src/cif/tokenizer.zig
const std = @import("std");
const testing = std.testing;

test "tokenize data block" {
    const source = "data_TEST\n_tag value\n";
    var tok = Tokenizer.init(source);
    const t1 = tok.next();
    try testing.expectEqual(TokenType.data, t1.type);
    try testing.expectEqualStrings("data_TEST", t1.text(source));
}

test "tokenize loop" {
    const source =
        \\loop_
        \\_tag1
        \\_tag2
        \\val1 val2
        \\val3 val4
    ;
    var tok = Tokenizer.init(source);
    try testing.expectEqual(TokenType.loop, tok.next().type);
    try testing.expectEqual(TokenType.tag, tok.next().type);
    try testing.expectEqual(TokenType.tag, tok.next().type);
    try testing.expectEqual(TokenType.value, tok.next().type);
    try testing.expectEqual(TokenType.value, tok.next().type);
    try testing.expectEqual(TokenType.value, tok.next().type);
    try testing.expectEqual(TokenType.value, tok.next().type);
    try testing.expectEqual(TokenType.eof, tok.next().type);
}
```

- [ ] **Step 3: Run tests to verify they fail**

Run: `cd ~/zreduce && zig build test 2>&1 | head -5`
Expected: Compilation error — `Tokenizer` not defined.

- [ ] **Step 4: Implement tokenizer**

```zig
// src/cif/tokenizer.zig
const std = @import("std");
const char_table = @import("char_table.zig");

pub const TokenType = enum {
    data,        // data_BLOCKNAME
    loop,        // loop_
    tag,         // _tag_name
    value,       // unquoted, single-quoted, double-quoted, or semicolon-delimited
    save_begin,  // save_FRAMENAME
    save_end,    // save_
    eof,
};

pub const Token = struct {
    type: TokenType,
    start: u32,
    end: u32,

    pub fn text(self: Token, source: []const u8) []const u8 {
        return source[self.start..self.end];
    }
};

pub const Tokenizer = struct {
    source: []const u8,
    pos: u32,

    pub fn init(source: []const u8) Tokenizer {
        return .{ .source = source, .pos = 0 };
    }

    pub fn next(self: *Tokenizer) Token {
        self.skipWhitespaceAndComments();

        if (self.pos >= self.source.len) {
            return .{ .type = .eof, .start = self.pos, .end = self.pos };
        }

        const c = self.source[self.pos];

        // Semicolon-delimited text field (must be at start of line)
        if (c == ';' and self.isStartOfLine()) {
            return self.readTextField();
        }

        // Quoted string
        if (c == '\'' or c == '"') {
            return self.readQuotedString(c);
        }

        // Tag
        if (c == '_') {
            return self.readTag();
        }

        // Keywords or unquoted value
        const start = self.pos;
        self.advanceToWhitespace();
        const word = self.source[start..self.pos];

        if (std.ascii.eqlIgnoreCase(word, "loop_")) {
            return .{ .type = .loop, .start = start, .end = self.pos };
        }

        if (word.len > 5 and std.ascii.eqlIgnoreCase(word[0..5], "data_")) {
            return .{ .type = .data, .start = start, .end = self.pos };
        }

        if (word.len > 5 and std.ascii.eqlIgnoreCase(word[0..5], "save_")) {
            if (word.len == 5) {
                return .{ .type = .save_end, .start = start, .end = self.pos };
            }
            return .{ .type = .save_begin, .start = start, .end = self.pos };
        }

        return .{ .type = .value, .start = start, .end = self.pos };
    }

    fn skipWhitespaceAndComments(self: *Tokenizer) void {
        while (self.pos < self.source.len) {
            const c = self.source[self.pos];
            if (char_table.isWhitespace(c)) {
                self.pos += 1;
                continue;
            }
            if (c == '#') {
                // Skip to end of line
                while (self.pos < self.source.len and self.source[self.pos] != '\n') {
                    self.pos += 1;
                }
                continue;
            }
            break;
        }
    }

    fn isStartOfLine(self: *const Tokenizer) bool {
        if (self.pos == 0) return true;
        return self.source[self.pos - 1] == '\n';
    }

    fn readTextField(self: *Tokenizer) Token {
        // Skip opening semicolon
        self.pos += 1;
        const start = self.pos;
        // Skip to end of line
        while (self.pos < self.source.len and self.source[self.pos] != '\n') {
            self.pos += 1;
        }
        if (self.pos < self.source.len) self.pos += 1; // skip newline
        // Read until semicolon at start of line
        while (self.pos < self.source.len) {
            if (self.source[self.pos] == ';' and self.isStartOfLine()) {
                const end = self.pos;
                self.pos += 1; // skip closing semicolon
                return .{ .type = .value, .start = @intCast(start), .end = @intCast(end) };
            }
            self.pos += 1;
        }
        return .{ .type = .value, .start = @intCast(start), .end = self.pos };
    }

    fn readQuotedString(self: *Tokenizer, quote: u8) Token {
        self.pos += 1; // skip opening quote
        const start = self.pos;
        while (self.pos < self.source.len) {
            if (self.source[self.pos] == quote) {
                // Quote ends only if followed by whitespace or EOF
                if (self.pos + 1 >= self.source.len or char_table.isWhitespace(self.source[self.pos + 1])) {
                    const end = self.pos;
                    self.pos += 1; // skip closing quote
                    return .{ .type = .value, .start = @intCast(start), .end = @intCast(end) };
                }
            }
            self.pos += 1;
        }
        return .{ .type = .value, .start = @intCast(start), .end = self.pos };
    }

    fn readTag(self: *Tokenizer) Token {
        const start = self.pos;
        self.advanceToWhitespace();
        return .{ .type = .tag, .start = start, .end = self.pos };
    }

    fn advanceToWhitespace(self: *Tokenizer) void {
        while (self.pos < self.source.len and !char_table.isWhitespace(self.source[self.pos])) {
            self.pos += 1;
        }
    }
};
```

- [ ] **Step 5: Run tests**

Run: `cd ~/zreduce && zig build test`
Expected: All tests pass.

- [ ] **Step 6: Add more edge case tests**

```zig
test "tokenize quoted strings" {
    const source = "'hello world' \"another string\"";
    var tok = Tokenizer.init(source);
    const t1 = tok.next();
    try testing.expectEqual(TokenType.value, t1.type);
    try testing.expectEqualStrings("hello world", t1.text(source));
    const t2 = tok.next();
    try testing.expectEqualStrings("another string", t2.text(source));
}

test "tokenize semicolon text field" {
    const source = "data_X\n_tag\n;multi\nline\nvalue\n;\n";
    var tok = Tokenizer.init(source);
    _ = tok.next(); // data_X
    _ = tok.next(); // _tag
    const t = tok.next();
    try testing.expectEqual(TokenType.value, t.type);
    try testing.expectEqualStrings("multi\nline\nvalue\n", t.text(source));
}

test "skip comments" {
    const source = "# comment\ndata_X\n# another\n_tag val";
    var tok = Tokenizer.init(source);
    try testing.expectEqual(TokenType.data, tok.next().type);
    try testing.expectEqual(TokenType.tag, tok.next().type);
    try testing.expectEqual(TokenType.value, tok.next().type);
}
```

Run: `cd ~/zreduce && zig build test`
Expected: All tests pass.

- [ ] **Step 7: Commit**

```bash
cd ~/zreduce && git add src/cif/ && git commit -m "feat: add CIF tokenizer with char table"
```

---

### Task 5: CIF Parser (Document/Block/Loop)

**Files:**
- Create: `~/zreduce/src/cif/types.zig`
- Create: `~/zreduce/src/cif/value.zig`
- Create: `~/zreduce/src/cif/parser.zig`
- Create: `~/zreduce/src/cif.zig`

- [ ] **Step 1: Create types.zig**

```zig
// src/cif/types.zig
const std = @import("std");

pub const ItemType = enum {
    pair,
    loop,
};

pub const Pair = struct {
    tag: []const u8,
    value: []const u8,
};

pub const Loop = struct {
    tags: std.ArrayList([]const u8),
    values: std.ArrayList([]const u8),
    allocator: std.mem.Allocator,

    pub fn init(allocator: std.mem.Allocator) Loop {
        return .{
            .tags = .empty,
            .values = .empty,
            .allocator = allocator,
        };
    }

    pub fn deinit(self: *Loop) void {
        self.tags.deinit(self.allocator);
        self.values.deinit(self.allocator);
    }

    pub fn width(self: *const Loop) usize {
        return self.tags.items.len;
    }

    pub fn length(self: *const Loop) usize {
        const w = self.width();
        if (w == 0) return 0;
        return self.values.items.len / w;
    }

    pub fn val(self: *const Loop, row: usize, col: usize) ?[]const u8 {
        const w = self.width();
        if (w == 0 or col >= w) return null;
        const idx = row * w + col;
        if (idx >= self.values.items.len) return null;
        return self.values.items[idx];
    }

    pub fn findTag(self: *const Loop, tag: []const u8) ?usize {
        for (self.tags.items, 0..) |t, i| {
            if (std.ascii.eqlIgnoreCase(t, tag)) return i;
        }
        return null;
    }
};

pub const Item = union(ItemType) {
    pair: Pair,
    loop: Loop,

    pub fn deinit(self: *Item) void {
        switch (self.*) {
            .loop => |*l| l.deinit(),
            .pair => {},
        }
    }
};

pub const Block = struct {
    name: []const u8,
    items: std.ArrayList(Item),
    allocator: std.mem.Allocator,

    pub fn init(allocator: std.mem.Allocator, name: []const u8) Block {
        return .{
            .name = name,
            .items = .empty,
            .allocator = allocator,
        };
    }

    pub fn deinit(self: *Block) void {
        for (self.items.items) |*item| item.deinit();
        self.items.deinit(self.allocator);
    }

    /// Find first loop containing the given tag name.
    pub fn findLoop(self: *const Block, tag: []const u8) ?*const Loop {
        for (self.items.items) |*item| {
            switch (item.*) {
                .loop => |*l| {
                    if (l.findTag(tag) != null) return l;
                },
                .pair => {},
            }
        }
        return null;
    }

    /// Find value for a tag-value pair.
    pub fn findValue(self: *const Block, tag: []const u8) ?[]const u8 {
        for (self.items.items) |item| {
            switch (item) {
                .pair => |p| {
                    if (std.ascii.eqlIgnoreCase(p.tag, tag)) return p.value;
                },
                .loop => {},
            }
        }
        return null;
    }
};

pub const Document = struct {
    blocks: std.ArrayList(Block),
    allocator: std.mem.Allocator,

    pub fn init(allocator: std.mem.Allocator) Document {
        return .{
            .blocks = .empty,
            .allocator = allocator,
        };
    }

    pub fn deinit(self: *Document) void {
        for (self.blocks.items) |*block| block.deinit();
        self.blocks.deinit(self.allocator);
    }

    /// Find block by name (case-insensitive).
    pub fn findBlock(self: *const Document, name: []const u8) ?*const Block {
        for (self.blocks.items) |*block| {
            if (std.ascii.eqlIgnoreCase(block.name, name)) return block;
        }
        return null;
    }
};
```

- [ ] **Step 2: Create value.zig**

```zig
// src/cif/value.zig
const std = @import("std");

/// Check if CIF value represents null/unknown (. or ?)
pub fn isNull(v: []const u8) bool {
    return v.len == 1 and (v[0] == '.' or v[0] == '?');
}

/// Extract string value, stripping quotes if present.
pub fn asString(v: []const u8) []const u8 {
    if (isNull(v)) return "";
    return v;
}

/// Parse as float, returning null for null values.
pub fn asFloat(v: []const u8) ?f32 {
    if (isNull(v)) return null;
    // CIF floats may have trailing uncertainty in parens: "1.234(5)"
    var end = v.len;
    for (v, 0..) |c, i| {
        if (c == '(') { end = i; break; }
    }
    return std.fmt.parseFloat(f32, v[0..end]) catch null;
}

/// Parse as float with default value for null.
pub fn asFloatOr(v: []const u8, default: f32) f32 {
    return asFloat(v) orelse default;
}

/// Parse as integer, returning null for null values.
pub fn asInt(comptime T: type, v: []const u8) ?T {
    if (isNull(v)) return null;
    return std.fmt.parseInt(T, v, 10) catch null;
}

/// Parse as integer with default value for null.
pub fn asIntOr(comptime T: type, v: []const u8, default: T) T {
    return asInt(T, v) orelse default;
}

const testing = std.testing;

test "isNull" {
    try testing.expect(isNull("."));
    try testing.expect(isNull("?"));
    try testing.expect(!isNull("hello"));
    try testing.expect(!isNull(".."));
}

test "asFloat with uncertainty" {
    const v = asFloat("1.234(5)");
    try testing.expect(v != null);
    try testing.expectApproxEqAbs(v.?, 1.234, 1e-4);
}

test "asFloat null" {
    try testing.expect(asFloat(".") == null);
    try testing.expect(asFloat("?") == null);
}
```

- [ ] **Step 3: Create parser.zig**

```zig
// src/cif/parser.zig
const std = @import("std");
const types = @import("types.zig");
const tokenizer_mod = @import("tokenizer.zig");
const Tokenizer = tokenizer_mod.Tokenizer;
const TokenType = tokenizer_mod.TokenType;
const Token = tokenizer_mod.Token;
const Document = types.Document;
const Block = types.Block;
const Loop = types.Loop;
const Item = types.Item;

pub const ParseError = error{
    UnexpectedToken,
    OutOfMemory,
};

/// Parse CIF source text into a Document.
pub fn readString(allocator: std.mem.Allocator, source: []const u8) ParseError!Document {
    var doc = Document.init(allocator);
    errdefer doc.deinit();

    var tok = Tokenizer.init(source);
    var current_block: ?*Block = null;

    while (true) {
        const token = tok.next();
        switch (token.type) {
            .eof => break,

            .data => {
                const name_text = token.text(source);
                const name = if (name_text.len > 5) name_text[5..] else "";
                doc.blocks.append(doc.allocator, Block.init(allocator, name)) catch return error.OutOfMemory;
                current_block = &doc.blocks.items[doc.blocks.items.len - 1];
            },

            .tag => {
                if (current_block) |block| {
                    const tag_name = token.text(source);
                    const val_token = tok.next();
                    const val = if (val_token.type == .value) val_token.text(source) else "";
                    block.items.append(block.allocator, Item{ .pair = .{
                        .tag = tag_name,
                        .value = val,
                    } }) catch return error.OutOfMemory;
                }
            },

            .loop => {
                if (current_block) |block| {
                    var loop = Loop.init(allocator);
                    errdefer loop.deinit();

                    // Read tags
                    while (true) {
                        const t = tok.next();
                        if (t.type == .tag) {
                            loop.tags.append(loop.allocator, t.text(source)) catch return error.OutOfMemory;
                        } else {
                            // First non-tag token is start of values
                            if (t.type == .value) {
                                loop.values.append(loop.allocator, t.text(source)) catch return error.OutOfMemory;
                            } else {
                                // loop with no values (or next keyword)
                                // Put back by re-parsing? No, we just break and the
                                // outer loop will handle it on next iteration.
                                // For simplicity, we'll just add the loop as-is.
                                block.items.append(block.allocator, Item{ .loop = loop }) catch return error.OutOfMemory;
                                // Need to handle the token we just consumed.
                                // Since our tokenizer doesn't support pushback,
                                // we handle this case specially.
                                switch (t.type) {
                                    .data => {
                                        const n = t.text(source);
                                        const nm = if (n.len > 5) n[5..] else "";
                                        doc.blocks.append(doc.allocator, Block.init(allocator, nm)) catch return error.OutOfMemory;
                                        current_block = &doc.blocks.items[doc.blocks.items.len - 1];
                                    },
                                    .eof => return doc,
                                    else => {},
                                }
                                continue;
                            }
                            break;
                        }
                    }

                    // Read remaining values
                    const w = loop.width();
                    if (w > 0) {
                        while (true) {
                            const t = tok.next();
                            if (t.type == .value) {
                                loop.values.append(loop.allocator, t.text(source)) catch return error.OutOfMemory;
                            } else {
                                block.items.append(block.allocator, Item{ .loop = loop }) catch return error.OutOfMemory;
                                switch (t.type) {
                                    .data => {
                                        const n = t.text(source);
                                        const nm = if (n.len > 5) n[5..] else "";
                                        doc.blocks.append(doc.allocator, Block.init(allocator, nm)) catch return error.OutOfMemory;
                                        current_block = &doc.blocks.items[doc.blocks.items.len - 1];
                                    },
                                    .loop => {
                                        // Another loop follows — continue outer loop
                                        // We need to parse this loop now.
                                        // Simplest: re-enter loop parsing via recursion or goto.
                                        // For now, just handle it inline.
                                        var loop2 = Loop.init(allocator);
                                        errdefer loop2.deinit();
                                        while (true) {
                                            const t2 = tok.next();
                                            if (t2.type == .tag) {
                                                loop2.tags.append(loop2.allocator, t2.text(source)) catch return error.OutOfMemory;
                                            } else {
                                                if (t2.type == .value) {
                                                    loop2.values.append(loop2.allocator, t2.text(source)) catch return error.OutOfMemory;
                                                }
                                                break;
                                            }
                                        }
                                        // Read loop2 values
                                        while (true) {
                                            const t2 = tok.next();
                                            if (t2.type == .value) {
                                                loop2.values.append(loop2.allocator, t2.text(source)) catch return error.OutOfMemory;
                                            } else {
                                                block.items.append(block.allocator, Item{ .loop = loop2 }) catch return error.OutOfMemory;
                                                // Handle leftover token...
                                                // This recursive problem needs a cleaner solution.
                                                // TODO: refactor with pending token pattern
                                                break;
                                            }
                                        }
                                    },
                                    .tag => {
                                        // Tag-value pair after loop
                                        const tag_name = t.text(source);
                                        const val_token = tok.next();
                                        const val = if (val_token.type == .value) val_token.text(source) else "";
                                        block.items.append(block.allocator, Item{ .pair = .{
                                            .tag = tag_name,
                                            .value = val,
                                        } }) catch return error.OutOfMemory;
                                    },
                                    .eof => return doc,
                                    else => {},
                                }
                                break;
                            }
                        }
                    }
                }
            },

            .value, .save_begin, .save_end => {
                // Skip unexpected tokens
            },
        }
    }

    return doc;
}

const testing = std.testing;

test "parse simple document" {
    const source =
        \\data_TEST
        \\_entry.id TEST
        \\loop_
        \\_atom.x
        \\_atom.y
        \\1.0 2.0
        \\3.0 4.0
    ;
    var doc = try readString(testing.allocator, source);
    defer doc.deinit();

    try testing.expectEqual(@as(usize, 1), doc.blocks.items.len);
    const block = &doc.blocks.items[0];
    try testing.expectEqualStrings("TEST", block.name);

    // Check pair
    const entry_id = block.findValue("_entry.id");
    try testing.expect(entry_id != null);
    try testing.expectEqualStrings("TEST", entry_id.?);

    // Check loop
    const loop = block.findLoop("_atom.x");
    try testing.expect(loop != null);
    try testing.expectEqual(@as(usize, 2), loop.?.width());
    try testing.expectEqual(@as(usize, 2), loop.?.length());
    try testing.expectEqualStrings("1.0", loop.?.val(0, 0).?);
    try testing.expectEqualStrings("4.0", loop.?.val(1, 1).?);
}
```

**Note:** The parser above uses an inline approach for handling the pending token problem. During implementation, this should be refactored to use the `pending` token pattern from zig-cif-graph-parser's `component_dict.zig` for cleaner control flow.

- [ ] **Step 4: Create cif.zig module re-export**

```zig
// src/cif.zig
pub const types = @import("cif/types.zig");
pub const char_table = @import("cif/char_table.zig");
pub const tokenizer = @import("cif/tokenizer.zig");
pub const parser = @import("cif/parser.zig");
pub const value = @import("cif/value.zig");

pub const Document = types.Document;
pub const Block = types.Block;
pub const Loop = types.Loop;
pub const readString = parser.readString;
pub const isNull = value.isNull;
pub const asFloat = value.asFloat;
pub const asFloatOr = value.asFloatOr;
pub const asInt = value.asInt;
pub const asString = value.asString;

test {
    @import("std").testing.refAllDecls(@This());
}
```

- [ ] **Step 5: Run all tests**

Run: `cd ~/zreduce && zig build test`
Expected: All tests pass.

- [ ] **Step 6: Register in root.zig and commit**

Add to root.zig:
```zig
pub const cif = @import("cif.zig");
```

```bash
cd ~/zreduce && git add src/cif/ src/cif.zig src/root.zig && git commit -m "feat: add CIF parser (tokenizer, types, parser, value helpers)"
```

---

### Task 6: Model Structs

**Files:**
- Create: `~/zreduce/src/model/atom.zig`
- Create: `~/zreduce/src/model/residue.zig`
- Create: `~/zreduce/src/model/chain.zig`
- Create: `~/zreduce/src/model/bond.zig`
- Create: `~/zreduce/src/model/model.zig`
- Create: `~/zreduce/src/model/neighbor.zig`
- Create: `~/zreduce/src/model.zig`

- [ ] **Step 1: Create atom.zig**

```zig
// src/model/atom.zig
const math = @import("../math.zig");
const element = @import("../element.zig");

pub const AtomFlags = element.AtomFlags;

pub const Atom = struct {
    pos: math.Vec3(f32),
    name: [4]u8 = .{ ' ', ' ', ' ', ' ' },
    name_len: u4 = 0,
    element_type: element.AtomType = .unknown,
    residue_idx: u32 = 0,
    altloc: u8 = ' ',
    occupancy: f32 = 1.0,
    b_factor: f32 = 0.0,
    is_hydrogen: bool = false,
    is_added: bool = false,      // true if added by zreduce
    vdw_radius: f32 = 1.70,
    flags: AtomFlags = .{},
    serial: u32 = 0,

    pub fn nameSlice(self: *const Atom) []const u8 {
        return self.name[0..self.name_len];
    }

    pub fn setName(self: *Atom, name: []const u8) void {
        const len: u4 = @intCast(@min(name.len, 4));
        self.name = .{ ' ', ' ', ' ', ' ' };
        for (0..len) |i| self.name[i] = name[i];
        self.name_len = len;
    }
};
```

- [ ] **Step 2: Create residue.zig, chain.zig, bond.zig**

```zig
// src/model/residue.zig
pub const EntityType = enum {
    polymer,
    non_polymer,
    water,
    unknown,
};

pub const Residue = struct {
    comp_id: [3]u8 = .{ ' ', ' ', ' ' },
    comp_id_len: u3 = 0,
    chain_idx: u16 = 0,
    seq_id: i32 = 0,
    ins_code: u8 = ' ',
    atom_start: u32 = 0,
    atom_end: u32 = 0,
    entity_type: EntityType = .unknown,

    pub fn compIdSlice(self: *const Residue) []const u8 {
        return self.comp_id[0..self.comp_id_len];
    }

    pub fn setCompId(self: *Residue, id: []const u8) void {
        const len: u3 = @intCast(@min(id.len, 3));
        self.comp_id = .{ ' ', ' ', ' ' };
        for (0..len) |i| self.comp_id[i] = id[i];
        self.comp_id_len = len;
    }

    pub fn atomCount(self: *const Residue) u32 {
        return self.atom_end - self.atom_start;
    }
};
```

```zig
// src/model/chain.zig
pub const Chain = struct {
    label_asym_id: [4]u8 = .{ ' ', ' ', ' ', ' ' },
    label_asym_id_len: u4 = 0,
    auth_asym_id: [4]u8 = .{ ' ', ' ', ' ', ' ' },
    auth_asym_id_len: u4 = 0,
    entity_id: u16 = 0,
    residue_start: u32 = 0,
    residue_end: u32 = 0,

    pub fn labelSlice(self: *const Chain) []const u8 {
        return self.label_asym_id[0..self.label_asym_id_len];
    }
};
```

```zig
// src/model/bond.zig
pub const BondOrder = enum(u3) {
    single,
    double,
    triple,
    aromatic,
    delocalized,
    unknown,

    pub fn fromString(s: []const u8) BondOrder {
        if (s.len == 0) return .unknown;
        const std = @import("std");
        var buf: [4]u8 = undefined;
        const len = @min(s.len, 4);
        for (0..len) |i| buf[i] = std.ascii.toUpper(s[i]);
        const upper = buf[0..len];
        if (std.mem.eql(u8, upper, "SING")) return .single;
        if (std.mem.eql(u8, upper, "DOUB")) return .double;
        if (std.mem.eql(u8, upper, "TRIP")) return .triple;
        if (std.mem.eql(u8, upper, "AROM")) return .aromatic;
        if (std.mem.eql(u8, upper, "DELO")) return .delocalized;
        return .unknown;
    }
};

pub const BondSource = enum(u3) {
    component_template,
    struct_conn,
    polymer_backbone,
    branch_link,
    inferred,
};

pub const Bond = struct {
    atom_1: u32,
    atom_2: u32,
    order: BondOrder = .single,
    source: BondSource = .inferred,
};
```

- [ ] **Step 3: Create model.zig (aggregate)**

```zig
// src/model/model.zig
const std = @import("std");
const atom_mod = @import("atom.zig");
const residue_mod = @import("residue.zig");
const chain_mod = @import("chain.zig");
const bond_mod = @import("bond.zig");

pub const Atom = atom_mod.Atom;
pub const Residue = residue_mod.Residue;
pub const Chain = chain_mod.Chain;
pub const Bond = bond_mod.Bond;
pub const BondOrder = bond_mod.BondOrder;
pub const BondSource = bond_mod.BondSource;
pub const EntityType = residue_mod.EntityType;

pub const Model = struct {
    atoms: std.ArrayList(Atom),
    residues: std.ArrayList(Residue),
    chains: std.ArrayList(Chain),
    bonds: std.ArrayList(Bond),
    allocator: std.mem.Allocator,

    pub fn init(allocator: std.mem.Allocator) Model {
        return .{
            .atoms = .empty,
            .residues = .empty,
            .chains = .empty,
            .bonds = .empty,
            .allocator = allocator,
        };
    }

    pub fn deinit(self: *Model) void {
        self.atoms.deinit(self.allocator);
        self.residues.deinit(self.allocator);
        self.chains.deinit(self.allocator);
        self.bonds.deinit(self.allocator);
    }

    /// Get atoms belonging to a residue.
    pub fn residueAtoms(self: *const Model, res: *const Residue) []const Atom {
        return self.atoms.items[res.atom_start..res.atom_end];
    }

    /// Get mutable atoms belonging to a residue.
    pub fn residueAtomsMut(self: *Model, res: *const Residue) []Atom {
        return self.atoms.items[res.atom_start..res.atom_end];
    }

    /// Find atom by name within a residue.
    pub fn findAtomInResidue(self: *const Model, res: *const Residue, name: []const u8) ?u32 {
        for (res.atom_start..res.atom_end) |i| {
            if (std.mem.eql(u8, self.atoms.items[i].nameSlice(), name)) {
                return @intCast(i);
            }
        }
        return null;
    }
};
```

- [ ] **Step 4: Create neighbor.zig (spatial cell list)**

```zig
// src/model/neighbor.zig
const std = @import("std");
const math = @import("../math.zig");

const Vec3 = math.Vec3(f32);
const Allocator = std.mem.Allocator;

pub const CellList = struct {
    atom_indices: []u32,
    cell_offsets: []u32,
    nx: u32,
    ny: u32,
    nz: u32,
    cell_size: f32,
    x_min: f32,
    y_min: f32,
    z_min: f32,
    allocator: Allocator,

    pub fn init(allocator: Allocator, positions: []const Vec3, cell_size: f32) !CellList {
        if (positions.len == 0) return error.NoAtoms;

        // Compute bounding box
        var x_min = positions[0].x;
        var x_max = positions[0].x;
        var y_min = positions[0].y;
        var y_max = positions[0].y;
        var z_min = positions[0].z;
        var z_max = positions[0].z;

        for (positions[1..]) |p| {
            x_min = @min(x_min, p.x); x_max = @max(x_max, p.x);
            y_min = @min(y_min, p.y); y_max = @max(y_max, p.y);
            z_min = @min(z_min, p.z); z_max = @max(z_max, p.z);
        }

        const nx: u32 = @intFromFloat(@max(1.0, (x_max - x_min) / cell_size + 1.0));
        const ny: u32 = @intFromFloat(@max(1.0, (y_max - y_min) / cell_size + 1.0));
        const nz: u32 = @intFromFloat(@max(1.0, (z_max - z_min) / cell_size + 1.0));
        const n_cells = nx * ny * nz;

        // Counting sort
        const counts = try allocator.alloc(u32, n_cells + 1);
        defer allocator.free(counts);
        @memset(counts, 0);

        for (positions, 0..) |p, i| {
            _ = i;
            const ci = cellIndex(p.x, p.y, p.z, x_min, y_min, z_min, cell_size, nx, ny, nz);
            counts[ci] += 1;
        }

        // Prefix sum
        const offsets = try allocator.alloc(u32, n_cells + 1);
        offsets[0] = 0;
        for (1..n_cells + 1) |i| offsets[i] = offsets[i - 1] + counts[i - 1];

        const indices = try allocator.alloc(u32, positions.len);
        const temp_offsets = try allocator.alloc(u32, n_cells);
        defer allocator.free(temp_offsets);
        @memcpy(temp_offsets, offsets[0..n_cells]);

        for (positions, 0..) |p, i| {
            const ci = cellIndex(p.x, p.y, p.z, x_min, y_min, z_min, cell_size, nx, ny, nz);
            indices[temp_offsets[ci]] = @intCast(i);
            temp_offsets[ci] += 1;
        }

        return .{
            .atom_indices = indices,
            .cell_offsets = offsets,
            .nx = nx, .ny = ny, .nz = nz,
            .cell_size = cell_size,
            .x_min = x_min, .y_min = y_min, .z_min = z_min,
            .allocator = allocator,
        };
    }

    pub fn deinit(self: *CellList) void {
        self.allocator.free(self.atom_indices);
        self.allocator.free(self.cell_offsets);
    }

    /// Iterate neighbor atom indices within radius of a query point.
    pub fn neighborsInRadius(self: *const CellList, query: Vec3, radius: f32, result: *std.ArrayList(u32), positions: []const Vec3) !void {
        const r2 = radius * radius;
        const r_cells: i32 = @intFromFloat(radius / self.cell_size + 1.0);

        const cx: i32 = @intFromFloat((query.x - self.x_min) / self.cell_size);
        const cy: i32 = @intFromFloat((query.y - self.y_min) / self.cell_size);
        const cz: i32 = @intFromFloat((query.z - self.z_min) / self.cell_size);

        var dz: i32 = -r_cells;
        while (dz <= r_cells) : (dz += 1) {
            const gz = cz + dz;
            if (gz < 0 or gz >= @as(i32, @intCast(self.nz))) continue;
            var dy: i32 = -r_cells;
            while (dy <= r_cells) : (dy += 1) {
                const gy = cy + dy;
                if (gy < 0 or gy >= @as(i32, @intCast(self.ny))) continue;
                var dx: i32 = -r_cells;
                while (dx <= r_cells) : (dx += 1) {
                    const gx = cx + dx;
                    if (gx < 0 or gx >= @as(i32, @intCast(self.nx))) continue;

                    const ci: u32 = @intCast(@as(u32, @intCast(gz)) * self.nx * self.ny + @as(u32, @intCast(gy)) * self.nx + @as(u32, @intCast(gx)));
                    const start = self.cell_offsets[ci];
                    const end = self.cell_offsets[ci + 1];

                    for (self.atom_indices[start..end]) |idx| {
                        const d = query.sub(positions[idx]);
                        if (d.dot(d) <= r2) {
                            try result.append(result.allocator, idx);
                        }
                    }
                }
            }
        }
    }
};

fn cellIndex(x: f32, y: f32, z: f32, x_min: f32, y_min: f32, z_min: f32, cs: f32, nx: u32, ny: u32, nz: u32) u32 {
    const ix: u32 = @intFromFloat(@max(0.0, (x - x_min) / cs));
    const iy: u32 = @intFromFloat(@max(0.0, (y - y_min) / cs));
    const iz: u32 = @intFromFloat(@max(0.0, (z - z_min) / cs));
    return @min(iz, nz - 1) * nx * ny + @min(iy, ny - 1) * nx + @min(ix, nx - 1);
}

const testing = std.testing;

test "CellList basic neighbor search" {
    const positions = [_]Vec3{
        .{ .x = 0.0, .y = 0.0, .z = 0.0 },
        .{ .x = 1.0, .y = 0.0, .z = 0.0 },
        .{ .x = 10.0, .y = 10.0, .z = 10.0 },
    };
    var cl = try CellList.init(testing.allocator, &positions, 3.0);
    defer cl.deinit();

    var result: std.ArrayList(u32) = .empty;
    defer result.deinit(testing.allocator);

    try cl.neighborsInRadius(.{ .x = 0.0, .y = 0.0, .z = 0.0 }, 2.0, &result, &positions);
    try testing.expectEqual(@as(usize, 2), result.items.len);
}
```

- [ ] **Step 5: Create model.zig re-export**

```zig
// src/model.zig
pub const atom = @import("model/atom.zig");
pub const residue = @import("model/residue.zig");
pub const chain = @import("model/chain.zig");
pub const bond = @import("model/bond.zig");
pub const model = @import("model/model.zig");
pub const neighbor = @import("model/neighbor.zig");

pub const Atom = model.Atom;
pub const Residue = model.Residue;
pub const Chain = model.Chain;
pub const Bond = model.Bond;
pub const Model = model.Model;
pub const CellList = neighbor.CellList;

test {
    @import("std").testing.refAllDecls(@This());
}
```

- [ ] **Step 6: Register in root.zig, run tests, commit**

Add to root.zig:
```zig
pub const model = @import("model.zig");
```

Run: `cd ~/zreduce && zig build test`
Expected: All tests pass.

```bash
cd ~/zreduce && git add src/model/ src/model.zig src/root.zig && git commit -m "feat: add molecular model structs with spatial neighbor list"
```

---

### Task 7: mmCIF Extraction (atom_site → Model)

**Files:**
- Create: `~/zreduce/src/mmcif.zig`
- Create: `~/zreduce/test_data/tiny.cif`

- [ ] **Step 1: Create minimal test mmCIF file**

```
data_TINY
#
_entry.id TINY
#
loop_
_atom_site.group_PDB
_atom_site.id
_atom_site.type_symbol
_atom_site.label_atom_id
_atom_site.label_comp_id
_atom_site.label_asym_id
_atom_site.label_seq_id
_atom_site.Cartn_x
_atom_site.Cartn_y
_atom_site.Cartn_z
_atom_site.occupancy
_atom_site.B_iso_or_equiv
_atom_site.label_alt_id
_atom_site.auth_asym_id
ATOM 1  N N   ALA A 1 1.000 2.000 3.000 1.00 10.0 . A
ATOM 2  C CA  ALA A 1 2.000 3.000 4.000 1.00 10.0 . A
ATOM 3  C C   ALA A 1 3.000 4.000 5.000 1.00 10.0 . A
ATOM 4  O O   ALA A 1 4.000 5.000 6.000 1.00 10.0 . A
ATOM 5  C CB  ALA A 1 2.500 2.500 3.500 1.00 10.0 . A
#
```

- [ ] **Step 2: Write failing test for mmCIF extraction**

```zig
test "parse tiny mmCIF" {
    const source = @embedFile("../test_data/tiny.cif");
    var mdl = try parseModel(testing.allocator, source);
    defer mdl.deinit();

    try testing.expectEqual(@as(usize, 5), mdl.atoms.items.len);
    try testing.expectEqual(@as(usize, 1), mdl.residues.items.len);
    try testing.expectEqual(@as(usize, 1), mdl.chains.items.len);

    // Check first atom coordinates
    try testing.expectApproxEqAbs(mdl.atoms.items[0].pos.x, 1.0, 1e-3);
    try testing.expectApproxEqAbs(mdl.atoms.items[0].pos.y, 2.0, 1e-3);

    // Check residue
    try testing.expectEqualStrings("ALA", mdl.residues.items[0].compIdSlice());
}
```

- [ ] **Step 3: Implement mmcif.zig**

```zig
// src/mmcif.zig
const std = @import("std");
const cif = @import("cif.zig");
const mdl = @import("model.zig");
const elem_mod = @import("element.zig");
const Model = mdl.Model;
const Atom = mdl.Atom;
const Residue = mdl.Residue;
const Chain = mdl.Chain;

pub const MmcifError = error{
    NoAtomSiteLoop,
    MissingCoordinateField,
    OutOfMemory,
};

/// Column indices for _atom_site fields.
const AtomSiteColumns = struct {
    cartn_x: ?usize = null,
    cartn_y: ?usize = null,
    cartn_z: ?usize = null,
    type_symbol: ?usize = null,
    label_atom_id: ?usize = null,
    label_comp_id: ?usize = null,
    label_asym_id: ?usize = null,
    label_seq_id: ?usize = null,
    auth_asym_id: ?usize = null,
    label_alt_id: ?usize = null,
    occupancy: ?usize = null,
    b_factor: ?usize = null,
    group_pdb: ?usize = null,
    id: ?usize = null,
};

/// Parse mmCIF source into a Model.
pub fn parseModel(allocator: std.mem.Allocator, source: []const u8) MmcifError!Model {
    var doc = cif.readString(allocator, source) catch return error.OutOfMemory;
    defer doc.deinit();

    var model = Model.init(allocator);
    errdefer model.deinit();

    // Find first block
    if (doc.blocks.items.len == 0) return model;
    const block = &doc.blocks.items[0];

    // Find _atom_site loop
    const loop = block.findLoop("_atom_site.Cartn_x") orelse return error.NoAtomSiteLoop;

    // Map columns
    var cols = AtomSiteColumns{};
    cols.cartn_x = loop.findTag("_atom_site.Cartn_x");
    cols.cartn_y = loop.findTag("_atom_site.Cartn_y");
    cols.cartn_z = loop.findTag("_atom_site.Cartn_z");
    cols.type_symbol = loop.findTag("_atom_site.type_symbol");
    cols.label_atom_id = loop.findTag("_atom_site.label_atom_id");
    cols.label_comp_id = loop.findTag("_atom_site.label_comp_id");
    cols.label_asym_id = loop.findTag("_atom_site.label_asym_id");
    cols.label_seq_id = loop.findTag("_atom_site.label_seq_id");
    cols.auth_asym_id = loop.findTag("_atom_site.auth_asym_id");
    cols.label_alt_id = loop.findTag("_atom_site.label_alt_id");
    cols.occupancy = loop.findTag("_atom_site.occupancy");
    cols.b_factor = loop.findTag("_atom_site.B_iso_or_equiv");
    cols.group_pdb = loop.findTag("_atom_site.group_PDB");
    cols.id = loop.findTag("_atom_site.id");

    if (cols.cartn_x == null or cols.cartn_y == null or cols.cartn_z == null) {
        return error.MissingCoordinateField;
    }

    // Parse rows
    var prev_chain_id: ?[]const u8 = null;
    var prev_seq_id: ?i32 = null;
    var prev_comp_id: ?[]const u8 = null;

    for (0..loop.length()) |row| {
        const x = cif.asFloatOr(loop.val(row, cols.cartn_x.?).?, 0.0);
        const y = cif.asFloatOr(loop.val(row, cols.cartn_y.?).?, 0.0);
        const z = cif.asFloatOr(loop.val(row, cols.cartn_z.?).?, 0.0);

        var atom = Atom{ .pos = .{ .x = x, .y = y, .z = z } };

        if (cols.type_symbol) |col| {
            const sym = cif.asString(loop.val(row, col) orelse "");
            atom.element_type = elem_mod.elementFromSymbol(sym);
            atom.is_hydrogen = (sym.len >= 1 and (sym[0] == 'H' or sym[0] == 'h'));
        }

        if (cols.label_atom_id) |col| {
            atom.setName(cif.asString(loop.val(row, col) orelse ""));
        }

        if (cols.label_alt_id) |col| {
            const alt = cif.asString(loop.val(row, col) orelse "");
            if (alt.len > 0 and alt[0] != '.') atom.altloc = alt[0];
        }

        if (cols.occupancy) |col| atom.occupancy = cif.asFloatOr(loop.val(row, col) orelse "1.0", 1.0);
        if (cols.b_factor) |col| atom.b_factor = cif.asFloatOr(loop.val(row, col) orelse "0.0", 0.0);
        if (cols.id) |col| atom.serial = @intCast(cif.asIntOr(u32, loop.val(row, col) orelse "0", 0));

        // Track residue/chain boundaries
        const chain_id = if (cols.label_asym_id) |col| (loop.val(row, col) orelse "") else "";
        const seq_id = if (cols.label_seq_id) |col| cif.asIntOr(i32, loop.val(row, col) orelse "", 0) else 0;
        const comp_id = if (cols.label_comp_id) |col| cif.asString(loop.val(row, col) orelse "") else "";

        // New chain?
        const new_chain = (prev_chain_id == null or !std.mem.eql(u8, chain_id, prev_chain_id.?));
        if (new_chain) {
            if (model.chains.items.len > 0) {
                model.chains.items[model.chains.items.len - 1].residue_end = @intCast(model.residues.items.len);
            }
            var chain = Chain{};
            const cid = cif.asString(chain_id);
            const len: u4 = @intCast(@min(cid.len, 4));
            for (0..len) |i| chain.label_asym_id[i] = cid[i];
            chain.label_asym_id_len = len;
            chain.residue_start = @intCast(model.residues.items.len);
            model.chains.append(model.allocator, chain) catch return error.OutOfMemory;
            prev_chain_id = chain_id;
            prev_seq_id = null; // force new residue
        }

        // New residue?
        const new_res = (prev_seq_id == null or seq_id != prev_seq_id.? or !std.mem.eql(u8, comp_id, prev_comp_id orelse ""));
        if (new_res) {
            if (model.residues.items.len > 0) {
                model.residues.items[model.residues.items.len - 1].atom_end = @intCast(model.atoms.items.len);
            }
            var res = Residue{};
            res.setCompId(comp_id);
            res.chain_idx = @intCast(model.chains.items.len - 1);
            res.seq_id = seq_id;
            res.atom_start = @intCast(model.atoms.items.len);
            model.residues.append(model.allocator, res) catch return error.OutOfMemory;
            prev_seq_id = seq_id;
            prev_comp_id = comp_id;
        }

        atom.residue_idx = @intCast(model.residues.items.len - 1);
        model.atoms.append(model.allocator, atom) catch return error.OutOfMemory;
    }

    // Close last residue and chain
    if (model.residues.items.len > 0) {
        model.residues.items[model.residues.items.len - 1].atom_end = @intCast(model.atoms.items.len);
    }
    if (model.chains.items.len > 0) {
        model.chains.items[model.chains.items.len - 1].residue_end = @intCast(model.residues.items.len);
    }

    return model;
}

const testing = std.testing;

test "parse tiny mmCIF" {
    const source = @embedFile("../test_data/tiny.cif");
    var mdl_result = try parseModel(testing.allocator, source);
    defer mdl_result.deinit();

    try testing.expectEqual(@as(usize, 5), mdl_result.atoms.items.len);
    try testing.expectEqual(@as(usize, 1), mdl_result.residues.items.len);
    try testing.expectEqual(@as(usize, 1), mdl_result.chains.items.len);
    try testing.expectApproxEqAbs(mdl_result.atoms.items[0].pos.x, 1.0, 1e-3);
    try testing.expectEqualStrings("ALA", mdl_result.residues.items[0].compIdSlice());
}
```

- [ ] **Step 4: Run tests**

Run: `cd ~/zreduce && zig build test`
Expected: All tests pass.

- [ ] **Step 5: Register and commit**

Add to root.zig:
```zig
pub const mmcif = @import("mmcif.zig");
```

```bash
cd ~/zreduce && git add src/mmcif.zig test_data/tiny.cif src/root.zig && git commit -m "feat: add mmCIF atom_site extraction into Model"
```

---

## Phase 2: Hydrogen Placement

### Task 8: Placement Geometry Functions (Type 1-6)

**Files:**
- Create: `~/zreduce/src/place/geometry.zig`

- [ ] **Step 1: Write failing tests for type 1 (tetrahedral) placement**

```zig
test "type1 HXR3 tetrahedral" {
    // Place H on sp3 carbon with 3 neighbors (like CA-HA)
    const center = Vec3(f64){ .x = 0.0, .y = 0.0, .z = 0.0 };
    const n1 = Vec3(f64){ .x = 1.0, .y = 0.0, .z = 0.0 };
    const n2 = Vec3(f64){ .x = 0.0, .y = 1.0, .z = 0.0 };
    const n3 = Vec3(f64){ .x = 0.0, .y = 0.0, .z = 1.0 };
    const h = placeHXR3(center, n1, n2, n3, 1.10);
    // H should be opposite to the 3 neighbors
    try testing.expect(h.x < 0);
    try testing.expect(h.y < 0);
    try testing.expect(h.z < 0);
    // Bond length should be 1.10
    try testing.expectApproxEqAbs(h.distance(center), 1.10, 1e-3);
}
```

- [ ] **Step 2: Implement all 6 placement types**

Direct port of AtomConn.cpp. See design spec for exact vector math. Each function takes reference atom positions and returns the hydrogen position.

```zig
// src/place/geometry.zig
const std = @import("std");
const math_mod = @import("../math.zig");
const Vec3 = math_mod.Vec3;

/// Type 1 (HXR3): Tetrahedral — H opposite to 3 neighbors.
pub fn placeHXR3(center: Vec3(f64), n1: Vec3(f64), n2: Vec3(f64), n3: Vec3(f64), bond_len: f64) Vec3(f64) {
    const v1 = n1.sub(center).normalize();
    const v2 = n2.sub(center).normalize();
    const v3 = n3.sub(center).normalize();
    const dir = v1.add(v2).add(v3).scaleTo(-bond_len);
    return center.add(dir);
}

/// Type 2 (H2XR2): Two H on sp2 atom. Returns one H; call twice with different fudge for both.
pub fn placeH2XR2(center: Vec3(f64), n1: Vec3(f64), n2: Vec3(f64), bond_len: f64, angle_deg: f64, dihedral_deg: f64) Vec3(f64) {
    const v1 = n1.sub(center).normalize();
    const v2 = n2.sub(center).normalize();
    const between = center.add(v1.add(v2).scale(0.5));
    return placeH3XR(center, between, n1, bond_len, angle_deg, dihedral_deg);
}

/// Type 3 (H3XR): Dihedral-controlled placement.
pub fn placeH3XR(a1: Vec3(f64), a2: Vec3(f64), a3: Vec3(f64), bond_len: f64, theta_deg: f64, phi_deg: f64) Vec3(f64) {
    const v21 = a1.sub(a2).normalize();
    const v23 = a3.sub(a2).normalize();
    const norm = v21.cross(v23).scaleTo(bond_len);

    const pos4 = math_mod.rotateAroundAxis(f64, a1.add(norm), a2, a1.sub(a2), phi_deg - 90.0);

    const v14 = pos4.sub(a1).normalize();
    const v12 = a2.sub(a1).normalize();
    const pos5 = v14.cross(v12).add(a1);

    return math_mod.rotateAroundAxis(f64, pos4, a1, pos5.sub(a1), 90.0 - theta_deg);
}

/// Type 4 (HXR2): Planar bisector with fudge factor.
pub fn placeHXR2Planar(center: Vec3(f64), n1: Vec3(f64), n2: Vec3(f64), bond_len: f64, fudge: f64) Vec3(f64) {
    const v1 = n1.sub(center).normalize();
    const v2 = n2.sub(center).normalize();
    // Interpolate with fudge: 0.5 + fudge biases toward n1 or n2
    const t = 0.5 + fudge;
    const interp = v1.scale(1.0 - t).add(v2.scale(t));
    return center.add(interp.scaleTo(-bond_len));
}

/// Type 5 (HXR2): Fractional angle distribution.
pub fn placeHXR2Frac(a1: Vec3(f64), a2: Vec3(f64), a3: Vec3(f64), bond_len: f64, fract: f64) Vec3(f64) {
    const v12 = a2.sub(a1).scaleTo(bond_len);
    const v13 = a3.sub(a1).scaleTo(bond_len);
    const pos4 = v12.cross(v13).add(a1);
    const cnca_angle = math_mod.angle(f64, a2, a1, a3);
    const hnca_angle = fract * (360.0 - cnca_angle);
    return math_mod.rotateAroundAxis(f64, a1.add(v12), pos4, a1.sub(pos4), hnca_angle);
}

/// Type 6 (HXY): Linear extension.
pub fn placeHXY(center: Vec3(f64), neighbor_atom: Vec3(f64), bond_len: f64) Vec3(f64) {
    const dir = center.sub(neighbor_atom).scaleTo(bond_len);
    return center.add(dir);
}

const testing = std.testing;

test "type1 HXR3 tetrahedral" {
    const center = Vec3(f64){ .x = 0.0, .y = 0.0, .z = 0.0 };
    const n1 = Vec3(f64){ .x = 1.0, .y = 0.0, .z = 0.0 };
    const n2 = Vec3(f64){ .x = 0.0, .y = 1.0, .z = 0.0 };
    const n3 = Vec3(f64){ .x = 0.0, .y = 0.0, .z = 1.0 };
    const h = placeHXR3(center, n1, n2, n3, 1.10);
    try testing.expect(h.x < 0);
    try testing.expect(h.y < 0);
    try testing.expect(h.z < 0);
    try testing.expectApproxEqAbs(h.distance(center), 1.10, 1e-3);
}

test "type6 HXY linear" {
    const center = Vec3(f64){ .x = 0.0, .y = 0.0, .z = 0.0 };
    const neighbor_atom = Vec3(f64){ .x = -1.5, .y = 0.0, .z = 0.0 };
    const h = placeHXY(center, neighbor_atom, 1.0);
    try testing.expectApproxEqAbs(h.x, 1.0, 1e-6);
    try testing.expectApproxEqAbs(h.y, 0.0, 1e-6);
}
```

- [ ] **Step 3: Run tests**

Run: `cd ~/zreduce && zig build test`
Expected: All tests pass.

- [ ] **Step 4: Commit**

```bash
cd ~/zreduce && git add src/place/geometry.zig && git commit -m "feat: add type 1-6 hydrogen placement geometry functions"
```

---

### Task 9: Standard Residue Placement Plans

**Files:**
- Create: `~/zreduce/src/place/standard.zig`

- [ ] **Step 1: Define PlacementPlan struct and MoverHint**

```zig
// src/place/standard.zig
const std = @import("std");
const element = @import("../element.zig");

pub const PlacementType = enum(u3) {
    hxr3,           // Type 1: tetrahedral (sp3, 3 neighbors)
    h2xr2,          // Type 2: two H on sp2
    h3xr,           // Type 3: dihedral-controlled
    hxr2_planar,    // Type 4: planar bisector
    hxr2_frac,      // Type 5: fractional angle
    hxy,            // Type 6: linear
};

pub const MoverHint = enum(u3) {
    none,
    rotate,          // OH, SH
    rotate_nh3,      // NH3+
    rotate_methyl,   // CH3
    flip_amide,      // Asn/Gln
    flip_his,        // His
};

pub const PlacementPlan = struct {
    h_name: [4]u8,
    placement_type: PlacementType,
    connected: [3][4]u8,      // reference atom names (up to 3)
    n_connected: u2,          // number of reference atoms used
    bond_len: f32,
    angle: f32 = 0.0,
    dihedral: f32 = 0.0,
    fudge: f32 = 0.0,
    atom_type: element.AtomType,
    mover_hint: MoverHint = .none,
};

fn name(comptime s: []const u8) [4]u8 {
    var buf: [4]u8 = .{ ' ', ' ', ' ', ' ' };
    for (s, 0..) |c, i| {
        if (i >= 4) break;
        buf[i] = c;
    }
    return buf;
}
```

- [ ] **Step 2: Add ALA plans as first residue**

```zig
const ala_plans = [_]PlacementPlan{
    // HA: tetrahedral on CA, neighbors N, C, CB
    .{ .h_name = name(" HA "), .placement_type = .hxr3,
       .connected = .{ name(" N  "), name(" C  "), name(" CB ") }, .n_connected = 3,
       .bond_len = 1.10, .atom_type = .H },
    // HB1, HB2, HB3: methyl on CB, dihedral-controlled
    .{ .h_name = name(" HB1"), .placement_type = .h3xr,
       .connected = .{ name(" CB "), name(" CA "), name(" N  ") }, .n_connected = 3,
       .bond_len = 1.10, .angle = 109.5, .dihedral = 180.0,
       .atom_type = .H, .mover_hint = .rotate_methyl },
    .{ .h_name = name(" HB2"), .placement_type = .h3xr,
       .connected = .{ name(" CB "), name(" CA "), name(" N  ") }, .n_connected = 3,
       .bond_len = 1.10, .angle = 109.5, .dihedral = 60.0,
       .atom_type = .H, .mover_hint = .rotate_methyl },
    .{ .h_name = name(" HB3"), .placement_type = .h3xr,
       .connected = .{ name(" CB "), name(" CA "), name(" N  ") }, .n_connected = 3,
       .bond_len = 1.10, .angle = 109.5, .dihedral = -60.0,
       .atom_type = .H, .mover_hint = .rotate_methyl },
    // Backbone H on N
    .{ .h_name = name(" H  "), .placement_type = .h3xr,
       .connected = .{ name(" N  "), name(" CA "), name(" C  ") }, .n_connected = 3,
       .bond_len = 1.02, .angle = 119.0, .dihedral = 180.0,
       .atom_type = .Hpol },
};
```

- [ ] **Step 3: Add remaining amino acids (GLY, VAL, LEU, ILE, etc.)**

This step adds the hardcoded plans for all 20 standard amino acids. Each amino acid follows the same pattern as ALA above. Key residues with special movers:

- **SER/THR/TYR**: OH hydrogen with `mover_hint = .rotate`
- **CYS**: SH hydrogen with `mover_hint = .rotate`
- **LYS**: NH3+ with `mover_hint = .rotate_nh3`
- **ASN/GLN**: amide H with `mover_hint = .flip_amide`
- **HIS**: ring H with `mover_hint = .flip_his`
- **ARG**: guanidinium H (no mover needed, rigid)
- **TRP**: ring NH with donor flag

Build up the complete table iteratively for each residue.

- [ ] **Step 4: Create lookup function**

```zig
/// Get placement plans for a standard residue. Returns null for non-standard.
pub fn getPlans(comp_id: []const u8) ?[]const PlacementPlan {
    const map = std.StaticStringMap([]const PlacementPlan).initComptime(.{
        .{ "ALA", &ala_plans },
        .{ "GLY", &gly_plans },
        .{ "VAL", &val_plans },
        .{ "LEU", &leu_plans },
        .{ "ILE", &ile_plans },
        .{ "PRO", &pro_plans },
        .{ "PHE", &phe_plans },
        .{ "TYR", &tyr_plans },
        .{ "TRP", &trp_plans },
        .{ "SER", &ser_plans },
        .{ "THR", &thr_plans },
        .{ "CYS", &cys_plans },
        .{ "MET", &met_plans },
        .{ "ASP", &asp_plans },
        .{ "GLU", &glu_plans },
        .{ "ASN", &asn_plans },
        .{ "GLN", &gln_plans },
        .{ "LYS", &lys_plans },
        .{ "ARG", &arg_plans },
        .{ "HIS", &his_plans },
    });
    return map.get(comp_id);
}

test "ALA plans" {
    const plans = getPlans("ALA");
    try std.testing.expect(plans != null);
    try std.testing.expectEqual(@as(usize, 5), plans.?.len);
}
```

- [ ] **Step 5: Run tests, commit**

Run: `cd ~/zreduce && zig build test`
Expected: All tests pass.

```bash
cd ~/zreduce && git add src/place/standard.zig && git commit -m "feat: add hardcoded H placement plans for standard amino acids"
```

---

### Task 10: CCD Loader

**Files:**
- Create: `~/zreduce/src/ccd.zig`

- [ ] **Step 1: Define Component/ComponentDict structs**

```zig
// src/ccd.zig
const std = @import("std");
const cif = @import("cif.zig");

pub const CompAtom = struct {
    name: [4]u8 = .{ ' ', ' ', ' ', ' ' },
    name_len: u4 = 0,
    element_symbol: [2]u8 = .{ ' ', ' ' },
    charge: i8 = 0,
    leaving: bool = false,
    aromatic: bool = false,
    ideal_x: f32 = 0.0,
    ideal_y: f32 = 0.0,
    ideal_z: f32 = 0.0,

    pub fn nameSlice(self: *const CompAtom) []const u8 {
        return self.name[0..self.name_len];
    }
};

pub const CompBond = struct {
    atom_idx_1: u16,
    atom_idx_2: u16,
    order: BondOrder,
    aromatic: bool = false,
};

pub const BondOrder = enum(u3) {
    single, double, triple, aromatic, delocalized, unknown,

    pub fn fromString(s: []const u8) BondOrder {
        if (s.len == 0) return .unknown;
        var buf: [4]u8 = undefined;
        const len = @min(s.len, 4);
        for (0..len) |i| buf[i] = std.ascii.toUpper(s[i]);
        const upper = buf[0..len];
        if (std.mem.eql(u8, upper, "SING")) return .single;
        if (std.mem.eql(u8, upper, "DOUB")) return .double;
        if (std.mem.eql(u8, upper, "TRIP")) return .triple;
        if (std.mem.eql(u8, upper, "AROM")) return .aromatic;
        if (std.mem.eql(u8, upper, "DELO")) return .delocalized;
        return .unknown;
    }
};

pub const Component = struct {
    comp_id: []const u8,
    comp_type: []const u8,
    atoms: []CompAtom,
    bonds: []CompBond,
};

pub const ComponentDict = struct {
    components: std.StringHashMap(Component),
    allocator: std.mem.Allocator,

    pub fn deinit(self: *ComponentDict) void {
        var it = self.components.iterator();
        while (it.next()) |entry| {
            self.allocator.free(entry.value_ptr.atoms);
            self.allocator.free(entry.value_ptr.bonds);
        }
        self.components.deinit();
    }

    pub fn get(self: *const ComponentDict, comp_id: []const u8) ?Component {
        return self.components.get(comp_id);
    }
};
```

- [ ] **Step 2: Implement streaming CCD parser**

Uses the streaming tokenizer pattern from zig-cif-graph-parser (does NOT build full Document since CCD has ~40K blocks).

```zig
/// Parse components.cif content (streaming, does not build full Document).
pub fn parseComponentDict(allocator: std.mem.Allocator, source: []const u8) !ComponentDict {
    var dict = ComponentDict{
        .components = std.StringHashMap(Component).init(allocator),
        .allocator = allocator,
    };
    errdefer dict.deinit();

    const Tok = cif.tokenizer.Tokenizer;
    var tok = Tok.init(source);

    var current_comp_id: ?[]const u8 = null;
    var current_comp_type: []const u8 = "";
    var current_atoms: std.ArrayList(CompAtom) = .empty;
    defer current_atoms.deinit(allocator);
    var current_bonds: std.ArrayList(CompBond) = .empty;
    defer current_bonds.deinit(allocator);
    var pending: ?cif.tokenizer.Token = null;

    while (true) {
        const token = pending orelse tok.next();
        pending = null;

        switch (token.type) {
            .eof => {
                try flushComponent(allocator, &dict.components, current_comp_id, current_comp_type, &current_atoms, &current_bonds);
                break;
            },
            .data => {
                try flushComponent(allocator, &dict.components, current_comp_id, current_comp_type, &current_atoms, &current_bonds);
                const text = token.text(source);
                current_comp_id = if (text.len > 5) text[5..] else "";
                current_comp_type = "";
                current_atoms.clearRetainingCapacity();
                current_bonds.clearRetainingCapacity();
            },
            .tag => {
                const tag = token.text(source);
                const val_token = tok.next();
                const val = if (val_token.type == .value) val_token.text(source) else "";
                if (std.ascii.eqlIgnoreCase(tag, "_chem_comp.type")) {
                    current_comp_type = val;
                }
            },
            .loop => {
                // Read loop tags, identify atom or bond loop, parse values
                // (Same pattern as zig-cif-graph-parser/component_dict.zig)
                pending = try parseLoop(allocator, source, &tok, &current_atoms, &current_bonds);
            },
            else => {},
        }
    }

    return dict;
}

fn flushComponent(
    allocator: std.mem.Allocator,
    components: *std.StringHashMap(Component),
    comp_id: ?[]const u8,
    comp_type: []const u8,
    atoms: *std.ArrayList(CompAtom),
    bonds: *std.ArrayList(CompBond),
) !void {
    if (comp_id) |id| {
        if (atoms.items.len > 0) {
            const a = try allocator.dupe(CompAtom, atoms.items);
            const b = try allocator.dupe(CompBond, bonds.items);
            try components.put(id, .{
                .comp_id = id,
                .comp_type = comp_type,
                .atoms = a,
                .bonds = b,
            });
        }
    }
}

fn parseLoop(
    allocator: std.mem.Allocator,
    source: []const u8,
    tok: *cif.tokenizer.Tokenizer,
    atoms: *std.ArrayList(CompAtom),
    bonds: *std.ArrayList(CompBond),
) !?cif.tokenizer.Token {
    // Read loop tags
    var tags: [32][]const u8 = undefined;
    var tag_count: usize = 0;

    while (true) {
        const t = tok.next();
        if (t.type == .tag) {
            if (tag_count < 32) {
                tags[tag_count] = t.text(source);
                tag_count += 1;
            }
        } else {
            // Identify loop type and parse values
            var col_atom_id: ?usize = null;
            var col_type_symbol: ?usize = null;
            var col_charge: ?usize = null;
            var col_leaving: ?usize = null;
            var col_atom1: ?usize = null;
            var col_atom2: ?usize = null;
            var col_order: ?usize = null;
            var col_arom: ?usize = null;
            var col_ideal_x: ?usize = null;
            var col_ideal_y: ?usize = null;
            var col_ideal_z: ?usize = null;

            for (tags[0..tag_count], 0..) |tg, i| {
                if (std.ascii.eqlIgnoreCase(tg, "_chem_comp_atom.atom_id")) col_atom_id = i;
                if (std.ascii.eqlIgnoreCase(tg, "_chem_comp_atom.type_symbol")) col_type_symbol = i;
                if (std.ascii.eqlIgnoreCase(tg, "_chem_comp_atom.charge")) col_charge = i;
                if (std.ascii.eqlIgnoreCase(tg, "_chem_comp_atom.pdbx_leaving_atom_flag")) col_leaving = i;
                if (std.ascii.eqlIgnoreCase(tg, "_chem_comp_atom.pdbx_model_Cartn_x_ideal")) col_ideal_x = i;
                if (std.ascii.eqlIgnoreCase(tg, "_chem_comp_atom.pdbx_model_Cartn_y_ideal")) col_ideal_y = i;
                if (std.ascii.eqlIgnoreCase(tg, "_chem_comp_atom.pdbx_model_Cartn_z_ideal")) col_ideal_z = i;
                if (std.ascii.eqlIgnoreCase(tg, "_chem_comp_bond.atom_id_1")) col_atom1 = i;
                if (std.ascii.eqlIgnoreCase(tg, "_chem_comp_bond.atom_id_2")) col_atom2 = i;
                if (std.ascii.eqlIgnoreCase(tg, "_chem_comp_bond.value_order")) col_order = i;
                if (std.ascii.eqlIgnoreCase(tg, "_chem_comp_bond.pdbx_aromatic_flag")) col_arom = i;
            }

            const is_atom_loop = col_atom_id != null and col_type_symbol != null;
            const is_bond_loop = col_atom1 != null and col_atom2 != null;

            // First value token already consumed
            if (t.type != .value) return t;

            var col: usize = 0;
            var current_atom = CompAtom{};
            var current_bond = CompBond{ .atom_idx_1 = 0, .atom_idx_2 = 0, .order = .unknown };
            var atom_name_for_bond: ?[]const u8 = null;
            _ = atom_name_for_bond;

            // Process first value
            if (is_atom_loop) processAtomValue(&current_atom, col, t.text(source), col_atom_id, col_type_symbol, col_charge, col_leaving, col_ideal_x, col_ideal_y, col_ideal_z);
            col += 1;

            while (true) {
                const vt = tok.next();
                if (vt.type != .value) {
                    // Flush last row
                    if (col > 0 and col >= tag_count) {
                        if (is_atom_loop) try atoms.append(allocator, current_atom);
                        if (is_bond_loop) try bonds.append(allocator, current_bond);
                    }
                    return vt;
                }

                if (is_atom_loop) processAtomValue(&current_atom, col, vt.text(source), col_atom_id, col_type_symbol, col_charge, col_leaving, col_ideal_x, col_ideal_y, col_ideal_z);

                col += 1;
                if (col >= tag_count) {
                    if (is_atom_loop) try atoms.append(allocator, current_atom);
                    if (is_bond_loop) try bonds.append(allocator, current_bond);
                    col = 0;
                    current_atom = CompAtom{};
                    current_bond = CompBond{ .atom_idx_1 = 0, .atom_idx_2 = 0, .order = .unknown };
                }
            }
        }
    }
}

fn processAtomValue(
    atom: *CompAtom,
    col: usize,
    val: []const u8,
    col_atom_id: ?usize,
    col_type_symbol: ?usize,
    col_charge: ?usize,
    col_leaving: ?usize,
    col_ideal_x: ?usize,
    col_ideal_y: ?usize,
    col_ideal_z: ?usize,
) void {
    const s = cif.asString(val);
    if (col_atom_id != null and col == col_atom_id.?) {
        const len: u4 = @intCast(@min(s.len, 4));
        for (0..len) |i| atom.name[i] = s[i];
        atom.name_len = len;
    }
    if (col_type_symbol != null and col == col_type_symbol.?) {
        if (s.len >= 1) atom.element_symbol[0] = s[0];
        if (s.len >= 2) atom.element_symbol[1] = s[1];
    }
    if (col_charge != null and col == col_charge.?) {
        atom.charge = cif.asIntOr(i8, val, 0);
    }
    if (col_leaving != null and col == col_leaving.?) {
        atom.leaving = s.len > 0 and (s[0] == 'Y' or s[0] == 'y');
    }
    if (col_ideal_x != null and col == col_ideal_x.?) atom.ideal_x = cif.asFloatOr(val, 0.0);
    if (col_ideal_y != null and col == col_ideal_y.?) atom.ideal_y = cif.asFloatOr(val, 0.0);
    if (col_ideal_z != null and col == col_ideal_z.?) atom.ideal_z = cif.asFloatOr(val, 0.0);
}
```

- [ ] **Step 3: Add test with inline CCD fragment**

```zig
test "parse CCD fragment" {
    const source =
        \\data_ALA
        \\_chem_comp.type 'L-peptide linking'
        \\loop_
        \\_chem_comp_atom.atom_id
        \\_chem_comp_atom.type_symbol
        \\_chem_comp_atom.charge
        \\_chem_comp_atom.pdbx_leaving_atom_flag
        \\N   N 0 N
        \\CA  C 0 N
        \\C   C 0 N
        \\O   O 0 N
        \\CB  C 0 N
        \\H   H 0 N
        \\HA  H 0 N
        \\HB1 H 0 N
        \\HB2 H 0 N
        \\HB3 H 0 N
        \\OXT O 0 Y
    ;
    var dict = try parseComponentDict(std.testing.allocator, source);
    defer dict.deinit();

    const ala = dict.get("ALA");
    try std.testing.expect(ala != null);
    try std.testing.expectEqual(@as(usize, 11), ala.?.atoms.len);
}
```

- [ ] **Step 4: Run tests, register, commit**

Run: `cd ~/zreduce && zig build test`
Expected: All tests pass.

Add to root.zig: `pub const ccd = @import("ccd.zig");`

```bash
cd ~/zreduce && git add src/ccd.zig src/root.zig && git commit -m "feat: add streaming CCD component dictionary parser"
```

---

### Task 11: Unified Placer + HET Placement

**Files:**
- Create: `~/zreduce/src/place/het.zig`
- Create: `~/zreduce/src/place/placer.zig`
- Create: `~/zreduce/src/place.zig`

- [ ] **Step 1: Implement het.zig (CCD-derived H placement)**

Derives placement plans from CCD bond topology by determining hybridization from bond orders and selecting appropriate geometry types.

- [ ] **Step 2: Implement placer.zig (unified entry point)**

```zig
/// Add hydrogens to the model.
/// Uses standard plans for known residues, CCD-derived plans for HET groups.
pub fn addHydrogens(
    model: *Model,
    ccd_dict: ?*const ComponentDict,
) !PlacementResult;
```

- [ ] **Step 3: Add integration test with tiny.cif**

Verify that HA is placed on ALA at correct position and bond length.

- [ ] **Step 4: Run tests, commit**

```bash
cd ~/zreduce && git add src/place/ src/place.zig src/root.zig && git commit -m "feat: add unified H placer with standard + CCD support"
```

---

## Phase 3: Optimization

### Task 12: Dot Sphere Generation

**Files:**
- Create: `~/zreduce/src/optimize/dot_sphere.zig`

- [ ] **Step 1: Write test for dot count at standard density**

```zig
test "dot sphere count at density 16" {
    const sphere = try DotSphere.generate(testing.allocator, 1.40, 16.0);
    defer sphere.deinit();
    // Expected: ~4π * 16 * 1.4² ≈ 394 dots
    try testing.expect(sphere.points.len > 350);
    try testing.expect(sphere.points.len < 450);
}
```

- [ ] **Step 2: Implement dot sphere generation**

Port the original reduce algorithm: concentric rings with alternating 5° offset, matching DotSph.cpp.

- [ ] **Step 3: Run tests, commit**

```bash
cd ~/zreduce && git add src/optimize/dot_sphere.zig && git commit -m "feat: add dot sphere generation matching original reduce"
```

---

### Task 13: Scorer (bump/H-bond)

**Files:**
- Create: `~/zreduce/src/optimize/scorer.zig`

- [ ] **Step 1: Define ScoringParams and ScoreResult**

All constants from original reduce (see design spec).

- [ ] **Step 2: Implement atomScore function**

Port the exact scoring logic from AtomPositions.cpp:
- Contact: `exp(-(gap/0.25)²)`
- Clash: `-10.0 * (-0.5 * gap)`
- H-bond: `+4.0 * (-0.5 * gap)`

- [ ] **Step 3: Test with known good/bad contacts**

- [ ] **Step 4: Commit**

```bash
cd ~/zreduce && git add src/optimize/scorer.zig && git commit -m "feat: add dot-sphere bump/H-bond scorer"
```

---

### Task 14: Mover Types (Rotators)

**Files:**
- Create: `~/zreduce/src/optimize/mover.zig`
- Create: `~/zreduce/src/optimize/rotator.zig`

- [ ] **Step 1: Define Mover interface**

```zig
pub const Mover = struct {
    kind: MoverKind,
    residue_idx: u32,
    atom_indices: []u32,
    n_orientations: u16,
    best_orientation: u16,
    current_orientation: u16,
    penalty: f32,

    pub fn applyOrientation(self: *Mover, model: *Model, idx: u16) void;
    pub fn orientationPenalty(self: *const Mover, idx: u16) f32;
};
```

- [ ] **Step 2: Implement SingleHRotator (OH, SH)**

12 orientations at 30° intervals around bond axis.

- [ ] **Step 3: Implement NH3Rotator**

3 orientations maintaining 120° H spacing.

- [ ] **Step 4: Test rotator with known geometry, commit**

```bash
cd ~/zreduce && git add src/optimize/mover.zig src/optimize/rotator.zig && git commit -m "feat: add OH/SH/NH3+ rotation movers"
```

---

### Task 15: Mover Types (Flippers)

**Files:**
- Create: `~/zreduce/src/optimize/flipper.zig`

- [ ] **Step 1: Implement AmideFlipper (Asn/Gln)**

2 orientations: original + O↔N swap. Penalties: 0.00, 0.50.

- [ ] **Step 2: Implement HisFlipper**

6 orientations (2 flip states × 3 protonation states). Penalties from design spec.

- [ ] **Step 3: Test flip coordinate swaps, commit**

```bash
cd ~/zreduce && git add src/optimize/flipper.zig && git commit -m "feat: add Asn/Gln/His flip movers"
```

---

### Task 16: Clique Detection

**Files:**
- Create: `~/zreduce/src/optimize/clique.zig`

- [ ] **Step 1: Build interaction graph**

Edge between two movers if any of their atoms are within scoring distance.

- [ ] **Step 2: Find connected components (cliques)**

Simple BFS/DFS on the interaction graph.

- [ ] **Step 3: Test with known mover layout, commit**

```bash
cd ~/zreduce && git add src/optimize/clique.zig && git commit -m "feat: add interaction graph and clique detection"
```

---

### Task 17: Optimizer (brute-force + vertex-cut)

**Files:**
- Create: `~/zreduce/src/optimize/optimizer.zig`
- Create: `~/zreduce/src/optimize.zig`

- [ ] **Step 1: Implement singleton optimization**

Score all orientations, pick best.

- [ ] **Step 2: Implement brute-force clique optimization**

Enumerate all combinations up to brute_force_limit (100K).

- [ ] **Step 3: Implement vertex-cut decomposition**

Port from cctbx reduce2's approach: find smallest vertex cut, recursively decompose.

- [ ] **Step 4: Create optimize.zig re-export, test with small clique, commit**

```bash
cd ~/zreduce && git add src/optimize/ src/optimize.zig src/root.zig && git commit -m "feat: add clique optimizer with brute-force and vertex-cut"
```

---

## Phase 4: Output + CLI

### Task 18: mmCIF Writer

**Files:**
- Create: `~/zreduce/src/writer/mmcif_writer.zig`

- [ ] **Step 1: Write atom_site loop with added H atoms**

Preserve original atoms, append H atoms, update atom IDs.

- [ ] **Step 2: Add _zreduce_log custom category**

- [ ] **Step 3: Test round-trip (read → add H → write → read), commit**

```bash
cd ~/zreduce && git add src/writer/mmcif_writer.zig && git commit -m "feat: add mmCIF writer with H atoms and zreduce log"
```

---

### Task 19: JSON Log Writer

**Files:**
- Create: `~/zreduce/src/writer/json_writer.zig`
- Create: `~/zreduce/src/writer.zig`

- [ ] **Step 1: Implement JSON output**

```json
{
  "version": "0.1.0",
  "input": "structure.cif",
  "hydrogens_added": 1234,
  "movers": [
    {"residue": "A.HIS.42", "type": "his_flip", "orientation": 2, "score": -3.45}
  ]
}
```

- [ ] **Step 2: Create writer.zig re-export, commit**

```bash
cd ~/zreduce && git add src/writer/ src/writer.zig src/root.zig && git commit -m "feat: add JSON log writer"
```

---

### Task 20: Full CLI Integration

**Files:**
- Modify: `~/zreduce/src/main.zig`

- [ ] **Step 1: Wire up complete pipeline in main.zig**

```zig
pub fn main() !void {
    // Parse args
    // Read mmCIF
    // Load CCD (if provided)
    // Build model
    // Place hydrogens
    // Build movers
    // Optimize (unless --no-opt)
    // Write output
}
```

- [ ] **Step 2: Add CLI flags**

`--dict`, `--output`, `--no-opt`, `--no-flip`, `--density`, `--limit`

- [ ] **Step 3: End-to-end test with real mmCIF, commit**

```bash
cd ~/zreduce && git add src/main.zig && git commit -m "feat: wire up full CLI pipeline"
```

---

### Task 21: Integration Tests + Regression

**Files:**
- Create: `~/zreduce/tests/integration_test.zig`

- [ ] **Step 1: Test against original reduce output**

Download a small PDB structure, run both original reduce and zreduce, compare H atom counts and positions (within tolerance).

- [ ] **Step 2: Test flip decisions**

Verify Asn/Gln/His flip decisions match original reduce on test structures.

- [ ] **Step 3: Commit**

```bash
cd ~/zreduce && git add tests/ && git commit -m "test: add integration tests comparing with original reduce"
```

---

## Dependencies Between Tasks

```
Task 1 (scaffold)
  ├── Task 2 (math)
  ├── Task 3 (element)
  └── Task 4 (tokenizer)
        └── Task 5 (parser)
              ├── Task 7 (mmcif → model)  ← Task 6 (model structs)
              └── Task 10 (CCD loader)
                    └── Task 11 (placer) ← Task 8 (geometry) + Task 9 (standard plans)
                          └── Task 12 (dot sphere)
                                └── Task 13 (scorer)
                                      └── Task 14 (rotators) + Task 15 (flippers)
                                            └── Task 16 (clique)
                                                  └── Task 17 (optimizer)
                                                        ├── Task 18 (mmcif writer)
                                                        ├── Task 19 (json writer)
                                                        └── Task 20 (CLI)
                                                              └── Task 21 (integration)
```

**Parallelizable tasks:**
- Tasks 2, 3, 4 can all be done in parallel
- Tasks 8, 9, 10 can all be done in parallel (after Task 5)
- Tasks 14, 15 can be done in parallel
- Tasks 18, 19 can be done in parallel
