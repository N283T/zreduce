# zreduce Design Specification

**Date:** 2026-03-26
**Status:** Draft
**Author:** nagaet + Claude

## Overview

zreduce is a from-scratch Zig implementation of the reduce hydrogen placement tool.
It reads mmCIF structures, adds hydrogen atoms using geometric rules and CCD bond
topology, then optimizes rotatable/flippable groups via clique-based search with
probe dot-sphere scoring.

### Goals

- Full feature parity with original C++ reduce (geometry, rotation, flip, scoring)
- mmCIF-native (input and output)
- CCD-driven hydrogen placement for HET groups
- High performance via specialized parsers, SIMD, and Zig comptime
- Clean, modular architecture informed by zsasa/zdssp/zig-cif-graph-parser patterns

### Non-Goals (for now)

- PDB format input/output
- Python bindings
- GUI

---

## Architecture

### Pipeline

```
mmCIF input
    │
    ▼
┌─────────────┐
│  CIF Parser  │  Specialized tokenizer + parser (atom_site, struct_conn)
└──────┬──────┘
       │
       ▼
┌─────────────┐
│  CCD Loader  │  Parse components.cif → ComponentDict (HET groups)
└──────┬──────┘
       │
       ▼
┌─────────────┐
│ Model Build  │  Atom[], Residue[], Chain[], Bond[] + neighbor list
└──────┬──────┘
       │
       ▼
┌─────────────────┐
│  H Placement     │  Geometric placement (type 1-6) using bond topology
│  (Standard+HET)  │  Standard: hardcoded tables, HET: CCD-derived
└──────┬──────────┘
       │
       ▼
┌─────────────────┐
│  Mover Registry  │  Identify rotatable/flippable groups → Mover objects
└──────┬──────────┘
       │
       ▼
┌─────────────────┐
│  Clique Builder  │  Build interaction graph → detect cliques
└──────┬──────────┘
       │
       ▼
┌─────────────────┐
│  Optimizer       │  Dot-sphere scoring + clique search (brute force + vertex-cut)
└──────┬──────────┘
       │
       ▼
┌─────────────────┐
│  Output Writer   │  mmCIF with H atoms + flip corrections
│                  │  Log: custom _zreduce category + JSON
└─────────────────┘
```

### Source Layout

```
~/zreduce/
├── build.zig
├── build.zig.zon
├── src/
│   ├── main.zig              # CLI entry point
│   ├── root.zig              # Library re-exports
│   │
│   ├── cif/                  # CIF parser subsystem (new, specialized)
│   │   ├── tokenizer.zig     # CIF token stream
│   │   ├── parser.zig        # Tokens → Document/Block/Loop
│   │   ├── types.zig         # Document, Block, Loop, Item
│   │   ├── value.zig         # Value extraction helpers
│   │   └── char_table.zig    # Character classification LUT
│   │
│   ├── mmcif.zig             # mmCIF atom_site extraction
│   ├── ccd.zig               # CCD components.cif parser → ComponentDict
│   │
│   ├── model/                # Molecular model
│   │   ├── atom.zig          # Atom struct
│   │   ├── residue.zig       # Residue struct (atom range, comp_id)
│   │   ├── chain.zig         # Chain struct (residue range)
│   │   ├── bond.zig          # Bond struct (order, source)
│   │   ├── model.zig         # Model aggregate (atoms, residues, chains, bonds)
│   │   └── neighbor.zig      # Spatial neighbor list (cell list)
│   │
│   ├── element.zig           # Element properties (VDW radii, comptime tables)
│   │
│   ├── place/                # Hydrogen placement
│   │   ├── geometry.zig      # Type 1-6 placement functions
│   │   ├── standard.zig      # Hardcoded plans for 20 AA + nucleic acids
│   │   ├── het.zig           # CCD-derived placement for HET groups
│   │   └── placer.zig        # Unified placement entry point
│   │
│   ├── optimize/             # Optimization engine
│   │   ├── mover.zig         # Mover interface + concrete types
│   │   ├── rotator.zig       # OH/SH/NH3+/methyl rotation movers
│   │   ├── flipper.zig       # Asn/Gln/His flip movers
│   │   ├── dot_sphere.zig    # Dot sphere generation (comptime + runtime)
│   │   ├── scorer.zig        # Bump/H-bond scoring
│   │   ├── clique.zig        # Clique detection (connected components)
│   │   └── optimizer.zig     # Brute-force + vertex-cut search
│   │
│   ├── writer/               # Output
│   │   ├── mmcif_writer.zig  # mmCIF output (atom_site + custom categories)
│   │   └── json_writer.zig   # JSON log output
│   │
│   └── math.zig              # Vec3, rotation, dihedral, cross product
│
├── test_data/                # Test mmCIF files
├── docs/
│   └── specs/                # This document
└── tests/                    # Integration tests
```

---

## Module Design

### 1. CIF Parser (`src/cif/`)

**Purpose:** Parse CIF/mmCIF files into a structured Document tree.

**Design:** New implementation referencing zig-cif-graph-parser's gemmi-port approach.
Specialized for reduce's needs but general enough to handle any CIF category.

```zig
// Core types
const Document = struct {
    blocks: []Block,
};

const Block = struct {
    name: []const u8,
    items: []Item,

    pub fn findLoop(self: *const Block, tag: []const u8) ?*const Loop;
    pub fn findValue(self: *const Block, tag: []const u8) ?[]const u8;
};

const Loop = struct {
    tags: [][]const u8,
    values: [][]const u8,  // row-major

    pub fn width(self: *const Loop) usize;
    pub fn length(self: *const Loop) usize;
    pub fn val(self: *const Loop, row: usize, col: usize) ?[]const u8;
    pub fn findTag(self: *const Loop, tag: []const u8) ?usize;
};
```

**Key optimizations:**
- Zero-copy: values are slices into mmap'd source buffer
- Case-insensitive tag lookup
- Gzip-aware file reading (via zlib linkage)

### 2. mmCIF Extraction (`src/mmcif.zig`)

**Purpose:** Extract atom coordinates and connectivity from mmCIF.

**Categories parsed:**
- `_atom_site` → Atom coordinates, element, residue identity, altloc
- `_struct_conn` → Disulfide bonds, metal coordination, etc.
- `_entity` → Entity type (polymer, non-polymer, water)

**Output:** Populated `Model` struct.

### 3. CCD Loader (`src/ccd.zig`)

**Purpose:** Parse components.cif (or individual CCD entries) for HET group topology.

**Categories parsed:**
- `_chem_comp` → Component metadata (type, name)
- `_chem_comp_atom` → Atom names, element, charge, leaving_atom_flag, aromatic_flag
- `_chem_comp_bond` → Bond pairs, bond order (SING/DOUB/TRIP/AROM), stereo

**Output:**
```zig
const ComponentDict = std.StringHashMap(Component);

const Component = struct {
    comp_id: [3]u8,
    comp_id_len: u3,
    type: CompType,            // peptide, RNA, DNA, non-polymer, etc.
    atoms: []CompAtom,
    bonds: []CompBond,
};

const CompAtom = struct {
    name: [4]u8,
    name_len: u4,
    element: Element,
    charge: i8,
    leaving: bool,             // leaving_atom_flag
    aromatic: bool,
    // Ideal coordinates (for validation/fallback)
    ideal_x: f32,
    ideal_y: f32,
    ideal_z: f32,
};

const CompBond = struct {
    atom_1: u16,               // index into atoms
    atom_2: u16,
    order: BondOrder,          // sing, doub, trip, arom
    stereo: StereoConfig,
};
```

**Hybridization derivation from CCD:**
```
sp3: atom has 4 bonds, all SING
sp2: atom has 3 bonds with at least one DOUB or AROM
sp:  atom has 2 bonds with a TRIP, or 2 DOUB
```

### 4. Model (`src/model/`)

**Purpose:** In-memory molecular model with flat array storage.

```zig
const Atom = struct {
    pos: Vec3(f32),
    name: [4]u8,
    name_len: u4,
    element: Element,
    residue_idx: u32,
    altloc: u8,
    occupancy: f32,
    b_factor: f32,
    is_hydrogen: bool,
    // Scoring properties (set during placement)
    vdw_radius: f32,
    flags: AtomFlags,          // donor, acceptor, aromatic, charged, etc.
};

const Residue = struct {
    comp_id: [3]u8,
    comp_id_len: u3,
    chain_idx: u16,
    seq_id: i32,
    ins_code: u8,
    atom_range: Range(u32),    // [start, end) into atom array
    entity_type: EntityType,
};

const Model = struct {
    atoms: []Atom,
    residues: []Residue,
    chains: []Chain,
    bonds: []Bond,
    neighbor_list: CellList,

    pub fn atomsOfResidue(self: *const Model, res: *const Residue) []Atom;
    pub fn neighborsOf(self: *const Model, atom_idx: u32, radius: f32) NeighborIterator;
};
```

### 5. Element Table (`src/element.zig`)

**Purpose:** Element properties for scoring and placement.

**Design:** Comptime-generated lookup tables matching original reduce's values.

```zig
// Comptime-generated from original reduce ElementInfo table
const AtomType = enum {
    H,          // 1.22 Å (non-polar H)
    Hpol,       // 1.05 Å (polar H, donor)
    Har,        // 1.05 Å (aromatic H)
    C,          // 1.70 Å
    Car,        // 1.75 Å (aromatic C, acceptor)
    C_eq_O,     // 1.65 Å (carbonyl C)
    N,          // 1.55 Å
    Nacc,       // 1.55 Å (N acceptor)
    O,          // 1.40 Å (acceptor)
    S,          // 1.80 Å (acceptor)
    P,          // 1.80 Å
    Se,         // 1.90 Å
    // ... metals, halogens
};

const AtomTypeInfo = struct {
    explicit_radius: f32,      // VDW radius with explicit H
    implicit_radius: f32,      // VDW radius with implicit H
    covalent_radius: f32,
    flags: AtomFlags,          // DONOR, ACCEPTOR, AROMATIC, CHARGED, etc.
};

// Comptime lookup
const atom_type_table: [AtomType.count]AtomTypeInfo = comptime blk: {
    // ... populate from original reduce constants
};
```

### 6. Hydrogen Placement (`src/place/`)

#### geometry.zig — Type 1-6 Placement Functions

Direct port of AtomConn.cpp's vector math:

```zig
/// Type 1 (HXR3): Tetrahedral — H opposite to 3 neighbors
pub fn placeHXR3(center: Vec3, n1: Vec3, n2: Vec3, n3: Vec3, bond_len: f32) Vec3 {
    const v1 = n1.sub(center).normalize();
    const v2 = n2.sub(center).normalize();
    const v3 = n3.sub(center).normalize();
    const dir = v1.add(v2).add(v3).scaleTo(-bond_len);
    return center.add(dir);
}

/// Type 2 (H2XR2): Two H on sp2
/// Type 3 (H3XR): Dihedral-controlled placement
/// Type 4 (HXR2): Planar bisector with fudge
/// Type 5 (HXR2): Fractional angle placement
/// Type 6 (HXY): Linear extension
```

#### standard.zig — Hardcoded Standard Residue Plans

```zig
const PlacementPlan = struct {
    h_name: [4]u8,
    placement_type: PlacementType,  // hxr3, h2xr2, h3xr, hxr2_planar, hxr2_frac, hxy
    connected: [3][4]u8,            // names of reference atoms
    bond_len: f32,
    angle: f32,
    dihedral: f32,
    atom_type: AtomType,            // for VDW radius and flags
    mover_hint: MoverHint,          // none, rotate, rotate_nh3, flip
};

// Example: ALA
const ala_plans = [_]PlacementPlan{
    .{ .h_name = " HA ", .placement_type = .hxr3,
       .connected = .{" N  ", " C  ", " CB "},
       .bond_len = 1.10, .angle = 0, .dihedral = 0,
       .atom_type = .H, .mover_hint = .none },
    // ... HB1, HB2, HB3 with rotate hint
};

// Lookup table
const standard_plans = std.StaticStringMap([]const PlacementPlan).initComptime(.{
    .{ "ALA", &ala_plans },
    .{ "GLY", &gly_plans },
    // ... 20 amino acids + nucleic acids
});
```

#### het.zig — CCD-Derived Placement

```zig
/// Derive placement plans from CCD component definition at runtime
pub fn derivePlans(
    allocator: Allocator,
    component: *const Component,
    existing_atoms: []const Atom,
) ![]PlacementPlan {
    // For each atom in component:
    //   1. Count bonds and determine hybridization (sp3/sp2/sp)
    //   2. Check if hydrogen (element == H)
    //   3. If H and not present in existing_atoms:
    //      - Find bonded heavy atom
    //      - Find heavy atom's other neighbors
    //      - Select placement type based on hybridization
    //      - Compute bond length from element types
    //      - Create PlacementPlan
}
```

### 7. Optimization Engine (`src/optimize/`)

#### dot_sphere.zig — Dot Sphere Generation

```zig
const DotSphere = struct {
    points: []Vec3(f32),
    density: f32,
    radius: f32,
};

/// Generate dot sphere matching original reduce algorithm
/// density: 16.0 dots/Å² (default)
/// Uses concentric ring placement with 5° alternating offset
pub fn generate(allocator: Allocator, radius: f32, density: f32) !DotSphere;

/// Comptime-generated unit sphere at standard density
/// Runtime: scale to desired radius
const unit_sphere: [N]Vec3(f32) = comptime generateUnitSphere(16.0);
```

#### scorer.zig — Dot-Sphere Scoring

Exact port of original reduce scoring constants and logic:

```zig
const ScoringParams = struct {
    gap_scale: f32 = 0.25,
    bump_weight: f32 = 10.0,
    hb_weight: f32 = 4.0,
    min_reg_hb_gap: f32 = 0.6,
    min_charged_hb_gap: f32 = 0.8,
    bad_bump_gap_cut: f32 = 0.4,
    dot_density: f32 = 16.0,
    probe_radius: f32 = 0.0,
};

const ScoreResult = struct {
    total: f32,
    bump_sub: f32,
    hb_sub: f32,
    has_bad_bump: bool,
};

/// Score a single atom against its environment using dot sphere probes
pub fn scoreAtom(
    atom: *const Atom,
    dot_sphere: *const DotSphere,
    neighbors: []const Atom,
    params: ScoringParams,
) ScoreResult {
    // For each dot on sphere (scaled to atom VDW + probe radius):
    //   1. Find minimum gap to any neighbor atom
    //   2. Classify: contact (gap > 0), clash (gap < 0), or H-bond
    //   3. Score:
    //      - Contact: exp(-(gap/0.25)²)
    //      - Clash: -10.0 * (-0.5 * gap)
    //      - H-bond: +4.0 * (-0.5 * gap)  (if valid H-bond geometry)
    //   4. Bad bump if gap <= -0.4
    // Normalize total by dot density
}
```

#### mover.zig — Mover Interface

```zig
const Mover = struct {
    kind: MoverKind,
    residue_idx: u32,
    atom_indices: []u32,       // atoms that move with this mover
    orientations: []Orientation,
    best_orientation: u16,
    current_orientation: u16,
    penalty: f32,              // bias toward original orientation

    const MoverKind = enum {
        single_h_rotator,      // OH, SH, SeH
        nh3_rotator,           // NH3+
        methyl_rotator,        // CH3
        aromatic_methyl,       // aromatic CH3
        amide_flip,            // Asn, Gln
        his_flip,              // His (6 orientations)
    };

    const Orientation = struct {
        positions: []Vec3(f32), // atom positions for this orientation
        penalty: f32,           // orientation-specific penalty
    };

    pub fn applyOrientation(self: *Mover, model: *Model, idx: u16) void;
    pub fn scoreOrientation(self: *const Mover, model: *const Model, idx: u16, scorer: *const Scorer) f32;
};
```

#### clique.zig — Clique Detection

```zig
/// Build interaction graph: edge between movers whose atoms can interact
pub fn buildInteractionGraph(
    movers: []const Mover,
    model: *const Model,
    cutoff: f32,
) InteractionGraph;

/// Find connected components (independent cliques)
pub fn findCliques(graph: *const InteractionGraph) [][]u32;
```

#### optimizer.zig — Search Algorithm

```zig
/// Optimize all movers
pub fn optimize(
    movers: []Mover,
    model: *Model,
    scorer: *const Scorer,
    config: OptConfig,
) void {
    const graph = clique.buildInteractionGraph(movers, model, config.cutoff);
    const cliques = clique.findCliques(&graph);

    for (cliques) |clq| {
        if (clq.len == 1) {
            optimizeSingleton(movers[clq[0]], model, scorer);
        } else if (totalStates(movers, clq) <= config.brute_force_limit) {
            optimizeBruteForce(movers, clq, model, scorer);
        } else {
            optimizeVertexCut(movers, clq, model, scorer, config);
        }
    }
}

const OptConfig = struct {
    brute_force_limit: u64 = 100_000,
    vertex_cut_max_depth: u32 = 3,
    penalty_magnitude: f32 = 0.01,
};
```

### 8. Flip Logic (`src/optimize/flipper.zig`)

Hardcoded flip tables for Asn/Gln/His:

```zig
const HisFlipState = enum(u3) {
    no_hd1 = 0,               // ND1 acceptor, NE2-HE2 donor
    no_he2 = 1,               // ND1-HD1 donor, NE2 acceptor
    both_h = 2,               // Both protonated (+1 charge)
    flip_no_hd1 = 3,          // Flipped ring, ND1 acceptor
    flip_no_he2 = 4,          // Flipped ring, NE2 acceptor
    flip_both_h = 5,          // Flipped ring, both protonated
};

const his_penalties = [6]f32{ 0.00, 0.00, 0.05, 0.50, 0.50, 0.55 };

const AsqFlipState = enum(u1) {
    original = 0,
    flipped = 1,               // O↔N swap
};

const asn_penalties = [2]f32{ 0.00, 0.50 };

/// Apply His flip: swap atom positions and update donor/acceptor flags
pub fn applyHisFlip(model: *Model, residue: *const Residue, state: HisFlipState) void {
    // Swap pairs: ND1↔CD2, CE1↔NE2 (for flipped states 3-5)
    // Recompute H positions using type 4 placement
    // Update atom flags (donor/acceptor)
}

/// Apply Asn/Gln flip: swap O↔N positions
pub fn applyAmideFlip(model: *Model, residue: *const Residue, state: AsqFlipState) void {
    // Swap OD1↔ND2 (Asn) or OE1↔NE2 (Gln)
    // Recompute H positions using type 3 placement
    // Update atom flags
}
```

### 9. Output (`src/writer/`)

#### mmcif_writer.zig

```zig
/// Write mmCIF with:
/// - Original atoms (with flip corrections applied)
/// - Added hydrogen atoms in _atom_site
/// - Custom _zreduce_log category with optimization details
pub fn write(
    writer: anytype,
    model: *const Model,
    movers: []const Mover,
    original_doc: *const cif.Document,
) !void;
```

Custom category:
```
loop_
_zreduce_log.residue_id
_zreduce_log.action           # add_h, rotate, flip
_zreduce_log.description
_zreduce_log.orientation
_zreduce_log.score_before
_zreduce_log.score_after
```

#### json_writer.zig

```zig
/// Write JSON log with detailed optimization results
pub fn writeLog(writer: anytype, movers: []const Mover) !void;
```

### 10. Math Utilities (`src/math.zig`)

```zig
pub fn Vec3(comptime T: type) type {
    return struct {
        x: T, y: T, z: T,

        pub fn add, sub, scale, scaleTo, dot, cross,
               length, distance, normalize,
               rotate, dihedral, angle: ...;
    };
}

/// Rotate point around axis by angle (degrees)
pub fn rotateAroundAxis(point: Vec3, origin: Vec3, axis: Vec3, degrees: f32) Vec3;

/// Compute dihedral angle between 4 points
pub fn dihedralAngle(a: Vec3, b: Vec3, c: Vec3, d: Vec3) f32;
```

---

## Implementation Phases

### Phase 1: Foundation
1. CIF parser (tokenizer → parser → Document)
2. Math library (Vec3, rotation, dihedral)
3. Element table (comptime VDW radii)
4. mmCIF extraction (atom_site → Model)

### Phase 2: Hydrogen Placement
5. Type 1-6 geometry functions
6. Standard residue placement plans (hardcoded 20 AA + nucleic acids)
7. CCD loader + HET placement
8. mmCIF writer (input + added H → output)

### Phase 3: Optimization
9. Dot sphere generation
10. Scorer (bump/H-bond scoring)
11. Mover types (rotators)
12. Mover types (flippers: Asn/Gln/His)
13. Neighbor list for spatial queries
14. Clique detection
15. Optimizer (singleton + brute force)
16. Optimizer (vertex-cut for large cliques)

### Phase 4: Polish
17. JSON log output
18. Custom mmCIF log category
19. CLI with full option support
20. Integration tests against original reduce output
21. Performance optimization (SIMD, threading)

---

## Key Design Decisions

### D1: Specialized CIF parser vs reuse

**Decision:** Build from scratch.
**Rationale:** Reduce needs very fast parsing of specific categories. A general-purpose
parser adds overhead. The specialized parser can skip irrelevant categories entirely
and extract only needed columns without building a full document tree for unused data.

### D2: Hardcoded standard residues vs all-CCD

**Decision:** Hardcode standard residues, CCD for HET only.
**Rationale:** Standard residues have well-known, stable H placement rules with
special-case handling (flips, rotations). Deriving these from CCD would lose the
domain-specific optimization hints (mover_hint, orientation penalties). HET groups
are too numerous and varied to hardcode.

### D3: Dot sphere scoring vs distance-based

**Decision:** Dot sphere (original reduce compatible).
**Rationale:** Maintains compatibility with established validation workflows.
The dot sphere approach captures surface geometry effects that simple distance
cutoffs miss.

### D4: Vertex-cut optimization

**Decision:** Implement both brute-force and vertex-cut.
**Rationale:** Most cliques are small (1-4 movers) and brute-force is optimal.
Vertex-cut handles the rare large cliques that original reduce would abandon.

### D5: f32 vs f64

**Decision:** f32 for storage, f64 for calculation where needed.
**Rationale:** Matches zsasa/zdssp convention. Coordinate precision of mmCIF
is ~0.001 Å, well within f32 range. Scoring accumulation may benefit from f64.

---

## Scoring Constants (from original reduce)

| Parameter | Value | Description |
|-----------|-------|-------------|
| `dot_density` | 16.0 dots/Å² | Dot sphere density |
| `probe_radius` | 0.0 Å | Probe sphere radius |
| `gap_scale` | 0.25 | Contact score falloff |
| `bump_weight` | 10.0 | Clash penalty multiplier |
| `hb_weight` | 4.0 | H-bond reward multiplier |
| `min_reg_hb_gap` | 0.6 Å | Regular H-bond gap threshold |
| `min_charged_hb_gap` | 0.8 Å | Charged H-bond gap threshold |
| `bad_bump_gap_cut` | 0.4 Å | Bad bump classification |

## VDW Radii (explicit H model)

| Type | Radius (Å) | Flags |
|------|------------|-------|
| H (non-polar) | 1.22 | — |
| Hpol/Har | 1.05 | donor |
| C | 1.70 | — |
| Car | 1.75 | acceptor |
| C=O | 1.65 | — |
| N | 1.55 | — |
| Nacc | 1.55 | acceptor |
| O | 1.40 | acceptor |
| S | 1.80 | acceptor |
| P | 1.80 | — |

## Flip Orientation Penalties

| Residue | State | Penalty |
|---------|-------|---------|
| His | no HD1 (default) | 0.00 |
| His | no HE2 | 0.00 |
| His | both protonated | 0.05 |
| His | flip no HD1 | 0.50 |
| His | flip no HE2 | 0.50 |
| His | flip both | 0.55 |
| Asn/Gln | original | 0.00 |
| Asn/Gln | flipped | 0.50 |
