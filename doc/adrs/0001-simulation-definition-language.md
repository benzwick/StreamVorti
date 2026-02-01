# ADR-0001: Simulation Definition Language (SDL)

**Status:** Proposed
**Date:** 2026-02-01
**Authors:** Benjamin Zwick

## Context

StreamVorti needs a domain-specific language for defining simulations that is:
- Expressive enough for complex multiphysics problems
- Simple enough for common cases
- Dimension-agnostic (1D, 2D, 3D, and beyond)
- Method-agnostic (FEM, FVM, meshless, particle methods, etc.)
- Compatible with existing mesh formats and workflows

Current commercial tools have various limitations:
- **COMSOL**: GUI-centric, expensive, proprietary format
- **ANSYS APDL**: Cryptic syntax, steep learning curve
- **Fluent**: Limited to CFD, TUI is arcane
- **STAR-CCM+**: Java macros are verbose
- **Abaqus**: Input files are rigid, keyword-based
- **OpenFOAM**: Multiple dictionaries, inconsistent syntax

We want a language that combines the best aspects of these tools while being more intuitive and mathematically precise.

## Decision

We adopt a Lisp-based SDL with the following core concepts:

### 1. Domain (Discrete Representation of Ω)

The `domain` unifies all spatial discretizations:

```lisp
;; Structured mesh on n-dimensional box
(domain (box (0 0) (1 1)) :mesh :n (40 40))
(domain (box (0 0 0) (1 1 1)) :mesh :n (20 20 20))

;; Point cloud (for meshless methods)
(domain (box (0 0) (1 1)) :points :n (40 40))
(domain (box (0 0) (1 1)) :points :h 0.025)

;; Particles (for SPH, DEM, MD)
(domain (box (0 0) (1 1)) :particles :n 10000)

;; Unstructured mesh via external generator
(domain (circle (0 0) 1) :mesh :generator :gmsh :h 0.05)

;; From file
(domain :file "cavity.mesh")
(domain :file "particles.csv" :type :points)
```

**Geometry primitives (dimension-agnostic):**
- `box` - n-dimensional hyperrectangle (interval, rectangle, cuboid, ...)
- `ball` - n-dimensional ball
- `simplex` - n-simplex (segment, triangle, tetrahedron, ...)
- `cylinder`, `cone`, `torus` - common 3D shapes
- CSG operations: `union`, `intersection`, `difference`

**Discretization types:**
- `:mesh` - nodes + connectivity (for FEM, FVM, FD)
- `:points` - nodes only (for meshless: DCPSE, RBF, MLS, RKPM, EFG)
- `:particles` - nodes + properties (for SPH, DEM, MD)

### 2. Boundaries (Named Subsets of ∂Ω)

Boundaries are purely geometric - no physics here:

```lisp
(boundaries
  ;; By coordinate predicate
  (lid    (= y 1))
  (bottom (= y 0))
  (walls  (or (= x 0) (= x 1)))

  ;; By mesh attribute (MFEM style)
  (inlet  (attribute 3))
  (outlet (attribute 4))

  ;; By named group (from mesh file)
  (wing   (group "wing-surface")))
```

**Predicate expressions:**
- Comparisons: `(= x 0)`, `(< y 0.5)`, `(<= z 1)`
- Boolean: `(and ...)`, `(or ...)`, `(not ...)`
- Functions: `(near x 0.5 :tol 1e-6)`, `(between y 0 1)`
- Geometric: `(on-surface "cad-face-3")`, `(distance-from point :< 0.1)`

### 3. Subdomains (Named Interior Regions)

For multi-material, multi-phase, or region-specific settings:

```lisp
(subdomains
  ;; By predicate
  (fluid  (< x 0.5))
  (solid  (>= x 0.5))

  ;; By mesh attribute
  (material-1 (attribute 1))
  (material-2 (attribute 2))

  ;; By geometry
  (inclusion (inside (ball (0.5 0.5) 0.1))))
```

### 4. Physics (Governing Equations + BCs)

Each physics block defines equations and boundary conditions:

```lisp
(physics :fluid
  :type :navier-stokes
  :formulation :velocity-pressure  ; or :stream-vorticity
  :Re 1000

  ;; BCs reference boundaries by name
  (bc lid    :velocity (1 0 0))
  (bc inlet  :velocity :parabolic :max 1.5)
  (bc outlet :pressure 0)
  (bc walls  :no-slip)

  ;; Initial conditions
  (ic :velocity (0 0 0))
  (ic :pressure 0))

(physics :thermal
  :type :heat-conduction
  :conductivity 1.0

  ;; Same boundaries, different conditions
  (bc lid    :temperature 100)
  (bc walls  :convection :h 10 :T-inf 20)
  (bc bottom :flux 0)  ; insulated

  (ic :temperature 20))
```

**Supported physics types:**
- Fluids: `:navier-stokes`, `:stokes`, `:euler`, `:potential-flow`
- Heat: `:heat-conduction`, `:heat-convection`, `:radiation`
- Structures: `:elasticity`, `:plasticity`, `:hyperelasticity`
- Electromagnetics: `:electrostatics`, `:magnetostatics`, `:maxwell`
- Multiphase: `:level-set`, `:vof`, `:phase-field`
- Particles: `:dem`, `:sph`, `:md`
- Custom: `:pde` with user-defined operators

### 5. Coupling (Multiphysics)

For coupled problems:

```lisp
;; Sequential coupling
(coupling :thermal-fluid
  :physics (thermal fluid)
  :type :sequential
  :order (thermal fluid)
  :iterations 10
  :tolerance 1e-6)

;; Strong coupling (monolithic)
(coupling :fsi
  :physics (fluid structure)
  :type :monolithic
  :interface moving-wall)

;; One-way coupling
(coupling :thermal-stress
  :physics (thermal structure)
  :type :one-way
  :from thermal :to structure
  :field :temperature)
```

### 6. Methods (Discretization Schemes)

Different numerical methods for different physics/regions:

```lisp
;; Global method
(method :dcpse
  :order 2
  :neighbors 25
  :support-radius 3.5)

;; Per-physics method
(physics :fluid
  :type :navier-stokes
  (method :fvm :scheme :upwind :limiter :van-leer))

(physics :structure
  :type :elasticity
  (method :fem :order 2 :integration :gauss))

;; Per-subdomain method
(subdomain solid
  (method :fem :order 1))
(subdomain fluid
  (method :fvm :scheme :central))
```

**Supported methods:**
- Mesh-based: `:fem`, `:fvm`, `:fdm`, `:dg`, `:spectral`
- Meshless: `:dcpse`, `:sph`, `:mls`, `:rbf`, `:rkpm`, `:efg`, `:mlpg`
- Isogeometric: `:iga`, `:nurbs`
- Particle: `:pic`, `:flip`, `:mpm`, `:dem`, `:md`
- Hybrid: `:coupled-fem-meshless`, `:coupled-particle-mesh`

### 7. Solvers

Time integration and linear/nonlinear solvers:

```lisp
;; Time integration
(time
  :method :bdf2          ; or :euler, :rk4, :newmark, :hht
  :dt 0.001
  :end 10.0
  :adaptive t
  :cfl 0.5
  :tolerance 1e-6)

;; Steady-state
(steady
  :method :newton
  :tolerance 1e-8
  :max-iterations 100
  :line-search t)

;; Linear solver
(linear-solver
  :method :gmres
  :preconditioner :ilu
  :tolerance 1e-10
  :max-iterations 1000)

;; Per-physics solver settings
(physics :fluid
  (solver
    :linear :gmres :preconditioner :amg
    :nonlinear :newton :damping 0.8))
```

### 8. Materials

Material properties for different subdomains:

```lisp
(materials
  (steel
    :density 7850
    :youngs-modulus 200e9
    :poissons-ratio 0.3
    :conductivity 50
    :specific-heat 500)

  (water
    :density 1000
    :viscosity 0.001
    :conductivity 0.6
    :specific-heat 4186)

  ;; Assign to subdomains
  (assign steel :to solid-region)
  (assign water :to fluid-region))

;; Or inline
(subdomain solid :material steel)
```

### 9. Output

Results and post-processing:

```lisp
(output
  :format :vtk           ; or :hdf5, :exodus, :cgns
  :directory "results/"
  :fields (velocity pressure temperature stress)

  ;; Time series
  :every 0.1             ; time interval
  :at (0 0.5 1.0 5.0)    ; specific times

  ;; Probes
  (probe :point (0.5 0.5) :fields (velocity pressure))
  (probe :line (0 0.5) (1 0.5) :n 100 :fields (velocity))

  ;; Derived quantities
  (monitor :drag :on wall :direction (1 0 0))
  (monitor :flux :on outlet :field temperature))
```

### 10. Complete Example

```lisp
(simulation "conjugate-heat-transfer" :dim 2

  ;; Domain with subdomains
  (domain (box (0 0) (2 1)) :mesh :n (80 40))

  (subdomains
    (fluid (< x 1))
    (solid (>= x 1)))

  (boundaries
    (inlet   (= x 0))
    (outlet  (= x 2))
    (walls   (or (= y 0) (= y 1)))
    (interface (= x 1)))

  ;; Materials
  (materials
    (water :density 1000 :viscosity 0.001
           :conductivity 0.6 :specific-heat 4186)
    (copper :density 8900 :conductivity 400
            :specific-heat 385))

  (subdomain fluid :material water)
  (subdomain solid :material copper)

  ;; Fluid physics (only in fluid subdomain)
  (physics :flow
    :type :navier-stokes
    :subdomain fluid

    (bc inlet  :velocity (0.1 0))
    (bc outlet :pressure 0)
    (bc walls  :no-slip)
    (bc interface :no-slip))

  ;; Thermal physics (both subdomains)
  (physics :thermal
    :type :heat-conduction

    (bc inlet  :temperature 20)
    (bc outlet :outflow)
    (bc walls  :insulated)

    ;; Conjugate interface - automatic continuity
    (interface interface :conjugate))

  ;; Coupling
  (coupling :flow-thermal
    :physics (flow thermal)
    :type :sequential
    :iterations 5)

  ;; Methods per subdomain
  (subdomain fluid
    (method :fvm :scheme :quick))
  (subdomain solid
    (method :fem :order 2))

  ;; Solver
  (time :method :bdf1 :dt 0.01 :end 100.0)

  ;; Output
  (output :format :vtk :every 1.0
    :fields (velocity pressure temperature)))
```

## Consequences

### Advantages

1. **Unified domain concept**: Mesh, points, and particles are all discrete representations of Ω, avoiding artificial distinctions.

2. **Separation of geometry and physics**: Boundaries and subdomains are defined once, then referenced by multiple physics. Same geometry can have different BCs for different physics.

3. **Dimension-agnostic**: The same syntax works for 1D, 2D, 3D (and higher). Mathematical terms like `box` generalize naturally.

4. **Method-agnostic**: Switch between FEM, FVM, meshless, particle methods by changing the `method` block, not rewriting the problem.

5. **Multiphysics-native**: Multiple `physics` blocks with coupling is first-class, not an afterthought.

6. **Readable**: Lisp syntax is minimal yet unambiguous. Keywords are self-documenting.

7. **Programmable**: Full Lisp power for parametric studies, loops, conditionals, user functions.

8. **Extensible**: New physics types, methods, and solvers can be added without syntax changes.

### Disadvantages

1. **Lisp syntax**: May be unfamiliar to users from ANSYS/STAR-CCM+ background. Parentheses can be intimidating.

2. **Verbosity for simple cases**: A basic lid-driven cavity still requires ~20 lines vs. a specialized CFD input format.

3. **Implementation complexity**: Supporting all advertised features (multiphysics, multi-method, etc.) is substantial engineering.

4. **Validation burden**: More flexibility means more combinations to test.

### Comparison to Existing Tools

| Feature | SDL | COMSOL | ANSYS | OpenFOAM | Abaqus |
|---------|-----|--------|-------|----------|--------|
| Geometry + mesh unified | ✓ | ✓ | partial | ✗ | ✗ |
| Dimension-agnostic | ✓ | ✓ | ✗ | ✗ | ✗ |
| Method-agnostic | ✓ | partial | ✗ | ✗ | ✗ |
| Multiphysics native | ✓ | ✓ | partial | ✗ | partial |
| Scriptable | ✓ | partial | ✓ | partial | ✓ |
| Open format | ✓ | ✗ | ✗ | ✓ | partial |
| Readable syntax | ✓ | N/A | ✗ | partial | ✗ |
| Meshless support | ✓ | ✗ | ✗ | ✗ | ✗ |
| Particle methods | ✓ | partial | partial | partial | partial |

### Migration Path

1. **Phase 1**: Implement core syntax for current StreamVorti capabilities (DCPSE, stream-vorticity, 2D)

2. **Phase 2**: Add mesh-based methods (FEM via MFEM), 3D support

3. **Phase 3**: Add multiphysics coupling, additional physics types

4. **Phase 4**: Add particle methods, isogeometric analysis

### Open Questions

1. **Coordinate names in higher dimensions**: Beyond (x, y, z), what names? (x₁, x₂, x₃, x₄) or (x, y, z, w)?

2. **Units**: Should we support units? `(velocity (1 m/s) (0 m/s))` vs. dimensionless?

3. **Parallel decomposition**: How to specify domain decomposition for MPI?

4. **Adaptive refinement**: Syntax for h/p/r refinement criteria?

5. **Optimization/sensitivity**: Syntax for defining objective functions, design variables?

## References

- COMSOL Multiphysics Reference Manual
- ANSYS APDL Command Reference
- OpenFOAM User Guide
- Abaqus Analysis User's Manual
- SfePy: Simple Finite Elements in Python (sfepy.org)
- MFEM: Modular Finite Element Methods Library (mfem.org)
- FEniCS Project Documentation
