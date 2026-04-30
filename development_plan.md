# JFVM Development Plan

Date: 2026-04-30

Goal: modernize the package internals, improve solver performance, and establish a robust testing/performance baseline while keeping the current user interface intact (or only very slightly changed with compatibility shims).

---

## Part I - Refactor + Efficiency Plan (implementation-oriented)

### 1) Current-state observations (from code review)

- Coordinate/geometry dispatch is encoded by floating-point tags in `MeshStructure.dimension` (examples: `1`, `1.5`, `2`, `2.5`, `2.8`, `3`, `3.2`).
- This numeric dispatch pattern is repeated in many files (`meshstructure`, `domainVariables`, `boundarycondition`, `diffusionterms`, `convectionTerms`, `calculusTerms`, `sourceTerms`, `transientTerms`, utilities), increasing maintenance cost and bug risk.
- Core solve API (`solveLinearPDE`) is hardcoded to `SparseMatrixCSC{Float64, Int64}` and uses direct `M\RHS` only.
- Time-stepping workflows rebuild many sparse structures in loops and allocate repeatedly.
- Type stability is limited by fields like `dimension::Real`, `dims::Array{Int,1}`, and value unions such as `Union{Array{<:Real}, DenseArray{Bool}}`.
- Existing `test/runtests.jl` is a smoke test only; there is no robust regression/performance test suite.
- `src/jfvm_test.jl` behaves like a demonstration script and does not enforce pass/fail correctness.

### 2) Non-breaking policy for this refactor

- Keep all currently exported function names and high-level workflows working.
- Keep all mesh constructor names (such as `createMesh1D`, `createMeshCylindrical2D`, etc.).
- Allow minor, controlled breakage only for internals and undocumented behavior.
- Any user-visible change must ship with:
	- deprecation warning (not immediate removal),
	- migration note,
	- regression tests proving old and new behavior equivalence.

### 3) Refactor architecture: replace floating-point dimension tags

#### 3.1 New coordinate/mesh taxonomy

Introduce explicit geometry/coordinate types:

```julia
abstract type AbstractCoordinateSystem end

struct Cartesian1D <: AbstractCoordinateSystem end
struct Cylindrical1D <: AbstractCoordinateSystem end
struct Cartesian2D <: AbstractCoordinateSystem end
struct Cylindrical2D <: AbstractCoordinateSystem end
struct Radial2D <: AbstractCoordinateSystem end
struct Cartesian3D <: AbstractCoordinateSystem end
struct Cylindrical3D <: AbstractCoordinateSystem end
```

Update mesh type to carry coordinate semantics directly:

```julia
struct MeshStructure{CS<:AbstractCoordinateSystem}
		coordinatesystem::CS
		dims::NTuple{3,Int}   # canonical storage, unused axes set to 1
		# existing geometry fields ...
end
```

Notes:
- Internal canonical dims should be fixed-size for type stability.
- Public constructors still accept 1D/2D/3D style arguments.
- Keep convenience helpers for active dimensions (`ndims_active(mesh)`), and geometry family checks (`is_cartesian(mesh)`, `is_cylindrical(mesh)`, `is_radial(mesh)`).

#### 3.2 Compatibility bridge

- Keep `dimension` available for one transition release as a computed compatibility property or mirrored field with deprecation warning.
- Maintain mapping:
	- `Cartesian1D -> 1.0`
	- `Cylindrical1D -> 1.5`
	- `Cartesian2D -> 2.0`
	- `Cylindrical2D -> 2.5`
	- `Radial2D -> 2.8`
	- `Cartesian3D -> 3.0`
	- `Cylindrical3D -> 3.2`
- Remove direct numeric comparisons from internal kernels incrementally and replace with dispatch/functions.

#### 3.3 Dispatch migration strategy

- Introduce internal dispatch entry points first, for example:
	- `diffusionTerm(::FaceValue, ::Cartesian2D)`
	- `diffusionTerm(::FaceValue, ::Cylindrical3D)`
- Temporary adapter:
	- `diffusionTerm(D::FaceValue) = diffusionTerm(D, D.domain.coordinatesystem)`
- Repeat for convection, divergence/gradient, source, transient, boundary, and averaging.
- Migrate one subsystem at a time, preserving behavior by regression tests before moving to next subsystem.

### 4) Data-structure and type-stability cleanup

#### 4.1 Prioritized type changes

- Change broad abstract fields to concrete parameterized forms where feasible.
- Replace `dims::Array{Int,1}` with tuple-backed representation internally.
- Rework value containers toward parameterized arrays:
	- `CellValue{T,N,A<:AbstractArray{T,N}}`
	- `FaceValue{T,AX,AY,AZ}` (or equivalent typed fields)
- Keep constructor front-end permissive (accept scalars, arrays, booleans), normalize internally.

#### 4.2 Allocation reduction targets

- Add `@views` and in-place broadcasting in hot paths.
- Eliminate avoidable `zeros(...)`/`repeat(...)` in term assembly.
- Cache frequently reused geometry vectors (cell widths, face widths, radial factors) per mesh.
- Introduce reusable sparse assembly workspaces for repeated transient solves.

### 5) Solver efficiency plan

#### 5.1 Linear solver abstraction layer

Add internal solver strategy object while preserving current API:

```julia
solveLinearPDE(mesh, A, b)                     # existing behavior preserved
solveLinearPDE(mesh, A, b; solver=DefaultLU())
solveLinearPDE!(mesh, A, b, phi; solver=...)
```

Recommended initial solver backends:
- Direct: built-in sparse `\` (default).
- Optional direct: MUMPS if available.
- Optional iterative: GMRES/BiCGStab/CG via `IterativeSolvers.jl` (for large systems).

#### 5.2 Factorization reuse

- Add optional reusable factorization path for repeated solves with fixed matrix pattern.
- In transient loops, separate matrix build and solve phases clearly.
- Add cache invalidation rules when coefficients change.

#### 5.3 Numeric generalization

- Relax strict `Float64` typing in solver signatures where safe.
- Preserve default Float64 path for stability and performance.

### 6) Testing strategy (must be built before aggressive refactor)

#### 6.1 Test suite structure

- Replace single smoke test with modules in `test/`:
	- mesh constructors and geometry invariants,
	- boundary condition assembly,
	- diffusion/convection/transient term assembly,
	- solver outputs for canonical PDE cases,
	- coordinate-system equivalence tests.

#### 6.2 Regression goldens

- Use small deterministic cases for each coordinate family.
- Validate:
	- matrix sparsity pattern,
	- residual norms,
	- solution norms/profiles,
	- mass/flux consistency checks.

#### 6.3 Compatibility tests

- Confirm old user-facing calls still run unchanged.
- Add deprecation tests for any intentionally adjusted behavior.

#### 6.4 Performance tests

- Add benchmark harness for representative problems (1D/2D/3D).
- Track:
	- assembly time,
	- solve time,
	- memory allocations,
	- iteration counts (iterative solvers).

### 7) Suggested implementation phases and milestones

#### Phase A - Stabilize baseline (1-2 weeks)

- Build robust tests from current behavior.
- Add benchmark scenarios and save baseline metrics.
- Fix obvious syntax/logic defects discovered by tests.

Exit criteria:
- tests pass consistently,
- benchmark baseline recorded,
- no API changes yet.

#### Phase B - Coordinate type introduction (1-2 weeks)

- Add coordinate structs and constructor wiring.
- Add compatibility mapping for old numeric dimensions.
- Migrate one low-risk subsystem first (for example source/transient terms).

Exit criteria:
- old and new pathways produce numerically equivalent outputs,
- no public constructor breakage.

#### Phase C - Full dispatch migration (2-4 weeks)

- Migrate diffusion, convection, boundary, calculus, averaging, utilities.
- Remove direct numeric comparisons from internals.
- Keep compatibility shim active.

Exit criteria:
- no internal `d==1.5`-style logic in core kernels,
- full test suite passes,
- no measurable regression in baseline benchmarks.

#### Phase D - Solver modernization + optimization (2-4 weeks)

- Add solver abstraction and optional iterative backends.
- Implement factorization reuse and key in-place optimizations.
- Benchmark and tune high-allocation hotspots.

Exit criteria:
- default solver behavior unchanged for existing users,
- improved performance on medium/large cases,
- documented solver options.

#### Phase E - Compatibility cleanup release (1 week)

- Keep compatibility warnings, prepare migration notes.
- Decide whether to retain or remove numeric `dimension` compatibility in next major version.

Exit criteria:
- clear migration guide,
- release notes with any slight breaks documented.

### 8) Main risks and mitigations

- Risk: hidden behavior differences across coordinate systems during migration.
	- Mitigation: matrix-level and solution-level regression tests per geometry.
- Risk: performance regressions from type refactor.
	- Mitigation: benchmark gate in CI, maintain fast path for Float64.
- Risk: optional solver dependencies increase maintenance burden.
	- Mitigation: keep optional backends behind feature flags and test matrix.
- Risk: over-refactor before tests are reliable.
	- Mitigation: enforce test-first sequence (Phase A before B/C/D).

---

## Part II - General package improvement roadmap (not yet implementation-ready)

This section is intentionally broader and exploratory.

### 1) Unstructured mesh support roadmap

#### 1.1 Target capabilities

- 2D polygonal and 3D polyhedral control volumes.
- Face-based connectivity with owner/neighbor topology.
- Native support for non-orthogonal corrections.
- Boundary patch groups with named tags.

#### 1.2 Architecture direction

- Introduce `AbstractMesh` with two concrete families:
	- `StructuredMesh` (current style),
	- `UnstructuredMesh` (new).
- Introduce geometry operator layer using mesh connectivity rather than index slicing.
- Separate discretization operators from storage layout.

#### 1.3 Ecosystem integration options

- Mesh import from Gmsh/VTK formats.
- Potential interoperability with Julia mesh ecosystems for preprocessing.

#### 1.4 Research items before implementation

- Data structure choice for sparse connectivity and cache efficiency.
- Non-orthogonal correction strategy and limiter design.
- Boundary condition API that works uniformly across structured/unstructured meshes.

### 2) Nonlinear solver roadmap

#### 2.1 Solver stack goals

- Picard (fixed-point) framework for mild nonlinearities.
- Newton or quasi-Newton framework for stronger nonlinear systems.
- Optional line search or trust-region globalisation.
- Residual/Jacobian monitoring and convergence diagnostics.

#### 2.2 API direction

- Add equation/problem object abstraction with callbacks:
	- residual assembly,
	- Jacobian assembly or Jacobian-free operator,
	- preconditioner assembly hooks.
- Keep existing linear workflows as a subset of the new formulation.

#### 2.3 Future accelerators

- Automatic differentiation support for Jacobians.
- Matrix-free Krylov methods for large coupled systems.
- Block preconditioners for multiphysics extensions.

### 3) Visualization and output roadmap

#### 3.1 Core output capabilities

- Built-in export to VTK/XDMF for external post-processing.
- Lightweight plotting recipes for 1D/2D in Julia plotting ecosystem.
- Optional interactive visualization backend for notebooks.

#### 3.2 UX goals

- Consistent plotting API for structured and future unstructured meshes.
- Field metadata (units, labels, variable names) carried with outputs.
- Standardized snapshot and animation utilities for transient runs.

### 4) Broader package improvements

- Documentation overhaul:
	- modern tutorials,
	- migration guide,
	- solver/mesh capability matrix.
- CI modernization:
	- multi-version Julia tests,
	- benchmark regression checks,
	- optional dependency matrix.
- Numerical quality:
	- conservation and boundedness test suite,
	- verification cases against analytical or manufactured solutions.
- Developer experience:
	- clearer module organization,
	- contribution guide,
	- coding standards for performance-critical kernels.

### 5) Long-term sequencing proposal

- Wave 1: finish Part I (stabilize/refactor/performance).
- Wave 2: prototype unstructured mesh core + minimal diffusion solver.
- Wave 3: nonlinear framework over structured meshes first, then unstructured.
- Wave 4: unified visualization/output pipeline across all mesh types.

### 6) Decision gates for Part II

Before committing to full implementation, define:

- target problem scales (small academic vs large industrial),
- preferred dependency strategy (minimal dependencies vs richer ecosystem),
- maintenance capacity for optional backends,
- performance targets for structured and unstructured modes.

---

## Immediate next actions

1. Start with Part I Phase A: create a real test suite from existing behavior.
2. Establish baseline benchmarks for at least one case in each geometry family.
3. Open a refactor branch focused only on coordinate-type introduction with compatibility shims.
