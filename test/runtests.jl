using JFVM
using Test
using SparseArrays

@testset "Coordinate system mapping" begin
	m1 = createMesh1D(5, 1.0)
	@test m1.coordinatesystem isa Cartesian1D
	@test m1.dimension == 1.0
	@test is_1d(m1)

	mc1 = createMeshCylindrical1D(5, 1.0)
	@test mc1.coordinatesystem isa Cylindrical1D
	@test mc1.dimension == 1.5
	@test is_1d(mc1)

	m2 = createMesh2D(4, 3, 1.0, 2.0)
	@test m2.coordinatesystem isa Cartesian2D
	@test is_2d(m2)

	mr2 = createMeshRadial2D(4, 5, 1.0, 2.0*pi)
	@test mr2.coordinatesystem isa Radial2D
	@test mr2.dimension == 2.8
	@test is_2d(mr2)

	m3 = createMesh3D(3, 2, 2, 1.0, 1.0, 1.0)
	@test m3.coordinatesystem isa Cartesian3D
	@test is_3d(m3)

	mc3 = createMeshCylindrical3D(3, 4, 2, 1.0, pi, 1.0)
	@test mc3.coordinatesystem isa Cylindrical3D
	@test is_3d(mc3)
end

@testset "Core PDE building blocks" begin
	m = createMesh2D(4, 3, 1.0, 1.0)
	BC = createBC(m)

	Mbc, RHSbc = boundaryConditionTerm(BC)
	n = prod(m.dims .+ 2)
	@test size(Mbc) == (n, n)
	@test length(RHSbc) == n

	phi = createCellVariable(m, 2.0, BC)
	D = harmonicMean(phi)
	Mdif = diffusionTerm(D)
	@test size(Mdif) == (n, n)

	u = createFaceVariable(m, 0.1)
	Mconv = convectionTerm(u)
	@test size(Mconv) == (n, n)

	divu = divergenceTerm(u)
	@test length(divu) == n

	gradphi = gradientTerm(phi)
	@test size(gradphi.xvalue) == (m.dims[1] + 1, m.dims[2])
end

@testset "Linear and explicit solvers" begin
	m = createMesh1D(6, 1.0)
	n = prod(m.dims .+ 2)
	A = spdiagm(0 => ones(Float64, n))
	b = collect(1.0:n)

	phi = solveLinearPDE(m, A, b)
	@test phi.value == reshape(b, tuple(m.dims .+ 2...))

	phi2 = createCellVariable(m, 0.0)
	solveLinearPDE!(m, A, b, phi2)
	@test phi2.value == reshape(b, tuple(m.dims .+ 2...))

	BC = createBC(m)
	rhs = zeros(Float64, n)
	phi_exp = solveExplicitPDE(phi, 0.1, rhs, BC)
	@test length(internalCells(phi_exp)) == m.dims[1]
end

@testset "Legacy smoke workflow" begin
	@test_nowarn JFVM_test()
end
