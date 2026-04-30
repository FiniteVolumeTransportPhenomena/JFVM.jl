# ===============================
# Written by AAE
# TU Delft, Spring 2014
# simulkade.com
# ===============================

# ================================ CREATE BOUNDARY CONDITION =======================
function createBC(m::MeshStructure)
  createBC(m, m.coordinatesystem)
end

function createBC(m::MeshStructure, ::Cartesian1D)
  BoundaryCondition(m,
      BorderValue([1.0], [0.0], [0.0], false),
      BorderValue([1.0], [0.0], [0.0], false),
      BorderValue([0.0], [0.0], [0.0], false),
      BorderValue([0.0], [0.0], [0.0], false),
      BorderValue([0.0], [0.0], [0.0], false),
      BorderValue([0.0], [0.0], [0.0], false))
end

createBC(m::MeshStructure, ::Cylindrical1D) = createBC(m, Cartesian1D())

function createBC(m::MeshStructure, ::Cartesian2D)
  Nx = m.dims[1]
  Ny = m.dims[2]
  BoundaryCondition(m,
      BorderValue(ones(1,Ny), zeros(1,Ny), zeros(1,Ny), false),
      BorderValue(ones(1,Ny), zeros(1,Ny), zeros(1,Ny), false),
      BorderValue(ones(Nx,1), zeros(Nx,1), zeros(Nx,1), false),
      BorderValue(ones(Nx,1), zeros(Nx,1), zeros(Nx,1), false),
      BorderValue([0.0], [0.0], [0.0], false),
      BorderValue([0.0], [0.0], [0.0], false))
    end

    createBC(m::MeshStructure, ::Cylindrical2D) = createBC(m, Cartesian2D())
    createBC(m::MeshStructure, ::Radial2D) = createBC(m, Cartesian2D())

    function createBC(m::MeshStructure, ::Cartesian3D)
  Nx = m.dims[1]
  Ny = m.dims[2]
  Nz = m.dims[3]
  BoundaryCondition(m,
      BorderValue(ones(1,Ny,Nz), zeros(1,Ny,Nz), zeros(1,Ny,Nz), false),
      BorderValue(ones(1,Ny,Nz), zeros(1,Ny,Nz), zeros(1,Ny,Nz), false),
      BorderValue(ones(Nx,1,Nz), zeros(Nx,1,Nz), zeros(Nx,1,Nz), false),
      BorderValue(ones(Nx,1,Nz), zeros(Nx,1,Nz), zeros(Nx,1,Nz), false),
      BorderValue(ones(Nx,Ny,1), zeros(Nx,Ny,1), zeros(Nx,Ny,1), false),
      BorderValue(ones(Nx,Ny,1), zeros(Nx,Ny,1), zeros(Nx,Ny,1), false))
end

createBC(m::MeshStructure, ::Cylindrical3D) = createBC(m, Cartesian3D())

function boundaryConditionTerm(BC::BoundaryCondition)
  boundaryConditionTerm(BC, BC.domain.coordinatesystem)
end

boundaryConditionTerm(BC::BoundaryCondition, ::Cartesian1D) = boundaryCondition1D(BC)
boundaryConditionTerm(BC::BoundaryCondition, ::Cylindrical1D) = boundaryCondition1D(BC)
boundaryConditionTerm(BC::BoundaryCondition, ::Cartesian2D) = boundaryCondition2D(BC)
boundaryConditionTerm(BC::BoundaryCondition, ::Cylindrical2D) = boundaryCondition2D(BC)
boundaryConditionTerm(BC::BoundaryCondition, ::Radial2D) = boundaryConditionRadial2D(BC)
boundaryConditionTerm(BC::BoundaryCondition, ::Cartesian3D) = boundaryCondition3D(BC)
boundaryConditionTerm(BC::BoundaryCondition, ::Cylindrical3D) = boundaryConditionCylindrical3D(BC)

# =============================== BOUNDARY CARTESIAN 1D ============================
function boundaryCondition1D(BC::BoundaryCondition)
# creates the matrix of coefficients and RHS for
# a boundary condition structure
Nx = BC.domain.dims[1]
G = [1:Nx+2;]
dx_1 = BC.domain.cellsize.x[1]
dx_end = BC.domain.cellsize.x[end]
# number of boundary nodes:
nb = 4

# define the vectors to be used for the creation of the sparse matrix
ii = zeros(Int64,nb)
jj = zeros(Int64,nb)
s = zeros(Float64, nb) # Float64 by default, but specify the type just in case

# define the RHS column vector
BCRHS = zeros(Nx+2)

q = 0
# Assign values to the boundary condition matrix and the RHS vector based
# on the BC structure
if !BC.right.periodic && !BC.left.periodic # non-periodic boundary condition
    # Right boundary
    i = Nx+2
    q=q+1
    ii[q] = G[i]
    jj[q] = G[i]
    s[q] = BC.right.b[1]/2.0 + BC.right.a[1]/dx_end
    q=q+1
    ii[q] = G[i]
    jj[q] = G[i-1]
    s[q] = BC.right.b[1]/2.0 - BC.right.a[1]/dx_end
    BCRHS[G[i]] = BC.right.c[1]

    # Left boundary
    i = 1
    q=q+1
    ii[q] = G[i]
    jj[q] = G[i+1]
    s[q] = -(BC.left.b[1]/2.0 + BC.left.a[1]/dx_1)
    q=q+1
    ii[q] = G[i]
    jj[q] = G[i]
    s[q] = -(BC.left.b[1]/2.0 - BC.left.a[1]/dx_1)
    BCRHS[G[i]] = -BC.left.c[1]
elseif BC.right.periodic || BC.left.periodic  # periodic boundary condition
    nb=8
    ii = zeros(Int64,nb)
    jj = zeros(Int64,nb)
    s = zeros(Float64, nb)
    # Right boundary
    i = Nx+2
    q=q+1
    ii[q] = G[Nx+2]
    jj[q] = G[Nx+2]
    s[q] .= 1.0
    q=q+1
    ii[q] = G[Nx+2]
    jj[q] = G[Nx+1]
    s[q] .= -1.0
    q=q+1
    ii[q] = G[Nx+2]
    jj[q] = G[1]
    s[q] = dx_end/dx_1
    q=q+1
    ii[q] = G[Nx+2]
    jj[q] = G[2]
    s[q] = -dx_end/dx_1
    BCRHS[G[i]] .= 0.0

    # Left boundary
    i = 1
    q=q+1
    ii[q] = G[1]
    jj[q] = G[1]
    s[q] .= 1.0
    q=q+1
    ii[q] = G[1]
    jj[q] = G[2]
    s[q] .= 1.0
    q=q+1
    ii[q] = G[1]
    jj[q] = G[Nx+1]
    s[q] .= -1.0;
    q=q+1
    ii[q] = G[1]
    jj[q] = G[Nx+2]
    s[q] .= -1.0;
    BCRHS[G[i]] .= 0.0
end

# Build the sparse matrix of the boundary conditions
BCMatrix = sparse(ii[1:q], jj[1:q], s[1:q], Nx+2, Nx+2)
(BCMatrix, BCRHS)

end



# =============================== BOUNDARY CARTESIAN 2D ============================
function boundaryCondition2D(BC::BoundaryCondition)
# creates the matrix of coefficients and RHS for
# a boundary condition structure
Nx, Ny = tuple(BC.domain.dims...)
G=reshape([1:(Nx+2)*(Ny+2);], Nx+2, Ny+2)
dx_1 = BC.domain.cellsize.x[1]
dx_end = BC.domain.cellsize.x[end]
dy_1 = BC.domain.cellsize.y[1]
dy_end = BC.domain.cellsize.y[end]

# number of boundary nodes:
nb = 8*(Nx+Ny+2)

# define the vectors to be used for the creation of the sparse matrix
ii = zeros(Int64,nb)
jj = zeros(Int64,nb)
s = zeros(Float64, nb) # Float64 by default, but specify the type just in case

# define the RHS column vector
BCRHS = zeros((Nx+2)*(Ny+2))

for q = 1:4
  ii[q] = BC.domain.corner[q]
  jj[q] = BC.domain.corner[q]
  s[q] = maximum(BC.top.b/2.0+BC.top.a/dy_end)
  BCRHS[BC.domain.corner[q]] = 0.0
end
q = 4
# Assign values to the boundary condition matrix and the RHS vector based
# on the BC structure
if !BC.top.periodic && !BC.bottom.periodic # non-periodic boundary condition
  # top boundary
  j=Ny+2
  for i=2:Nx+1
    q+=1
    ii[q] = G[i,j]
    jj[q] = G[i,j]
    s[q] = BC.top.b[i-1]/2.0 + BC.top.a[i-1]/dy_end
    q+=1
    ii[q] = G[i,j]
    jj[q] = G[i,j-1]
    s[q] = BC.top.b[i-1]/2.0 - BC.top.a[i-1]/dy_end
    BCRHS[G[i,j]] = BC.top.c[i-1]
  end
  # bottom boundary
  j=1
  for i=2:Nx+1
    q+=1
    ii[q] = G[i,j]
    jj[q] = G[i,j+1]
    s[q] = -(BC.bottom.b[i-1]/2.0 + BC.bottom.a[i-1]/dy_1)
    q+=1
    ii[q] = G[i,j]
    jj[q] = G[i,j]
    s[q] = -(BC.bottom.b[i-1]/2.0 - BC.bottom.a[i-1]/dy_1)
    BCRHS[G[i,j]] = -(BC.bottom.c[i-1])
  end
elseif BC.top.periodic || BC.bottom.periodic  # periodic boundary condition
  # top boundary
  j=Ny+2
  for i=2:Nx+1
    q+=1
    ii[q] = G[i,j]
    jj[q] = G[i,j]
    s[q] .= 1.0
    q+=1
    ii[q] = G[i,j]
    jj[q] = G[i,j-1]
    s[q] .= -1.0
    q+=1
    ii[q] = G[i,j]
    jj[q] = G[i,1]
    s[q] = dy_end/dy_1
    q+=1
    ii[q] = G[i,j]
    jj[q] = G[i,2]
    s[q] = -dy_end/dy_1
    BCRHS[G[i,j]] .= 0.0
  end
  # bottom boundary
  j=1
  for i=2:Nx+1
    q+=1
    ii[q] = G[i,j]
    jj[q] = G[i,j]
    s[q] .= 1.0
    q+=1
    ii[q] = G[i,j]
    jj[q] = G[i,j+1]
    s[q] .= 1.0
    q+=1
    ii[q] = G[i,j]
    jj[q] = G[i,Ny+1]
    s[q] .= -1.0
    q+=1
    ii[q] = G[i,j]
    jj[q] = G[i,Ny+2]
    s[q] .= -1.0
    BCRHS[G[i,j]] .= 0.0
  end
end

if !BC.right.periodic && !BC.left.periodic # non-periodic boundary condition
  # Right boundary
  i=Nx+2
  for j=2:Ny+1
    q+=1
    ii[q] = G[i,j]
    jj[q] = G[i,j]
    s[q] = BC.right.b[j-1]/2.0 + BC.right.a[j-1]/dx_end
    q+=1
    ii[q] = G[i,j]
    jj[q] = G[i-1,j]
    s[q] = BC.right.b[j-1]/2.0 - BC.right.a[j-1]/dx_end
    BCRHS[G[i,j]] = BC.right.c[j-1]
  end
  # Left boundary
  i = 1
  for j=2:Ny+1
    q+=1
    ii[q] = G[i,j]
    jj[q] = G[i+1,j]
    s[q] = -(BC.left.b[j-1]/2.0 + BC.left.a[j-1]/dx_1)
    q+=1
    ii[q] = G[i,j]
    jj[q] = G[i,j]
    s[q] = -(BC.left.b[j-1]/2.0 - BC.left.a[j-1]/dx_1)
    BCRHS[G[i,j]] = -(BC.left.c[j-1])
  end
elseif BC.right.periodic || BC.left.periodic  # periodic boundary condition
  # Right boundary
  i=Nx+2
  for j=2:Ny+1
    q+=1
    ii[q] = G[i,j]
    jj[q] = G[i,j]
    s[q] .= 1.0
    q+=1
    ii[q] = G[i,j]
    jj[q] = G[i-1,j]
    s[q] .= -1.0
    q+=1
    ii[q] = G[i,j]
    jj[q] = G[1,j]
    s[q] = dx_end/dx_1
    q+=1
    ii[q] = G[i,j]
    jj[q] = G[2,j]
    s[q] = -dx_end/dx_1
    BCRHS[G[i,j]] .= 0.0
  end
  # Left boundary
  i = 1;
  for j=2:Ny+1
    q+=1
    ii[q] = G[i,j]
    jj[q] = G[i,j]
    s[q] .= 1.0
    q+=1
    ii[q] = G[i,j]
    jj[q] = G[i+1,j]
    s[q] .= 1.0
    q+=1
    ii[q] = G[i,j]
    jj[q] = G[Nx+1,j]
    s[q] .= -1.0
    q+=1
    ii[q] = G[i,j]
    jj[q] = G[Nx+2,j]
    s[q] .= -1.0
    BCRHS[G[i,j]] .= 0.0
  end
end

# Build the sparse matrix of the boundary conditions
BCMatrix = sparse(ii[1:q], jj[1:q], s[1:q], (Nx+2)*(Ny+2), (Nx+2)*(Ny+2))
(BCMatrix, BCRHS)
end

# =============================== BOUNDARY Radial 2D ============================
function boundaryConditionRadial2D(BC::BoundaryCondition)
# creates the matrix of coefficients and RHS for
# a boundary condition structure
Nx, Ntheta = tuple(BC.domain.dims...)
G=reshape([1:(Nx+2)*(Ntheta+2);], Nx+2, Ntheta+2)
dx_1 = BC.domain.cellsize.x[1]
dx_end = BC.domain.cellsize.x[end]
dtheta_1 = BC.domain.cellsize.y[1]
dtheta_end = BC.domain.cellsize.y[end]
rp = BC.domain.cellcenters.x
# number of boundary nodes:
nb = 8*(Nx+Ntheta+2)

# define the vectors to be used for the creation of the sparse matrix
ii = zeros(Int64,nb)
jj = zeros(Int64,nb)
s = zeros(Float64, nb) # Float64 by default, but specify the type just in case

# define the RHS column vector
BCRHS = zeros((Nx+2)*(Ntheta+2))

for q = 1:4
  ii[q] = BC.domain.corner[q]
  jj[q] = BC.domain.corner[q]
  s[q] = maximum(BC.top.b/2.0+BC.top.a./(dtheta_end*rp))
  BCRHS[BC.domain.corner[q]] = 0.0
end
q = 4
# Assign values to the boundary condition matrix and the RHS vector based
# on the BC structure
if !BC.top.periodic && !BC.bottom.periodic # non-periodic boundary condition
  # top boundary
  j=Ntheta+2
  for i=2:Nx+1
    q+=1
    ii[q] = G[i,j]
    jj[q] = G[i,j]
    s[q] = BC.top.b[i-1]/2.0 + BC.top.a[i-1]/(dtheta_end*rp[i-1])
    q+=1
    ii[q] = G[i,j]
    jj[q] = G[i,j-1]
    s[q] = BC.top.b[i-1]/2.0 - BC.top.a[i-1]/(dtheta_end*rp[i-1])
    BCRHS[G[i,j]] = BC.top.c[i-1]
  end
  # bottom boundary
  j=1
  for i=2:Nx+1
    q+=1
    ii[q] = G[i,j]
    jj[q] = G[i,j+1]
    s[q] = -(BC.bottom.b[i-1]/2.0 + BC.bottom.a[i-1]/(dtheta_1*rp[i-1]))
    q+=1
    ii[q] = G[i,j]
    jj[q] = G[i,j]
    s[q] = -(BC.bottom.b[i-1]/2.0 - BC.bottom.a[i-1]/(dtheta_1*rp[i-1]))
    BCRHS[G[i,j]] = -(BC.bottom.c[i-1])
  end
elseif BC.top.periodic || BC.bottom.periodic  # periodic boundary condition
  # top boundary
  j=Ntheta+2
  for i=2:Nx+1
    q+=1
    ii[q] = G[i,j]
    jj[q] = G[i,j]
    s[q] = 1.0
    q+=1
    ii[q] = G[i,j]
    jj[q] = G[i,j-1]
    s[q] = -1.0
    q+=1
    ii[q] = G[i,j]
    jj[q] = G[i,1]
    s[q] = dtheta_end/dtheta_1
    q+=1
    ii[q] = G[i,j]
    jj[q] = G[i,2]
    s[q] = -dtheta_end/dtheta_1
    BCRHS[G[i,j]] = 0.0
  end
  # bottom boundary
  j=1
  for i=2:Nx+1
    q+=1
    ii[q] = G[i,j]
    jj[q] = G[i,j]
    s[q] = 1.0
    q+=1
    ii[q] = G[i,j]
    jj[q] = G[i,j+1]
    s[q] = 1.0
    q+=1
    ii[q] = G[i,j]
    jj[q] = G[i,Ntheta+1]
    s[q] = -1.0
    q+=1
    ii[q] = G[i,j]
    jj[q] = G[i,Ntheta+2]
    s[q] = -1.0
    BCRHS[G[i,j]] = 0.0
  end
end

if !BC.right.periodic && !BC.left.periodic # non-periodic boundary condition
  # Right boundary
  i=Nx+2
  for j=2:Ntheta+1
    q+=1
    ii[q] = G[i,j]
    jj[q] = G[i,j]
    s[q] = BC.right.b[j-1]/2.0 + BC.right.a[j-1]/dx_end
    q+=1
    ii[q] = G[i,j]
    jj[q] = G[i-1,j]
    s[q] = BC.right.b[j-1]/2.0 - BC.right.a[j-1]/dx_end
    BCRHS[G[i,j]] = BC.right.c[j-1]
  end
  # Left boundary
  i = 1
  for j=2:Ntheta+1
    q+=1
    ii[q] = G[i,j]
    jj[q] = G[i+1,j]
    s[q] = -(BC.left.b[j-1]/2.0 + BC.left.a[j-1]/dx_1)
    q+=1
    ii[q] = G[i,j]
    jj[q] = G[i,j]
    s[q] = -(BC.left.b[j-1]/2.0 - BC.left.a[j-1]/dx_1)
    BCRHS[G[i,j]] = -(BC.left.c[j-1])
  end
elseif BC.right.periodic || BC.left.periodic  # periodic boundary condition
  # Right boundary
  i=Nx+2
  for j=2:Ntheta+1
    q+=1
    ii[q] = G[i,j]
    jj[q] = G[i,j]
    s[q] .= 1.0
    q+=1
    ii[q] = G[i,j]
    jj[q] = G[i-1,j]
    s[q] .= -1.0
    q+=1
    ii[q] = G[i,j]
    jj[q] = G[1,j]
    s[q] = dx_end/dx_1
    q+=1
    ii[q] = G[i,j]
    jj[q] = G[2,j]
    s[q] = -dx_end/dx_1
    BCRHS[G[i,j]] .= 0.0
  end
  # Left boundary
  i = 1;
  for j=2:Ntheta+1
    q+=1
    ii[q] = G[i,j]
    jj[q] = G[i,j]
    s[q] .= 1.0
    q+=1
    ii[q] = G[i,j]
    jj[q] = G[i+1,j]
    s[q] .= 1.0
    q+=1
    ii[q] = G[i,j]
    jj[q] = G[Nx+1,j]
    s[q] .= -1.0
    q+=1
    ii[q] = G[i,j]
    jj[q] = G[Nx+2,j]
    s[q] .= -1.0
    BCRHS[G[i,j]] .= 0.0
  end
end

# Build the sparse matrix of the boundary conditions
BCMatrix = sparse(ii[1:q], jj[1:q], s[1:q], (Nx+2)*(Ntheta+2), (Nx+2)*(Ntheta+2))
(BCMatrix, BCRHS)
end


# ===================================== BOUNDARY 3D =====================================
function boundaryCondition3D(BC::BoundaryCondition)
# creates the matrix of coefficients and RHS for
# a boundary condition structure
Nx, Ny, Nz = tuple(BC.domain.dims...)
G=reshape([1:(Nx+2)*(Ny+2)*(Nz+2);], Nx+2, Ny+2, Nz+2)
dx_1 = BC.domain.cellsize.x[1]
dx_end = BC.domain.cellsize.x[end]
dy_1 = BC.domain.cellsize.y[1]
dy_end = BC.domain.cellsize.y[end]
dz_1 = BC.domain.cellsize.z[1]
dz_end = BC.domain.cellsize.z[end]

# number of boundary nodes (exact number is 2[(m+1)(n+1)*(n+1)*(p+1)+(m+1)*p+1]:
nb = 8*((Nx+1)*(Ny+1)+(Nx+1)*(Nz+1)+(Ny+1)*(Nz+1))

# define the vectors to be used for the creation of the sparse matrix
ii = zeros(Int64, nb)
jj = zeros(Int64, nb)
s = zeros(Float64, nb)

# define the RHS column vector
BCRHS = zeros(Float64, (Nx+2)*(Ny+2)*(Nz+2))

# assign value to the corner nodes (useless cells)
q = 1:8
ii[q] = BC.domain.corner
jj[q] = BC.domain.corner
s[q] .= 1.0
BCRHS[BC.domain.corner] .= 0.0

# assign values to the edges (useless cells)
q = q[end].+[1:length(BC.domain.edge);]
ii[q] = BC.domain.edge
jj[q] = BC.domain.edge
s[q] .= 1.0
BCRHS[BC.domain.edge] .= 0.0

# Assign values to the boundary condition matrix and the RHS vector based
# on the BC structure
if !BC.top.periodic && !BC.bottom.periodic
    # top boundary
    j=Ny+2
    i=2:Nx+1
    k=2:Nz+1
    q = q[end].+[1:Nx*Nz;]
    ii[q] = G[i,j,k]
    jj[q] = G[i,j,k]
    s[q] = BC.top.b/2.0 + BC.top.a/dy_end
    q = q[end].+[1:Nx*Nz;]
    ii[q] = G[i,j,k]
    jj[q] = G[i,j-1,k]
    s[q] = BC.top.b/2.0 - BC.top.a/dy_end
    BCRHS[G[i,j,k]] = BC.top.c

    # Bottom boundary
    j=1
    i=2:Nx+1
    k=2:Nz+1
    q = q[end].+[1:Nx*Nz;]
    ii[q] = G[i,j,k]
    jj[q] = G[i,j+1,k]
    s[q] = -(BC.bottom.b/2.0 + BC.bottom.a/dy_1) # consider the reverse direction of normal
    q = q[end].+[1:Nx*Nz;]
    ii[q] = G[i,j,k]
    jj[q] = G[i,j,k]
    s[q] = -(BC.bottom.b/2.0 - BC.bottom.a/dy_1) # consider the reverse direction of normal
    BCRHS[G[i,j,k]] = -(BC.bottom.c)
elseif BC.top.periodic || BC.bottom.periodic # periodic
    # top boundary
    j=Ny+2
    i=2:Nx+1
    k=2:Nz+1
    q = q[end].+[1:Nx*Nz;]
    ii[q] = G[i,j,k]
    jj[q] = G[i,j,k]
    s[q] = 1
    q = q[end].+[1:Nx*Nz;]
    ii[q] = G[i,j,k]
    jj[q] = G[i,j-1,k]
    s[q] = -1
    q = q[end].+[1:Nx*Nz;]
    ii[q] = G[i,j,k]
    jj[q] = G[i,1,k]
    s[q] = dy_end/dy_1
    q = q[end].+[1:Nx*Nz;]
    ii[q] = G[i,j,k]
    jj[q] = G[i,2,k]
    s[q] = -dy_end/dy_1
    BCRHS[G[i,j,k]] .= 0.0

    # Bottom boundary
    j=1
    i=2:Nx+1
    k=2:Nz+1
    q = q[end].+[1:Nx*Nz;]
    ii[q] = G[i,j,k]
    jj[q] = G[i,j,k]
    s[q] .= 1.0
    q = q[end].+[1:Nx*Nz;]
    ii[q] = G[i,j,k]
    jj[q] = G[i,j+1,k]
    s[q] .= 1.0
    q = q[end].+[1:Nx*Nz;]
    ii[q] = G[i,j,k]
    jj[q] = G[i,Ny+1,k]
    s[q] .= -1.0
    q = q[end].+[1:Nx*Nz;]
    ii[q] = G[i,j,k]
    jj[q] = G[i,Ny+2,k]
    s[q] .= -1.0
    BCRHS[G[i,j,k]] .= 0.0
end

if !BC.right.periodic && !BC.left.periodic
    # Right boundary
    i=Nx+2
    j=2:Ny+1
    k=2:Nz+1
    q = q[end].+[1:Ny*Nz;]
    ii[q] = G[i,j,k]
    jj[q] = G[i,j,k]
    s[q] = BC.right.b/2.0 + BC.right.a/dx_end
    q = q[end].+[1:Ny*Nz;]
    ii[q] = G[i,j,k]
    jj[q] = G[i-1,j,k]
    s[q] = BC.right.b/2.0 - BC.right.a/dx_end
    BCRHS[G[i,j,k]] = BC.right.c

    # Left boundary
    i = 1
    j=2:Ny+1
    k=2:Nz+1
    q = q[end].+[1:Ny*Nz;]
    ii[q] = G[i,j,k]
    jj[q] = G[i+1,j,k]
    s[q] = -(BC.left.b/2.0 + BC.left.a/dx_1)
    # consider the reverse direction of normal
    q = q[end].+[1:Ny*Nz;]
    ii[q] = G[i,j,k]
    jj[q] = G[i,j,k]
    s[q] = -(BC.left.b/2.0 - BC.left.a/dx_1)  # consider the reverse direction of normal
    BCRHS[G[i,j,k]] = -(BC.left.c)
elseif BC.right.periodic || BC.left.periodic # periodic
    # Right boundary
    i=Nx+2
    j=2:Ny+1
    k=2:Nz+1
    q = q[end].+[1:Ny*Nz;]
    ii[q] = G[i,j,k]
    jj[q] = G[i,j,k]
    s[q] .= 1.0
    q = q[end].+[1:Ny*Nz;]
    ii[q] = G[i,j,k]
    jj[q] = G[i-1,j,k]
    s[q] .= -1.0
    q = q[end].+[1:Ny*Nz;]
    ii[q] = G[i,j,k]
    jj[q] = G[1,j,k]
    s[q] = dx_end/dx_1
    q = q[end].+[1:Ny*Nz;]
    ii[q] = G[i,j,k]
    jj[q] = G[2,j,k]
    s[q] = -dx_end/dx_1
    BCRHS[G[i,j,k]] .= 0.0

    # Left boundary
    i = 1
    j=2:Ny+1
    k=2:Nz+1
    q = q[end].+[1:Ny*Nz;]
    ii[q] = G[i,j,k]
    jj[q] = G[i,j,k]
    s[q] .= 1.0
    q = q[end].+[1:Ny*Nz;]
    ii[q] = G[i,j,k]
    jj[q] = G[i+1,j,k]
    s[q] .= 1.0
    q = q[end].+[1:Ny*Nz;]
    ii[q] = G[i,j,k]
    jj[q] = G[Nx+1,j,k]
    s[q] .= -1.0
    q = q[end].+[1:Ny*Nz;]
    ii[q] = G[i,j,k]
    jj[q] = G[Nx+2,j,k]
    s[q] .= -1.0
    BCRHS[G[i,j,k]] .= 0.0
end

if !BC.front.periodic && !BC.back.periodic
    # Back boundary
    k=1
    i = 2:Nx+1
    j=2:Ny+1
    q = q[end].+[1:Nx*Ny;]
    ii[q] = G[i,j,k]
    jj[q] = G[i,j,k+1]
    s[q] = -(BC.back.b/2.0 + BC.back.a/dz_1)  # consider the reverse direction of normal
    q = q[end].+[1:Nx*Ny;]
    ii[q] = G[i,j,k]
    jj[q] = G[i,j,k]
    s[q] = -(BC.back.b/2.0 - BC.back.a/dz_1)  # consider the reverse direction of normal
    BCRHS[G[i,j,k]] = -(BC.back.c)

    # Front boundary
    k=Nz+2
    i = 2:Nx+1
    j=2:Ny+1
    q = q[end].+[1:Nx*Ny;]
    ii[q] = G[i,j,k]
    jj[q] = G[i,j,k]
    s[q] = BC.front.b/2.0 + BC.front.a/dz_end
    q = q[end].+[1:Nx*Ny;]
    ii[q] = G[i,j,k]
    jj[q] = G[i,j,k-1]
    s[q] = BC.front.b/2.0 - BC.front.a/dz_end
    BCRHS[G[i,j,k]] = BC.front.c
elseif BC.front.periodic || BC.back.periodic  # periodic
    # Back boundary
    k=1
    i = 2:Nx+1
    j=2:Ny+1
    q = q[end].+[1:Nx*Ny;]
    ii[q] = G[i,j,k]
    jj[q] = G[i,j,k]
    s[q] .= 1.0
    q = q[end].+[1:Nx*Ny;]
    ii[q] = G[i,j,k]
    jj[q] = G[i,j,k+1]
    s[q] .= 1.0
    q = q[end].+[1:Nx*Ny;]
    ii[q] = G[i,j,k]
    jj[q] = G[i,j,Nz+1]
    s[q] .= -1.0
    q = q[end].+[1:Nx*Ny;]
    ii[q] = G[i,j,k]
    jj[q] = G[i,j,Nz+2]
    s[q] .= -1.0
    BCRHS[G[i,j,k]] .= 0.0

    # Front boundary
    k=Nz+2
    i = 2:Nx+1
    j=2:Ny+1
    q = q[end].+[1:Nx*Ny;]
    ii[q] = G[i,j,k]
    jj[q] = G[i,j,k]
    s[q] .= 1.0
    q = q[end].+[1:Nx*Ny;]
    ii[q] = G[i,j,k]
    jj[q] = G[i,j,k-1]
    s[q] .= -1.0
    q = q[end].+[1:Nx*Ny;]
    ii[q] = G[i,j,k]
    jj[q] = G[i,j,1]
    s[q] = dz_end/dz_1
    q = q[end].+[1:Nx*Ny;]
    ii[q] = G[i,j,k]
    jj[q] = G[i,j,2]
    s[q] = -dz_end/dz_1
    BCRHS[G[i,j,k]] .= 0.0
end

# Build the sparse matrix of the boundary conditions
BCMatrix = sparse(ii[1:q[end]], jj[1:q[end]], s[1:q[end]],
    (Nx+2)*(Ny+2)*(Nz+2), (Nx+2)*(Ny+2)*(Nz+2))

(BCMatrix, BCRHS)
end

# ===================================== BOUNDARY Cylindrical 3D =====================================
function boundaryConditionCylindrical3D(BC::BoundaryCondition)
# creates the matrix of coefficients and RHS for
# a boundary condition structure
Nx, Ntheta, Nz = tuple(BC.domain.dims...)
G=reshape([1:(Nx+2)*(Ntheta+2)*(Nz+2);], Nx+2, Ntheta+2, Nz+2)
dx_1 = BC.domain.cellsize.x[1]
dx_end = BC.domain.cellsize.x[end]
dtheta_1 = BC.domain.cellsize.y[1]
dtheta_end = BC.domain.cellsize.y[end]
dz_1 = BC.domain.cellsize.z[1]
dz_end = BC.domain.cellsize.z[end]

rp = BC.domain.cellcenters.x
#rp = repmat(BC.domain.cellcenters.x', 1, Nz)

# number of boundary nodes (exact number is 2[(m+1)(n+1)*(n+1)*(p+1)+(m+1)*p+1]:
nb = 8*((Nx+1)*(Ntheta+1)+(Nx+1)*(Nz+1)+(Ntheta+1)*(Nz+1))

# define the vectors to be used for the creation of the sparse matrix
ii = zeros(Int64, nb)
jj = zeros(Int64, nb)
s = zeros(Float64, nb)

# define the RHS column vector
BCRHS = zeros(Float64, (Nx+2)*(Ntheta+2)*(Nz+2))

# assign value to the corner nodes (useless cells)
q = 1:8
ii[q] = BC.domain.corner
jj[q] = BC.domain.corner
s[q] .= 1.0
BCRHS[BC.domain.corner] .= 0.0

# assign values to the edges (useless cells)
q = q[end].+[1:length(BC.domain.edge);]
ii[q] = BC.domain.edge
jj[q] = BC.domain.edge
s[q] .= 1.0
BCRHS[BC.domain.edge] .= 0.0

# Assign values to the boundary condition matrix and the RHS vector based
# on the BC structure
if !BC.top.periodic && !BC.bottom.periodic
    # top boundary
    j=Ntheta+2
    i=2:Nx+1
    k=2:Nz+1
    q = q[end].+[1:Nx*Nz;]
    ii[q] = G[i,j,k]
    jj[q] = G[i,j,k]
    s[q] = BC.top.b/2.0 + BC.top.a./(dtheta_end*rp)
    q = q[end].+[1:Nx*Nz;]
    ii[q] = G[i,j,k]
    jj[q] = G[i,j-1,k]
    s[q] = BC.top.b/2.0 - BC.top.a./(dtheta_end*rp)
    BCRHS[G[i,j,k]] = BC.top.c

    # Bottom boundary
    j=1
    i=2:Nx+1
    k=2:Nz+1
    q = q[end].+[1:Nx*Nz;]
    ii[q] = G[i,j,k]
    jj[q] = G[i,j+1,k]
    s[q] = -(BC.bottom.b/2.0 + BC.bottom.a./(dtheta_1*rp)) # consider the reverse direction of normal
    q = q[end].+[1:Nx*Nz;]
    ii[q] = G[i,j,k]
    jj[q] = G[i,j,k]
    s[q] = -(BC.bottom.b/2.0 - BC.bottom.a./(dtheta_1*rp)) # consider the reverse direction of normal
    BCRHS[G[i,j,k]] = -(BC.bottom.c)
elseif BC.top.periodic || BC.bottom.periodic # periodic
    # top boundary
    j=Ntheta+2
    i=2:Nx+1
    k=2:Nz+1
    q = q[end].+[1:Nx*Nz;]
    ii[q] = G[i,j,k]
    jj[q] = G[i,j,k]
    s[q] = 1
    q = q[end].+[1:Nx*Nz;]
    ii[q] = G[i,j,k]
    jj[q] = G[i,j,k]
    s[q] = 1
    q = q[end].+[1:Nx*Nz;]
    ii[q] = G[i,j,k]
    jj[q] = G[i,2,k]
    s[q] = -1
    q = q[end].+[1:Nx*Nz;]
    ii[q] = G[i,j,k]
    jj[q] = G[i,2,k]
    s[q] = -1
    BCRHS[G[i,j,k]] .= 0.0

    # Bottom boundary
    j=1
    i=2:Nx+1
    k=2:Nz+1
    q = q[end].+[1:Nx*Nz;]
    ii[q] = G[i,j,k]
    jj[q] = G[i,j,k]
    s[q] .= 1.0
    q = q[end].+[1:Nx*Nz;]
    ii[q] = G[i,j,k]
    jj[q] = G[i,j,k]
    s[q] .= 1.0
    q = q[end].+[1:Nx*Nz;]
    ii[q] = G[i,j,k]
    jj[q] = G[i,Ntheta+1,k]
    s[q] .= -1.0
    q = q[end].+[1:Nx*Nz;]
    ii[q] = G[i,j,k]
    jj[q] = G[i,Ntheta+1,k]
    s[q] .= -1.0
    BCRHS[G[i,j,k]] .= 0.0
end

if !BC.right.periodic && !BC.left.periodic
    # Right boundary
    i=Nx+2
    j=2:Ntheta+1
    k=2:Nz+1
    q = q[end].+[1:Ntheta*Nz;]
    ii[q] = G[i,j,k]
    jj[q] = G[i,j,k]
    s[q] = BC.right.b/2.0 + BC.right.a/dx_end
    q = q[end].+[1:Ntheta*Nz;]
    ii[q] = G[i,j,k]
    jj[q] = G[i-1,j,k]
    s[q] = BC.right.b/2.0 - BC.right.a/dx_end
    BCRHS[G[i,j,k]] = BC.right.c

    # Left boundary
    i = 1
    j=2:Ntheta+1
    k=2:Nz+1
    q = q[end].+[1:Ntheta*Nz;]
    ii[q] = G[i,j,k]
    jj[q] = G[i+1,j,k]
    s[q] = -(BC.left.b/2.0 + BC.left.a/dx_1)
    # consider the reverse direction of normal
    q = q[end].+[1:Ntheta*Nz;]
    ii[q] = G[i,j,k]
    jj[q] = G[i,j,k]
    s[q] = -(BC.left.b/2.0 - BC.left.a/dx_1)  # consider the reverse direction of normal
    BCRHS[G[i,j,k]] = -(BC.left.c)
elseif BC.right.periodic || BC.left.periodic # periodic
    # Right boundary
    i=Nx+2
    j=2:Ntheta+1
    k=2:Nz+1
    q = q[end].+[1:Ntheta*Nz;]
    ii[q] = G[i,j,k]
    jj[q] = G[i,j,k]
    s[q] .= 1.0
    q = q[end].+[1:Ntheta*Nz;]
    ii[q] = G[i,j,k]
    jj[q] = G[i,j,k]
    s[q] .= 1.0
    q = q[end].+[1:Ntheta*Nz;]
    ii[q] = G[i,j,k]
    jj[q] = G[2,j,k]
    s[q] .= -1.0
    q = q[end].+[1:Ntheta*Nz;]
    ii[q] = G[i,j,k]
    jj[q] = G[2,j,k]
    s[q] .= -1.0
    BCRHS[G[i,j,k]] .= 0.0

    # Left boundary
    i = 1
    j=2:Ntheta+1
    k=2:Nz+1
    q = q[end].+[1:Ntheta*Nz;]
    ii[q] = G[i,j,k]
    jj[q] = G[i,j,k]
    s[q] .= 1.0
    q = q[end].+[1:Ntheta*Nz;]
    ii[q] = G[i,j,k]
    jj[q] = G[i,j,k]
    s[q] .= 1.0
    q = q[end].+[1:Ntheta*Nz;]
    ii[q] = G[i,j,k]
    jj[q] = G[Nx+1,j,k]
    s[q] .= -1.0
    q = q[end].+[1:Ntheta*Nz;]
    ii[q] = G[i,j,k]
    jj[q] = G[Nx+1,j,k]
    s[q] .= -1.0
    BCRHS[G[i,j,k]] .= 0.0
end

if !BC.front.periodic && !BC.back.periodic
    # Back boundary
    k=1
    i = 2:Nx+1
    j=2:Ntheta+1
    q = q[end].+[1:Nx*Ntheta;]
    ii[q] = G[i,j,k]
    jj[q] = G[i,j,k+1]
    s[q] = -(BC.back.b/2.0 + BC.back.a/dz_1)  # consider the reverse direction of normal
    q = q[end].+[1:Nx*Ntheta;]
    ii[q] = G[i,j,k]
    jj[q] = G[i,j,k]
    s[q] = -(BC.back.b/2.0 - BC.back.a/dz_1)  # consider the reverse direction of normal
    BCRHS[G[i,j,k]] = -(BC.back.c)

    # Front boundary
    k=Nz+2
    i = 2:Nx+1
    j=2:Ntheta+1
    q = q[end].+[1:Nx*Ntheta;]
    ii[q] = G[i,j,k]
    jj[q] = G[i,j,k]
    s[q] = BC.front.b/2.0 + BC.front.a/dz_end
    q = q[end].+[1:Nx*Ntheta;]
    ii[q] = G[i,j,k]
    jj[q] = G[i,j,k-1]
    s[q] = BC.front.b/2.0 - BC.front.a/dz_end
    BCRHS[G[i,j,k]] = BC.front.c
elseif BC.front.periodic || BC.back.periodic  # periodic
    # Back boundary
    k=1
    i = 2:Nx+1
    j=2:Ntheta+1
    q = q[end].+[1:Nx*Ntheta;]
    ii[q] = G[i,j,k]
    jj[q] = G[i,j,k]
    s[q] .= 1.0
    q = q[end].+[1:Nx*Ntheta;]
    ii[q] = G[i,j,k]
    jj[q] = G[i,j,k]
    s[q] .= 1.0
    q = q[end].+[1:Nx*Ntheta;]
    ii[q] = G[i,j,k]
    jj[q] = G[i,j,Nz+1]
    s[q] .= -1.0
    q = q[end].+[1:Nx*Ntheta;]
    ii[q] = G[i,j,k]
    jj[q] = G[i,j,Nz+1]
    s[q] .= -1.0
    BCRHS[G[i,j,k]] .= 0.0

    # Front boundary
    k=Nz+2
    i = 2:Nx+1
    j=2:Ntheta+1
    q = q[end].+[1:Nx*Ntheta;]
    ii[q] = G[i,j,k]
    jj[q] = G[i,j,k]
    s[q] .= 1.0
    q = q[end].+[1:Nx*Ntheta;]
    ii[q] = G[i,j,k]
    jj[q] = G[i,j,k]
    s[q] .= 1.0
    q = q[end].+[1:Nx*Ntheta;]
    ii[q] = G[i,j,k]
    jj[q] = G[i,j,2]
    s[q] .= -1.0
    q = q[end].+[1:Nx*Ntheta;]
    ii[q] = G[i,j,k]
    jj[q] = G[i,j,2]
    s[q] .= -1.0
    BCRHS[G[i,j,k]] .= 0.0
end

# Build the sparse matrix of the boundary conditions
BCMatrix = sparse(ii[1:q[end]], jj[1:q[end]], s[1:q[end]],
    (Nx+2)*(Ntheta+2)*(Nz+2), (Nx+2)*(Ntheta+2)*(Nz+2))

(BCMatrix, BCRHS)
end


# =============================== Cell Boundary ====================================
# Assign values to ghost cells
function cellBoundary!(phi::CellValue, BC::BoundaryCondition)
d=phi.domain.dimension
if d ==1 || d==1.5
  cellBoundary1D!(phi, BC);
elseif d==2 || d==2.5
  cellBoundary2D!(phi, BC)
elseif d==2.8
  cellBoundaryRadial2D!(phi, BC)
elseif d==3
  cellBoundary3D!(phi, BC)
elseif d==3.2
  cellBoundaryCylindrical3D!(phi, BC)
end
end

# =========================== 1D Cartesian ================================
function cellBoundary1D!(phi::CellValue, BC::BoundaryCondition)
dx_1 = phi.domain.cellsize.x[1]
dx_end = phi.domain.cellsize.x[end]
# boundary condition (a d\phi/dx + b \phi = c, a column vector of [d a])
# a (phi(i)-phi(i-1))/dx + b (phi(i)+phi(i-1))/2 = c
# phi(i) (a/dx+b/2) + phi(i-1) (-a/dx+b/2) = c
# Right boundary, i=m+2
# phi(i) (a/dx+b/2) = c- phi(i-1) (-a/dx+b/2)
# Left boundary, i=2
#  phi(i-1) (-a/dx+b/2) = c - phi(i) (a/dx+b/2)
# define the new phi
if !BC.left.periodic && !BC.right.periodic
    phiBC = [(BC.left.c[1]-phi.value[2]*(BC.left.a[1]/dx_1+BC.left.b[1]/2.0))/(-BC.left.a[1]/dx_1+BC.left.b[1]/2.0);
             phi.value[2:end-1];
             (BC.right.c[1]-phi.value[end-1]*(-BC.right.a[1]/dx_end+BC.right.b[1]/2.0))/(BC.right.a[1]/dx_end+BC.right.b[1]/2.0)]
else
    phiBC = [phi.value[end]; phi.value[2:end-1]; phi.value[1]]
end
phi.value[:] = phiBC[:]
phi
end

# =========================== 2D Cartesian ================================
function cellBoundary2D!(phi::CellValue, BC::BoundaryCondition)
# extract data from the mesh structure
Nx, Ny = tuple(BC.domain.dims...)
G=reshape([1:(Nx+2)*(Ny+2);], Nx+2, Ny+2)
dx_1 = BC.domain.cellsize.x[1]
dx_end = BC.domain.cellsize.x[end]
dy_1 = BC.domain.cellsize.y[1]
dy_end = BC.domain.cellsize.y[end]

# define the output matrix
phiBC = zeros(Nx+2, Ny+2)
phiBC[:] = phi.value[:]

# Assign values to the boundary values
if !BC.top.periodic && !BC.bottom.periodic
    # top boundary
    j=Ny+2
    for i = 2:Nx+1
      phiBC[i,j]= (BC.top.c[i-1]-phi.value[i,end-1].*(-BC.top.a[i-1]/dy_end+BC.top.b[i-1]/2.0))/(BC.top.a[i-1]/dy_end+BC.top.b[i-1]/2.0)
    end
    # Bottom boundary
    j=1;
    for i = 2:Nx+1
      phiBC[i,j]= (BC.bottom.c[i-1]-phi.value[i,2].*(BC.bottom.a[i-1]/dy_1+BC.bottom.b[i-1]/2.0))/(-BC.bottom.a[i-1]/dy_1+BC.bottom.b[i-1]/2.0)
    end
else
    # top boundary
    j=Ny+2
    for i = 2:Nx+1
      phiBC[i,j]= phi.value[i,2]
    end
    # Bottom boundary
    j=1
    for i = 2:Nx+1
      phiBC[i,j]= phi.value[i,end-1]
    end
end

if !BC.left.periodic && !BC.right.periodic
    # Right boundary
    i = Nx+2
    for j = 2:Ny+1
      phiBC[i,j]= (BC.right.c[j-1]-phi.value[end-1,j]*(-BC.right.a[j-1]/dx_end+BC.right.b[j-1]/2.0))/(BC.right.a[j-1]/dx_end+BC.right.b[j-1]/2.0)
    end

    # Left boundary
    i = 1
    for j = 2:Ny+1
      phiBC[i,j]= (BC.left.c[j-1]-phi.value[2,j]*(BC.left.a[j-1]/dx_1+BC.left.b[j-1]/2.0))/(-BC.left.a[j-1]/dx_1+BC.left.b[j-1]/2)
    end
else
    # Right boundary
    i = Nx+2
    for j = 2:Ny+1
      phiBC[i,j]= phi.value[2,j]
    end
    # Left boundary
    i = 1
    for j = 2:Ny+1
      phiBC[i,j]= phi.value[end-1,j]
    end
end
phi.value[:] = phiBC[:]
phi
end

# ======================================== Cell Boundary 3D ===================================
function cellBoundary3D!(phi::CellValue, BC::BoundaryCondition)
# extract data from the mesh structure
Nx, Ny, Nz = tuple(BC.domain.dims...)
G=reshape([1:(Nx+2)*(Ny+2)*(Nz+2);], Nx+2, Ny+2, Nz+2)
dx_1 = BC.domain.cellsize.x[1]
dx_end = BC.domain.cellsize.x[end]
dy_1 = BC.domain.cellsize.y[1]
dy_end = BC.domain.cellsize.y[end]
dz_1 = BC.domain.cellsize.z[1]
dz_end = BC.domain.cellsize.z[end]

# define the output matrix
phiBC = zeros(Nx+2, Ny+2, Nz+2)
phiBC[:] = phi.value[:]

# Assign values to the boundary values
if !BC.top.periodic && !BC.bottom.periodic
    # top boundary
    j=Ny+2
    i = 2:Nx+1
    k = 2:Nz+1
    phiBC[i,j,k]= (BC.top.c-(phi.value[i,end-1:end-1,k]).*(-BC.top.a/dy_end+BC.top.b/2.0))./(BC.top.a/dy_end+BC.top.b/2.0)

    # Bottom boundary
    j=1
    i = 2:Nx+1
    k = 2:Nz+1
    phiBC[i,j,k]= (BC.bottom.c-(phi.value[i,2:2,k]).*(BC.bottom.a/dy_1+BC.bottom.b/2.0))./(-BC.bottom.a/dy_1+BC.bottom.b/2.0)
else
    # top boundary
    j=Ny+2
    i = 2:Nx+1
    k = 2:Nz+1
    phiBC[i,j,k]= phi.value[i,2,k]

    # Bottom boundary
    j=1
    i = 2:Nx+1
    k = 2:Nz+1
    phiBC[i,j,k]= phi.value[i,end-1,k]
end

if !BC.left.periodic && !BC.right.periodic
    # Right boundary
    i = Nx+2
    j = 2:Ny+1
    k = 2:Nz+1
    phiBC[i,j,k]= (BC.right.c-(phi.value[end-1:end-1,j,k]).*(-BC.right.a/dx_end+BC.right.b/2.0))./(BC.right.a/dx_end+BC.right.b/2.0)

    # Left boundary
    i = 1;
    j = 2:Ny+1
    k = 2:Nz+1
    phiBC[i,j,k]= (BC.left.c-(phi.value[2:2,j,k]).*(BC.left.a/dx_1+BC.left.b/2.0))./(-BC.left.a/dx_1+BC.left.b/2.0)
else
    # Right boundary
    i = Nx+2
    j = 2:Ny+1
    k = 2:Nz+1
    phiBC[i,j,k]= phi.value[2,j,k]

    # Left boundary
    i = 1
    j = 2:Ny+1
    k = 2:Nz+1
    phiBC[i,j,k]= phi.value[end-1,j,k]
end

if !BC.bottom.periodic && !BC.top.periodic
    # front boundary
    i = 2:Nx+1
    j = 2:Ny+1
    k = Nz+2
    phiBC[i,j,k]= (BC.front.c-(phi.value[i,j,end-1:end-1]).*(-BC.front.a/dz_end+BC.front.b/2.0))./(BC.front.a/dz_end+BC.front.b/2.0)

    # back boundary
    i = 2:Nx+1
    j = 2:Ny+1
    k = 1
    phiBC[i,j,k]= (BC.back.c-(phi.value[i,j,2:2]).*(BC.back.a/dz_1+BC.back.b/2.0))./(-BC.back.a/dz_1+BC.back.b/2.0)
else
    # front boundary
    i = 2:Nx+1
    j = 2:Ny+1
    k = Nz+2
    phiBC[i,j,k]= phi.value[i,j,2]

    # back boundary
    i = 2:Nx+1
    j = 2:Ny+1
    k = 1
    phiBC[i,j,k]= phi.value[j,j,end-1]
end
phi.value[:] = phiBC[:]
phi
end


# ===================================== Cell Boundary Radial 2D ===========================
function cellBoundaryRadial2D!(phi::CellValue, BC::BoundaryCondition)
# extract data from the mesh structure
Nr, Ntheta = tuple(BC.domain.dims...)
G=reshape([1:(Nr+2)*(Ntheta+2);], Nr+2, Ntheta+2)
dr_1 = BC.domain.cellsize.x[1]
dr_end = BC.domain.cellsize.x[end]
dtheta_1 = BC.domain.cellsize.y[1]
dtheta_end = BC.domain.cellsize.y[end]
rp = BC.domain.cellcenters.x

# define the output matrix
phiBC = zeros(Nr+2, Ntheta+2)
phiBC[:] = phi.value[:]

# Assign values to the boundary values
if !BC.top.periodic && !BC.bottom.periodic
    # top boundary
    j=Ntheta+2
    for i = 2:Nr+1
      phiBC[i,j]= (BC.top.c[i-1]-phi.value[i,end-1].*(-BC.top.a[i-1]/(dtheta_end*rp[i-1])+BC.top.b[i-1]/2.0))/(BC.top.a[i-1]/(dtheta_end*rp[i-1])+BC.top.b[i-1]/2.0)
    end
    # Bottom boundary
    j=1;
    for i = 2:Nr+1
      phiBC[i,j]= (BC.bottom.c[i-1]-phi.value[i,2].*(BC.bottom.a[i-1]/(dtheta_1*rp[i-1])+BC.bottom.b[i-1]/2.0))/(-BC.bottom.a[i-1]/(dtheta_1*rp[i-1])+BC.bottom.b[i-1]/2.0)
    end
else
    # top boundary
    j=Ntheta+2
    for i = 2:Nr+1
      phiBC[i,j]= phi.value[i,2]
    end
    # Bottom boundary
    j=1
    for i = 2:Nr+1
      phiBC[i,j]= phi.value[i,end-1]
    end
end

if !BC.left.periodic && !BC.right.periodic
    # Right boundary
    i = Nr+2
    for j = 2:Ntheta+1
      phiBC[i,j]= (BC.right.c[j-1]-phi.value[end-1,j]*(-BC.right.a[j-1]/dr_end+BC.right.b[j-1]/2.0))/(BC.right.a[j-1]/dr_end+BC.right.b[j-1]/2.0)
    end

    # Left boundary
    i = 1
    for j = 2:Ntheta+1
      phiBC[i,j]= (BC.left.c[j-1]-phi.value[2,j]*(BC.left.a[j-1]/dr_1+BC.left.b[j-1]/2.0))/(-BC.left.a[j-1]/dr_1+BC.left.b[j-1]/2.0)
    end
else
    # Right boundary
    i = Nr+2
    for j = 2:Ntheta+1
      phiBC[i,j]= phi.value[2,j]
    end
    # Left boundary
    i = 1
    for j = 2:Ntheta+1
      phiBC[i,j]= phi.value[end-1,j]
    end
end
phi.value[:] = phiBC[:]
phi
end


# ===================================== Cell Boundary Cylindrical 3D ===========================
function cellBoundaryCylindrical3D!(phi::CellValue, BC::BoundaryCondition)
# extract data from the mesh structure
Nr, Ntheta, Nz = tuple(BC.domain.dims...)
G=reshape([1:(Nr+2)*(Ntheta+2)*(Nz+2);], Nr+2, Ntheta+2, Nz+2)
dr_1 = BC.domain.cellsize.x[1]
dr_end = BC.domain.cellsize.x[end]
dtheta_1 = BC.domain.cellsize.y[1]
dtheta_end = BC.domain.cellsize.y[end]
dz_1 = BC.domain.cellsize.z[1]
dz_end = BC.domain.cellsize.z[end]
#rp = zeros(Nx,1,Nz)
#rp[:,1,:] = repmat(phi.domain.cellcenters.x, 1, Nz)
rp = phi.domain.cellcenters.x

# define the output matrix
phiBC = zeros(Nr+2, Ntheta+2, Nz+2)
phiBC[:] = phi.value[:]

# Assign values to the boundary values
if !BC.top.periodic && !BC.bottom.periodic
    # top boundary
    j=Ntheta+2
    i = 2:Nr+1
    k = 2:Nz+1
    phiBC[i,j,k]= (BC.top.c-(phi.value[i,end-1:end-1,k]).*(-BC.top.a./(dtheta_end*rp)+BC.top.b/2.0))./(BC.top.a./(dtheta_end*rp)+BC.top.b/2.0)

    # Bottom boundary
    j=1
    i = 2:Nr+1
    k = 2:Nz+1
    phiBC[i,j,k]= (BC.bottom.c-(phi.value[i,2:2,k]).*(BC.bottom.a./(dtheta_1*rp)+BC.bottom.b/2.0))./(-BC.bottom.a./(dtheta_1*rp)+BC.bottom.b/2.0)
else
    # top boundary
    j=Ntheta+2
    i = 2:Nr+1
    k = 2:Nz+1
    phiBC[i,j,k]= phi.value[i,2,k]

    # Bottom boundary
    j=1
    i = 2:Nr+1
    k = 2:Nz+1
    phiBC[i,j,k]= phi.value[i,end-1,k]
end

if !BC.left.periodic && !BC.right.periodic
    # Right boundary
    i = Nr+2
    j = 2:Ntheta+1
    k = 2:Nz+1
    phiBC[i,j,k]= (BC.right.c-(phi.value[end-1:end-1,j,k]).*(-BC.right.a/dr_end+BC.right.b/2.0))./(BC.right.a/dr_end+BC.right.b/2.0)

    # Left boundary
    i = 1;
    j = 2:Ntheta+1
    k = 2:Nz+1
    phiBC[i,j,k]= (BC.left.c-(phi.value[2:2,j,k]).*(BC.left.a/dr_1+BC.left.b/2.0))./(-BC.left.a/dr_1+BC.left.b/2.0)
else
    # Right boundary
    i = Nr+2
    j = 2:Ntheta+1
    k = 2:Nz+1
    phiBC[i,j,k]= phi.value[2,j,k]

    # Left boundary
    i = 1
    j = 2:Ntheta+1
    k = 2:Nz+1
    phiBC[i,j,k]= phi.value[end-1,j,k]
end

if !BC.bottom.periodic && !BC.top.periodic
    # front boundary
    i = 2:Nr+1
    j = 2:Ntheta+1
    k = Nz+2
    phiBC[i,j,k]= (BC.front.c-(phi.value[i,j,end-1:end-1]).*(-BC.front.a/dz_end+BC.front.b/2.0))./(BC.front.a/dz_end+BC.front.b/2.0)

    # back boundary
    i = 2:Nr+1
    j = 2:Ntheta+1
    k = 1
    phiBC[i,j,k]= (BC.back.c-(phi.value[i,j,2:2]).*(BC.back.a/dz_1+BC.back.b/2.0))./(-BC.back.a/dz_1+BC.back.b/2.0)
else
    # front boundary
    i = 2:Nr+1
    j = 2:Ntheta+1
    k = Nz+2
    phiBC[i,j,k]= phi.value[i,j,2]

    # back boundary
    i = 2:Nr+1
    j = 2:Ntheta+1
    k = 1
    phiBC[i,j,k]= phi.value[j,j,end-1]
end
phi.value[:] = phiBC[:]
phi
end
