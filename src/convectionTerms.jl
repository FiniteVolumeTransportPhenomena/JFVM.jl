# ===============================
# Written by AAE
# TU Delft, Spring 2014
# simulkade.com
# edited on multiple occasions at
# Technical University of Denmark
# ===============================


function convectionTerm(u::FaceValue)
  convectionTerm(u, u.domain.coordinatesystem)
end

convectionTerm(u::FaceValue, ::Cartesian1D) = convectionTerm1D(u)
convectionTerm(u::FaceValue, ::Cylindrical1D) = convectionTermCylindrical1D(u)
convectionTerm(u::FaceValue, ::Cartesian2D) = first(convectionTerm2D(u))
convectionTerm(u::FaceValue, ::Cylindrical2D) = first(convectionTermCylindrical2D(u))
convectionTerm(u::FaceValue, ::Radial2D) = first(convectionTermRadial2D(u))
convectionTerm(u::FaceValue, ::Cartesian3D) = first(convectionTerm3D(u))
convectionTerm(u::FaceValue, ::Cylindrical3D) = first(convectionTermCylindrical3D(u))

function convectionUpwindTerm(u::FaceValue)
  convectionUpwindTerm(u, u.domain.coordinatesystem)
end

convectionUpwindTerm(u::FaceValue, ::Cartesian1D) = convectionUpwindTerm1D(u)
convectionUpwindTerm(u::FaceValue, ::Cylindrical1D) = convectionUpwindTermCylindrical1D(u)
convectionUpwindTerm(u::FaceValue, ::Cartesian2D) = first(convectionUpwindTerm2D(u))
convectionUpwindTerm(u::FaceValue, ::Cylindrical2D) = first(convectionUpwindTermCylindrical2D(u))
convectionUpwindTerm(u::FaceValue, ::Radial2D) = first(convectionUpwindTermRadial2D(u))
convectionUpwindTerm(u::FaceValue, ::Cartesian3D) = first(convectionUpwindTerm3D(u))
convectionUpwindTerm(u::FaceValue, ::Cylindrical3D) = first(convectionUpwindTermCylindrical3D(u))

function convectionUpwindTerm(u::FaceValue, u_upwind::FaceValue)
  convectionUpwindTerm(u, u_upwind, u.domain.coordinatesystem)
end

convectionUpwindTerm(u::FaceValue, u_upwind::FaceValue, ::Cartesian1D) = convectionUpwindTerm1D(u, u_upwind)
convectionUpwindTerm(u::FaceValue, u_upwind::FaceValue, ::Cylindrical1D) = convectionUpwindTermCylindrical1D(u, u_upwind)
convectionUpwindTerm(u::FaceValue, u_upwind::FaceValue, ::Cartesian2D) = first(convectionUpwindTerm2D(u, u_upwind))
convectionUpwindTerm(u::FaceValue, u_upwind::FaceValue, ::Cylindrical2D) = first(convectionUpwindTermCylindrical2D(u, u_upwind))
convectionUpwindTerm(u::FaceValue, u_upwind::FaceValue, ::Radial2D) = first(convectionUpwindTermRadial2D(u, u_upwind))
convectionUpwindTerm(u::FaceValue, u_upwind::FaceValue, ::Cartesian3D) = first(convectionUpwindTerm3D(u, u_upwind))
convectionUpwindTerm(u::FaceValue, u_upwind::FaceValue, ::Cylindrical3D) = first(convectionUpwindTermCylindrical3D(u, u_upwind))

function convectionTvdTerm(u::FaceValue, phi::CellValue, FL::Function)
  convectionTvdTerm(u, phi, FL, u.domain.coordinatesystem)
end

function convectionTvdTerm(u::FaceValue, phi::CellValue, FL::Function, ::Cartesian1D)
  M, RHS = convectionTvdTerm1D(u::FaceValue, phi::CellValue, FL::Function)
  return (M, RHS)
end

function convectionTvdTerm(u::FaceValue, phi::CellValue, FL::Function, ::Cylindrical1D)
  M, RHS = convectionTvdTermCylindrical1D(u::FaceValue, phi::CellValue, FL::Function)
  return (M, RHS)
end

function convectionTvdTerm(u::FaceValue, phi::CellValue, FL::Function, ::Cartesian2D)
  M, RHS, Mx, My, RHSx, RHSy = convectionTvdTerm2D(u::FaceValue, phi::CellValue, FL::Function)
  return (M, RHS)
end

function convectionTvdTerm(u::FaceValue, phi::CellValue, FL::Function, ::Cylindrical2D)
  M, RHS, Mx, My, RHSx, RHSy = convectionTvdTermCylindrical2D(u::FaceValue, phi::CellValue, FL::Function)
  return (M, RHS)
end

function convectionTvdTerm(u::FaceValue, phi::CellValue, FL::Function, ::Radial2D)
  M, RHS, Mx, My, RHSx, RHSy = convectionTvdTermRadial2D(u::FaceValue, phi::CellValue, FL::Function)
  return (M, RHS)
end

function convectionTvdTerm(u::FaceValue, phi::CellValue, FL::Function, ::Cartesian3D)
  M, RHS, Mx, My, Mz, RHSx, RHSy, RHSz = convectionTvdTerm3D(u::FaceValue, phi::CellValue, FL::Function)
  return (M, RHS)
end

function convectionTvdTerm(u::FaceValue, phi::CellValue, FL::Function, ::Cylindrical3D)
  M, RHS, Mx, My, Mz, RHSx, RHSy, RHSz = convectionTvdTermCylindrical3D(u::FaceValue, phi::CellValue, FL::Function)
  return (M, RHS)
end

function convectionTvdTerm(u::FaceValue, phi::CellValue, FL::Function, u_upwind)
    M   = convectionUpwindTerm(u, u_upwind)
    RHS = convectionTvdRHS(u, phi, FL, u_upwind)
    (M, RHS)
end

function convectionTvdRHS(u::FaceValue, phi::CellValue, FL::Function)
  convectionTvdRHS(u, phi, FL, u.domain.coordinatesystem)
end

convectionTvdRHS(u::FaceValue, phi::CellValue, FL::Function, ::Cartesian1D) = convectionTvdRHS1D(u::FaceValue, phi::CellValue, FL::Function)
convectionTvdRHS(u::FaceValue, phi::CellValue, FL::Function, ::Cylindrical1D) = convectionTvdRHSCylindrical1D(u::FaceValue, phi::CellValue, FL::Function)
convectionTvdRHS(u::FaceValue, phi::CellValue, FL::Function, ::Cartesian2D) = first(convectionTvdRHS2D(u::FaceValue, phi::CellValue, FL::Function))
convectionTvdRHS(u::FaceValue, phi::CellValue, FL::Function, ::Cylindrical2D) = first(convectionTvdRHSCylindrical2D(u::FaceValue, phi::CellValue, FL::Function))
convectionTvdRHS(u::FaceValue, phi::CellValue, FL::Function, ::Radial2D) = first(convectionTvdRHSRadial2D(u::FaceValue, phi::CellValue, FL::Function))
convectionTvdRHS(u::FaceValue, phi::CellValue, FL::Function, ::Cartesian3D) = first(convectionTvdRHS3D(u::FaceValue, phi::CellValue, FL::Function))
convectionTvdRHS(u::FaceValue, phi::CellValue, FL::Function, ::Cylindrical3D) = first(convectionTvdRHSCylindrical3D(u::FaceValue, phi::CellValue, FL::Function))

function convectionTvdRHS(u::FaceValue, phi::CellValue, FL::Function, u_upwind::FaceValue)
  convectionTvdRHS(u, phi, FL, u_upwind, u.domain.coordinatesystem)
end

convectionTvdRHS(u::FaceValue, phi::CellValue, FL::Function, u_upwind::FaceValue, ::Cartesian1D) = convectionTvdRHS1D(u::FaceValue, phi::CellValue, FL::Function, u_upwind::FaceValue)
convectionTvdRHS(u::FaceValue, phi::CellValue, FL::Function, u_upwind::FaceValue, ::Cylindrical1D) = convectionTvdRHSCylindrical1D(u::FaceValue, phi::CellValue, FL::Function, u_upwind::FaceValue)
convectionTvdRHS(u::FaceValue, phi::CellValue, FL::Function, u_upwind::FaceValue, ::Cartesian2D) = first(convectionTvdRHS2D(u::FaceValue, phi::CellValue, FL::Function, u_upwind::FaceValue))
convectionTvdRHS(u::FaceValue, phi::CellValue, FL::Function, u_upwind::FaceValue, ::Cylindrical2D) = first(convectionTvdRHSCylindrical2D(u::FaceValue, phi::CellValue, FL::Function, u_upwind::FaceValue))
convectionTvdRHS(u::FaceValue, phi::CellValue, FL::Function, u_upwind::FaceValue, ::Radial2D) = first(convectionTvdRHSRadial2D(u::FaceValue, phi::CellValue, FL::Function, u_upwind::FaceValue))
convectionTvdRHS(u::FaceValue, phi::CellValue, FL::Function, u_upwind::FaceValue, ::Cartesian3D) = first(convectionTvdRHS3D(u::FaceValue, phi::CellValue, FL::Function, u_upwind::FaceValue))
convectionTvdRHS(u::FaceValue, phi::CellValue, FL::Function, u_upwind::FaceValue, ::Cylindrical3D) = first(convectionTvdRHSCylindrical3D(u::FaceValue, phi::CellValue, FL::Function, u_upwind::FaceValue))

# =================== 1D Convection Terms =======================
# =============== 1D Convection Terms Cartesian Central =================
function convectionTerm1D(u::FaceValue)
# u is a face variable

# extract data from the mesh structure
Nx = u.domain.dims[1]
G = [1:Nx+2;]
#DX = u.domain.cellsize.x

DXe = u.domain.cellsize.x[3:end]
DXw = u.domain.cellsize.x[1:end-2]
DXp = u.domain.cellsize.x[2:end-1]

# define the vectors to store the sparse matrix data
iix = zeros(Int64, 3*(Nx+2))
jjx = zeros(Int64, 3*(Nx+2))
sx = zeros(Float64, 3*(Nx+2))

# reassign the east, west for code readability
ue = u.xvalue[2:Nx+1]./(DXp.+DXe)
uw = u.xvalue[1:Nx]./(DXp.+DXw)

# calculate the coefficients for the internal cells
AE = reshape(ue,Nx)
AW = reshape(-uw,Nx)
APx = reshape((ue.*DXe.-uw.*DXw)./DXp,Nx)

# build the sparse matrix based on the numbering system
rowx_index = reshape(G[2:Nx+1],Nx) # main diagonal x
iix[1:3*Nx] .= repeat(rowx_index,3)
jjx[1:3*Nx] .= [reshape(G[1:Nx],Nx); reshape(G[2:Nx+1],Nx); reshape(G[3:Nx+2],Nx)]
sx[1:3*Nx] .= [AW; APx; AE]

# build the sparse matrix
kx = 3*Nx
M = sparse(iix[1:kx], jjx[1:kx], sx[1:kx], Nx+2, Nx+2)

end


# =============== 1D Convection Terms Cylindrical Central =================
function convectionTermCylindrical1D(u::FaceValue)
# u is a face variable

# extract data from the mesh structure
Nx = u.domain.dims[1]
G = [1:Nx+2;]
DXe = u.domain.cellsize.x[3:end]
DXw = u.domain.cellsize.x[1:end-2]
DXp = u.domain.cellsize.x[2:end-1]
rp = u.domain.cellcenters.x
rf = u.domain.facecenters.x

# define the vectors to store the sparse matrix data
iix = zeros(Int64, 3*(Nx+2))
jjx = zeros(Int64, 3*(Nx+2))
sx = zeros(Float64, 3*(Nx+2))

# reassign the east, west for code readability
ue = rf[2:Nx+1].*u.xvalue[2:Nx+1]./(rp.*(DXp.+DXe))
uw = rf[1:Nx].*u.xvalue[1:Nx]./(rp.*(DXp.+DXw))

# calculate the coefficients for the internal cells
AE = reshape(ue,Nx)
AW = reshape(-uw,Nx)
APx = reshape((ue.*DXe-uw.*DXw)./DXp,Nx)

# build the sparse matrix based on the numbering system
rowx_index = reshape(G[2:Nx+1],Nx) # main diagonal x
iix[1:3*Nx] .= repeat(rowx_index,3)
jjx[1:3*Nx] .= [reshape(G[1:Nx],Nx); reshape(G[2:Nx+1],Nx); reshape(G[3:Nx+2],Nx)]
sx[1:3*Nx] .= [AW; APx; AE]

# build the sparse matrix
kx = 3*Nx
M = sparse(iix[1:kx], jjx[1:kx], sx[1:kx], Nx+2, Nx+2)

end


# =============== 1D Convection Term Cartesian Upwind =================
function convectionUpwindTerm1D(u::FaceValue)
# u is a face variable

# extract data from the mesh structure
Nx = u.domain.dims[1]
G = [1:Nx+2;]
DXp = u.domain.cellsize.x[2:end-1]

# define the vectors to store the sparse matrix data
iix = zeros(Int64, 3*(Nx+2))
jjx = zeros(Int64, 3*(Nx+2))
sx = zeros(Float64, 3*(Nx+2))

# find the velocity direction for the upwind scheme
ue_min = min.(u.xvalue[2:Nx+1],0.0)
ue_max = max.(u.xvalue[2:Nx+1],0.0)
uw_min = min.(u.xvalue[1:Nx],0.0)
uw_max = max.(u.xvalue[1:Nx],0.0)

# calculate the coefficients for the internal cells
AE = reshape(ue_min./DXp,Nx)
AW = reshape(-uw_max./DXp,Nx)
APx = reshape((ue_max.-uw_min)./DXp,Nx)

# correct for the cells next to the boundary
# Left boundary:
APx[1] = APx[1].-uw_max[1]/(2.0*DXp[1])
AW[1] = AW[1]/2.0
# Right boundary:
AE[end] = AE[end]/2.0
APx[end] = APx[end] + ue_min[end]/(2.0*DXp[end])

# build the sparse matrix based on the numbering system
rowx_index = reshape(G[2:Nx+1],Nx) # main diagonal x
iix[1:3*Nx] .= repeat(rowx_index,3)
jjx[1:3*Nx] .= [reshape(G[1:Nx],Nx); reshape(G[2:Nx+1],Nx); reshape(G[3:Nx+2],Nx)]
sx[1:3*Nx] .= [AW; APx; AE]

# build the sparse matrix
kx = 3*Nx
M = sparse(iix[1:kx], jjx[1:kx], sx[1:kx], Nx+2, Nx+2)

end

"""
convectionUpwindTerm1D(u::FaceValue, u_upwind::FaceValue)
To be used for the discretization of the nonlinear equations.
"""
function convectionUpwindTerm1D(u::FaceValue, u_upwind::FaceValue)
# u is a face variable

# extract data from the mesh structure
Nx = u.domain.dims[1]
G = [1:Nx+2;]
DXp = u.domain.cellsize.x[2:end-1]

# define the vectors to store the sparse matrix data
iix = zeros(Int64, 3*(Nx+2))
jjx = zeros(Int64, 3*(Nx+2))
sx = zeros(Float64, 3*(Nx+2))

# find the velocity direction for the upwind scheme
ue_min = u.xvalue[2:Nx+1]
ue_max = u.xvalue[2:Nx+1]
uw_min = u.xvalue[1:Nx]
uw_max = u.xvalue[1:Nx]

ue_min[u_upwind.xvalue[2:Nx+1].>0.0] .= 0.0
ue_max[u_upwind.xvalue[2:Nx+1].<0.0] .= 0.0
uw_min[u_upwind.xvalue[1:Nx].>0.0] .= 0.0
uw_max[u_upwind.xvalue[1:Nx].<0.0] .= 0.0

# calculate the coefficients for the internal cells
AE = reshape(ue_min./DXp,Nx)
AW = reshape(.-uw_max./DXp,Nx)
APx = reshape((ue_max.-uw_min)./DXp,Nx)

# correct for the cells next to the boundary
# Left boundary:
APx[1] = APx[1].-uw_max[1]/(2.0*DXp[1])
AW[1] = AW[1]/2.0
# Right boundary:
AE[end] = AE[end]/2.0
APx[end] = APx[end] + ue_min[end]/(2.0*DXp[end])

# build the sparse matrix based on the numbering system
rowx_index = reshape(G[2:Nx+1],Nx) # main diagonal x
iix[1:3*Nx] .= repeat(rowx_index,3)
jjx[1:3*Nx] .= [reshape(G[1:Nx],Nx); reshape(G[2:Nx+1],Nx); reshape(G[3:Nx+2],Nx)]
sx[1:3*Nx] .= [AW; APx; AE]

# build the sparse matrix
kx = 3*Nx
M = sparse(iix[1:kx], jjx[1:kx], sx[1:kx], Nx+2, Nx+2)

end

# =================== 1D Convection Terms Cylindrical Upwind =====================
function convectionUpwindTermCylindrical1D(u::FaceValue)
# u is a face variable

# extract data from the mesh structure
Nx = u.domain.dims[1]
G = [1:Nx+2;]
DXp = u.domain.cellsize.x[2:end-1]
rp = u.domain.cellcenters.x
rf = u.domain.facecenters.x

# define the vectors to store the sparse matrix data
iix = zeros(Int64, 3*(Nx+2))
jjx = zeros(Int64, 3*(Nx+2))
sx = zeros(Float64, 3*(Nx+2))

# reassign the east, west for code readability
re = rf[2:Nx+1]
rw = rf[1:Nx]

# find the velocity direction for the upwind scheme
ue_min = min.(u.xvalue[2:Nx+1],0.0)
ue_max = max.(u.xvalue[2:Nx+1],0.0)
uw_min = min.(u.xvalue[1:Nx],0.0)
uw_max = max.(u.xvalue[1:Nx],0.0)

# calculate the coefficients for the internal cells
AE = reshape(re.*ue_min./(DXp.*rp),Nx)
AW = reshape(-rw.*uw_max./(DXp.*rp),Nx)
APx = reshape((re.*ue_max.-rw.*uw_min)./(DXp.*rp),Nx)

# correct for the cells next to the boundary
# Left boundary:
APx[1] = APx[1].-rw[1]*uw_max[1]/(2.0*DXp[1]*rp[1])
AW[1] = AW[1]/2.0
# Right boundary:
AE[end] = AE[end]/2.0
APx[end] = APx[end] + re[end]*ue_min[end]/(2.0*DXp[end]*rp[end])

# build the sparse matrix based on the numbering system
rowx_index = reshape(G[2:Nx+1],Nx) # main diagonal x
iix[1:3*Nx] .= repeat(rowx_index,3)
jjx[1:3*Nx] .= [reshape(G[1:Nx],Nx); reshape(G[2:Nx+1],Nx); reshape(G[3:Nx+2],Nx)]
sx[1:3*Nx] .= [AW; APx; AE]

# build the sparse matrix
kx = 3*Nx
M = sparse(iix[1:kx], jjx[1:kx], sx[1:kx], Nx+2, Nx+2)

end

"""
function convectionUpwindTermCylindrical1D(u::FaceValue, u_upwind::FaceValue)
for the linearization of the nonlinear terms
"""
function convectionUpwindTermCylindrical1D(u::FaceValue, u_upwind::FaceValue)
# u is a face variable

# extract data from the mesh structure
Nx = u.domain.dims[1]
G = [1:Nx+2;]
DXp = u.domain.cellsize.x[2:end-1]
rp = u.domain.cellcenters.x
rf = u.domain.facecenters.x

# define the vectors to store the sparse matrix data
iix = zeros(Int64, 3*(Nx+2))
jjx = zeros(Int64, 3*(Nx+2))
sx = zeros(Float64, 3*(Nx+2))

# reassign the east, west for code readability
re = rf[2:Nx+1]
rw = rf[1:Nx]

# find the velocity direction for the upwind scheme
ue_min = u.xvalue[2:Nx+1]
ue_max = u.xvalue[2:Nx+1]
uw_min = u.xvalue[1:Nx]
uw_max = u.xvalue[1:Nx]

ue_min[u_upwind.xvalue[2:Nx+1].>0.0] .= 0.0
ue_max[u_upwind.xvalue[2:Nx+1].<0.0] .= 0.0
uw_min[u_upwind.xvalue[1:Nx].>0.0] .= 0.0
uw_max[u_upwind.xvalue[1:Nx].<0.0] .= 0.0


# calculate the coefficients for the internal cells
AE = reshape(re.*ue_min./(DXp.*rp),Nx)
AW = reshape(.-rw.*uw_max./(DXp.*rp),Nx)
APx = reshape((re.*ue_max.-rw.*uw_min)./(DXp.*rp),Nx)

# correct for the cells next to the boundary
# Left boundary:
APx[1] = APx[1].-rw[1]*uw_max[1]/(2.0*DXp[1]*rp[1])
AW[1] = AW[1]/2.0
# Right boundary:
AE[end] = AE[end]/2.0
APx[end] = APx[end] + re[end]*ue_min[end]/(2.0*DXp[end]*rp[end])

# build the sparse matrix based on the numbering system
rowx_index = reshape(G[2:Nx+1],Nx) # main diagonal x
iix[1:3*Nx] .= repeat(rowx_index,3)
jjx[1:3*Nx] .= [reshape(G[1:Nx],Nx); reshape(G[2:Nx+1],Nx); reshape(G[3:Nx+2],Nx)]
sx[1:3*Nx] .= [AW; APx; AE]

# build the sparse matrix
kx = 3*Nx
M = sparse(iix[1:kx], jjx[1:kx], sx[1:kx], Nx+2, Nx+2)

end

# ================= 1D Convection Terms Cylindrical TVD =====================
function convectionTvdTermCylindrical1D(u::FaceValue, phi::CellValue, FL::Function)
# u is a face variable
# phi is a cell variable

# a function to avoid division by zero
eps1 = 1.0e-20
fsign(phi_in) = (abs.(phi_in).>=eps1).*phi_in+eps1*(phi_in.==0.0)+eps1*(abs.(phi_in).<eps1).*sign.(phi_in)

# extract data from the mesh structure
Nx = u.domain.dims[1]
G = [1:Nx+2;]
DXp = u.domain.cellsize.x[2:end-1]
dx = 0.5.*(u.domain.cellsize.x[1:end-1].+u.domain.cellsize.x[2:end])
RHS = zeros(Float64, Nx+2)
psi_p = zeros(Float64, Nx+1)
psi_m = zeros(Float64, Nx+1)
r = u.domain.cellcenters.x
rf = u.domain.facecenters.x

# define the vectors to store the sparse matrix data
iix = zeros(Int64, 3*(Nx+2))
jjx = zeros(Int64, 3*(Nx+2))
sx = zeros(Float64, 3*(Nx+2))

# calculate the upstream to downstream gradient ratios for u>0 (+ ratio)
dphi_p = (phi.value[2:Nx+2].-phi.value[1:Nx+1])./dx
rp = dphi_p[1:end-1]./fsign(dphi_p[2:end])
psi_p[2:Nx+1] .= 0.5.*FL.(rp).*(phi.value[3:Nx+2].-phi.value[2:Nx+1])
psi_p[1] = 0.0 # left boundary will be handled explicitly

# calculate the upstream to downstream gradient ratios for u<0 (- ratio)
rm = dphi_p[2:end]./fsign(dphi_p[1:end-1])
psi_m[1:Nx] .= 0.5.*FL.(rm).*(phi.value[1:Nx].-phi.value[2:Nx+1])
psi_m[Nx+1] = 0.0 # right boundary will be handled explicitly

# reassign the east, west for code readability
re = rf[2:Nx+1]
rw = rf[1:Nx]

# find the velocity direction for the upwind scheme
ue_min = min.(u.xvalue[2:Nx+1],0.0)
ue_max = max.(u.xvalue[2:Nx+1],0.0)
uw_min = min.(u.xvalue[1:Nx],0.0)
uw_max = max.(u.xvalue[1:Nx],0.0)

# calculate the TVD correction term
RHS[2:Nx+1] .= .-(1.0./(DXp.*r)).*(re.*(ue_max.*psi_p[2:Nx+1].+ue_min.*psi_m[2:Nx+1])-
              rw.*(uw_max.*psi_p[1:Nx].+uw_min.*psi_m[1:Nx]))

# calculate the coefficients for the internal cells
AE = reshape(re.*ue_min./(DXp.*r),Nx)
AW = reshape(-rw.*uw_max./(DXp.*r),Nx)
APx = reshape((re.*ue_max.-rw.*uw_min)./(DXp.*r),Nx)

# correct for the cells next to the boundary
# Left boundary:
APx[1] = APx[1].-rw[1]*uw_max[1]/(2.0*DXp[1]*r[1])
AW[1] = AW[1]/2.0
# Right boundary:
AE[end] = AE[end]/2.0
APx[end] = APx[end] + re[end]*ue_min[end]/(2.0*DXp[end]*r[end])

# build the sparse matrix based on the numbering system
rowx_index = reshape(G[2:Nx+1],Nx) # main diagonal x
iix[1:3*Nx] .= repeat(rowx_index,3)
jjx[1:3*Nx] .= [reshape(G[1:Nx],Nx); reshape(G[2:Nx+1],Nx); reshape(G[3:Nx+2],Nx)]
sx[1:3*Nx] .= [AW; APx; AE]

# build the sparse matrix
kx = 3*Nx
M = sparse(iix[1:kx], jjx[1:kx], sx[1:kx], Nx+2, Nx+2)

(M, RHS)
end

######################### Convection TVD only RHS ############################
function convectionTvdRHSCylindrical1D(u::FaceValue, phi::CellValue, FL::Function)
# u is a face variable
# phi is a cell variable

# a function to avoid division by zero
eps1 = 1.0e-20
fsign(phi_in) = (abs.(phi_in).>=eps1).*phi_in+eps1*(phi_in.==0.0)+eps1*(abs.(phi_in).<eps1).*sign.(phi_in)

# extract data from the mesh structure
Nx = u.domain.dims[1]
G = [1:Nx+2;]
DXp = u.domain.cellsize.x[2:end-1]
dx = 0.5.*(u.domain.cellsize.x[1:end-1].+u.domain.cellsize.x[2:end])
RHS = zeros(Float64, Nx+2)
psi_p = zeros(Float64, Nx+1)
psi_m = zeros(Float64, Nx+1)
r = u.domain.cellcenters.x
rf = u.domain.facecenters.x

# calculate the upstream to downstream gradient ratios for u>0 (+ ratio)
dphi_p = (phi.value[2:Nx+2].-phi.value[1:Nx+1])./dx
rp = dphi_p[1:end-1]./fsign(dphi_p[2:end])
psi_p[2:Nx+1] .= 0.5.*FL.(rp).*(phi.value[3:Nx+2].-phi.value[2:Nx+1])
psi_p[1] = 0.0 # left boundary will be handled explicitly

# calculate the upstream to downstream gradient ratios for u<0 (- ratio)
rm = dphi_p[2:end]./fsign(dphi_p[1:end-1])
psi_m[1:Nx] .= 0.5.*FL.(rm).*(phi.value[1:Nx].-phi.value[2:Nx+1])
psi_m[Nx+1] = 0.0 # right boundary will be handled explicitly

# reassign the east, west for code readability
re = rf[2:Nx+1]
rw = rf[1:Nx]

# find the velocity direction for the upwind scheme
ue_min = min.(u.xvalue[2:Nx+1],0.0)
ue_max = max.(u.xvalue[2:Nx+1],0.0)
uw_min = min.(u.xvalue[1:Nx],0.0)
uw_max = max.(u.xvalue[1:Nx],0.0)

# calculate the TVD correction term
RHS[2:Nx+1] .= .-(1.0./(DXp.*r)).*(re.*(ue_max.*psi_p[2:Nx+1].+ue_min.*psi_m[2:Nx+1])-
              rw.*(uw_max.*psi_p[1:Nx].+uw_min.*psi_m[1:Nx]))
return RHS
end

function convectionTvdRHSCylindrical1D(u::FaceValue, phi::CellValue, FL::Function, u_upwind::FaceValue)
# u is a face variable
# phi is a cell variable

# a function to avoid division by zero
eps1 = 1.0e-20
fsign(phi_in) = (abs.(phi_in).>=eps1).*phi_in+eps1*(phi_in.==0.0)+eps1*(abs.(phi_in).<eps1).*sign.(phi_in)

# extract data from the mesh structure
Nx = u.domain.dims[1]
G = [1:Nx+2;]
DXp = u.domain.cellsize.x[2:end-1]
dx = 0.5.*(u.domain.cellsize.x[1:end-1].+u.domain.cellsize.x[2:end])
RHS = zeros(Float64, Nx+2)
psi_p = zeros(Float64, Nx+1)
psi_m = zeros(Float64, Nx+1)
r = u.domain.cellcenters.x
rf = u.domain.facecenters.x

# calculate the upstream to downstream gradient ratios for u>0 (+ ratio)
dphi_p = (phi.value[2:Nx+2].-phi.value[1:Nx+1])./dx
rp = dphi_p[1:end-1]./fsign(dphi_p[2:end])
psi_p[2:Nx+1] .= 0.5.*FL.(rp).*(phi.value[3:Nx+2].-phi.value[2:Nx+1])
psi_p[1] = 0.0 # left boundary will be handled explicitly

# calculate the upstream to downstream gradient ratios for u<0 (- ratio)
rm = dphi_p[2:end]./fsign(dphi_p[1:end-1])
psi_m[1:Nx] .= 0.5.*FL.(rm).*(phi.value[1:Nx].-phi.value[2:Nx+1])
psi_m[Nx+1] = 0.0 # right boundary will be handled explicitly

# reassign the east, west for code readability
re = rf[2:Nx+1]
rw = rf[1:Nx]

# find the velocity direction for the upwind scheme
ue_min = u.xvalue[2:Nx+1]
ue_max = u.xvalue[2:Nx+1]
uw_min = u.xvalue[1:Nx]
uw_max = u.xvalue[1:Nx]

ue_min[u_upwind.xvalue[2:Nx+1].>0.0] .= 0.0
ue_max[u_upwind.xvalue[2:Nx+1].<0.0] .= 0.0
uw_min[u_upwind.xvalue[1:Nx].>0.0] .= 0.0
uw_max[u_upwind.xvalue[1:Nx].<0.0] .= 0.0

# calculate the TVD correction term
RHS[2:Nx+1] .= -(1.0./(DXp.*r)).*(re.*(ue_max.*psi_p[2:Nx+1].+ue_min.*psi_m[2:Nx+1])-
              rw.*(uw_max.*psi_p[1:Nx].+uw_min.*psi_m[1:Nx]))
return RHS
end

# ================= 1D Convection Terms Cartesian TVD =====================
function convectionTvdTerm1D(u::FaceValue, phi::CellValue, FL::Function)
# u is a face variable
# phi is a cell variable

# a function to avoid division by zero
eps1 = 1.0e-20
fsign(phi_in) = (abs.(phi_in).>=eps1).*phi_in+eps1*(phi_in.==0.0)+eps1*(abs.(phi_in).<eps1).*sign.(phi_in)

# extract data from the mesh structure
Nx = u.domain.dims[1]
G = [1:Nx+2;]
DXp = u.domain.cellsize.x[2:end-1]
dx = 0.5.*(u.domain.cellsize.x[1:end-1].+u.domain.cellsize.x[2:end])
RHS = zeros(Float64, Nx+2)
psi_p = zeros(Float64, Nx+1)
psi_m = zeros(Float64, Nx+1)

# define the vectors to store the sparse matrix data
iix = zeros(Int64, 3*(Nx+2))
jjx = zeros(Int64, 3*(Nx+2))
sx = zeros(Float64, 3*(Nx+2))

# calculate the upstream to downstream gradient ratios for u>0 (+ ratio)
dphi_p = (phi.value[2:Nx+2].-phi.value[1:Nx+1])./dx
rp = dphi_p[1:end-1]./fsign(dphi_p[2:end])
psi_p[2:Nx+1] .= 0.5.*FL.(rp).*(phi.value[3:Nx+2].-phi.value[2:Nx+1])
psi_p[1] = 0.0 # left boundary will be handled explicitly

# calculate the upstream to downstream gradient ratios for u<0 (- ratio)
rm = dphi_p[2:end]./fsign(dphi_p[1:end-1])
psi_m[1:Nx] .= 0.5.*FL.(rm).*(phi.value[1:Nx].-phi.value[2:Nx+1])
psi_m[Nx+1] = 0.0 # right boundary will be handled explicitly

# find the velocity direction for the upwind scheme
ue_min = min.(u.xvalue[2:Nx+1],0.0)
ue_max = max.(u.xvalue[2:Nx+1],0.0)
uw_min = min.(u.xvalue[1:Nx],0.0)
uw_max = max.(u.xvalue[1:Nx],0.0)

# calculate the TVD correction term
RHS[2:Nx+1] .= -(1.0./DXp).*((ue_max.*psi_p[2:Nx+1].+ue_min.*psi_m[2:Nx+1])-
              (uw_max.*psi_p[1:Nx].+uw_min.*psi_m[1:Nx]))

# calculate the coefficients for the internal cells
AE = reshape(ue_min./DXp,Nx)
AW = reshape(-uw_max./DXp,Nx)
APx = reshape((ue_max.-uw_min)./DXp,Nx)

# correct for the cells next to the boundary
# Left boundary:
APx[1] = APx[1].-uw_max[1]/(2.0*DXp[1])
AW[1] = AW[1]/2.0
# Right boundary:
AE[end] = AE[end]/2.0
APx[end] = APx[end] + ue_min[end]/(2.0*DXp[end])

# build the sparse matrix based on the numbering system
rowx_index = reshape(G[2:Nx+1],Nx) # main diagonal x
iix[1:3*Nx] .= repeat(rowx_index,3)
jjx[1:3*Nx] .= [reshape(G[1:Nx],Nx); reshape(G[2:Nx+1],Nx); reshape(G[3:Nx+2],Nx)]
sx[1:3*Nx] .= [AW; APx; AE]

# build the sparse matrix
kx = 3*Nx
M = sparse(iix[1:kx], jjx[1:kx], sx[1:kx], Nx+2, Nx+2)

(M, RHS)
end

# ================= 1D Convection RHS Cartesian TVD =====================
function convectionTvdRHS1D(u::FaceValue, phi::CellValue, FL::Function)
# u is a face variable
# phi is a cell variable

# a function to avoid division by zero
eps1 = 1.0e-20
fsign(phi_in) = (abs.(phi_in).>=eps1).*phi_in+eps1*(phi_in.==0.0)+eps1*(abs.(phi_in).<eps1).*sign.(phi_in)

# extract data from the mesh structure
Nx = u.domain.dims[1]
G = [1:Nx+2;]
DXp = u.domain.cellsize.x[2:end-1]
dx = 0.5.*(u.domain.cellsize.x[1:end-1].+u.domain.cellsize.x[2:end])
RHS = zeros(Float64, Nx+2)
psi_p = zeros(Float64, Nx+1)
psi_m = zeros(Float64, Nx+1)

# calculate the upstream to downstream gradient ratios for u>0 (+ ratio)
dphi_p = (phi.value[2:Nx+2].-phi.value[1:Nx+1])./dx
rp = dphi_p[1:end-1]./fsign(dphi_p[2:end])
psi_p[2:Nx+1] .= 0.5.*FL.(rp).*(phi.value[3:Nx+2].-phi.value[2:Nx+1])
psi_p[1] = 0.0 # left boundary will be handled explicitly

# calculate the upstream to downstream gradient ratios for u<0 (- ratio)
rm = dphi_p[2:end]./fsign(dphi_p[1:end-1])
psi_m[1:Nx] .= 0.5.*FL.(rm).*(phi.value[1:Nx].-phi.value[2:Nx+1])
psi_m[Nx+1] = 0.0 # right boundary will be handled explicitly

# find the velocity direction for the upwind scheme
ue_min = min.(u.xvalue[2:Nx+1],0.0)
ue_max = max.(u.xvalue[2:Nx+1],0.0)
uw_min = min.(u.xvalue[1:Nx],0.0)
uw_max = max.(u.xvalue[1:Nx],0.0)

# calculate the TVD correction term
RHS[2:Nx+1] .= -(1.0./DXp).*((ue_max.*psi_p[2:Nx+1].+ue_min.*psi_m[2:Nx+1])-
              (uw_max.*psi_p[1:Nx].+uw_min.*psi_m[1:Nx]))
return RHS
end

function convectionTvdRHS1D(u::FaceValue, phi::CellValue, FL::Function, u_upwind::FaceValue)
# u is a face variable
# phi is a cell variable

# a function to avoid division by zero
eps1 = 1.0e-20
fsign(phi_in) = (abs.(phi_in).>=eps1).*phi_in+eps1*(phi_in.==0.0)+eps1*(abs.(phi_in).<eps1).*sign.(phi_in)

# extract data from the mesh structure
Nx = u.domain.dims[1]
G = [1:Nx+2;]
DXp = u.domain.cellsize.x[2:end-1]
dx = 0.5.*(u.domain.cellsize.x[1:end-1].+u.domain.cellsize.x[2:end])
RHS = zeros(Float64, Nx+2)
psi_p = zeros(Float64, Nx+1)
psi_m = zeros(Float64, Nx+1)

# calculate the upstream to downstream gradient ratios for u>0 (+ ratio)
dphi_p = (phi.value[2:Nx+2].-phi.value[1:Nx+1])./dx
rp = dphi_p[1:end-1]./fsign(dphi_p[2:end])
psi_p[2:Nx+1] .= 0.5.*FL.(rp).*(phi.value[3:Nx+2].-phi.value[2:Nx+1])
psi_p[1] = 0.0 # left boundary will be handled explicitly

# calculate the upstream to downstream gradient ratios for u<0 (- ratio)
rm = dphi_p[2:end]./fsign(dphi_p[1:end-1])
psi_m[1:Nx] .= 0.5.*FL.(rm).*(phi.value[1:Nx].-phi.value[2:Nx+1])
psi_m[Nx+1] = 0.0 # right boundary will be handled explicitly

# find the velocity direction for the upwind scheme
ue_min = u.xvalue[2:Nx+1]
ue_max = u.xvalue[2:Nx+1]
uw_min = u.xvalue[1:Nx]
uw_max = u.xvalue[1:Nx]

ue_min[u_upwind.xvalue[2:Nx+1].>0.0] .= 0.0
ue_max[u_upwind.xvalue[2:Nx+1].<0.0] .= 0.0
uw_min[u_upwind.xvalue[1:Nx].>0.0] .= 0.0
uw_max[u_upwind.xvalue[1:Nx].<0.0] .= 0.0

# calculate the TVD correction term
RHS[2:Nx+1] .= -(1.0./DXp).*((ue_max.*psi_p[2:Nx+1].+ue_min.*psi_m[2:Nx+1])-
              (uw_max.*psi_p[1:Nx].+uw_min.*psi_m[1:Nx]))
return RHS
end


# =================== 2D Convection Terms =======================
# ================ 2D Convection Term Central ===================
function convectionTerm2D(u::FaceValue)
# u is a face variable
# extract data from the mesh structure
Nx = u.domain.dims[1]
Ny = u.domain.dims[2]
G=reshape([1:(Nx+2)*(Ny+2);], Nx+2, Ny+2)
DXe = u.domain.cellsize.x[3:end]
DXw = u.domain.cellsize.x[1:end-2]
DXp = u.domain.cellsize.x[2:end-1]
DY = zeros( 1, Ny+2)
DY[:] .= u.domain.cellsize.y
DYn = DY[:,3:end]
DYs = DY[:,1:end-2]
DYp = DY[:,2:end-1]

# define the vectors to store the sparse matrix data
iix = zeros(Int64, 3*(Nx+2)*(Ny+2))
iiy = zeros(Int64, 3*(Nx+2)*(Ny+2))
jjx = zeros(Int64, 3*(Nx+2)*(Ny+2))
jjy = zeros(Int64, 3*(Nx+2)*(Ny+2))
sx = zeros(Float64, 3*(Nx+2)*(Ny+2))
sy = zeros(Float64, 3*(Nx+2)*(Ny+2))
mnx = Nx*Ny
mny = Nx*Ny

# reassign the east, west for code readability
ue = u.xvalue[2:Nx+1,:]./(DXp+DXe)
uw = u.xvalue[1:Nx,:]./(DXp+DXw)
vn = u.yvalue[:,2:Ny+1]./(DYp+DYn)
vs = u.yvalue[:,1:Ny]./(DYp+DYs)

# calculate the coefficients for the internal cells
AE = reshape(ue,mnx)
AW = reshape(-uw,mnx)
AN = reshape(vn,mny)
AS = reshape(-vs,mny)
APx = reshape((ue.*DXe-uw.*DXw)./DXp,mnx)
APy = reshape((vn.*DYn-vs.*DYs)./DYp,mny)

# build the sparse matrix based on the numbering system
rowx_index = reshape(G[2:Nx+1,2:Ny+1],mnx) # main diagonal x
iix[1:3*mnx] .= repeat(rowx_index,3)
rowy_index = reshape(G[2:Nx+1,2:Ny+1],mny) # main diagonal y
iiy[1:3*mny] .= repeat(rowy_index,3)
jjx[1:3*mnx] .= [reshape(G[1:Nx,2:Ny+1],mnx); reshape(G[2:Nx+1,2:Ny+1],mnx); reshape(G[3:Nx+2,2:Ny+1],mnx)]
jjy[1:3*mny] .= [reshape(G[2:Nx+1,1:Ny],mny); reshape(G[2:Nx+1,2:Ny+1],mny); reshape(G[2:Nx+1,3:Ny+2],mny)]
sx[1:3*mnx] .= [AW; APx; AE]
sy[1:3*mny] .= [AS; APy; AN]

# build the sparse matrix
kx = 3*mnx
ky = 3*mny
Mx = sparse(iix[1:kx], jjx[1:kx], sx[1:kx], (Nx+2)*(Ny+2), (Nx+2)*(Ny+2))
My = sparse(iiy[1:ky], jjy[1:ky], sy[1:ky], (Nx+2)*(Ny+2), (Nx+2)*(Ny+2))
M = Mx + My
(M, Mx, My)
end

# ================ 2D Convection Term Upwind ===================
function convectionUpwindTerm2D(u::FaceValue)
# u is a face variable
# extract data from the mesh structure
Nx = u.domain.dims[1]
Ny = u.domain.dims[2]
G=reshape([1:(Nx+2)*(Ny+2);], Nx+2, Ny+2)
DXp = u.domain.cellsize.x[2:end-1]
DYp = zeros( 1, Ny)
DYp[:] .= u.domain.cellsize.y[2:end-1]

# define the vectors to store the sparse matrix data
iix = zeros(Int64, 3*(Nx+2)*(Ny+2))
iiy = zeros(Int64, 3*(Nx+2)*(Ny+2))
jjx = zeros(Int64, 3*(Nx+2)*(Ny+2))
jjy = zeros(Int64, 3*(Nx+2)*(Ny+2))
sx = zeros(Float64, 3*(Nx+2)*(Ny+2))
sy = zeros(Float64, 3*(Nx+2)*(Ny+2))
mnx = Nx*Ny
mny = Nx*Ny

# find the velocity direction for the upwind scheme
ue_min = min.(u.xvalue[2:Nx+1,:],0.0)
ue_max = max.(u.xvalue[2:Nx+1,:],0.0)
uw_min = min.(u.xvalue[1:Nx,:],0.0)
uw_max = max.(u.xvalue[1:Nx,:],0.0)
vn_min = min.(u.yvalue[:,2:Ny+1],0.0)
vn_max = max.(u.yvalue[:,2:Ny+1],0.0)
vs_min = min.(u.yvalue[:,1:Ny],0.0)
vs_max = max.(u.yvalue[:,1:Ny],0.0)

# calculate the coefficients for the internal cells, not reshape
AE = ue_min./DXp
AW = -uw_max./DXp
AN = vn_min./DYp
AS = -vs_max./DYp
APx = (ue_max.-uw_min)./DXp
APy = (vn_max.-vs_min)./DYp

# Also correct for the boundary cells (not the ghost cells)
# Left boundary:
APx[1,:] .= APx[1,:].-uw_max[1,:]/(2.0*DXp[1])
AW[1,:] .= AW[1,:]/2.0
# Right boundary:
AE[end,:] .= AE[end,:]/2.0
APx[end,:] .= APx[end,:].+ue_min[end,:]/(2.0*DXp[end])
# Bottom boundary:
APy[:,1] .= APy[:,1].-vs_max[:,1]/(2.0*DYp[1])
AS[:,1] .= AS[:,1]/2.0
# Top boundary:
AN[:,end] .= AN[:,end]/2.0
APy[:,end] .= APy[:,end].+vn_min[:,end]/(2.0*DYp[end])

# now reshape
AE = reshape(AE,mnx)
AW = reshape(AW,mnx)
AN = reshape(AN,mny)
AS = reshape(AS,mny)
APx = reshape(APx,mnx)
APy = reshape(APy,mny)

# build the sparse matrix based on the numbering system
rowx_index = reshape(G[2:Nx+1,2:Ny+1],mnx) # main diagonal x
iix[1:3*mnx] .= repeat(rowx_index,3)
rowy_index = reshape(G[2:Nx+1,2:Ny+1],mny) # main diagonal y
iiy[1:3*mny] .= repeat(rowy_index,3)
jjx[1:3*mnx] .= [reshape(G[1:Nx,2:Ny+1],mnx); reshape(G[2:Nx+1,2:Ny+1],mnx); reshape(G[3:Nx+2,2:Ny+1],mnx)]
jjy[1:3*mny] .= [reshape(G[2:Nx+1,1:Ny],mny); reshape(G[2:Nx+1,2:Ny+1],mny); reshape(G[2:Nx+1,3:Ny+2],mny)]
sx[1:3*mnx] .= [AW; APx; AE]
sy[1:3*mny] .= [AS; APy; AN]

# build the sparse matrix
kx = 3*mnx
ky = 3*mny
Mx = sparse(iix[1:kx], jjx[1:kx], sx[1:kx], (Nx+2)*(Ny+2), (Nx+2)*(Ny+2))
My = sparse(iiy[1:ky], jjy[1:ky], sy[1:ky], (Nx+2)*(Ny+2), (Nx+2)*(Ny+2))
M = Mx + My
(M, Mx, My)
end

function convectionUpwindTerm2D(u::FaceValue, u_upwind::FaceValue)
# u is a face variable
# extract data from the mesh structure
Nx = u.domain.dims[1]
Ny = u.domain.dims[2]
G=reshape([1:(Nx+2)*(Ny+2);], Nx+2, Ny+2)
DXp = u.domain.cellsize.x[2:end-1]
DYp = zeros( 1, Ny)
DYp[:] .= u.domain.cellsize.y[2:end-1]

# define the vectors to store the sparse matrix data
iix = zeros(Int64, 3*(Nx+2)*(Ny+2))
iiy = zeros(Int64, 3*(Nx+2)*(Ny+2))
jjx = zeros(Int64, 3*(Nx+2)*(Ny+2))
jjy = zeros(Int64, 3*(Nx+2)*(Ny+2))
sx = zeros(Float64, 3*(Nx+2)*(Ny+2))
sy = zeros(Float64, 3*(Nx+2)*(Ny+2))
mnx = Nx*Ny
mny = Nx*Ny

# find the velocity direction for the upwind scheme
ue_min = u.xvalue[2:Nx+1,:]
ue_max = u.xvalue[2:Nx+1,:]
uw_min = u.xvalue[1:Nx,:]
uw_max = u.xvalue[1:Nx,:]
vn_min = u.yvalue[:,2:Ny+1]
vn_max = u.yvalue[:,2:Ny+1]
vs_min = u.yvalue[:,1:Ny]
vs_max = u.yvalue[:,1:Ny]

ue_min[u_upwind.xvalue[2:Nx+1,:].>0.0] .= 0.0
ue_max[u_upwind.xvalue[2:Nx+1,:].<0.0] .= 0.0
uw_min[u_upwind.xvalue[1:Nx,:].>0.0] .= 0.0
uw_max[u_upwind.xvalue[1:Nx,:].<0.0] .= 0.0
vn_min[u_upwind.yvalue[:,2:Ny+1].>0.0] .= 0.0
vn_max[u_upwind.yvalue[:,2:Ny+1].<0.0] .= 0.0
vs_min[u_upwind.yvalue[:,1:Ny].>0.0] .= 0.0
vs_max[u_upwind.yvalue[:,1:Ny].<0.0] .= 0.0

# calculate the coefficients for the internal cells, not reshape
AE = ue_min./DXp
AW = -uw_max./DXp
AN = vn_min./DYp
AS = -vs_max./DYp
APx = (ue_max.-uw_min)./DXp
APy = (vn_max.-vs_min)./DYp

# Also correct for the boundary cells (not the ghost cells)
# Left boundary:
APx[1,:] .= APx[1,:].-uw_max[1,:]/(2.0*DXp[1])
AW[1,:] .= AW[1,:]/2.0
# Right boundary:
AE[end,:] .= AE[end,:]/2.0
APx[end,:] .= APx[end,:].+ue_min[end,:]/(2.0*DXp[end])
# Bottom boundary:
APy[:,1] .= APy[:,1].-vs_max[:,1]/(2.0*DYp[1])
AS[:,1] .= AS[:,1]/2.0
# Top boundary:
AN[:,end] .= AN[:,end]/2.0
APy[:,end] .= APy[:,end].+vn_min[:,end]/(2.0*DYp[end])

# now reshape
AE = reshape(AE,mnx)
AW = reshape(AW,mnx)
AN = reshape(AN,mny)
AS = reshape(AS,mny)
APx = reshape(APx,mnx)
APy = reshape(APy,mny)

# build the sparse matrix based on the numbering system
rowx_index = reshape(G[2:Nx+1,2:Ny+1],mnx) # main diagonal x
iix[1:3*mnx] .= repeat(rowx_index,3)
rowy_index = reshape(G[2:Nx+1,2:Ny+1],mny) # main diagonal y
iiy[1:3*mny] .= repeat(rowy_index,3)
jjx[1:3*mnx] .= [reshape(G[1:Nx,2:Ny+1],mnx); reshape(G[2:Nx+1,2:Ny+1],mnx); reshape(G[3:Nx+2,2:Ny+1],mnx)]
jjy[1:3*mny] .= [reshape(G[2:Nx+1,1:Ny],mny); reshape(G[2:Nx+1,2:Ny+1],mny); reshape(G[2:Nx+1,3:Ny+2],mny)]
sx[1:3*mnx] .= [AW; APx; AE]
sy[1:3*mny] .= [AS; APy; AN]

# build the sparse matrix
kx = 3*mnx
ky = 3*mny
Mx = sparse(iix[1:kx], jjx[1:kx], sx[1:kx], (Nx+2)*(Ny+2), (Nx+2)*(Ny+2))
My = sparse(iiy[1:ky], jjy[1:ky], sy[1:ky], (Nx+2)*(Ny+2), (Nx+2)*(Ny+2))
M = Mx + My
(M, Mx, My)
end

# ================ 2D Convection Term TVD =======================
function convectionTvdTerm2D(u::FaceValue, phi::CellValue, FL::Function)
# u is a face variable
# phi is a cell variable

# a function to avoid division by zero
eps1 = 1.0e-20
fsign(phi_in) = (abs.(phi_in).>=eps1).*phi_in+eps1*(phi_in.==0.0)+eps1*(abs.(phi_in).<eps1).*sign.(phi_in)

# extract data from the mesh structure
Nx = u.domain.dims[1]
Ny = u.domain.dims[2]
G=reshape([1:(Nx+2)*(Ny+2);], Nx+2, Ny+2)
DXp = u.domain.cellsize.x[2:end-1]
DYp = zeros( 1, Ny)
DYp[:] .= u.domain.cellsize.y[2:end-1]
dx=0.5.*(u.domain.cellsize.x[1:end-1].+u.domain.cellsize.x[2:end])
dy=zeros( 1, Ny+1)
dy[:] .= 0.5.*(u.domain.cellsize.y[1:end-1].+u.domain.cellsize.y[2:end])

psiX_p = zeros(Nx+1,Ny)
psiX_m = zeros(Nx+1,Ny)
psiY_p = zeros(Nx,Ny+1)
psiY_m = zeros(Nx,Ny+1)

# define the vectors to store the sparse matrix data
iix = zeros(Int64, 3*(Nx+2)*(Ny+2))
iiy = zeros(Int64, 3*(Nx+2)*(Ny+2))
jjx = zeros(Int64, 3*(Nx+2)*(Ny+2))
jjy = zeros(Int64, 3*(Nx+2)*(Ny+2))
sx = zeros(Float64, 3*(Nx+2)*(Ny+2))
sy = zeros(Float64, 3*(Nx+2)*(Ny+2))
mnx = Nx*Ny
mny = Nx*Ny

# calculate the upstream to downstream gradient ratios for u>0 (+ ratio)
# x direction
dphiX_p = (phi.value[2:Nx+2, 2:Ny+1].-phi.value[1:Nx+1, 2:Ny+1])./dx
rX_p = dphiX_p[1:end-1,:]./fsign(dphiX_p[2:end,:])
psiX_p[2:Nx+1,:] .= 0.5.*FL.(rX_p).*(phi.value[3:Nx+2,2:Ny+1].-
		    phi.value[2:Nx+1, 2:Ny+1])
psiX_p[1, :] .= 0.0 # left boundary will be handled in the main matrix
# y direction
dphiY_p = (phi.value[2:Nx+1, 2:Ny+2].-phi.value[2:Nx+1, 1:Ny+1])./dy
rY_p = dphiY_p[:,1:end-1]./fsign(dphiY_p[:,2:end])
psiY_p[:,2:Ny+1] .= 0.5.*FL.(rY_p).*(phi.value[2:Nx+1,3:Ny+2].-
		  phi.value[2:Nx+1, 2:Ny+1])
psiY_p[:,1] .= 0.0 # Bottom boundary will be handled in the main matrix

# calculate the upstream to downstream gradient ratios for u<0 (- ratio)
# x direction
rX_m = dphiX_p[2:end,:]./fsign(dphiX_p[1:end-1,:])
psiX_m[1:Nx,:] .= 0.5.*FL.(rX_m).*(phi.value[1:Nx, 2:Ny+1].-
		phi.value[2:Nx+1, 2:Ny+1])
psiX_m[Nx+1,:] .= 0.0 # right boundary
# y direction
rY_m = dphiY_p[:,2:end]./fsign(dphiY_p[:,1:end-1])
psiY_m[:,1:Ny] .= 0.5.*FL.(rY_m).*(phi.value[2:Nx+1, 1:Ny].-
	      phi.value[2:Nx+1, 2:Ny+1])
psiY_m[:, Ny+1] .= 0.0 # top boundary will be handled in the main matrix

# find the velocity direction for the upwind scheme
ue_min = min.(u.xvalue[2:Nx+1,:],0.0)
ue_max = max.(u.xvalue[2:Nx+1,:],0.0)
uw_min = min.(u.xvalue[1:Nx,:],0.0)
uw_max = max.(u.xvalue[1:Nx,:],0.0)
vn_min = min.(u.yvalue[:,2:Ny+1],0.0)
vn_max = max.(u.yvalue[:,2:Ny+1],0.0)
vs_min = min.(u.yvalue[:,1:Ny],0.0)
vs_max = max.(u.yvalue[:,1:Ny],0.0)

# calculate the coefficients for the internal cells, not reshape
AE = ue_min./DXp
AW = -uw_max./DXp
AN = vn_min./DYp
AS = -vs_max./DYp
APx = (ue_max.-uw_min)./DXp
APy = (vn_max.-vs_min)./DYp

# Also correct for the boundary cells (not the ghost cells)
# Left boundary:
APx[1,:] .= APx[1,:].-uw_max[1,:]/(2.0*DXp[1])
AW[1,:] .= AW[1,:]/2.0
# Right boundary:
AE[end,:] .= AE[end,:]/2.0
APx[end,:] .= APx[end,:].+ue_min[end,:]/(2.0*DXp[end])
# Bottom boundary:
APy[:,1] .= APy[:,1].-vs_max[:,1]/(2.0*DYp[1])
AS[:,1] .= AS[:,1]/2.0
# Top boundary:
AN[:,end] .= AN[:,end]/2.0
APy[:,end] .= APy[:,end].+vn_min[:,end]/(2.0*DYp[end])

# now reshape
AE = reshape(AE,mnx)
AW = reshape(AW,mnx)
AN = reshape(AN,mny)
AS = reshape(AS,mny)
APx = reshape(APx,mnx)
APy = reshape(APy,mny)

# build the sparse matrix based on the numbering system
rowx_index = reshape(G[2:Nx+1,2:Ny+1],mnx) # main diagonal x
iix[1:3*mnx] .= repeat(rowx_index,3)
rowy_index = reshape(G[2:Nx+1,2:Ny+1],mny) # main diagonal y
iiy[1:3*mny] .= repeat(rowy_index,3)
jjx[1:3*mnx] .= [reshape(G[1:Nx,2:Ny+1],mnx); reshape(G[2:Nx+1,2:Ny+1],mnx); reshape(G[3:Nx+2,2:Ny+1],mnx)]
jjy[1:3*mny] .= [reshape(G[2:Nx+1,1:Ny],mny); reshape(G[2:Nx+1,2:Ny+1],mny); reshape(G[2:Nx+1,3:Ny+2],mny)]
sx[1:3*mnx] .= [AW; APx; AE]
sy[1:3*mny] .= [AS; APy; AN]

# calculate the TVD correction term
div_x = -(1.0./DXp).*((ue_max.*psiX_p[2:Nx+1,:].+ue_min.*psiX_m[2:Nx+1,:])-
              (uw_max.*psiX_p[1:Nx,:].+uw_min.*psiX_m[1:Nx,:]))
div_y = -(1.0./DYp).*((vn_max.*psiY_p[:,2:Ny+1].+vn_min.*psiY_m[:,2:Ny+1])-
              (vs_max.*psiY_p[:,1:Ny].+vs_min.*psiY_m[:,1:Ny]))

# define the RHS Vector
RHS = zeros((Nx+2)*(Ny+2))
RHSx = zeros((Nx+2)*(Ny+2))
RHSy = zeros((Nx+2)*(Ny+2))

# assign the values of the RHS vector
RHS[rowx_index] .= reshape(div_x+div_y,Nx*Ny)
RHSx[rowx_index] .= reshape(div_x,Nx*Ny)
RHSy[rowy_index] .= reshape(div_y,Nx*Ny)

# build the sparse matrix
kx = 3*mnx
ky = 3*mny
Mx = sparse(iix[1:kx], jjx[1:kx], sx[1:kx], (Nx+2)*(Ny+2), (Nx+2)*(Ny+2))
My = sparse(iiy[1:ky], jjy[1:ky], sy[1:ky], (Nx+2)*(Ny+2), (Nx+2)*(Ny+2))
M = Mx + My
(M, RHS, Mx, My, RHSx, RHSy)
end

# ================ 2D Convection Term TVD =======================
function convectionTvdRHS2D(u::FaceValue, phi::CellValue, FL::Function)
# u is a face variable
# phi is a cell variable

# a function to avoid division by zero
eps1 = 1.0e-20
fsign(phi_in) = (abs.(phi_in).>=eps1).*phi_in+eps1*(phi_in.==0.0)+eps1*(abs.(phi_in).<eps1).*sign.(phi_in)

# extract data from the mesh structure
Nx = u.domain.dims[1]
Ny = u.domain.dims[2]
G=reshape([1:(Nx+2)*(Ny+2);], Nx+2, Ny+2)
DXp = u.domain.cellsize.x[2:end-1]
DYp = zeros( 1, Ny)
DYp[:] .= u.domain.cellsize.y[2:end-1]
dx=0.5.*(u.domain.cellsize.x[1:end-1].+u.domain.cellsize.x[2:end])
dy=zeros( 1, Ny+1)
dy[:] .=0.5.*(u.domain.cellsize.y[1:end-1].+u.domain.cellsize.y[2:end])

psiX_p = zeros(Nx+1,Ny)
psiX_m = zeros(Nx+1,Ny)
psiY_p = zeros(Nx,Ny+1)
psiY_m = zeros(Nx,Ny+1)

# calculate the upstream to downstream gradient ratios for u>0 (+ ratio)
# x direction
dphiX_p = (phi.value[2:Nx+2, 2:Ny+1].-phi.value[1:Nx+1, 2:Ny+1])./dx
rX_p = dphiX_p[1:end-1,:]./fsign(dphiX_p[2:end,:])
psiX_p[2:Nx+1,:] .= 0.5.*FL.(rX_p).*(phi.value[3:Nx+2,2:Ny+1].-
		    phi.value[2:Nx+1, 2:Ny+1])
psiX_p[1, :] .= 0.0 # left boundary will be handled in the main matrix
# y direction
dphiY_p = (phi.value[2:Nx+1, 2:Ny+2].-phi.value[2:Nx+1, 1:Ny+1])./dy
rY_p = dphiY_p[:,1:end-1]./fsign(dphiY_p[:,2:end])
psiY_p[:,2:Ny+1] .= 0.5.*FL.(rY_p).*(phi.value[2:Nx+1,3:Ny+2].-
		  phi.value[2:Nx+1, 2:Ny+1])
psiY_p[:,1] .= 0.0 # Bottom boundary will be handled in the main matrix

# calculate the upstream to downstream gradient ratios for u<0 (- ratio)
# x direction
rX_m = dphiX_p[2:end,:]./fsign(dphiX_p[1:end-1,:])
psiX_m[1:Nx,:] .= 0.5.*FL.(rX_m).*(phi.value[1:Nx, 2:Ny+1].-
		phi.value[2:Nx+1, 2:Ny+1])
psiX_m[Nx+1,:] .= 0.0 # right boundary
# y direction
rY_m = dphiY_p[:,2:end]./fsign(dphiY_p[:,1:end-1])
psiY_m[:,1:Ny] .= 0.5.*FL.(rY_m).*(phi.value[2:Nx+1, 1:Ny].-
	      phi.value[2:Nx+1, 2:Ny+1])
psiY_m[:, Ny+1] .= 0.0 # top boundary will be handled in the main matrix

# find the velocity direction for the upwind scheme
ue_min = min.(u.xvalue[2:Nx+1,:],0.0)
ue_max = max.(u.xvalue[2:Nx+1,:],0.0)
uw_min = min.(u.xvalue[1:Nx,:],0.0)
uw_max = max.(u.xvalue[1:Nx,:],0.0)
vn_min = min.(u.yvalue[:,2:Ny+1],0.0)
vn_max = max.(u.yvalue[:,2:Ny+1],0.0)
vs_min = min.(u.yvalue[:,1:Ny],0.0)
vs_max = max.(u.yvalue[:,1:Ny],0.0)

# calculate the TVD correction term
div_x = -(1.0./DXp).*((ue_max.*psiX_p[2:Nx+1,:].+ue_min.*psiX_m[2:Nx+1,:])-
              (uw_max.*psiX_p[1:Nx,:].+uw_min.*psiX_m[1:Nx,:]))
div_y = -(1.0./DYp).*((vn_max.*psiY_p[:,2:Ny+1].+vn_min.*psiY_m[:,2:Ny+1])-
              (vs_max.*psiY_p[:,1:Ny].+vs_min.*psiY_m[:,1:Ny]))

# define the RHS Vector
RHS = zeros((Nx+2)*(Ny+2))
RHSx = zeros((Nx+2)*(Ny+2))
RHSy = zeros((Nx+2)*(Ny+2))

# assign the values of the RHS vector
mnx = Nx*Ny
mny = Nx*Ny
rowx_index = reshape(G[2:Nx+1,2:Ny+1],mnx) # main diagonal x
rowy_index = reshape(G[2:Nx+1,2:Ny+1],mny) # main diagonal y
RHS[rowx_index] .= reshape(div_x+div_y,Nx*Ny)
RHSx[rowx_index] .= reshape(div_x,Nx*Ny)
RHSy[rowy_index] .= reshape(div_y,Nx*Ny)

(RHS, RHSx, RHSy)
end

function convectionTvdRHS2D(u::FaceValue, phi::CellValue, FL::Function, u_upwind::FaceValue)
# u is a face variable
# phi is a cell variable

# a function to avoid division by zero
eps1 = 1.0e-20
fsign(phi_in) = (abs.(phi_in).>=eps1).*phi_in+eps1*(phi_in.==0.0)+eps1*(abs.(phi_in).<eps1).*sign.(phi_in)

# extract data from the mesh structure
Nx = u.domain.dims[1]
Ny = u.domain.dims[2]
G=reshape([1:(Nx+2)*(Ny+2);], Nx+2, Ny+2)
DXp = u.domain.cellsize.x[2:end-1]
DYp = zeros( 1, Ny)
DYp[:] .= u.domain.cellsize.y[2:end-1]
dx=0.5.*(u.domain.cellsize.x[1:end-1].+u.domain.cellsize.x[2:end])
dy=zeros( 1, Ny+1)
dy[:] .=0.5.*(u.domain.cellsize.y[1:end-1].+u.domain.cellsize.y[2:end])

psiX_p = zeros(Nx+1,Ny)
psiX_m = zeros(Nx+1,Ny)
psiY_p = zeros(Nx,Ny+1)
psiY_m = zeros(Nx,Ny+1)

# calculate the upstream to downstream gradient ratios for u>0 (+ ratio)
# x direction
dphiX_p = (phi.value[2:Nx+2, 2:Ny+1].-phi.value[1:Nx+1, 2:Ny+1])./dx
rX_p = dphiX_p[1:end-1,:]./fsign(dphiX_p[2:end,:])
psiX_p[2:Nx+1,:] .= 0.5.*FL.(rX_p).*(phi.value[3:Nx+2,2:Ny+1].-
		    phi.value[2:Nx+1, 2:Ny+1])
psiX_p[1, :] .= 0.0 # left boundary will be handled in the main matrix
# y direction
dphiY_p = (phi.value[2:Nx+1, 2:Ny+2].-phi.value[2:Nx+1, 1:Ny+1])./dy
rY_p = dphiY_p[:,1:end-1]./fsign(dphiY_p[:,2:end])
psiY_p[:,2:Ny+1] .= 0.5.*FL.(rY_p).*(phi.value[2:Nx+1,3:Ny+2].-
		  phi.value[2:Nx+1, 2:Ny+1])
psiY_p[:,1] .= 0.0 # Bottom boundary will be handled in the main matrix

# calculate the upstream to downstream gradient ratios for u<0 (- ratio)
# x direction
rX_m = dphiX_p[2:end,:]./fsign(dphiX_p[1:end-1,:])
psiX_m[1:Nx,:] .= 0.5.*FL.(rX_m).*(phi.value[1:Nx, 2:Ny+1].-
		phi.value[2:Nx+1, 2:Ny+1])
psiX_m[Nx+1,:] .= 0.0 # right boundary
# y direction
rY_m = dphiY_p[:,2:end]./fsign(dphiY_p[:,1:end-1])
psiY_m[:,1:Ny] .= 0.5.*FL.(rY_m).*(phi.value[2:Nx+1, 1:Ny].-
	      phi.value[2:Nx+1, 2:Ny+1])
psiY_m[:, Ny+1] .= 0.0 # top boundary will be handled in the main matrix

# find the velocity direction for the upwind scheme
ue_min = u.xvalue[2:Nx+1,:]
ue_max = u.xvalue[2:Nx+1,:]
uw_min = u.xvalue[1:Nx,:]
uw_max = u.xvalue[1:Nx,:]
vn_min = u.yvalue[:,2:Ny+1]
vn_max = u.yvalue[:,2:Ny+1]
vs_min = u.yvalue[:,1:Ny]
vs_max = u.yvalue[:,1:Ny]

ue_min[u_upwind.xvalue[2:Nx+1,:].>0.0] .= 0.0
ue_max[u_upwind.xvalue[2:Nx+1,:].<0.0] .= 0.0
uw_min[u_upwind.xvalue[1:Nx,:].>0.0] .= 0.0
uw_max[u_upwind.xvalue[1:Nx,:].<0.0] .= 0.0
vn_min[u_upwind.yvalue[:,2:Ny+1].>0.0] .= 0.0
vn_max[u_upwind.yvalue[:,2:Ny+1].<0.0] .= 0.0
vs_min[u_upwind.yvalue[:,1:Ny].>0.0] .= 0.0
vs_max[u_upwind.yvalue[:,1:Ny].<0.0] .= 0.0

# calculate the TVD correction term
div_x = -(1.0./DXp).*((ue_max.*psiX_p[2:Nx+1,:].+ue_min.*psiX_m[2:Nx+1,:])-
              (uw_max.*psiX_p[1:Nx,:].+uw_min.*psiX_m[1:Nx,:]))
div_y = -(1.0./DYp).*((vn_max.*psiY_p[:,2:Ny+1].+vn_min.*psiY_m[:,2:Ny+1])-
              (vs_max.*psiY_p[:,1:Ny].+vs_min.*psiY_m[:,1:Ny]))

# define the RHS Vector
RHS = zeros((Nx+2)*(Ny+2))
RHSx = zeros((Nx+2)*(Ny+2))
RHSy = zeros((Nx+2)*(Ny+2))

# assign the values of the RHS vector
mnx = Nx*Ny
mny = Nx*Ny
rowx_index = reshape(G[2:Nx+1,2:Ny+1],mnx) # main diagonal x
rowy_index = reshape(G[2:Nx+1,2:Ny+1],mny) # main diagonal y
RHS[rowx_index] .= reshape(div_x+div_y,Nx*Ny)
RHSx[rowx_index] .= reshape(div_x,Nx*Ny)
RHSy[rowy_index] .= reshape(div_y,Nx*Ny)

(RHS, RHSx, RHSy)
end

# ================ 2D Convection Term Cylindrical Central ===================
function convectionTermCylindrical2D(u::FaceValue)
# u is a face variable
# extract data from the mesh structure
Nr = u.domain.dims[1]
Nz = u.domain.dims[2]
G=reshape([1:(Nr+2)*(Nz+2);], Nr+2, Nz+2)
DXe = u.domain.cellsize.x[3:end]
DXw = u.domain.cellsize.x[1:end-2]
DXp = u.domain.cellsize.x[2:end-1]
DZ = zeros( 1, Nz+2)
DZ[:] .= u.domain.cellsize.y
DYn = DZ[:,3:end]
DYs = DZ[:,1:end-2]
DYp = DZ[:,2:end-1]
rp = repeat(u.domain.cellcenters.x, 1, Nz)
rf = repeat(u.domain.facecenters.x, 1, Nz)

# define the vectors to store the sparse matrix data
iix = zeros(Int64, 3*(Nr+2)*(Nz+2))
iiy = zeros(Int64, 3*(Nr+2)*(Nz+2))
jjx = zeros(Int64, 3*(Nr+2)*(Nz+2))
jjy = zeros(Int64, 3*(Nr+2)*(Nz+2))
sx = zeros(Float64, 3*(Nr+2)*(Nz+2))
sy = zeros(Float64, 3*(Nr+2)*(Nz+2))
mnx = Nr*Nz
mny = Nr*Nz

# reassign the east, west for code readability
ue = rf[2:Nr+1,:].*u.xvalue[2:Nr+1,:]./(rp.*(DXp+DXe))
uw = rf[1:Nr,:].*u.xvalue[1:Nr,:]./(rp.*(DXp+DXw))
vn = u.yvalue[:,2:Nz+1]./(DYp+DYn)
vs = u.yvalue[:,1:Nz]./(DYp+DYs)
re = rf[2:Nr+1,:]
rw = rf[1:Nr,:]

# calculate the coefficients for the internal cells
AE = reshape(ue,mnx)
AW = reshape(-uw,mnx)
AN = reshape(vn,mny)
AS = reshape(-vs,mny)
APx = reshape((ue.*DXe-uw.*DXw)./DXp,mnx)
APy = reshape((vn.*DYn-vs.*DYs)./DYp,mny)

# build the sparse matrix based on the numbering system
rowx_index = reshape(G[2:Nr+1,2:Nz+1],mnx) # main diagonal x
iix[1:3*mnx] .= repeat(rowx_index,3)
rowy_index = reshape(G[2:Nr+1,2:Nz+1],mny) # main diagonal y
iiy[1:3*mny] .= repeat(rowy_index,3)
jjx[1:3*mnx] .= [reshape(G[1:Nr,2:Nz+1],mnx); reshape(G[2:Nr+1,2:Nz+1],mnx); reshape(G[3:Nr+2,2:Nz+1],mnx)]
jjy[1:3*mny] .= [reshape(G[2:Nr+1,1:Nz],mny); reshape(G[2:Nr+1,2:Nz+1],mny); reshape(G[2:Nr+1,3:Nz+2],mny)]
sx[1:3*mnx] .= [AW; APx; AE]
sy[1:3*mny] .= [AS; APy; AN]

# build the sparse matrix
kx = 3*mnx
ky = 3*mny
Mx = sparse(iix[1:kx], jjx[1:kx], sx[1:kx], (Nr+2)*(Nz+2), (Nr+2)*(Nz+2))
My = sparse(iiy[1:ky], jjy[1:ky], sy[1:ky], (Nr+2)*(Nz+2), (Nr+2)*(Nz+2))
M = Mx + My
(M, Mx, My)
end



# ================ 2D Convection Term Cylindrical Upwind ==============================
function convectionUpwindTermCylindrical2D(u::FaceValue)
# u is a face variable
# extract data from the mesh structure
Nr = u.domain.dims[1]
Nz = u.domain.dims[2]
G=reshape([1:(Nr+2)*(Nz+2);], Nr+2, Nz+2)
DRp = u.domain.cellsize.x[2:end-1]
DZp = zeros( 1, Nz)
DZp[:] .= u.domain.cellsize.y[2:end-1]
rp = repeat(u.domain.cellcenters.x, 1, Nz)
rf = repeat(u.domain.facecenters.x, 1, Nz)

# define the vectors to store the sparse matrix data
iix = zeros(Int64, 3*(Nr+2)*(Nz+2))
iiy = zeros(Int64, 3*(Nr+2)*(Nz+2))
jjx = zeros(Int64, 3*(Nr+2)*(Nz+2))
jjy = zeros(Int64, 3*(Nr+2)*(Nz+2))
sx = zeros(Float64, 3*(Nr+2)*(Nz+2))
sy = zeros(Float64, 3*(Nr+2)*(Nz+2))
mnx = Nr*Nz
mny = Nr*Nz

re = rf[2:Nr+1,:]
rw = rf[1:Nr,:]

# find the velocity direction for the upwind scheme
ue_min = min.(u.xvalue[2:Nr+1,:],0.0)
ue_max = max.(u.xvalue[2:Nr+1,:],0.0)
uw_min = min.(u.xvalue[1:Nr,:],0.0)
uw_max = max.(u.xvalue[1:Nr,:],0.0)
vn_min = min.(u.yvalue[:,2:Nz+1],0.0)
vn_max = max.(u.yvalue[:,2:Nz+1],0.0)
vs_min = min.(u.yvalue[:,1:Nz],0.0)
vs_max = max.(u.yvalue[:,1:Nz],0.0)

# calculate the coefficients for the internal cells, do not reshape yet
AE = re.*ue_min./(DRp.*rp)
AW = -rw.*uw_max./(DRp.*rp)
AN = vn_min./DZp
AS = -vs_max./DZp
APx = (re.*ue_max.-rw.*uw_min)./(DRp.*rp)
APy = (vn_max.-vs_min)./DZp

# Also correct for the boundary cells (not the ghost cells)
# Left boundary:
APx[1,:] .= APx[1,:].-rw[1,:].*uw_max[1,:]./(2.0*DRp[1]*rp[1,:])
AW[1,:] .= AW[1,:]/2.0
# Right boundary:
AE[end,:] .= AE[end,:]/2.0
APx[end,:] .= APx[end,:].+re[end,:].*ue_min[end,:]./(2.0*DRp[end]*rp[end,:])
# Bottom boundary:
APy[:,1] .= APy[:,1].-vs_max[:,1]/(2.0*DZp[1])
AS[:,1] .= AS[:,1]/2.0
# Top boundary:
AN[:,end] .= AN[:,end]/2.0
APy[:,end] .= APy[:,end].+vn_min[:,end]/(2.0*DZp[end])

# now reshape
AE = reshape(AE,mnx)
AW = reshape(AW,mnx)
AN = reshape(AN,mny)
AS = reshape(AS,mny)
APx = reshape(APx,mnx)
APy = reshape(APy,mny)

# build the sparse matrix based on the numbering system
rowx_index = reshape(G[2:Nr+1,2:Nz+1],mnx) # main diagonal x
iix[1:3*mnx] .= repeat(rowx_index,3)
rowy_index = reshape(G[2:Nr+1,2:Nz+1],mny) # main diagonal y
iiy[1:3*mny] .= repeat(rowy_index,3)
jjx[1:3*mnx] .= [reshape(G[1:Nr,2:Nz+1],mnx); reshape(G[2:Nr+1,2:Nz+1],mnx); reshape(G[3:Nr+2,2:Nz+1],mnx)]
jjy[1:3*mny] .= [reshape(G[2:Nr+1,1:Nz],mny); reshape(G[2:Nr+1,2:Nz+1],mny); reshape(G[2:Nr+1,3:Nz+2],mny)]
sx[1:3*mnx] .= [AW; APx; AE]
sy[1:3*mny] .= [AS; APy; AN]

# build the sparse matrix
kx = 3*mnx
ky = 3*mny
Mx = sparse(iix[1:kx], jjx[1:kx], sx[1:kx], (Nr+2)*(Nz+2), (Nr+2)*(Nz+2))
My = sparse(iiy[1:ky], jjy[1:ky], sy[1:ky], (Nr+2)*(Nz+2), (Nr+2)*(Nz+2))
M = Mx + My
(M, Mx, My)
end

function convectionUpwindTermCylindrical2D(u::FaceValue, u_upwind::FaceValue)
# u is a face variable
# extract data from the mesh structure
Nr = u.domain.dims[1]
Nz = u.domain.dims[2]
G=reshape([1:(Nr+2)*(Nz+2);], Nr+2, Nz+2)
DRp = u.domain.cellsize.x[2:end-1]
DZp = zeros( 1, Nz)
DZp[:] .= u.domain.cellsize.y[2:end-1]
rp = repeat(u.domain.cellcenters.x, 1, Nz)
rf = repeat(u.domain.facecenters.x, 1, Nz)

# define the vectors to store the sparse matrix data
iix = zeros(Int64, 3*(Nr+2)*(Nz+2))
iiy = zeros(Int64, 3*(Nr+2)*(Nz+2))
jjx = zeros(Int64, 3*(Nr+2)*(Nz+2))
jjy = zeros(Int64, 3*(Nr+2)*(Nz+2))
sx = zeros(Float64, 3*(Nr+2)*(Nz+2))
sy = zeros(Float64, 3*(Nr+2)*(Nz+2))
mnx = Nr*Nz
mny = Nr*Nz

re = rf[2:Nr+1,:]
rw = rf[1:Nr,:]

# find the velocity direction for the upwind scheme
ue_min = u.xvalue[2:Nr+1,:]
ue_max = u.xvalue[2:Nr+1,:]
uw_min = u.xvalue[1:Nr,:]
uw_max = u.xvalue[1:Nr,:]
vn_min = u.yvalue[:,2:Nz+1]
vn_max = u.yvalue[:,2:Nz+1]
vs_min = u.yvalue[:,1:Nz]
vs_max = u.yvalue[:,1:Nz]

ue_min[u_upwind.xvalue[2:Nr+1,:].>0.0] .= 0.0
ue_max[u_upwind.xvalue[2:Nr+1,:].<0.0] .= 0.0
uw_min[u_upwind.xvalue[1:Nr,:].>0.0] .= 0.0
uw_max[u_upwind.xvalue[1:Nr,:].<0.0] .= 0.0
vn_min[u_upwind.yvalue[:,2:Nz+1].>0.0] .= 0.0
vn_max[u_upwind.yvalue[:,2:Nz+1].<0.0] .= 0.0
vs_min[u_upwind.yvalue[:,1:Nz].>0.0] .= 0.0
vs_max[u_upwind.yvalue[:,1:Nz].<0.0] .= 0.0

# calculate the coefficients for the internal cells, do not reshape yet
AE = re.*ue_min./(DRp.*rp)
AW = -rw.*uw_max./(DRp.*rp)
AN = vn_min./DZp
AS = -vs_max./DZp
APx = (re.*ue_max.-rw.*uw_min)./(DRp.*rp)
APy = (vn_max.-vs_min)./DZp

# Also correct for the boundary cells (not the ghost cells)
# Left boundary:
APx[1,:] .= APx[1,:].-rw[1,:].*uw_max[1,:]./(2.0*DRp[1]*rp[1,:])
AW[1,:] .= AW[1,:]/2.0
# Right boundary:
AE[end,:] .= AE[end,:]/2.0
APx[end,:] .= APx[end,:].+re[end,:].*ue_min[end,:]./(2.0*DRp[end]*rp[end,:])
# Bottom boundary:
APy[:,1] .= APy[:,1].-vs_max[:,1]/(2.0*DZp[1])
AS[:,1] .= AS[:,1]/2.0
# Top boundary:
AN[:,end] .= AN[:,end]/2.0
APy[:,end] .= APy[:,end].+vn_min[:,end]/(2.0*DZp[end])

# now reshape
AE = reshape(AE,mnx)
AW = reshape(AW,mnx)
AN = reshape(AN,mny)
AS = reshape(AS,mny)
APx = reshape(APx,mnx)
APy = reshape(APy,mny)

# build the sparse matrix based on the numbering system
rowx_index = reshape(G[2:Nr+1,2:Nz+1],mnx) # main diagonal x
iix[1:3*mnx] .= repeat(rowx_index,3)
rowy_index = reshape(G[2:Nr+1,2:Nz+1],mny) # main diagonal y
iiy[1:3*mny] .= repeat(rowy_index,3)
jjx[1:3*mnx] .= [reshape(G[1:Nr,2:Nz+1],mnx); reshape(G[2:Nr+1,2:Nz+1],mnx); reshape(G[3:Nr+2,2:Nz+1],mnx)]
jjy[1:3*mny] .= [reshape(G[2:Nr+1,1:Nz],mny); reshape(G[2:Nr+1,2:Nz+1],mny); reshape(G[2:Nr+1,3:Nz+2],mny)]
sx[1:3*mnx] .= [AW; APx; AE]
sy[1:3*mny] .= [AS; APy; AN]

# build the sparse matrix
kx = 3*mnx
ky = 3*mny
Mx = sparse(iix[1:kx], jjx[1:kx], sx[1:kx], (Nr+2)*(Nz+2), (Nr+2)*(Nz+2))
My = sparse(iiy[1:ky], jjy[1:ky], sy[1:ky], (Nr+2)*(Nz+2), (Nr+2)*(Nz+2))
M = Mx + My
(M, Mx, My)
end

# ==================== 2D Convection Term Cylindrical TVD ===========================
function convectionTvdTermCylindrical2D(u::FaceValue, phi::CellValue, FL::Function)
# u is a face variable
# phi is a cell variable

# a function to avoid division by zero
eps1 = 1.0e-20
fsign(phi_in) = (abs.(phi_in).>=eps1).*phi_in+eps1*(phi_in.==0.0)+eps1*(abs.(phi_in).<eps1).*sign.(phi_in)

# extract data from the mesh structure
Nr = u.domain.dims[1]
Nz = u.domain.dims[2]
G=reshape([1:(Nr+2)*(Nz+2);], Nr+2, Nz+2)
DRp = u.domain.cellsize.x[2:end-1]
DZp = zeros( 1, Nz)
DZp[:] .= u.domain.cellsize.y[2:end-1]
dr=0.5.*(u.domain.cellsize.x[1:end-1].+u.domain.cellsize.x[2:end])
dz=zeros( 1, Nz+1)
dz[:] .=0.5.*(u.domain.cellsize.y[1:end-1].+u.domain.cellsize.y[2:end])
rp = repeat(u.domain.cellcenters.x, 1, Nz)
rf = repeat(u.domain.facecenters.x, 1, Nz)
psiX_p = zeros(Nr+1,Nz)
psiX_m = zeros(Nr+1,Nz)
psiY_p = zeros(Nr,Nz+1)
psiY_m = zeros(Nr,Nz+1)

# define the vectors to store the sparse matrix data
iix = zeros(Int64, 3*(Nr+2)*(Nz+2))
iiy = zeros(Int64, 3*(Nr+2)*(Nz+2))
jjx = zeros(Int64, 3*(Nr+2)*(Nz+2))
jjy = zeros(Int64, 3*(Nr+2)*(Nz+2))
sx = zeros(Float64, 3*(Nr+2)*(Nz+2))
sy = zeros(Float64, 3*(Nr+2)*(Nz+2))
mnx = Nr*Nz
mny = Nr*Nz

# calculate the upstream to downstream gradient ratios for u>0 (+ ratio)
# x direction
dphiX_p = (phi.value[2:Nr+2, 2:Nz+1].-phi.value[1:Nr+1, 2:Nz+1])./dr
rX_p = dphiX_p[1:end-1,:]./fsign(dphiX_p[2:end,:])
psiX_p[2:Nr+1,:] .= 0.5.*FL.(rX_p).*(phi.value[3:Nr+2,2:Nz+1].-
		    phi.value[2:Nr+1, 2:Nz+1])
psiX_p[1, :] .= 0.0 # left boundary will be handled in the main matrix
# y direction
dphiY_p = (phi.value[2:Nr+1, 2:Nz+2].-phi.value[2:Nr+1, 1:Nz+1])./dz
rY_p = dphiY_p[:,1:end-1]./fsign(dphiY_p[:,2:end])
psiY_p[:,2:Nz+1] .= 0.5.*FL.(rY_p).*(phi.value[2:Nr+1,3:Nz+2].-
		  phi.value[2:Nr+1, 2:Nz+1])
psiY_p[:,1] .= 0.0 # Bottom boundary will be handled in the main matrix

# calculate the upstream to downstream gradient ratios for u<0 (- ratio)
# x direction
rX_m = dphiX_p[2:end,:]./fsign(dphiX_p[1:end-1,:])
psiX_m[1:Nr,:] .= 0.5.*FL.(rX_m).*(phi.value[1:Nr, 2:Nz+1].-
		phi.value[2:Nr+1, 2:Nz+1])
psiX_m[Nr+1,:] .= 0.0 # right boundary
# y direction
rY_m = dphiY_p[:,2:end]./fsign(dphiY_p[:,1:end-1])
psiY_m[:,1:Nz] .= 0.5.*FL.(rY_m).*(phi.value[2:Nr+1, 1:Nz].-
	      phi.value[2:Nr+1, 2:Nz+1])
psiY_m[:, Nz+1] .= 0.0 # top boundary will be handled in the main matrix

re = rf[2:Nr+1,:]
rw = rf[1:Nr,:]

# find the velocity direction for the upwind scheme
ue_min = min.(u.xvalue[2:Nr+1,:],0.0)
ue_max = max.(u.xvalue[2:Nr+1,:],0.0)
uw_min = min.(u.xvalue[1:Nr,:],0.0)
uw_max = max.(u.xvalue[1:Nr,:],0.0)
vn_min = min.(u.yvalue[:,2:Nz+1],0.0)
vn_max = max.(u.yvalue[:,2:Nz+1],0.0)
vs_min = min.(u.yvalue[:,1:Nz],0.0)
vs_max = max.(u.yvalue[:,1:Nz],0.0)

# calculate the coefficients for the internal cells, do not reshape yet
AE = re.*ue_min./(DRp.*rp)
AW = -rw.*uw_max./(DRp.*rp)
AN = vn_min./DZp
AS = -vs_max./DZp
APx = (re.*ue_max.-rw.*uw_min)./(DRp.*rp)
APy = (vn_max.-vs_min)./DZp

# Also correct for the boundary cells (not the ghost cells)
# Left boundary:
APx[1,:] .= APx[1,:].-rw[1,:].*uw_max[1,:]./(2.0*DRp[1]*rp[1,:])
AW[1,:] .= AW[1,:]/2.0
# Right boundary:
AE[end,:] .= AE[end,:]/2.0
APx[end,:] .= APx[end,:].+re[end,:].*ue_min[end,:]./(2.0*DRp[end]*rp[end,:])
# Bottom boundary:
APy[:,1] .= APy[:,1].-vs_max[:,1]/(2.0*DZp[1])
AS[:,1] .= AS[:,1]/2.0
# Top boundary:
AN[:,end] .= AN[:,end]/2.0
APy[:,end] .= APy[:,end].+vn_min[:,end]/(2.0*DZp[end])

# now reshape
AE = reshape(AE,mnx)
AW = reshape(AW,mnx)
AN = reshape(AN,mny)
AS = reshape(AS,mny)
APx = reshape(APx,mnx)
APy = reshape(APy,mny)

# build the sparse matrix based on the numbering system
rowx_index = reshape(G[2:Nr+1,2:Nz+1],mnx) # main diagonal x
iix[1:3*mnx] .= repeat(rowx_index,3)
rowy_index = reshape(G[2:Nr+1,2:Nz+1],mny) # main diagonal y
iiy[1:3*mny] .= repeat(rowy_index,3)
jjx[1:3*mnx] .= [reshape(G[1:Nr,2:Nz+1],mnx); reshape(G[2:Nr+1,2:Nz+1],mnx); reshape(G[3:Nr+2,2:Nz+1],mnx)]
jjy[1:3*mny] .= [reshape(G[2:Nr+1,1:Nz],mny); reshape(G[2:Nr+1,2:Nz+1],mny); reshape(G[2:Nr+1,3:Nz+2],mny)]
sx[1:3*mnx] .= [AW; APx; AE]
sy[1:3*mny] .= [AS; APy; AN]

# calculate the TVD correction term
div_x = -(1.0./(DRp.*rp)).*(re.*(ue_max.*psiX_p[2:Nr+1,:].+ue_min.*psiX_m[2:Nr+1,:])-
              rw.*(uw_max.*psiX_p[1:Nr,:].+uw_min.*psiX_m[1:Nr,:]))
div_y = -(1.0./DZp).*((vn_max.*psiY_p[:,2:Nz+1].+vn_min.*psiY_m[:,2:Nz+1])-
              (vs_max.*psiY_p[:,1:Nz].+vs_min.*psiY_m[:,1:Nz]))

# define the RHS Vector
RHS = zeros((Nr+2)*(Nz+2))
RHSx = zeros((Nr+2)*(Nz+2))
RHSy = zeros((Nr+2)*(Nz+2))

# assign the values of the RHS vector
RHS[rowx_index] .= reshape(div_x+div_y,Nr*Nz)
RHSx[rowx_index] .= reshape(div_x,Nr*Nz)
RHSy[rowy_index] .= reshape(div_y,Nr*Nz)

# build the sparse matrix
kx = 3*mnx
ky = 3*mny
Mx = sparse(iix[1:kx], jjx[1:kx], sx[1:kx], (Nr+2)*(Nz+2), (Nr+2)*(Nz+2))
My = sparse(iiy[1:ky], jjy[1:ky], sy[1:ky], (Nr+2)*(Nz+2), (Nr+2)*(Nz+2))
M = Mx + My
(M, RHS, Mx, My, RHSx, RHSy)
end


# ==================== 2D Convection RHS Cylindrical TVD ===========================
function convectionTvdRHSCylindrical2D(u::FaceValue, phi::CellValue, FL::Function)
# u is a face variable
# phi is a cell variable

# a function to avoid division by zero
eps1 = 1.0e-20
fsign(phi_in) = (abs.(phi_in).>=eps1).*phi_in+eps1*(phi_in.==0.0)+eps1*(abs.(phi_in).<eps1).*sign.(phi_in)

# extract data from the mesh structure
Nr = u.domain.dims[1]
Nz = u.domain.dims[2]
G=reshape([1:(Nr+2)*(Nz+2);], Nr+2, Nz+2)
DRp = u.domain.cellsize.x[2:end-1]
DZp = zeros( 1, Nz)
DZp[:] .= u.domain.cellsize.y[2:end-1]
dr=0.5.*(u.domain.cellsize.x[1:end-1].+u.domain.cellsize.x[2:end])
dz=zeros( 1, Nz+1)
dz[:] .=0.5.*(u.domain.cellsize.y[1:end-1].+u.domain.cellsize.y[2:end])
rp = repeat(u.domain.cellcenters.x, 1, Nz)
rf = repeat(u.domain.facecenters.x, 1, Nz)
psiX_p = zeros(Nr+1,Nz)
psiX_m = zeros(Nr+1,Nz)
psiY_p = zeros(Nr,Nz+1)
psiY_m = zeros(Nr,Nz+1)

# calculate the upstream to downstream gradient ratios for u>0 (+ ratio)
# x direction
dphiX_p = (phi.value[2:Nr+2, 2:Nz+1].-phi.value[1:Nr+1, 2:Nz+1])./dr
rX_p = dphiX_p[1:end-1,:]./fsign(dphiX_p[2:end,:])
psiX_p[2:Nr+1,:] .= 0.5.*FL.(rX_p).*(phi.value[3:Nr+2,2:Nz+1].-
		    phi.value[2:Nr+1, 2:Nz+1])
psiX_p[1, :] .= 0.0 # left boundary will be handled in the main matrix
# y direction
dphiY_p = (phi.value[2:Nr+1, 2:Nz+2].-phi.value[2:Nr+1, 1:Nz+1])./dz
rY_p = dphiY_p[:,1:end-1]./fsign(dphiY_p[:,2:end])
psiY_p[:,2:Nz+1] .= 0.5.*FL.(rY_p).*(phi.value[2:Nr+1,3:Nz+2].-
		  phi.value[2:Nr+1, 2:Nz+1])
psiY_p[:,1] .= 0.0 # Bottom boundary will be handled in the main matrix

# calculate the upstream to downstream gradient ratios for u<0 (- ratio)
# x direction
rX_m = dphiX_p[2:end,:]./fsign(dphiX_p[1:end-1,:])
psiX_m[1:Nr,:] .= 0.5.*FL.(rX_m).*(phi.value[1:Nr, 2:Nz+1].-
		phi.value[2:Nr+1, 2:Nz+1])
psiX_m[Nr+1,:] .= 0.0 # right boundary
# y direction
rY_m = dphiY_p[:,2:end]./fsign(dphiY_p[:,1:end-1])
psiY_m[:,1:Nz] .= 0.5.*FL.(rY_m).*(phi.value[2:Nr+1, 1:Nz].-
	      phi.value[2:Nr+1, 2:Nz+1])
psiY_m[:, Nz+1] .= 0.0 # top boundary will be handled in the main matrix

re = rf[2:Nr+1,:]
rw = rf[1:Nr,:]

# find the velocity direction for the upwind scheme
ue_min = min.(u.xvalue[2:Nr+1,:],0.0)
ue_max = max.(u.xvalue[2:Nr+1,:],0.0)
uw_min = min.(u.xvalue[1:Nr,:],0.0)
uw_max = max.(u.xvalue[1:Nr,:],0.0)
vn_min = min.(u.yvalue[:,2:Nz+1],0.0)
vn_max = max.(u.yvalue[:,2:Nz+1],0.0)
vs_min = min.(u.yvalue[:,1:Nz],0.0)
vs_max = max.(u.yvalue[:,1:Nz],0.0)

# calculate the TVD correction term
div_x = -(1.0./(DRp.*rp)).*(re.*(ue_max.*psiX_p[2:Nr+1,:].+ue_min.*psiX_m[2:Nr+1,:])-
              rw.*(uw_max.*psiX_p[1:Nr,:].+uw_min.*psiX_m[1:Nr,:]))
div_y = -(1.0./DZp).*((vn_max.*psiY_p[:,2:Nz+1].+vn_min.*psiY_m[:,2:Nz+1])-
              (vs_max.*psiY_p[:,1:Nz].+vs_min.*psiY_m[:,1:Nz]))

# define the RHS Vector
RHS = zeros((Nr+2)*(Nz+2))
RHSx = zeros((Nr+2)*(Nz+2))
RHSy = zeros((Nr+2)*(Nz+2))

# assign the values of the RHS vector
mnx = Nr*Nz
mny = Nr*Nz
rowx_index = reshape(G[2:Nr+1,2:Nz+1],mnx) # main diagonal x
rowy_index = reshape(G[2:Nr+1,2:Nz+1],mny) # main diagonal y
RHS[rowx_index] .= reshape(div_x+div_y,Nr*Nz)
RHSx[rowx_index] .= reshape(div_x,Nr*Nz)
RHSy[rowy_index] .= reshape(div_y,Nr*Nz)

(RHS, RHSx, RHSy)
end

function convectionTvdRHSCylindrical2D(u::FaceValue, phi::CellValue, FL::Function, u_upwind::FaceValue)
# u is a face variable
# phi is a cell variable

# a function to avoid division by zero
eps1 = 1.0e-20
fsign(phi_in) = (abs.(phi_in).>=eps1).*phi_in+eps1*(phi_in.==0.0)+eps1*(abs.(phi_in).<eps1).*sign.(phi_in)

# extract data from the mesh structure
Nr = u.domain.dims[1]
Nz = u.domain.dims[2]
G=reshape([1:(Nr+2)*(Nz+2);], Nr+2, Nz+2)
DRp = u.domain.cellsize.x[2:end-1]
DZp = zeros( 1, Nz)
DZp[:] .= u.domain.cellsize.y[2:end-1]
dr=0.5.*(u.domain.cellsize.x[1:end-1].+u.domain.cellsize.x[2:end])
dz=zeros( 1, Nz+1)
dz[:] .=0.5.*(u.domain.cellsize.y[1:end-1].+u.domain.cellsize.y[2:end])
rp = repeat(u.domain.cellcenters.x, 1, Nz)
rf = repeat(u.domain.facecenters.x, 1, Nz)
psiX_p = zeros(Nr+1,Nz)
psiX_m = zeros(Nr+1,Nz)
psiY_p = zeros(Nr,Nz+1)
psiY_m = zeros(Nr,Nz+1)

# calculate the upstream to downstream gradient ratios for u>0 (+ ratio)
# x direction
dphiX_p = (phi.value[2:Nr+2, 2:Nz+1].-phi.value[1:Nr+1, 2:Nz+1])./dr
rX_p = dphiX_p[1:end-1,:]./fsign(dphiX_p[2:end,:])
psiX_p[2:Nr+1,:] .= 0.5.*FL.(rX_p).*(phi.value[3:Nr+2,2:Nz+1].-
		    phi.value[2:Nr+1, 2:Nz+1])
psiX_p[1, :] .= 0.0 # left boundary will be handled in the main matrix
# y direction
dphiY_p = (phi.value[2:Nr+1, 2:Nz+2].-phi.value[2:Nr+1, 1:Nz+1])./dz
rY_p = dphiY_p[:,1:end-1]./fsign(dphiY_p[:,2:end])
psiY_p[:,2:Nz+1] .= 0.5.*FL.(rY_p).*(phi.value[2:Nr+1,3:Nz+2].-
		  phi.value[2:Nr+1, 2:Nz+1])
psiY_p[:,1] .= 0.0 # Bottom boundary will be handled in the main matrix

# calculate the upstream to downstream gradient ratios for u<0 (- ratio)
# x direction
rX_m = dphiX_p[2:end,:]./fsign(dphiX_p[1:end-1,:])
psiX_m[1:Nr,:] .= 0.5.*FL.(rX_m).*(phi.value[1:Nr, 2:Nz+1].-
		phi.value[2:Nr+1, 2:Nz+1])
psiX_m[Nr+1,:] .= 0.0 # right boundary
# y direction
rY_m = dphiY_p[:,2:end]./fsign(dphiY_p[:,1:end-1])
psiY_m[:,1:Nz] .= 0.5.*FL.(rY_m).*(phi.value[2:Nr+1, 1:Nz].-
	      phi.value[2:Nr+1, 2:Nz+1])
psiY_m[:, Nz+1] .= 0.0 # top boundary will be handled in the main matrix

re = rf[2:Nr+1,:]
rw = rf[1:Nr,:]

# find the velocity direction for the upwind scheme
ue_min = u.xvalue[2:Nr+1,:]
ue_max = u.xvalue[2:Nr+1,:]
uw_min = u.xvalue[1:Nr,:]
uw_max = u.xvalue[1:Nr,:]
vn_min = u.yvalue[:,2:Nz+1]
vn_max = u.yvalue[:,2:Nz+1]
vs_min = u.yvalue[:,1:Nz]
vs_max = u.yvalue[:,1:Nz]

ue_min[u_upwind.xvalue[2:Nr+1,:].>0.0] .= 0.0
ue_max[u_upwind.xvalue[2:Nr+1,:].<0.0] .= 0.0
uw_min[u_upwind.xvalue[1:Nr,:].>0.0] .= 0.0
uw_max[u_upwind.xvalue[1:Nr,:].<0.0] .= 0.0
vn_min[u_upwind.yvalue[:,2:Nz+1].>0.0] .= 0.0
vn_max[u_upwind.yvalue[:,2:Nz+1].<0.0] .= 0.0
vs_min[u_upwind.yvalue[:,1:Nz].>0.0] .= 0.0
vs_max[u_upwind.yvalue[:,1:Nz].<0.0] .= 0.0

# calculate the TVD correction term
div_x = -(1.0./(DRp.*rp)).*(re.*(ue_max.*psiX_p[2:Nr+1,:].+ue_min.*psiX_m[2:Nr+1,:])-
              rw.*(uw_max.*psiX_p[1:Nr,:].+uw_min.*psiX_m[1:Nr,:]))
div_y = -(1.0./DZp).*((vn_max.*psiY_p[:,2:Nz+1].+vn_min.*psiY_m[:,2:Nz+1])-
              (vs_max.*psiY_p[:,1:Nz].+vs_min.*psiY_m[:,1:Nz]))

# define the RHS Vector
RHS = zeros((Nr+2)*(Nz+2))
RHSx = zeros((Nr+2)*(Nz+2))
RHSy = zeros((Nr+2)*(Nz+2))

# assign the values of the RHS vector
mnx = Nr*Nz
mny = Nr*Nz
rowx_index = reshape(G[2:Nr+1,2:Nz+1],mnx) # main diagonal x
rowy_index = reshape(G[2:Nr+1,2:Nz+1],mny) # main diagonal y
RHS[rowx_index] .= reshape(div_x+div_y,Nr*Nz)
RHSx[rowx_index] .= reshape(div_x,Nr*Nz)
RHSy[rowy_index] .= reshape(div_y,Nr*Nz)

(RHS, RHSx, RHSy)
end

# ================ 2D Convection Term Radial Central ===================
function convectionTermRadial2D(u::FaceValue)
# u is a face variable
# extract data from the mesh structure
Nr = u.domain.dims[1]
Ntheta = u.domain.dims[2]
G=reshape([1:(Nr+2)*(Ntheta+2);], Nr+2, Ntheta+2)
DRe = u.domain.cellsize.x[3:end]
DRw = u.domain.cellsize.x[1:end-2]
DRp = u.domain.cellsize.x[2:end-1]
DTHETA = zeros( 1, Ntheta+2)
DTHETA[:] .= u.domain.cellsize.y
DTHETAn = DTHETA[:,3:end]
DTHETAs = DTHETA[:,1:end-2]
DTHETAp = DTHETA[:,2:end-1]
rp = u.domain.cellcenters.x
rf = u.domain.facecenters.x

# define the vectors to store the sparse matrix data
iix = zeros(Int64, 3*(Nr+2)*(Ntheta+2))
iiy = zeros(Int64, 3*(Nr+2)*(Ntheta+2))
jjx = zeros(Int64, 3*(Nr+2)*(Ntheta+2))
jjy = zeros(Int64, 3*(Nr+2)*(Ntheta+2))
sx = zeros(Float64, 3*(Nr+2)*(Ntheta+2))
sy = zeros(Float64, 3*(Nr+2)*(Ntheta+2))
mnx = Nr*Ntheta
mny = Nr*Ntheta

re = rf[2:Nr+1,:]
rw = rf[1:Nr,:]

# reassign the east, west for code readability
ue = re.*u.xvalue[2:Nr+1,:]./(rp.*(DRp+DRe))
uw = rw.*u.xvalue[1:Nr,:]./(rp.*(DRp+DRw))
vn = u.yvalue[:,2:Ntheta+1]./(rp.*(DTHETAp+DTHETAn))
vs = u.yvalue[:,1:Ntheta]./(rp.*(DTHETAp+DTHETAs))

# calculate the coefficients for the internal cells
AE = reshape(ue,mnx)
AW = reshape(-uw,mnx)
AN = reshape(vn,mny)
AS = reshape(-vs,mny)
APx = reshape((ue.*DRe-uw.*DRw)./DRp,mnx)
APy = reshape((vn.*DTHETAn-vs.*DTHETAs)./DTHETAp,mny)

# build the sparse matrix based on the numbering system
rowx_index = reshape(G[2:Nr+1,2:Ntheta+1],mnx) # main diagonal x
iix[1:3*mnx] .= repeat(rowx_index,3)
rowy_index = reshape(G[2:Nr+1,2:Ntheta+1],mny) # main diagonal y
iiy[1:3*mny] .= repeat(rowy_index,3)
jjx[1:3*mnx] .= [reshape(G[1:Nr,2:Ntheta+1],mnx); reshape(G[2:Nr+1,2:Ntheta+1],mnx); reshape(G[3:Nr+2,2:Ntheta+1],mnx)]
jjy[1:3*mny] .= [reshape(G[2:Nr+1,1:Ntheta],mny); reshape(G[2:Nr+1,2:Ntheta+1],mny); reshape(G[2:Nr+1,3:Ntheta+2],mny)]
sx[1:3*mnx] .= [AW; APx; AE]
sy[1:3*mny] .= [AS; APy; AN]

# build the sparse matrix
kx = 3*mnx
ky = 3*mny
Mx = sparse(iix[1:kx], jjx[1:kx], sx[1:kx], (Nr+2)*(Ntheta+2), (Nr+2)*(Ntheta+2))
My = sparse(iiy[1:ky], jjy[1:ky], sy[1:ky], (Nr+2)*(Ntheta+2), (Nr+2)*(Ntheta+2))
M = Mx + My
(M, Mx, My)
end



# ================ 2D Convection Term Radial Upwind ==============================
function convectionUpwindTermRadial2D(u::FaceValue)
# u is a face variable
# extract data from the mesh structure
Nr = u.domain.dims[1]
Ntheta = u.domain.dims[2]
G=reshape([1:(Nr+2)*(Ntheta+2);], Nr+2, Ntheta+2)
DRp = u.domain.cellsize.x[2:end-1]
DTHETAp = zeros( 1, Ntheta)
DTHETAp[:] .= u.domain.cellsize.y[2:end-1]
rp = u.domain.cellcenters.x
rf = u.domain.facecenters.x

# define the vectors to store the sparse matrix data
iix = zeros(Int64, 3*(Nr+2)*(Ntheta+2))
iiy = zeros(Int64, 3*(Nr+2)*(Ntheta+2))
jjx = zeros(Int64, 3*(Nr+2)*(Ntheta+2))
jjy = zeros(Int64, 3*(Nr+2)*(Ntheta+2))
sx = zeros(Float64, 3*(Nr+2)*(Ntheta+2))
sy = zeros(Float64, 3*(Nr+2)*(Ntheta+2))
mnx = Nr*Ntheta
mny = Nr*Ntheta

re = rf[2:Nr+1,:]
rw = rf[1:Nr,:]

# find the velocity direction for the upwind scheme
ue_min = min.(u.xvalue[2:Nr+1,:],0.0)
ue_max = max.(u.xvalue[2:Nr+1,:],0.0)
uw_min = min.(u.xvalue[1:Nr,:],0.0)
uw_max = max.(u.xvalue[1:Nr,:],0.0)
vn_min = min.(u.yvalue[:,2:Ntheta+1],0.0)
vn_max = max.(u.yvalue[:,2:Ntheta+1],0.0)
vs_min = min.(u.yvalue[:,1:Ntheta],0.0)
vs_max = max.(u.yvalue[:,1:Ntheta],0.0)

# calculate the coefficients for the internal cells, do not reshape yet
AE = re.*ue_min./(DRp.*rp)
AW = -rw.*uw_max./(DRp.*rp)
AN = vn_min./(DTHETAp.*rp)
AS = -vs_max./(DTHETAp.*rp)
APx = (re.*ue_max.-rw.*uw_min)./(DRp.*rp)
APy = (vn_max.-vs_min)./(DTHETAp.*rp)

# Also correct for the boundary cells (not the ghost cells)
# Left boundary:
APx[1,:] .= APx[1,:].-rw[1,:].*uw_max[1,:]./(2.0*DRp[1]*rp[1,:])
AW[1,:] .= AW[1,:]/2.0
# Right boundary:
AE[end,:] .= AE[end,:]/2.0
APx[end,:] .= APx[end,:].+re[end,:].*ue_min[end,:]./(2.0*DRp[end]*rp[end,:])
# Bottom boundary:
APy[:,1] .= APy[:,1].-vs_max[:,1]./(2.0*DTHETAp[1]*rp[:,1])
AS[:,1] .= AS[:,1]/2.0
# Top boundary:
AN[:,end] .= AN[:,end]/2.0
APy[:,end] .= APy[:,end].+vn_min[:,end]./(2.0*DTHETAp[end]*rp[:,end])

# now reshape
AE = reshape(AE,mnx)
AW = reshape(AW,mnx)
AN = reshape(AN,mny)
AS = reshape(AS,mny)
APx = reshape(APx,mnx)
APy = reshape(APy,mny)

# build the sparse matrix based on the numbering system
rowx_index = reshape(G[2:Nr+1,2:Ntheta+1],mnx) # main diagonal x
iix[1:3*mnx] .= repeat(rowx_index,3)
rowy_index = reshape(G[2:Nr+1,2:Ntheta+1],mny) # main diagonal y
iiy[1:3*mny] .= repeat(rowy_index,3)
jjx[1:3*mnx] .= [reshape(G[1:Nr,2:Ntheta+1],mnx); reshape(G[2:Nr+1,2:Ntheta+1],mnx); reshape(G[3:Nr+2,2:Ntheta+1],mnx)]
jjy[1:3*mny] .= [reshape(G[2:Nr+1,1:Ntheta],mny); reshape(G[2:Nr+1,2:Ntheta+1],mny); reshape(G[2:Nr+1,3:Ntheta+2],mny)]
sx[1:3*mnx] .= [AW; APx; AE]
sy[1:3*mny] .= [AS; APy; AN]

# build the sparse matrix
kx = 3*mnx
ky = 3*mny
Mx = sparse(iix[1:kx], jjx[1:kx], sx[1:kx], (Nr+2)*(Ntheta+2), (Nr+2)*(Ntheta+2))
My = sparse(iiy[1:ky], jjy[1:ky], sy[1:ky], (Nr+2)*(Ntheta+2), (Nr+2)*(Ntheta+2))
M = Mx + My
(M, Mx, My)
end

function convectionUpwindTermRadial2D(u::FaceValue, u_upwind::FaceValue)
# u is a face variable
# extract data from the mesh structure
Nr = u.domain.dims[1]
Ntheta = u.domain.dims[2]
G=reshape([1:(Nr+2)*(Ntheta+2);], Nr+2, Ntheta+2)
DRp = u.domain.cellsize.x[2:end-1]
DTHETAp = zeros( 1, Ntheta)
DTHETAp[:] .= u.domain.cellsize.y[2:end-1]
rp = u.domain.cellcenters.x
rf = u.domain.facecenters.x

# define the vectors to store the sparse matrix data
iix = zeros(Int64, 3*(Nr+2)*(Ntheta+2))
iiy = zeros(Int64, 3*(Nr+2)*(Ntheta+2))
jjx = zeros(Int64, 3*(Nr+2)*(Ntheta+2))
jjy = zeros(Int64, 3*(Nr+2)*(Ntheta+2))
sx = zeros(Float64, 3*(Nr+2)*(Ntheta+2))
sy = zeros(Float64, 3*(Nr+2)*(Ntheta+2))
mnx = Nr*Ntheta
mny = Nr*Ntheta

re = rf[2:Nr+1,:]
rw = rf[1:Nr,:]

# find the velocity direction for the upwind scheme
ue_min = u.xvalue[2:Nr+1,:]
ue_max = u.xvalue[2:Nr+1,:]
uw_min = u.xvalue[1:Nr,:]
uw_max = u.xvalue[1:Nr,:]
vn_min = u.yvalue[:,2:Ntheta+1]
vn_max = u.yvalue[:,2:Ntheta+1]
vs_min = u.yvalue[:,1:Ntheta]
vs_max = u.yvalue[:,1:Ntheta]

ue_min[u_upwind.xvalue[2:Nr+1,:].>0.0] .= 0.0
ue_max[u_upwind.xvalue[2:Nr+1,:].<0.0] .= 0.0
uw_min[u_upwind.xvalue[1:Nr,:].>0.0] .= 0.0
uw_max[u_upwind.xvalue[1:Nr,:].<0.0] .= 0.0
vn_min[u_upwind.yvalue[:,2:Ntheta+1].>0.0] .= 0.0
vn_max[u_upwind.yvalue[:,2:Ntheta+1].<0.0] .= 0.0
vs_min[u_upwind.yvalue[:,1:Ntheta].>0.0] .= 0.0
vs_max[u_upwind.yvalue[:,1:Ntheta].<0.0] .= 0.0

# calculate the coefficients for the internal cells, do not reshape yet
AE = re.*ue_min./(DRp.*rp)
AW = -rw.*uw_max./(DRp.*rp)
AN = vn_min./(DTHETAp.*rp)
AS = -vs_max./(DTHETAp.*rp)
APx = (re.*ue_max.-rw.*uw_min)./(DRp.*rp)
APy = (vn_max.-vs_min)./(DTHETAp.*rp)

# Also correct for the boundary cells (not the ghost cells)
# Left boundary:
APx[1,:] .= APx[1,:].-rw[1,:].*uw_max[1,:]./(2.0*DRp[1]*rp[1,:])
AW[1,:] .= AW[1,:]/2.0
# Right boundary:
AE[end,:] .= AE[end,:]/2.0
APx[end,:] .= APx[end,:].+re[end,:].*ue_min[end,:]./(2.0*DRp[end]*rp[end,:])
# Bottom boundary:
APy[:,1] .= APy[:,1].-vs_max[:,1]./(2.0*DTHETAp[1]*rp[:,1])
AS[:,1] .= AS[:,1]/2.0
# Top boundary:
AN[:,end] .= AN[:,end]/2.0
APy[:,end] .= APy[:,end].+vn_min[:,end]./(2.0*DTHETAp[end]*rp[:,end])

# now reshape
AE = reshape(AE,mnx)
AW = reshape(AW,mnx)
AN = reshape(AN,mny)
AS = reshape(AS,mny)
APx = reshape(APx,mnx)
APy = reshape(APy,mny)

# build the sparse matrix based on the numbering system
rowx_index = reshape(G[2:Nr+1,2:Ntheta+1],mnx) # main diagonal x
iix[1:3*mnx] .= repeat(rowx_index,3)
rowy_index = reshape(G[2:Nr+1,2:Ntheta+1],mny) # main diagonal y
iiy[1:3*mny] .= repeat(rowy_index,3)
jjx[1:3*mnx] .= [reshape(G[1:Nr,2:Ntheta+1],mnx); reshape(G[2:Nr+1,2:Ntheta+1],mnx); reshape(G[3:Nr+2,2:Ntheta+1],mnx)]
jjy[1:3*mny] .= [reshape(G[2:Nr+1,1:Ntheta],mny); reshape(G[2:Nr+1,2:Ntheta+1],mny); reshape(G[2:Nr+1,3:Ntheta+2],mny)]
sx[1:3*mnx] .= [AW; APx; AE]
sy[1:3*mny] .= [AS; APy; AN]

# build the sparse matrix
kx = 3*mnx
ky = 3*mny
Mx = sparse(iix[1:kx], jjx[1:kx], sx[1:kx], (Nr+2)*(Ntheta+2), (Nr+2)*(Ntheta+2))
My = sparse(iiy[1:ky], jjy[1:ky], sy[1:ky], (Nr+2)*(Ntheta+2), (Nr+2)*(Ntheta+2))
M = Mx + My
(M, Mx, My)
end


# ================ 2D Convection Term Radial TVD =======================
function convectionTvdTermRadial2D(u::FaceValue, phi::CellValue, FL::Function)
# u is a face variable
# phi is a cell variable

# a function to avoid division by zero
eps1 = 1.0e-20
fsign(phi_in) = (abs.(phi_in).>=eps1).*phi_in+eps1*(phi_in.==0.0)+eps1*(abs.(phi_in).<eps1).*sign.(phi_in)

# extract data from the mesh structure
Nr = u.domain.dims[1]
Ntheta = u.domain.dims[2]
G=reshape([1:(Nr+2)*(Ntheta+2);], Nr+2, Ntheta+2)
DRp = u.domain.cellsize.x[2:end-1]
DTHETAp = zeros( 1, Ntheta)
DTHETAp[:] .= u.domain.cellsize.y[2:end-1]
dr = 0.5.*(u.domain.cellsize.x[1:end-1].+u.domain.cellsize.x[2:end])
dtheta = zeros( 1, 1+Ntheta)
dtheta[:] .= 0.5.*(u.domain.cellsize.y[1:end-1].+u.domain.cellsize.y[2:end])
rp = u.domain.cellcenters.x
rf = u.domain.facecenters.x
psiX_p = zeros(Nr+1,Ntheta)
psiX_m = zeros(Nr+1,Ntheta)
psiY_p = zeros(Nr,Ntheta+1)
psiY_m = zeros(Nr,Ntheta+1)

# define the vectors to store the sparse matrix data
iix = zeros(Int64, 3*(Nr+2)*(Ntheta+2))
iiy = zeros(Int64, 3*(Nr+2)*(Ntheta+2))
jjx = zeros(Int64, 3*(Nr+2)*(Ntheta+2))
jjy = zeros(Int64, 3*(Nr+2)*(Ntheta+2))
sx = zeros(Float64, 3*(Nr+2)*(Ntheta+2))
sy = zeros(Float64, 3*(Nr+2)*(Ntheta+2))
mnx = Nr*Ntheta
mny = Nr*Ntheta

# calculate the upstream to downstream gradient ratios for u>0 (+ ratio)
# x direction
dphiX_p = (phi.value[2:Nr+2, 2:Ntheta+1].-phi.value[1:Nr+1, 2:Ntheta+1])./dr
rX_p = dphiX_p[1:end-1,:]./fsign(dphiX_p[2:end,:])
psiX_p[2:Nr+1,:] .= 0.5.*FL.(rX_p).*(phi.value[3:Nr+2,2:Ntheta+1].-
		    phi.value[2:Nr+1, 2:Ntheta+1])
psiX_p[1, :] .= 0.0 # left boundary will be handled in the main matrix
# y direction
dphiY_p = (phi.value[2:Nr+1, 2:Ntheta+2].-phi.value[2:Nr+1, 1:Ntheta+1])./dtheta
rY_p = dphiY_p[:,1:end-1]./fsign(dphiY_p[:,2:end])
psiY_p[:,2:Ntheta+1] .= 0.5.*FL.(rY_p).*(phi.value[2:Nr+1,3:Ntheta+2].-
		  phi.value[2:Nr+1, 2:Ntheta+1])
psiY_p[:,1] .= 0.0 # Bottom boundary will be handled in the main matrix

# calculate the upstream to downstream gradient ratios for u<0 (- ratio)
# x direction
rX_m = dphiX_p[2:end,:]./fsign(dphiX_p[1:end-1,:])
psiX_m[1:Nr,:] .= 0.5.*FL.(rX_m).*(phi.value[1:Nr, 2:Ntheta+1].-
		phi.value[2:Nr+1, 2:Ntheta+1])
psiX_m[Nr+1,:] .= 0.0 # right boundary
# y direction
rY_m = dphiY_p[:,2:end]./fsign(dphiY_p[:,1:end-1])
psiY_m[:,1:Ntheta] .= 0.5.*FL.(rY_m).*(phi.value[2:Nr+1, 1:Ntheta].-
	      phi.value[2:Nr+1, 2:Ntheta+1])
psiY_m[:, Ntheta+1] .= 0.0 # top boundary will be handled in the main matrix

re = rf[2:Nr+1,:]
rw = rf[1:Nr,:]

# find the velocity direction for the upwind scheme
ue_min = min.(u.xvalue[2:Nr+1,:],0.0)
ue_max = max.(u.xvalue[2:Nr+1,:],0.0)
uw_min = min.(u.xvalue[1:Nr,:],0.0)
uw_max = max.(u.xvalue[1:Nr,:],0.0)
vn_min = min.(u.yvalue[:,2:Ntheta+1],0.0)
vn_max = max.(u.yvalue[:,2:Ntheta+1],0.0)
vs_min = min.(u.yvalue[:,1:Ntheta],0.0)
vs_max = max.(u.yvalue[:,1:Ntheta],0.0)

# calculate the coefficients for the internal cells, do not reshape yet
AE = re.*ue_min./(DRp.*rp)
AW = -rw.*uw_max./(DRp.*rp)
AN = vn_min./(DTHETAp.*rp)
AS = -vs_max./(DTHETAp.*rp)
APx = (re.*ue_max.-rw.*uw_min)./(DRp.*rp)
APy = (vn_max.-vs_min)./(DTHETAp.*rp)

# Also correct for the boundary cells (not the ghost cells)
# Left boundary:
APx[1,:] .= APx[1,:].-rw[1,:].*uw_max[1,:]./(2.0*DRp[1]*rp[1,:])
AW[1,:] .= AW[1,:]/2.0
# Right boundary:
AE[end,:] .= AE[end,:]/2.0
APx[end,:] .= APx[end,:].+re[end,:].*ue_min[end,:]./(2.0*DRp[end]*rp[end,:])
# Bottom boundary:
APy[:,1] .= APy[:,1].-vs_max[:,1]./(2.0*DTHETAp[1]*rp[:,1])
AS[:,1] .= AS[:,1]/2.0
# Top boundary:
AN[:,end] .= AN[:,end]/2.0
APy[:,end] .= APy[:,end].+vn_min[:,end]./(2.0*DTHETAp[end]*rp[:,end])

# now reshape
AE = reshape(AE,mnx)
AW = reshape(AW,mnx)
AN = reshape(AN,mny)
AS = reshape(AS,mny)
APx = reshape(APx,mnx)
APy = reshape(APy,mny)

# build the sparse matrix based on the numbering system
rowx_index = reshape(G[2:Nr+1,2:Ntheta+1],mnx) # main diagonal x
iix[1:3*mnx] .= repeat(rowx_index,3)
rowy_index = reshape(G[2:Nr+1,2:Ntheta+1],mny) # main diagonal y
iiy[1:3*mny] .= repeat(rowy_index,3)
jjx[1:3*mnx] .= [reshape(G[1:Nr,2:Ntheta+1],mnx); reshape(G[2:Nr+1,2:Ntheta+1],mnx); reshape(G[3:Nr+2,2:Ntheta+1],mnx)]
jjy[1:3*mny] .= [reshape(G[2:Nr+1,1:Ntheta],mny); reshape(G[2:Nr+1,2:Ntheta+1],mny); reshape(G[2:Nr+1,3:Ntheta+2],mny)]
sx[1:3*mnx] .= [AW; APx; AE]
sy[1:3*mny] .= [AS; APy; AN]

# calculate the TVD correction term
div_x = -(1.0./(DRp.*rp)).*(re.*(ue_max.*psiX_p[2:Nr+1,:].+ue_min.*psiX_m[2:Nr+1,:])-
              rw.*(uw_max.*psiX_p[1:Nr,:].+uw_min.*psiX_m[1:Nr,:]))
div_y = -(1.0./(DTHETAp.*rp)).*((vn_max.*psiY_p[:,2:Ntheta+1].+vn_min.*psiY_m[:,2:Ntheta+1])-
              (vs_max.*psiY_p[:,1:Ntheta].+vs_min.*psiY_m[:,1:Ntheta]))

# define the RHS Vector
RHS = zeros((Nr+2)*(Ntheta+2))
RHSx = zeros((Nr+2)*(Ntheta+2))
RHSy = zeros((Nr+2)*(Ntheta+2))

# assign the values of the RHS vector
RHS[rowx_index] .= reshape(div_x+div_y,Nr*Ntheta)
RHSx[rowx_index] .= reshape(div_x,Nr*Ntheta)
RHSy[rowy_index] .= reshape(div_y,Nr*Ntheta)

# build the sparse matrix
kx = 3*mnx
ky = 3*mny
Mx = sparse(iix[1:kx], jjx[1:kx], sx[1:kx], (Nr+2)*(Ntheta+2), (Nr+2)*(Ntheta+2))
My = sparse(iiy[1:ky], jjy[1:ky], sy[1:ky], (Nr+2)*(Ntheta+2), (Nr+2)*(Ntheta+2))
M = Mx + My
(M, RHS, Mx, My, RHSx, RHSy)
end


# ================ 2D Convection RHS Radial TVD =======================
function convectionTvdRHSRadial2D(u::FaceValue, phi::CellValue, FL::Function)
# u is a face variable
# phi is a cell variable

# a function to avoid division by zero
eps1 = 1.0e-20
fsign(phi_in) = (abs.(phi_in).>=eps1).*phi_in+eps1*(phi_in.==0.0)+eps1*(abs.(phi_in).<eps1).*sign.(phi_in)

# extract data from the mesh structure
Nr = u.domain.dims[1]
Ntheta = u.domain.dims[2]
G=reshape([1:(Nr+2)*(Ntheta+2);], Nr+2, Ntheta+2)
DRp = u.domain.cellsize.x[2:end-1]
DTHETAp = zeros( 1, Ntheta)
DTHETAp[:] .= u.domain.cellsize.y[2:end-1]
dr = 0.5.*(u.domain.cellsize.x[1:end-1].+u.domain.cellsize.x[2:end])
dtheta = zeros( 1, 1+Ntheta)
dtheta[:] .= 0.5.*(u.domain.cellsize.y[1:end-1].+u.domain.cellsize.y[2:end])
rp = u.domain.cellcenters.x
rf = u.domain.facecenters.x
psiX_p = zeros(Nr+1,Ntheta)
psiX_m = zeros(Nr+1,Ntheta)
psiY_p = zeros(Nr,Ntheta+1)
psiY_m = zeros(Nr,Ntheta+1)

# calculate the upstream to downstream gradient ratios for u>0 (+ ratio)
# x direction
dphiX_p = (phi.value[2:Nr+2, 2:Ntheta+1].-phi.value[1:Nr+1, 2:Ntheta+1])./dr
rX_p = dphiX_p[1:end-1,:]./fsign(dphiX_p[2:end,:])
psiX_p[2:Nr+1,:] .= 0.5.*FL.(rX_p).*(phi.value[3:Nr+2,2:Ntheta+1].-
		    phi.value[2:Nr+1, 2:Ntheta+1])
psiX_p[1, :] .= 0.0 # left boundary will be handled in the main matrix
# y direction
dphiY_p = (phi.value[2:Nr+1, 2:Ntheta+2].-phi.value[2:Nr+1, 1:Ntheta+1])./dtheta
rY_p = dphiY_p[:,1:end-1]./fsign(dphiY_p[:,2:end])
psiY_p[:,2:Ntheta+1] .= 0.5.*FL.(rY_p).*(phi.value[2:Nr+1,3:Ntheta+2].-
		  phi.value[2:Nr+1, 2:Ntheta+1])
psiY_p[:,1] .= 0.0 # Bottom boundary will be handled in the main matrix

# calculate the upstream to downstream gradient ratios for u<0 (- ratio)
# x direction
rX_m = dphiX_p[2:end,:]./fsign(dphiX_p[1:end-1,:])
psiX_m[1:Nr,:] .= 0.5.*FL.(rX_m).*(phi.value[1:Nr, 2:Ntheta+1].-
		phi.value[2:Nr+1, 2:Ntheta+1])
psiX_m[Nr+1,:] .= 0.0 # right boundary
# y direction
rY_m = dphiY_p[:,2:end]./fsign(dphiY_p[:,1:end-1])
psiY_m[:,1:Ntheta] .= 0.5.*FL.(rY_m).*(phi.value[2:Nr+1, 1:Ntheta].-
	      phi.value[2:Nr+1, 2:Ntheta+1])
psiY_m[:, Ntheta+1] .= 0.0 # top boundary will be handled in the main matrix

re = rf[2:Nr+1,:]
rw = rf[1:Nr,:]

# find the velocity direction for the upwind scheme
ue_min = min.(u.xvalue[2:Nr+1,:],0.0)
ue_max = max.(u.xvalue[2:Nr+1,:],0.0)
uw_min = min.(u.xvalue[1:Nr,:],0.0)
uw_max = max.(u.xvalue[1:Nr,:],0.0)
vn_min = min.(u.yvalue[:,2:Ntheta+1],0.0)
vn_max = max.(u.yvalue[:,2:Ntheta+1],0.0)
vs_min = min.(u.yvalue[:,1:Ntheta],0.0)
vs_max = max.(u.yvalue[:,1:Ntheta],0.0)

# calculate the TVD correction term
div_x = -(1.0./(DRp.*rp)).*(re.*(ue_max.*psiX_p[2:Nr+1,:].+ue_min.*psiX_m[2:Nr+1,:])-
              rw.*(uw_max.*psiX_p[1:Nr,:].+uw_min.*psiX_m[1:Nr,:]))
div_y = -(1.0./(DTHETAp.*rp)).*((vn_max.*psiY_p[:,2:Ntheta+1].+vn_min.*psiY_m[:,2:Ntheta+1])-
              (vs_max.*psiY_p[:,1:Ntheta].+vs_min.*psiY_m[:,1:Ntheta]))

# define the RHS Vector
RHS = zeros((Nr+2)*(Ntheta+2))
RHSx = zeros((Nr+2)*(Ntheta+2))
RHSy = zeros((Nr+2)*(Ntheta+2))

# assign the values of the RHS vector
mnx = Nr*Ntheta
mny = Nr*Ntheta
rowx_index = reshape(G[2:Nr+1,2:Ntheta+1],mnx) # main diagonal x
rowy_index = reshape(G[2:Nr+1,2:Ntheta+1],mny) # main diagonal y
RHS[rowx_index] .= reshape(div_x+div_y,Nr*Ntheta)
RHSx[rowx_index] .= reshape(div_x,Nr*Ntheta)
RHSy[rowy_index] .= reshape(div_y,Nr*Ntheta)

(RHS, RHSx, RHSy)
end

function convectionTvdRHSRadial2D(u::FaceValue, phi::CellValue, FL::Function, u_upwind::FaceValue)
# u is a face variable
# phi is a cell variable

# a function to avoid division by zero
eps1 = 1.0e-20
fsign(phi_in) = (abs.(phi_in).>=eps1).*phi_in+eps1*(phi_in.==0.0)+eps1*(abs.(phi_in).<eps1).*sign.(phi_in)

# extract data from the mesh structure
Nr = u.domain.dims[1]
Ntheta = u.domain.dims[2]
G=reshape([1:(Nr+2)*(Ntheta+2);], Nr+2, Ntheta+2)
DRp = u.domain.cellsize.x[2:end-1]
DTHETAp = zeros( 1, Ntheta)
DTHETAp[:] .= u.domain.cellsize.y[2:end-1]
dr = 0.5.*(u.domain.cellsize.x[1:end-1].+u.domain.cellsize.x[2:end])
dtheta = zeros( 1, 1+Ntheta)
dtheta[:] .= 0.5.*(u.domain.cellsize.y[1:end-1].+u.domain.cellsize.y[2:end])
rp = u.domain.cellcenters.x
rf = u.domain.facecenters.x
psiX_p = zeros(Nr+1,Ntheta)
psiX_m = zeros(Nr+1,Ntheta)
psiY_p = zeros(Nr,Ntheta+1)
psiY_m = zeros(Nr,Ntheta+1)

# calculate the upstream to downstream gradient ratios for u>0 (+ ratio)
# x direction
dphiX_p = (phi.value[2:Nr+2, 2:Ntheta+1].-phi.value[1:Nr+1, 2:Ntheta+1])./dr
rX_p = dphiX_p[1:end-1,:]./fsign(dphiX_p[2:end,:])
psiX_p[2:Nr+1,:] .= 0.5.*FL.(rX_p).*(phi.value[3:Nr+2,2:Ntheta+1].-
		    phi.value[2:Nr+1, 2:Ntheta+1])
psiX_p[1, :] .= 0.0 # left boundary will be handled in the main matrix
# y direction
dphiY_p = (phi.value[2:Nr+1, 2:Ntheta+2].-phi.value[2:Nr+1, 1:Ntheta+1])./dtheta
rY_p = dphiY_p[:,1:end-1]./fsign(dphiY_p[:,2:end])
psiY_p[:,2:Ntheta+1] .= 0.5.*FL.(rY_p).*(phi.value[2:Nr+1,3:Ntheta+2].-
		  phi.value[2:Nr+1, 2:Ntheta+1])
psiY_p[:,1] .= 0.0 # Bottom boundary will be handled in the main matrix

# calculate the upstream to downstream gradient ratios for u<0 (- ratio)
# x direction
rX_m = dphiX_p[2:end,:]./fsign(dphiX_p[1:end-1,:])
psiX_m[1:Nr,:] .= 0.5.*FL.(rX_m).*(phi.value[1:Nr, 2:Ntheta+1].-
		phi.value[2:Nr+1, 2:Ntheta+1])
psiX_m[Nr+1,:] .= 0.0 # right boundary
# y direction
rY_m = dphiY_p[:,2:end]./fsign(dphiY_p[:,1:end-1])
psiY_m[:,1:Ntheta] .= 0.5.*FL.(rY_m).*(phi.value[2:Nr+1, 1:Ntheta].-
	      phi.value[2:Nr+1, 2:Ntheta+1])
psiY_m[:, Ntheta+1] .= 0.0 # top boundary will be handled in the main matrix

re = rf[2:Nr+1,:]
rw = rf[1:Nr,:]

# find the velocity direction for the upwind scheme
ue_min = u.xvalue[2:Nr+1,:]
ue_max = u.xvalue[2:Nr+1,:]
uw_min = u.xvalue[1:Nr,:]
uw_max = u.xvalue[1:Nr,:]
vn_min = u.yvalue[:,2:Ntheta+1]
vn_max = u.yvalue[:,2:Ntheta+1]
vs_min = u.yvalue[:,1:Ntheta]
vs_max = u.yvalue[:,1:Ntheta]

ue_min[u_upwind.xvalue[2:Nr+1,:].>0.0] .= 0.0
ue_max[u_upwind.xvalue[2:Nr+1,:].<0.0] .= 0.0
uw_min[u_upwind.xvalue[1:Nr,:].>0.0] .= 0.0
uw_max[u_upwind.xvalue[1:Nr,:].<0.0] .= 0.0
vn_min[u_upwind.yvalue[:,2:Ntheta+1].>0.0] .= 0.0
vn_max[u_upwind.yvalue[:,2:Ntheta+1].<0.0] .= 0.0
vs_min[u_upwind.yvalue[:,1:Ntheta].>0.0] .= 0.0
vs_max[u_upwind.yvalue[:,1:Ntheta].<0.0] .= 0.0

# calculate the TVD correction term
div_x = -(1.0./(DRp.*rp)).*(re.*(ue_max.*psiX_p[2:Nr+1,:].+ue_min.*psiX_m[2:Nr+1,:])-
              rw.*(uw_max.*psiX_p[1:Nr,:].+uw_min.*psiX_m[1:Nr,:]))
div_y = -(1.0./(DTHETAp.*rp)).*((vn_max.*psiY_p[:,2:Ntheta+1].+vn_min.*psiY_m[:,2:Ntheta+1])-
              (vs_max.*psiY_p[:,1:Ntheta].+vs_min.*psiY_m[:,1:Ntheta]))

# define the RHS Vector
RHS = zeros((Nr+2)*(Ntheta+2))
RHSx = zeros((Nr+2)*(Ntheta+2))
RHSy = zeros((Nr+2)*(Ntheta+2))

# assign the values of the RHS vector
mnx = Nr*Ntheta
mny = Nr*Ntheta
rowx_index = reshape(G[2:Nr+1,2:Ntheta+1],mnx) # main diagonal x
rowy_index = reshape(G[2:Nr+1,2:Ntheta+1],mny) # main diagonal y
RHS[rowx_index] .= reshape(div_x+div_y,Nr*Ntheta)
RHSx[rowx_index] .= reshape(div_x,Nr*Ntheta)
RHSy[rowy_index] .= reshape(div_y,Nr*Ntheta)

(RHS, RHSx, RHSy)
end

# ============================= 3D Convection Term ===============================
function convectionTerm3D(u::FaceValue)
# u is a face variable
# extract data from the mesh structure
Nx = u.domain.dims[1]
Ny = u.domain.dims[2]
Nz = u.domain.dims[3]
G=reshape([1:(Nx+2)*(Ny+2)*(Nz+2);], Nx+2, Ny+2, Nz+2)
DXe = u.domain.cellsize.x[3:end]
DXw = u.domain.cellsize.x[1:end-2]
DXp = u.domain.cellsize.x[2:end-1]
DY = zeros( 1, Ny+2)
DY[:] .= u.domain.cellsize.y
DYn = DY[:,3:end]
DYs = DY[:,1:end-2]
DYp = DY[:,2:end-1]
DZ = zeros( 1, 1, Nz+2)
DZ[:] .= u.domain.cellsize.z
DZf = DZ[:,:,3:end]
DZb = DZ[:,:,1:end-2]
DZp = DZ[:,:,2:end-1]

# define the vectors to stores the sparse matrix data
iix = zeros(Int64, 3*(Nx+2)*(Ny+2)*(Nz+2))
jjx = zeros(Int64, 3*(Nx+2)*(Ny+2)*(Nz+2))
sx = zeros(Float64, 3*(Nx+2)*(Ny+2)*(Nz+2))
iiy = zeros(Int64, 3*(Nx+2)*(Ny+2)*(Nz+2))
jjy = zeros(Int64, 3*(Nx+2)*(Ny+2)*(Nz+2))
sy = zeros(Float64, 3*(Nx+2)*(Ny+2)*(Nz+2))
iiz = zeros(Int64, 3*(Nx+2)*(Ny+2)*(Nz+2))
jjz = zeros(Int64, 3*(Nx+2)*(Ny+2)*(Nz+2))
sz = zeros(Float64, 3*(Nx+2)*(Ny+2)*(Nz+2))
mnx = Nx*Ny*Nz
mny = Nx*Ny*Nz
mnz = Nx*Ny*Nz

# reassign the east, west, north, and south velocity vectors for the
# code readability
ue = u.xvalue[2:Nx+1,:,:]./(DXp+DXe)
uw = u.xvalue[1:Nx,:,:]./(DXp+DXw)
vn = u.yvalue[:,2:Ny+1,:]./(DYp+DYn)
vs = u.yvalue[:,1:Ny,:]./(DYp+DYs)
wf = u.zvalue[:,:,2:Nz+1]./(DZp+DZf)
wb = u.zvalue[:,:,1:Nz]./(DZp+DZb)

# calculate the coefficients for the internal cells
AE = reshape(ue,mnx)
AW = reshape(-uw,mnx)
AN = reshape(vn,mny)
AS = reshape(-vs,mny)
AF = reshape(wf,mnz)
AB = reshape(-wb,mnz)
APx = reshape((ue.*DXe-uw.*DXw)./DXp,mnx)
APy = reshape((vn.*DYn-vs.*DYs)./DYp,mny)
APz = reshape((wf.*DZf-wb.*DZb)./DZp,mnz)

# build the sparse matrix based on the numbering system
rowx_index = reshape(G[2:Nx+1,2:Ny+1,2:Nz+1],mnx)  # main diagonal x
iix[1:3*mnx] .= repeat(rowx_index,3)
rowy_index = reshape(G[2:Nx+1,2:Ny+1,2:Nz+1],mny)  # main diagonal y
iiy[1:3*mny] .= repeat(rowy_index,3)
rowz_index = reshape(G[2:Nx+1,2:Ny+1,2:Nz+1],mnz)  # main diagonal z
iiz[1:3*mnz] .= repeat(rowz_index,3)
jjx[1:3*mnx] .= [reshape(G[1:Nx,2:Ny+1,2:Nz+1],mnx); reshape(G[2:Nx+1,2:Ny+1,2:Nz+1],mnx); reshape(G[3:Nx+2,2:Ny+1,2:Nz+1],mnx)]
jjy[1:3*mny] .= [reshape(G[2:Nx+1,1:Ny,2:Nz+1],mny); reshape(G[2:Nx+1,2:Ny+1,2:Nz+1],mny); reshape(G[2:Nx+1,3:Ny+2,2:Nz+1],mny)]
jjz[1:3*mnz] .= [reshape(G[2:Nx+1,2:Ny+1,1:Nz],mnz); reshape(G[2:Nx+1,2:Ny+1,2:Nz+1],mnz); reshape(G[2:Nx+1,2:Ny+1,3:Nz+2],mnz)]
sx[1:3*mnx] .= [AW; APx; AE]
sy[1:3*mny] .= [AS; APy; AN]
sz[1:3*mnz] .= [AB; APz; AF]

# build the sparse matrix
kx = 3*mnx
ky = 3*mny
kz = 3*mnz
Mx = sparse(iix[1:kx], jjx[1:kx], sx[1:kx], (Nx+2)*(Ny+2)*(Nz+2), (Nx+2)*(Ny+2)*(Nz+2))
My = sparse(iiy[1:ky], jjy[1:ky], sy[1:ky], (Nx+2)*(Ny+2)*(Nz+2), (Nx+2)*(Ny+2)*(Nz+2))
Mz = sparse(iiz[1:kz], jjz[1:kz], sz[1:kz], (Nx+2)*(Ny+2)*(Nz+2), (Nx+2)*(Ny+2)*(Nz+2))
M = Mx + My + Mz

(M, Mx, My, Mz)
end


# ============================= 3D Convection Upwind Term ===============================
function convectionUpwindTerm3D(u::FaceValue)
# u is a face variable
# extract data from the mesh structure
Nx = u.domain.dims[1]
Ny = u.domain.dims[2]
Nz = u.domain.dims[3]
G=reshape([1:(Nx+2)*(Ny+2)*(Nz+2);], Nx+2, Ny+2, Nz+2)
DXp = u.domain.cellsize.x[2:end-1]
DY = zeros( 1, Ny+2)
DY[:] .= u.domain.cellsize.y
DYp = DY[:,2:end-1]
DZ = zeros( 1, 1, Nz+2)
DZ[:] .= u.domain.cellsize.z
DZp = DZ[:,:,2:end-1]

# define the vectors to stores the sparse matrix data
iix = zeros(Int64, 3*(Nx+2)*(Ny+2)*(Nz+2))
jjx = zeros(Int64, 3*(Nx+2)*(Ny+2)*(Nz+2))
sx = zeros(Float64, 3*(Nx+2)*(Ny+2)*(Nz+2))
iiy = zeros(Int64, 3*(Nx+2)*(Ny+2)*(Nz+2))
jjy = zeros(Int64, 3*(Nx+2)*(Ny+2)*(Nz+2))
sy = zeros(Float64, 3*(Nx+2)*(Ny+2)*(Nz+2))
iiz = zeros(Int64, 3*(Nx+2)*(Ny+2)*(Nz+2))
jjz = zeros(Int64, 3*(Nx+2)*(Ny+2)*(Nz+2))
sz = zeros(Float64, 3*(Nx+2)*(Ny+2)*(Nz+2))
mnx = Nx*Ny*Nz
mny = Nx*Ny*Nz
mnz = Nx*Ny*Nz

# find the velocity direction for the upwind scheme
ue_min = min.(u.xvalue[2:Nx+1,:,:],0)
ue_max = max.(u.xvalue[2:Nx+1,:,:],0)
uw_min = min.(u.xvalue[1:Nx,:,:],0)
uw_max = max.(u.xvalue[1:Nx,:,:],0)
vn_min = min.(u.yvalue[:,2:Ny+1,:],0)
vn_max = max.(u.yvalue[:,2:Ny+1,:],0)
vs_min = min.(u.yvalue[:,1:Ny,:],0)
vs_max = max.(u.yvalue[:,1:Ny,:],0)
wf_min = min.(u.zvalue[:,:,2:Nz+1],0)
wf_max = max.(u.zvalue[:,:,2:Nz+1],0)
wb_min = min.(u.zvalue[:,:,1:Nz],0)
wb_max = max.(u.zvalue[:,:,1:Nz],0)

# calculate the coefficients for the internal cells
AE = ue_min./DXp
AW = -uw_max./DXp
AN = vn_min./DYp
AS = -vs_max./DYp
AF = wf_min./DZp
AB = -wb_max./DZp
APx = (ue_max.-uw_min)./DXp
APy = (vn_max.-vs_min)./DYp
APz = (wf_max.-wb_min)./DZp

# Also correct for the boundary cells (not the ghost cells)
# Left boundary:
APx[1,:,:] .= APx[1,:,:].-uw_max[1,:,:]/(2.0*DXp[1])
AW[1,:,:] .= AW[1,:,:]/2.0
# Right boundary:
AE[end,:,:] .= AE[end,:,:]/2.0
APx[end,:,:] .= APx[end,:,:].+ue_min[end,:,:]/(2.0*DXp[end])
# Bottom boundary:
APy[:,1,:] .= APy[:,1,:].-vs_max[:,1,:]/(2.0*DYp[1])
AS[:,1,:] .= AS[:,1,:]/2.0
# Top boundary:
AN[:,end,:] .= AN[:,end,:]/2.0
APy[:,end,:] .= APy[:,end,:].+vn_min[:,end,:]/(2.0*DYp[end])
# Back boundary:
APz[:,:,1] .= APz[:,:,1].-wb_max[:,:,1]/(2.0*DZp[1])
AB[:,:,1] .= AB[:,:,1]/2.0
# Front boundary:
AF[:,:,end] .= AF[:,:,end]/2.0
APz[:,:,end] .= APz[:,:,end].+wf_min[:,:,end]/(2.0*DZp[end])

AE = reshape(AE,mnx)
AW = reshape(AW,mnx)
AN = reshape(AN,mny)
AS = reshape(AS,mny)
AF = reshape(AF,mnz)
AB = reshape(AB,mnz)
APx = reshape(APx,mnx)
APy = reshape(APy,mny)
APz = reshape(APz,mnz)

# build the sparse matrix based on the numbering system
rowx_index = reshape(G[2:Nx+1,2:Ny+1,2:Nz+1],mnx)  # main diagonal x
iix[1:3*mnx] .= repeat(rowx_index,3)
rowy_index = reshape(G[2:Nx+1,2:Ny+1,2:Nz+1],mny)  # main diagonal y
iiy[1:3*mny] .= repeat(rowy_index,3)
rowz_index = reshape(G[2:Nx+1,2:Ny+1,2:Nz+1],mnz)  # main diagonal z
iiz[1:3*mnz] .= repeat(rowz_index,3)
jjx[1:3*mnx] .= [reshape(G[1:Nx,2:Ny+1,2:Nz+1],mnx); reshape(G[2:Nx+1,2:Ny+1,2:Nz+1],mnx); reshape(G[3:Nx+2,2:Ny+1,2:Nz+1],mnx)]
jjy[1:3*mny] .= [reshape(G[2:Nx+1,1:Ny,2:Nz+1],mny); reshape(G[2:Nx+1,2:Ny+1,2:Nz+1],mny); reshape(G[2:Nx+1,3:Ny+2,2:Nz+1],mny)]
jjz[1:3*mnz] .= [reshape(G[2:Nx+1,2:Ny+1,1:Nz],mnz); reshape(G[2:Nx+1,2:Ny+1,2:Nz+1],mnz); reshape(G[2:Nx+1,2:Ny+1,3:Nz+2],mnz)]
sx[1:3*mnx] .= [AW; APx; AE]
sy[1:3*mny] .= [AS; APy; AN]
sz[1:3*mnz] .= [AB; APz; AF]

# build the sparse matrix
kx = 3*mnx
ky = 3*mny
kz = 3*mnz
Mx = sparse(iix[1:kx], jjx[1:kx], sx[1:kx], (Nx+2)*(Ny+2)*(Nz+2), (Nx+2)*(Ny+2)*(Nz+2))
My = sparse(iiy[1:ky], jjy[1:ky], sy[1:ky], (Nx+2)*(Ny+2)*(Nz+2), (Nx+2)*(Ny+2)*(Nz+2))
Mz = sparse(iiz[1:kz], jjz[1:kz], sz[1:kz], (Nx+2)*(Ny+2)*(Nz+2), (Nx+2)*(Ny+2)*(Nz+2))
M = Mx + My + Mz;

(M, Mx, My, Mz)
end

function convectionUpwindTerm3D(u::FaceValue, u_upwind::FaceValue)
# u is a face variable
# extract data from the mesh structure
Nx = u.domain.dims[1]
Ny = u.domain.dims[2]
Nz = u.domain.dims[3]
G=reshape([1:(Nx+2)*(Ny+2)*(Nz+2);], Nx+2, Ny+2, Nz+2)
DXp = u.domain.cellsize.x[2:end-1]
DY = zeros( 1, Ny+2)
DY[:] .= u.domain.cellsize.y
DYp = DY[:,2:end-1]
DZ = zeros( 1, 1, Nz+2)
DZ[:] .= u.domain.cellsize.z
DZp = DZ[:,:,2:end-1]

# define the vectors to stores the sparse matrix data
iix = zeros(Int64, 3*(Nx+2)*(Ny+2)*(Nz+2))
jjx = zeros(Int64, 3*(Nx+2)*(Ny+2)*(Nz+2))
sx = zeros(Float64, 3*(Nx+2)*(Ny+2)*(Nz+2))
iiy = zeros(Int64, 3*(Nx+2)*(Ny+2)*(Nz+2))
jjy = zeros(Int64, 3*(Nx+2)*(Ny+2)*(Nz+2))
sy = zeros(Float64, 3*(Nx+2)*(Ny+2)*(Nz+2))
iiz = zeros(Int64, 3*(Nx+2)*(Ny+2)*(Nz+2))
jjz = zeros(Int64, 3*(Nx+2)*(Ny+2)*(Nz+2))
sz = zeros(Float64, 3*(Nx+2)*(Ny+2)*(Nz+2))
mnx = Nx*Ny*Nz
mny = Nx*Ny*Nz
mnz = Nx*Ny*Nz

# find the velocity direction for the upwind scheme
ue_min = u.xvalue[2:Nx+1,:,:]
ue_max = u.xvalue[2:Nx+1,:,:]
uw_min = u.xvalue[1:Nx,:,:]
uw_max = u.xvalue[1:Nx,:,:]
vn_min = u.yvalue[:,2:Ny+1,:]
vn_max = u.yvalue[:,2:Ny+1,:]
vs_min = u.yvalue[:,1:Ny,:]
vs_max = u.yvalue[:,1:Ny,:]
wf_min = u.zvalue[:,:,2:Nz+1]
wf_max = u.zvalue[:,:,2:Nz+1]
wb_min = u.zvalue[:,:,1:Nz]
wb_max = u.zvalue[:,:,1:Nz]

ue_min[u_upwind.xvalue[2:Nx+1,:,:].>0.0] .= 0.0
ue_max[u_upwind.xvalue[2:Nx+1,:,:].<0.0] .= 0.0
uw_min[u_upwind.xvalue[1:Nx,:,:].>0.0] .= 0.0
uw_max[u_upwind.xvalue[1:Nx,:,:].<0.0] .= 0.0
vn_min[u_upwind.yvalue[:,2:Ny+1,:].>0.0] .= 0.0
vn_max[u_upwind.yvalue[:,2:Ny+1,:].<0.0] .= 0.0
vs_min[u_upwind.yvalue[:,1:Ny,:].>0.0] .= 0.0
vs_max[u_upwind.yvalue[:,1:Ny,:].<0.0] .= 0.0
wf_min[u_upwind.zvalue[:,:,2:Nz+1].>0.0] .= 0.0
wf_max[u_upwind.zvalue[:,:,2:Nz+1].<0.0] .= 0.0
wb_min[u_upwind.zvalue[:,:,1:Nz].>0.0] .= 0.0
wb_max[u_upwind.zvalue[:,:,1:Nz].<0.0] .= 0.0

# calculate the coefficients for the internal cells
AE = ue_min./DXp
AW = -uw_max./DXp
AN = vn_min./DYp
AS = -vs_max./DYp
AF = wf_min./DZp
AB = -wb_max./DZp
APx = (ue_max.-uw_min)./DXp
APy = (vn_max.-vs_min)./DYp
APz = (wf_max.-wb_min)./DZp

# Also correct for the boundary cells (not the ghost cells)
# Left boundary:
APx[1,:,:] .= APx[1,:,:].-uw_max[1,:,:]/(2.0*DXp[1])
AW[1,:,:] .= AW[1,:,:]/2.0
# Right boundary:
AE[end,:,:] .= AE[end,:,:]/2.0
APx[end,:,:] .= APx[end,:,:].+ue_min[end,:,:]/(2.0*DXp[end])
# Bottom boundary:
APy[:,1,:] .= APy[:,1,:].-vs_max[:,1,:]/(2.0*DYp[1])
AS[:,1,:] .= AS[:,1,:]/2.0
# Top boundary:
AN[:,end,:] .= AN[:,end,:]/2.0
APy[:,end,:] .= APy[:,end,:].+vn_min[:,end,:]/(2.0*DYp[end])
# Back boundary:
APz[:,:,1] .= APz[:,:,1].-wb_max[:,:,1]/(2.0*DZp[1])
AB[:,:,1] .= AB[:,:,1]/2.0
# Front boundary:
AF[:,:,end] .= AF[:,:,end]/2.0
APz[:,:,end] .= APz[:,:,end].+wf_min[:,:,end]/(2.0*DZp[end])

AE = reshape(AE,mnx)
AW = reshape(AW,mnx)
AN = reshape(AN,mny)
AS = reshape(AS,mny)
AF = reshape(AF,mnz)
AB = reshape(AB,mnz)
APx = reshape(APx,mnx)
APy = reshape(APy,mny)
APz = reshape(APz,mnz)

# build the sparse matrix based on the numbering system
rowx_index = reshape(G[2:Nx+1,2:Ny+1,2:Nz+1],mnx)  # main diagonal x
iix[1:3*mnx] .= repeat(rowx_index,3)
rowy_index = reshape(G[2:Nx+1,2:Ny+1,2:Nz+1],mny)  # main diagonal y
iiy[1:3*mny] .= repeat(rowy_index,3)
rowz_index = reshape(G[2:Nx+1,2:Ny+1,2:Nz+1],mnz)  # main diagonal z
iiz[1:3*mnz] .= repeat(rowz_index,3)
jjx[1:3*mnx] .= [reshape(G[1:Nx,2:Ny+1,2:Nz+1],mnx); reshape(G[2:Nx+1,2:Ny+1,2:Nz+1],mnx); reshape(G[3:Nx+2,2:Ny+1,2:Nz+1],mnx)]
jjy[1:3*mny] .= [reshape(G[2:Nx+1,1:Ny,2:Nz+1],mny); reshape(G[2:Nx+1,2:Ny+1,2:Nz+1],mny); reshape(G[2:Nx+1,3:Ny+2,2:Nz+1],mny)]
jjz[1:3*mnz] .= [reshape(G[2:Nx+1,2:Ny+1,1:Nz],mnz); reshape(G[2:Nx+1,2:Ny+1,2:Nz+1],mnz); reshape(G[2:Nx+1,2:Ny+1,3:Nz+2],mnz)]
sx[1:3*mnx] .= [AW; APx; AE]
sy[1:3*mny] .= [AS; APy; AN]
sz[1:3*mnz] .= [AB; APz; AF]

# build the sparse matrix
kx = 3*mnx
ky = 3*mny
kz = 3*mnz
Mx = sparse(iix[1:kx], jjx[1:kx], sx[1:kx], (Nx+2)*(Ny+2)*(Nz+2), (Nx+2)*(Ny+2)*(Nz+2))
My = sparse(iiy[1:ky], jjy[1:ky], sy[1:ky], (Nx+2)*(Ny+2)*(Nz+2), (Nx+2)*(Ny+2)*(Nz+2))
Mz = sparse(iiz[1:kz], jjz[1:kz], sz[1:kz], (Nx+2)*(Ny+2)*(Nz+2), (Nx+2)*(Ny+2)*(Nz+2))
M = Mx + My + Mz;

(M, Mx, My, Mz)
end


# ============================= 3D Convection TVD Term ===============================
function convectionTvdTerm3D(u::FaceValue, phi::CellValue, FL::Function)
# u is a face variable
# a function to avoid division by zero
eps1 = 1.0e-20
fsign(phi_in) = (abs.(phi_in).>=eps1).*phi_in+eps1*(phi_in.==0.0)+eps1*(abs.(phi_in).<eps1).*sign.(phi_in)
# extract data from the mesh structure
Nx = u.domain.dims[1]
Ny = u.domain.dims[2]
Nz = u.domain.dims[3]
G=reshape([1:(Nx+2)*(Ny+2)*(Nz+2);], Nx+2, Ny+2, Nz+2)
DXp = u.domain.cellsize.x[2:end-1]
DY = zeros( 1, Ny+2)
DY[:] .= u.domain.cellsize.y
DYp = DY[:,2:end-1]
DZ = zeros( 1, 1, Nz+2)
DZ[:] .= u.domain.cellsize.z
DZp = DZ[:,:,2:end-1]
dx=0.5.*(u.domain.cellsize.x[1:end-1].+u.domain.cellsize.x[2:end])
dy=zeros( 1, Ny+1)
dy[:] .=0.5.*(u.domain.cellsize.y[1:end-1].+u.domain.cellsize.y[2:end])
dz=zeros( 1, 1, Nz+1)
dz[:] .=0.5.*(u.domain.cellsize.z[1:end-1].+u.domain.cellsize.z[2:end])
psiX_p = zeros(Nx+1,Ny,Nz)
psiX_m = zeros(Nx+1,Ny,Nz)
psiY_p = zeros(Nx,Ny+1,Nz)
psiY_m = zeros(Nx,Ny+1,Nz)
psiZ_p = zeros(Nx,Ny,Nz+1)
psiZ_m = zeros(Nx,Ny,Nz+1)

# define the vectors to stores the sparse matrix data
iix = zeros(Int64, 3*(Nx+2)*(Ny+2)*(Nz+2))
jjx = zeros(Int64, 3*(Nx+2)*(Ny+2)*(Nz+2))
sx = zeros(Float64, 3*(Nx+2)*(Ny+2)*(Nz+2))
iiy = zeros(Int64, 3*(Nx+2)*(Ny+2)*(Nz+2))
jjy = zeros(Int64, 3*(Nx+2)*(Ny+2)*(Nz+2))
sy = zeros(Float64, 3*(Nx+2)*(Ny+2)*(Nz+2))
iiz = zeros(Int64, 3*(Nx+2)*(Ny+2)*(Nz+2))
jjz = zeros(Int64, 3*(Nx+2)*(Ny+2)*(Nz+2))
sz = zeros(Float64, 3*(Nx+2)*(Ny+2)*(Nz+2))
mnx = Nx*Ny*Nz
mny = Nx*Ny*Nz
mnz = Nx*Ny*Nz

# calculate the upstream to downstream gradient ratios for u>0 (+ ratio)
# x direction
dphiX_p = (phi.value[2:Nx+2, 2:Ny+1, 2:Nz+1].-phi.value[1:Nx+1, 2:Ny+1, 2:Nz+1])./dx
rX_p = dphiX_p[1:end-1,:,:]./fsign(dphiX_p[2:end,:,:])
psiX_p[2:Nx+1,:,:] .= 0.5.*FL.(rX_p).*(phi.value[3:Nx+2,2:Ny+1,2:Nz+1].-phi.value[2:Nx+1,2:Ny+1,2:Nz+1])
psiX_p[1,:,:] .= 0.0  # left boundary
# y direction
dphiY_p = (phi.value[2:Nx+1, 2:Ny+2, 2:Nz+1].-phi.value[2:Nx+1, 1:Ny+1, 2:Nz+1])./dy
rY_p = dphiY_p[:,1:end-1,:]./fsign(dphiY_p[:,2:end,:])
psiY_p[:,2:Ny+1,:] .= 0.5.*FL.(rY_p).*(phi.value[2:Nx+1,3:Ny+2,2:Nz+1].-phi.value[2:Nx+1, 2:Ny+1,2:Nz+1])
psiY_p[:,1,:] .= 0.0  # Bottom boundary
# z direction
dphiZ_p = (phi.value[2:Nx+1, 2:Ny+1, 2:Nz+2].-phi.value[2:Nx+1, 2:Ny+1, 1:Nz+1])./dz
rZ_p = dphiZ_p[:,:,1:end-1]./fsign(dphiZ_p[:,:,2:end])
psiZ_p[:,:,2:Nz+1] .= 0.5.*FL.(rZ_p).*(phi.value[2:Nx+1,2:Ny+1,3:Nz+2].-phi.value[2:Nx+1,2:Ny+1,2:Nz+1])
psiZ_p[:,:,1] .= 0.0  # Back boundary

# calculate the upstream to downstream gradient ratios for u<0 (- ratio)
# x direction
rX_m = dphiX_p[2:end,:,:]./fsign(dphiX_p[1:end-1,:,:])
psiX_m[1:Nx,:,:] .= 0.5.*FL.(rX_m).*(phi.value[1:Nx, 2:Ny+1, 2:Nz+1].-phi.value[2:Nx+1, 2:Ny+1, 2:Nz+1])
psiX_m[Nx+1,:,:] .= 0.0  # right boundary
# y direction
rY_m = dphiY_p[:,2:end,:]./fsign(dphiY_p[:,1:end-1,:])
psiY_m[:,1:Ny,:] .= 0.5.*FL.(rY_m).*(phi.value[2:Nx+1,1:Ny,2:Nz+1].-phi.value[2:Nx+1,2:Ny+1,2:Nz+1])
psiY_m[:,Ny+1,:] .= 0.0  # top boundary
# z direction
rZ_m = dphiZ_p[:,:,2:end]./fsign(dphiZ_p[:,:,1:end-1])
psiZ_m[:,:,1:Nz] .= 0.5.*FL.(rZ_m).*(phi.value[2:Nx+1,2:Ny+1,1:Nz].-phi.value[2:Nx+1,2:Ny+1,2:Nz+1])
psiZ_m[:,:,Nz+1] .= 0.0  # front boundary

# find the velocity direction for the upwind scheme
ue_min = min.(u.xvalue[2:Nx+1,:,:],0)
ue_max = max.(u.xvalue[2:Nx+1,:,:],0)
uw_min = min.(u.xvalue[1:Nx,:,:],0)
uw_max = max.(u.xvalue[1:Nx,:,:],0)
vn_min = min.(u.yvalue[:,2:Ny+1,:],0)
vn_max = max.(u.yvalue[:,2:Ny+1,:],0)
vs_min = min.(u.yvalue[:,1:Ny,:],0)
vs_max = max.(u.yvalue[:,1:Ny,:],0)
wf_min = min.(u.zvalue[:,:,2:Nz+1],0)
wf_max = max.(u.zvalue[:,:,2:Nz+1],0)
wb_min = min.(u.zvalue[:,:,1:Nz],0)
wb_max = max.(u.zvalue[:,:,1:Nz],0)

# calculate the coefficients for the internal cells
AE = ue_min./DXp
AW = -uw_max./DXp
AN = vn_min./DYp
AS = -vs_max./DYp
AF = wf_min./DZp
AB = -wb_max./DZp
APx = (ue_max.-uw_min)./DXp
APy = (vn_max.-vs_min)./DYp
APz = (wf_max.-wb_min)./DZp

# Also correct for the boundary cells (not the ghost cells)
# Left boundary:
APx[1,:,:] .= APx[1,:,:].-uw_max[1,:,:]/(2.0*DXp[1])
AW[1,:,:] .= AW[1,:,:]/2.0
# Right boundary:
AE[end,:,:] .= AE[end,:,:]/2.0
APx[end,:,:] .= APx[end,:,:].+ue_min[end,:,:]/(2.0*DXp[end])
# Bottom boundary:
APy[:,1,:] .= APy[:,1,:].-vs_max[:,1,:]/(2.0*DYp[1])
AS[:,1,:] .= AS[:,1,:]/2.0
# Top boundary:
AN[:,end,:] .= AN[:,end,:]/2.0
APy[:,end,:] .= APy[:,end,:].+vn_min[:,end,:]/(2.0*DYp[end])
# Back boundary:
APz[:,:,1] .= APz[:,:,1].-wb_max[:,:,1]/(2.0*DZp[1])
AB[:,:,1] .= AB[:,:,1]/2.0
# Front boundary:
AF[:,:,end] .= AF[:,:,end]/2.0
APz[:,:,end] .= APz[:,:,end].+wf_min[:,:,end]/(2.0*DZp[end])

AE = reshape(AE,mnx)
AW = reshape(AW,mnx)
AN = reshape(AN,mny)
AS = reshape(AS,mny)
AF = reshape(AF,mnz)
AB = reshape(AB,mnz)
APx = reshape(APx,mnx)
APy = reshape(APy,mny)
APz = reshape(APz,mnz)

# build the sparse matrix based on the numbering system
rowx_index = reshape(G[2:Nx+1,2:Ny+1,2:Nz+1],mnx)  # main diagonal x
iix[1:3*mnx] .= repeat(rowx_index,3)
rowy_index = reshape(G[2:Nx+1,2:Ny+1,2:Nz+1],mny)  # main diagonal y
iiy[1:3*mny] .= repeat(rowy_index,3)
rowz_index = reshape(G[2:Nx+1,2:Ny+1,2:Nz+1],mnz)  # main diagonal z
iiz[1:3*mnz] .= repeat(rowz_index,3)
jjx[1:3*mnx] .= [reshape(G[1:Nx,2:Ny+1,2:Nz+1],mnx); reshape(G[2:Nx+1,2:Ny+1,2:Nz+1],mnx); reshape(G[3:Nx+2,2:Ny+1,2:Nz+1],mnx)]
jjy[1:3*mny] .= [reshape(G[2:Nx+1,1:Ny,2:Nz+1],mny); reshape(G[2:Nx+1,2:Ny+1,2:Nz+1],mny); reshape(G[2:Nx+1,3:Ny+2,2:Nz+1],mny)]
jjz[1:3*mnz] .= [reshape(G[2:Nx+1,2:Ny+1,1:Nz],mnz); reshape(G[2:Nx+1,2:Ny+1,2:Nz+1],mnz); reshape(G[2:Nx+1,2:Ny+1,3:Nz+2],mnz)]
sx[1:3*mnx] .= [AW; APx; AE]
sy[1:3*mny] .= [AS; APy; AN]
sz[1:3*mnz] .= [AB; APz; AF]

# calculate the TVD correction term
div_x = -(1.0./DXp).*((ue_max.*psiX_p[2:Nx+1,:,:].+ue_min.*psiX_m[2:Nx+1,:,:])-
              (uw_max.*psiX_p[1:Nx,:,:].+uw_min.*psiX_m[1:Nx,:,:]))
div_y = -(1.0./DYp).*((vn_max.*psiY_p[:,2:Ny+1,:].+vn_min.*psiY_m[:,2:Ny+1,:])-
              (vs_max.*psiY_p[:,1:Ny,:].+vs_min.*psiY_m[:,1:Ny,:]))
div_z = -(1.0./DZp).*((wf_max.*psiZ_p[:,:,2:Nz+1].+wf_min.*psiZ_m[:,:,2:Nz+1])-
              (wb_max.*psiZ_p[:,:,1:Nz].+wb_min.*psiZ_m[:,:,1:Nz]))

# define the RHS Vector
RHS = zeros((Nx+2)*(Ny+2)*(Nz+2))
RHSx = zeros((Nx+2)*(Ny+2)*(Nz+2))
RHSy = zeros((Nx+2)*(Ny+2)*(Nz+2))
RHSz = zeros((Nx+2)*(Ny+2)*(Nz+2))

# assign the values of the RHS vector
row_index = rowx_index
RHS[row_index] .= reshape(div_x+div_y+div_z,Nx*Ny*Nz)
RHSx[rowx_index] .= reshape(div_x,Nx*Ny*Nz)
RHSy[rowy_index] .= reshape(div_y,Nx*Ny*Nz)
RHSz[rowz_index] .= reshape(div_z,Nx*Ny*Nz)

# build the sparse matrix
kx = 3*mnx
ky = 3*mny
kz = 3*mnz
Mx = sparse(iix[1:kx], jjx[1:kx], sx[1:kx], (Nx+2)*(Ny+2)*(Nz+2), (Nx+2)*(Ny+2)*(Nz+2))
My = sparse(iiy[1:ky], jjy[1:ky], sy[1:ky], (Nx+2)*(Ny+2)*(Nz+2), (Nx+2)*(Ny+2)*(Nz+2))
Mz = sparse(iiz[1:kz], jjz[1:kz], sz[1:kz], (Nx+2)*(Ny+2)*(Nz+2), (Nx+2)*(Ny+2)*(Nz+2))
M = Mx + My + Mz;

(M, RHS, Mx, My, Mz, RHSx, RHSy, RHSz)
end


# ============================= 3D Convection TVD RHS ===============================
function convectionTvdRHS3D(u::FaceValue, phi::CellValue, FL::Function)
# u is a face variable
# a function to avoid division by zero
eps1 = 1.0e-20
fsign(phi_in) = (abs.(phi_in).>=eps1).*phi_in+eps1*(phi_in.==0.0)+eps1*(abs.(phi_in).<eps1).*sign.(phi_in)
# extract data from the mesh structure
Nx = u.domain.dims[1]
Ny = u.domain.dims[2]
Nz = u.domain.dims[3]
G=reshape([1:(Nx+2)*(Ny+2)*(Nz+2);], Nx+2, Ny+2, Nz+2)
DXp = u.domain.cellsize.x[2:end-1]
DY = zeros( 1, Ny+2)
DY[:] .= u.domain.cellsize.y
DYp = DY[:,2:end-1]
DZ = zeros( 1, 1, Nz+2)
DZ[:] .= u.domain.cellsize.z
DZp = DZ[:,:,2:end-1]
dx=0.5.*(u.domain.cellsize.x[1:end-1].+u.domain.cellsize.x[2:end])
dy=zeros( 1, Ny+1)
dy[:] .=0.5.*(u.domain.cellsize.y[1:end-1].+u.domain.cellsize.y[2:end])
dz=zeros( 1, 1, Nz+1)
dz[:] .=0.5.*(u.domain.cellsize.z[1:end-1].+u.domain.cellsize.z[2:end])
psiX_p = zeros(Nx+1,Ny,Nz)
psiX_m = zeros(Nx+1,Ny,Nz)
psiY_p = zeros(Nx,Ny+1,Nz)
psiY_m = zeros(Nx,Ny+1,Nz)
psiZ_p = zeros(Nx,Ny,Nz+1)
psiZ_m = zeros(Nx,Ny,Nz+1)

# calculate the upstream to downstream gradient ratios for u>0 (+ ratio)
# x direction
dphiX_p = (phi.value[2:Nx+2, 2:Ny+1, 2:Nz+1].-phi.value[1:Nx+1, 2:Ny+1, 2:Nz+1])./dx
rX_p = dphiX_p[1:end-1,:,:]./fsign(dphiX_p[2:end,:,:])
psiX_p[2:Nx+1,:,:] .= 0.5.*FL.(rX_p).*(phi.value[3:Nx+2,2:Ny+1,2:Nz+1].-phi.value[2:Nx+1,2:Ny+1,2:Nz+1])
psiX_p[1,:,:] .= 0.0  # left boundary
# y direction
dphiY_p = (phi.value[2:Nx+1, 2:Ny+2, 2:Nz+1].-phi.value[2:Nx+1, 1:Ny+1, 2:Nz+1])./dy
rY_p = dphiY_p[:,1:end-1,:]./fsign(dphiY_p[:,2:end,:])
psiY_p[:,2:Ny+1,:] .= 0.5.*FL.(rY_p).*(phi.value[2:Nx+1,3:Ny+2,2:Nz+1].-phi.value[2:Nx+1, 2:Ny+1,2:Nz+1])
psiY_p[:,1,:] .= 0.0  # Bottom boundary
# z direction
dphiZ_p = (phi.value[2:Nx+1, 2:Ny+1, 2:Nz+2].-phi.value[2:Nx+1, 2:Ny+1, 1:Nz+1])./dz
rZ_p = dphiZ_p[:,:,1:end-1]./fsign(dphiZ_p[:,:,2:end])
psiZ_p[:,:,2:Nz+1] .= 0.5.*FL.(rZ_p).*(phi.value[2:Nx+1,2:Ny+1,3:Nz+2].-phi.value[2:Nx+1,2:Ny+1,2:Nz+1])
psiZ_p[:,:,1] .= 0.0  # Back boundary

# calculate the upstream to downstream gradient ratios for u<0 (- ratio)
# x direction
rX_m = dphiX_p[2:end,:,:]./fsign(dphiX_p[1:end-1,:,:])
psiX_m[1:Nx,:,:] .= 0.5.*FL.(rX_m).*(phi.value[1:Nx, 2:Ny+1, 2:Nz+1].-phi.value[2:Nx+1, 2:Ny+1, 2:Nz+1])
psiX_m[Nx+1,:,:] .= 0.0  # right boundary
# y direction
rY_m = dphiY_p[:,2:end,:]./fsign(dphiY_p[:,1:end-1,:])
psiY_m[:,1:Ny,:] .= 0.5.*FL.(rY_m).*(phi.value[2:Nx+1,1:Ny,2:Nz+1].-phi.value[2:Nx+1,2:Ny+1,2:Nz+1])
psiY_m[:,Ny+1,:] .= 0.0  # top boundary
# z direction
rZ_m = dphiZ_p[:,:,2:end]./fsign(dphiZ_p[:,:,1:end-1])
psiZ_m[:,:,1:Nz] .= 0.5.*FL.(rZ_m).*(phi.value[2:Nx+1,2:Ny+1,1:Nz].-phi.value[2:Nx+1,2:Ny+1,2:Nz+1])
psiZ_m[:,:,Nz+1] .= 0.0  # front boundary

# find the velocity direction for the upwind scheme
ue_min = min.(u.xvalue[2:Nx+1,:,:],0)
ue_max = max.(u.xvalue[2:Nx+1,:,:],0)
uw_min = min.(u.xvalue[1:Nx,:,:],0)
uw_max = max.(u.xvalue[1:Nx,:,:],0)
vn_min = min.(u.yvalue[:,2:Ny+1,:],0)
vn_max = max.(u.yvalue[:,2:Ny+1,:],0)
vs_min = min.(u.yvalue[:,1:Ny,:],0)
vs_max = max.(u.yvalue[:,1:Ny,:],0)
wf_min = min.(u.zvalue[:,:,2:Nz+1],0)
wf_max = max.(u.zvalue[:,:,2:Nz+1],0)
wb_min = min.(u.zvalue[:,:,1:Nz],0)
wb_max = max.(u.zvalue[:,:,1:Nz],0)

# calculate the TVD correction term
div_x = -(1.0./DXp).*((ue_max.*psiX_p[2:Nx+1,:,:].+ue_min.*psiX_m[2:Nx+1,:,:])-
              (uw_max.*psiX_p[1:Nx,:,:].+uw_min.*psiX_m[1:Nx,:,:]))
div_y = -(1.0./DYp).*((vn_max.*psiY_p[:,2:Ny+1,:].+vn_min.*psiY_m[:,2:Ny+1,:])-
              (vs_max.*psiY_p[:,1:Ny,:].+vs_min.*psiY_m[:,1:Ny,:]))
div_z = -(1.0./DZp).*((wf_max.*psiZ_p[:,:,2:Nz+1].+wf_min.*psiZ_m[:,:,2:Nz+1])-
              (wb_max.*psiZ_p[:,:,1:Nz].+wb_min.*psiZ_m[:,:,1:Nz]))

# define the RHS Vector
RHS = zeros((Nx+2)*(Ny+2)*(Nz+2))
RHSx = zeros((Nx+2)*(Ny+2)*(Nz+2))
RHSy = zeros((Nx+2)*(Ny+2)*(Nz+2))
RHSz = zeros((Nx+2)*(Ny+2)*(Nz+2))

# assign the values of the RHS vector
mnx = Nx*Ny*Nz
mny = Nx*Ny*Nz
mnz = Nx*Ny*Nz
rowx_index = reshape(G[2:Nx+1,2:Ny+1,2:Nz+1],mnx)  # main diagonal x
rowy_index = reshape(G[2:Nx+1,2:Ny+1,2:Nz+1],mny)  # main diagonal y
rowz_index = reshape(G[2:Nx+1,2:Ny+1,2:Nz+1],mnz)  # main diagonal z
row_index = rowx_index
RHS[row_index] .= reshape(div_x+div_y+div_z,Nx*Ny*Nz)
RHSx[rowx_index] .= reshape(div_x,Nx*Ny*Nz)
RHSy[rowy_index] .= reshape(div_y,Nx*Ny*Nz)
RHSz[rowz_index] .= reshape(div_z,Nx*Ny*Nz)

(RHS, RHSx, RHSy, RHSz)
end

function convectionTvdRHS3D(u::FaceValue, phi::CellValue, FL::Function, u_upwind::FaceValue)
# u is a face variable
# a function to avoid division by zero
eps1 = 1.0e-20
fsign(phi_in) = (abs.(phi_in).>=eps1).*phi_in+eps1*(phi_in.==0.0)+eps1*(abs.(phi_in).<eps1).*sign.(phi_in)
# extract data from the mesh structure
Nx = u.domain.dims[1]
Ny = u.domain.dims[2]
Nz = u.domain.dims[3]
G=reshape([1:(Nx+2)*(Ny+2)*(Nz+2);], Nx+2, Ny+2, Nz+2)
DXp = u.domain.cellsize.x[2:end-1]
DY = zeros( 1, Ny+2)
DY[:] .= u.domain.cellsize.y
DYp = DY[:,2:end-1]
DZ = zeros( 1, 1, Nz+2)
DZ[:] .= u.domain.cellsize.z
DZp = DZ[:,:,2:end-1]
dx=0.5.*(u.domain.cellsize.x[1:end-1].+u.domain.cellsize.x[2:end])
dy=zeros( 1, Ny+1)
dy[:] .=0.5.*(u.domain.cellsize.y[1:end-1].+u.domain.cellsize.y[2:end])
dz=zeros( 1, 1, Nz+1)
dz[:] .=0.5.*(u.domain.cellsize.z[1:end-1].+u.domain.cellsize.z[2:end])
psiX_p = zeros(Nx+1,Ny,Nz)
psiX_m = zeros(Nx+1,Ny,Nz)
psiY_p = zeros(Nx,Ny+1,Nz)
psiY_m = zeros(Nx,Ny+1,Nz)
psiZ_p = zeros(Nx,Ny,Nz+1)
psiZ_m = zeros(Nx,Ny,Nz+1)

# calculate the upstream to downstream gradient ratios for u>0 (+ ratio)
# x direction
dphiX_p = (phi.value[2:Nx+2, 2:Ny+1, 2:Nz+1].-phi.value[1:Nx+1, 2:Ny+1, 2:Nz+1])./dx
rX_p = dphiX_p[1:end-1,:,:]./fsign(dphiX_p[2:end,:,:])
psiX_p[2:Nx+1,:,:] .= 0.5.*FL.(rX_p).*(phi.value[3:Nx+2,2:Ny+1,2:Nz+1].-phi.value[2:Nx+1,2:Ny+1,2:Nz+1])
psiX_p[1,:,:] .= 0.0  # left boundary
# y direction
dphiY_p = (phi.value[2:Nx+1, 2:Ny+2, 2:Nz+1].-phi.value[2:Nx+1, 1:Ny+1, 2:Nz+1])./dy
rY_p = dphiY_p[:,1:end-1,:]./fsign(dphiY_p[:,2:end,:])
psiY_p[:,2:Ny+1,:] .= 0.5.*FL.(rY_p).*(phi.value[2:Nx+1,3:Ny+2,2:Nz+1].-phi.value[2:Nx+1, 2:Ny+1,2:Nz+1])
psiY_p[:,1,:] .= 0.0  # Bottom boundary
# z direction
dphiZ_p = (phi.value[2:Nx+1, 2:Ny+1, 2:Nz+2].-phi.value[2:Nx+1, 2:Ny+1, 1:Nz+1])./dz
rZ_p = dphiZ_p[:,:,1:end-1]./fsign(dphiZ_p[:,:,2:end])
psiZ_p[:,:,2:Nz+1] .= 0.5.*FL.(rZ_p).*(phi.value[2:Nx+1,2:Ny+1,3:Nz+2].-phi.value[2:Nx+1,2:Ny+1,2:Nz+1])
psiZ_p[:,:,1] .= 0.0  # Back boundary

# calculate the upstream to downstream gradient ratios for u<0 (- ratio)
# x direction
rX_m = dphiX_p[2:end,:,:]./fsign(dphiX_p[1:end-1,:,:])
psiX_m[1:Nx,:,:] .= 0.5.*FL.(rX_m).*(phi.value[1:Nx, 2:Ny+1, 2:Nz+1].-phi.value[2:Nx+1, 2:Ny+1, 2:Nz+1])
psiX_m[Nx+1,:,:] .= 0.0  # right boundary
# y direction
rY_m = dphiY_p[:,2:end,:]./fsign(dphiY_p[:,1:end-1,:])
psiY_m[:,1:Ny,:] .= 0.5.*FL.(rY_m).*(phi.value[2:Nx+1,1:Ny,2:Nz+1].-phi.value[2:Nx+1,2:Ny+1,2:Nz+1])
psiY_m[:,Ny+1,:] .= 0.0  # top boundary
# z direction
rZ_m = dphiZ_p[:,:,2:end]./fsign(dphiZ_p[:,:,1:end-1])
psiZ_m[:,:,1:Nz] .= 0.5.*FL.(rZ_m).*(phi.value[2:Nx+1,2:Ny+1,1:Nz].-phi.value[2:Nx+1,2:Ny+1,2:Nz+1])
psiZ_m[:,:,Nz+1] .= 0.0  # front boundary

# find the velocity direction for the upwind scheme
ue_min = u.xvalue[2:Nx+1,:,:]
ue_max = u.xvalue[2:Nx+1,:,:]
uw_min = u.xvalue[1:Nx,:,:]
uw_max = u.xvalue[1:Nx,:,:]
vn_min = u.yvalue[:,2:Ny+1,:]
vn_max = u.yvalue[:,2:Ny+1,:]
vs_min = u.yvalue[:,1:Ny,:]
vs_max = u.yvalue[:,1:Ny,:]
wf_min = u.zvalue[:,:,2:Nz+1]
wf_max = u.zvalue[:,:,2:Nz+1]
wb_min = u.zvalue[:,:,1:Nz]
wb_max = u.zvalue[:,:,1:Nz]

ue_min[u_upwind.xvalue[2:Nx+1,:,:].>0.0] .= 0.0
ue_max[u_upwind.xvalue[2:Nx+1,:,:].<0.0] .= 0.0
uw_min[u_upwind.xvalue[1:Nx,:,:].>0.0] .= 0.0
uw_max[u_upwind.xvalue[1:Nx,:,:].<0.0] .= 0.0
vn_min[u_upwind.yvalue[:,2:Ny+1,:].>0.0] .= 0.0
vn_max[u_upwind.yvalue[:,2:Ny+1,:].<0.0] .= 0.0
vs_min[u_upwind.yvalue[:,1:Ny,:].>0.0] .= 0.0
vs_max[u_upwind.yvalue[:,1:Ny,:].<0.0] .= 0.0
wf_min[u_upwind.zvalue[:,:,2:Nz+1].>0.0] .= 0.0
wf_max[u_upwind.zvalue[:,:,2:Nz+1].<0.0] .= 0.0
wb_min[u_upwind.zvalue[:,:,1:Nz].>0.0] .= 0.0
wb_max[u_upwind.zvalue[:,:,1:Nz].<0.0] .= 0.0

# calculate the TVD correction term
div_x = -(1.0./DXp).*((ue_max.*psiX_p[2:Nx+1,:,:].+ue_min.*psiX_m[2:Nx+1,:,:])-
              (uw_max.*psiX_p[1:Nx,:,:].+uw_min.*psiX_m[1:Nx,:,:]))
div_y = -(1.0./DYp).*((vn_max.*psiY_p[:,2:Ny+1,:].+vn_min.*psiY_m[:,2:Ny+1,:])-
              (vs_max.*psiY_p[:,1:Ny,:].+vs_min.*psiY_m[:,1:Ny,:]))
div_z = -(1.0./DZp).*((wf_max.*psiZ_p[:,:,2:Nz+1].+wf_min.*psiZ_m[:,:,2:Nz+1])-
              (wb_max.*psiZ_p[:,:,1:Nz].+wb_min.*psiZ_m[:,:,1:Nz]))

# define the RHS Vector
RHS = zeros((Nx+2)*(Ny+2)*(Nz+2))
RHSx = zeros((Nx+2)*(Ny+2)*(Nz+2))
RHSy = zeros((Nx+2)*(Ny+2)*(Nz+2))
RHSz = zeros((Nx+2)*(Ny+2)*(Nz+2))

# assign the values of the RHS vector
mnx = Nx*Ny*Nz
mny = Nx*Ny*Nz
mnz = Nx*Ny*Nz
rowx_index = reshape(G[2:Nx+1,2:Ny+1,2:Nz+1],mnx)  # main diagonal x
rowy_index = reshape(G[2:Nx+1,2:Ny+1,2:Nz+1],mny)  # main diagonal y
rowz_index = reshape(G[2:Nx+1,2:Ny+1,2:Nz+1],mnz)  # main diagonal z
row_index = rowx_index
RHS[row_index] .= reshape(div_x+div_y+div_z,Nx*Ny*Nz)
RHSx[rowx_index] .= reshape(div_x,Nx*Ny*Nz)
RHSy[rowy_index] .= reshape(div_y,Nx*Ny*Nz)
RHSz[rowz_index] .= reshape(div_z,Nx*Ny*Nz)

(RHS, RHSx, RHSy, RHSz)
end

# ============================= 3D Cylindrical Convection Term ===============================
function convectionTermCylindrical3D(u::FaceValue)
# u is a face variable
# extract data from the mesh structure
Nr = u.domain.dims[1]
Ntheta = u.domain.dims[2]
Nz = u.domain.dims[3]
G=reshape([1:(Nr+2)*(Ntheta+2)*(Nz+2);], Nr+2, Ntheta+2, Nz+2)
DRe = u.domain.cellsize.x[3:end]
DRw = u.domain.cellsize.x[1:end-2]
DRp = u.domain.cellsize.x[2:end-1]
DTHETA = zeros( 1, Ntheta+2)
DTHETA[:] .= u.domain.cellsize.y
DTHETAn = DTHETA[:,3:end]
DTHETAs = DTHETA[:,1:end-2]
DTHETAp = DTHETA[:,2:end-1]
DZ = zeros( 1, 1, Nz+2)
DZ[:] .= u.domain.cellsize.z
DZf = DZ[:,:,3:end]
DZb = DZ[:,:,1:end-2]
DZp = DZ[:,:,2:end-1]
rp = u.domain.cellcenters.x
rf = u.domain.facecenters.x

# define the vectors to stores the sparse matrix data
iix = zeros(Int64, 3*(Nr+2)*(Ntheta+2)*(Nz+2))
jjx = zeros(Int64, 3*(Nr+2)*(Ntheta+2)*(Nz+2))
sx = zeros(Float64, 3*(Nr+2)*(Ntheta+2)*(Nz+2))
iiy = zeros(Int64, 3*(Nr+2)*(Ntheta+2)*(Nz+2))
jjy = zeros(Int64, 3*(Nr+2)*(Ntheta+2)*(Nz+2))
sy = zeros(Float64, 3*(Nr+2)*(Ntheta+2)*(Nz+2))
iiz = zeros(Int64, 3*(Nr+2)*(Ntheta+2)*(Nz+2))
jjz = zeros(Int64, 3*(Nr+2)*(Ntheta+2)*(Nz+2))
sz = zeros(Float64, 3*(Nr+2)*(Ntheta+2)*(Nz+2))
mnx = Nr*Ntheta*Nz
mny = Nr*Ntheta*Nz
mnz = Nr*Ntheta*Nz

re = rf[2:Nr+1]
rw = rf[1:Nr]

# reassign the east, west, north, and south velocity vectors for the
# code readability
ue = rf[2:Nr+1].*u.xvalue[2:Nr+1,:,:]./(rp.*(DRp+DRe))
uw = rf[1:Nr].*u.xvalue[1:Nr,:,:]./(rp.*(DRp+DRw))
vn = u.yvalue[:,2:Ntheta+1,:]./(rp.*(DTHETAp+DTHETAn))
vs = u.yvalue[:,1:Ntheta,:]./(rp.*(DTHETAp+DTHETAs))
wf = u.zvalue[:,:,2:Nz+1]./(DZp+DZf)
wb = u.zvalue[:,:,1:Nz]./(DZp+DZb)

# calculate the coefficients for the internal cells
AE = reshape(ue,mnx)
AW = reshape(-uw,mnx)
AN = reshape(vn,mny)
AS = reshape(-vs,mny)
AF = reshape(wf,mnz)
AB = reshape(-wb,mnz)
APx = reshape((DRe.*ue-DRw.*uw)./DRp,mnx)
APy = reshape((DTHETAn.*vn-DTHETAs.*vs)./DTHETAp,mny)
APz = reshape((DZf.*wf-DZb.*wb)./DZp,mnz)

# build the sparse matrix based on the numbering system
rowx_index = reshape(G[2:Nr+1,2:Ntheta+1,2:Nz+1],mnx)  # main diagonal x
iix[1:3*mnx] .= repeat(rowx_index,3)
rowy_index = reshape(G[2:Nr+1,2:Ntheta+1,2:Nz+1],mny)  # main diagonal y
iiy[1:3*mny] .= repeat(rowy_index,3)
rowz_index = reshape(G[2:Nr+1,2:Ntheta+1,2:Nz+1],mnz)  # main diagonal z
iiz[1:3*mnz] .= repeat(rowz_index,3)
jjx[1:3*mnx] .= [reshape(G[1:Nr,2:Ntheta+1,2:Nz+1],mnx); reshape(G[2:Nr+1,2:Ntheta+1,2:Nz+1],mnx); reshape(G[3:Nr+2,2:Ntheta+1,2:Nz+1],mnx)]
jjy[1:3*mny] .= [reshape(G[2:Nr+1,1:Ntheta,2:Nz+1],mny); reshape(G[2:Nr+1,2:Ntheta+1,2:Nz+1],mny); reshape(G[2:Nr+1,3:Ntheta+2,2:Nz+1],mny)]
jjz[1:3*mnz] .= [reshape(G[2:Nr+1,2:Ntheta+1,1:Nz],mnz); reshape(G[2:Nr+1,2:Ntheta+1,2:Nz+1],mnz); reshape(G[2:Nr+1,2:Ntheta+1,3:Nz+2],mnz)]
sx[1:3*mnx] .= [AW; APx; AE]
sy[1:3*mny] .= [AS; APy; AN]
sz[1:3*mnz] .= [AB; APz; AF]

# build the sparse matrix
kx = 3*mnx
ky = 3*mny
kz = 3*mnz
Mx = sparse(iix[1:kx], jjx[1:kx], sx[1:kx], (Nr+2)*(Ntheta+2)*(Nz+2), (Nr+2)*(Ntheta+2)*(Nz+2))
My = sparse(iiy[1:ky], jjy[1:ky], sy[1:ky], (Nr+2)*(Ntheta+2)*(Nz+2), (Nr+2)*(Ntheta+2)*(Nz+2))
Mz = sparse(iiz[1:kz], jjz[1:kz], sz[1:kz], (Nr+2)*(Ntheta+2)*(Nz+2), (Nr+2)*(Ntheta+2)*(Nz+2))
M = Mx + My + Mz

(M, Mx, My, Mz)
end


# ============================= 3D Cylindrical Convection Upwind Term ===============================
function convectionUpwindTermCylindrical3D(u::FaceValue)
# u is a face variable
# extract data from the mesh structure
Nr = u.domain.dims[1]
Ntheta = u.domain.dims[2]
Nz = u.domain.dims[3]
G=reshape([1:(Nr+2)*(Ntheta+2)*(Nz+2);], Nr+2, Ntheta+2, Nz+2)
DRe = u.domain.cellsize.x[3:end]
DRw = u.domain.cellsize.x[1:end-2]
DRp = u.domain.cellsize.x[2:end-1]
DTHETA = zeros( 1, Ntheta+2)
DTHETA[:] .= u.domain.cellsize.y
DTHETAn = DTHETA[:,3:end]
DTHETAs = DTHETA[:,1:end-2]
DTHETAp = DTHETA[:,2:end-1]
DZ = zeros( 1, 1, Nz+2)
DZ[:] .= u.domain.cellsize.z
DZf = DZ[:,:,3:end]
DZb = DZ[:,:,1:end-2]
DZp = DZ[:,:,2:end-1]
rp = u.domain.cellcenters.x
rf = u.domain.facecenters.x

# define the vectors to stores the sparse matrix data
iix = zeros(Int64, 3*(Nr+2)*(Ntheta+2)*(Nz+2))
jjx = zeros(Int64, 3*(Nr+2)*(Ntheta+2)*(Nz+2))
sx = zeros(Float64, 3*(Nr+2)*(Ntheta+2)*(Nz+2))
iiy = zeros(Int64, 3*(Nr+2)*(Ntheta+2)*(Nz+2))
jjy = zeros(Int64, 3*(Nr+2)*(Ntheta+2)*(Nz+2))
sy = zeros(Float64, 3*(Nr+2)*(Ntheta+2)*(Nz+2))
iiz = zeros(Int64, 3*(Nr+2)*(Ntheta+2)*(Nz+2))
jjz = zeros(Int64, 3*(Nr+2)*(Ntheta+2)*(Nz+2))
sz = zeros(Float64, 3*(Nr+2)*(Ntheta+2)*(Nz+2))
mnx = Nr*Ntheta*Nz
mny = Nr*Ntheta*Nz
mnz = Nr*Ntheta*Nz

re = rf[2:Nr+1,:]
rw = rf[1:Nr,:]

# find the velocity direction for the upwind scheme
ue_min = min.(u.xvalue[2:Nr+1,:,:],0)
ue_max = max.(u.xvalue[2:Nr+1,:,:],0)
uw_min = min.(u.xvalue[1:Nr,:,:],0)
uw_max = max.(u.xvalue[1:Nr,:,:],0)
vn_min = min.(u.yvalue[:,2:Ntheta+1,:],0)
vn_max = max.(u.yvalue[:,2:Ntheta+1,:],0)
vs_min = min.(u.yvalue[:,1:Ntheta,:],0)
vs_max = max.(u.yvalue[:,1:Ntheta,:],0)
wf_min = min.(u.zvalue[:,:,2:Nz+1],0)
wf_max = max.(u.zvalue[:,:,2:Nz+1],0)
wb_min = min.(u.zvalue[:,:,1:Nz],0)
wb_max = max.(u.zvalue[:,:,1:Nz],0)

# calculate the coefficients for the internal cells
AE = re.*ue_min./(DRp.*rp)
AW = -rw.*uw_max./(DRp.*rp)
AN = vn_min./(DTHETAp.*rp)
AS = -vs_max./(DTHETAp.*rp)
AF = wf_min./DZp
AB = -wb_max./DZp
APx = (re.*ue_max.-rw.*uw_min)./(DRp.*rp)
APy = (vn_max.-vs_min)./(DTHETAp.*rp)
APz = (wf_max.-wb_min)./DZp

# Also correct for the boundary cells (not the ghost cells)
# Left boundary:
APx[1,:,:] .= APx[1,:,:].-rw[1,:,:].*uw_max[1,:,:]./(2.0*DRp[1]*rp[1,:,:])
AW[1,:,:] .= AW[1,:,:]/2.0
# Right boundary:
AE[end,:,:] .= AE[end,:,:]/2.0
APx[end,:,:] .= APx[end,:,:].+re[end,:,:].*ue_min[end,:,:]./(2.0*DRp[end]*rp[end,:,:])
# Bottom boundary:
APy[:,1,:] .= APy[:,1,:].-vs_max[:,1,:]./(2.0*DTHETAp[1]*rp[:,1,:])
AS[:,1,:] .= AS[:,1,:]/2.0
# Top boundary:
AN[:,end,:] .= AN[:,end,:]/2.0
APy[:,end,:] .= APy[:,end,:].+vn_min[:,end,:]./(2.0*DTHETAp[end]*rp[:,end,:])
# Back boundary:
APz[:,:,1] .= APz[:,:,1].-wb_max[:,:,1]/(2.0*DZp[1])
AB[:,:,1] .= AB[:,:,1]/2.0
# Front boundary:
AF[:,:,end] .= AF[:,:,end]/2.0
APz[:,:,end] .= APz[:,:,end].+wf_min[:,:,end]/(2.0*DZp[end])

AE = reshape(AE,mnx)
AW = reshape(AW,mnx)
AN = reshape(AN,mny)
AS = reshape(AS,mny)
AF = reshape(AF,mnz)
AB = reshape(AB,mnz)
APx = reshape(APx,mnx)
APy = reshape(APy,mny)
APz = reshape(APz,mnz)

# build the sparse matrix based on the numbering system
rowx_index = reshape(G[2:Nr+1,2:Ntheta+1,2:Nz+1],mnx)  # main diagonal x
iix[1:3*mnx] .= repeat(rowx_index,3)
rowy_index = reshape(G[2:Nr+1,2:Ntheta+1,2:Nz+1],mny)  # main diagonal y
iiy[1:3*mny] .= repeat(rowy_index,3)
rowz_index = reshape(G[2:Nr+1,2:Ntheta+1,2:Nz+1],mnz)  # main diagonal z
iiz[1:3*mnz] .= repeat(rowz_index,3)
jjx[1:3*mnx] .= [reshape(G[1:Nr,2:Ntheta+1,2:Nz+1],mnx); reshape(G[2:Nr+1,2:Ntheta+1,2:Nz+1],mnx); reshape(G[3:Nr+2,2:Ntheta+1,2:Nz+1],mnx)]
jjy[1:3*mny] .= [reshape(G[2:Nr+1,1:Ntheta,2:Nz+1],mny); reshape(G[2:Nr+1,2:Ntheta+1,2:Nz+1],mny); reshape(G[2:Nr+1,3:Ntheta+2,2:Nz+1],mny)]
jjz[1:3*mnz] .= [reshape(G[2:Nr+1,2:Ntheta+1,1:Nz],mnz); reshape(G[2:Nr+1,2:Ntheta+1,2:Nz+1],mnz); reshape(G[2:Nr+1,2:Ntheta+1,3:Nz+2],mnz)]
sx[1:3*mnx] .= [AW; APx; AE]
sy[1:3*mny] .= [AS; APy; AN]
sz[1:3*mnz] .= [AB; APz; AF]

# build the sparse matrix
kx = 3*mnx
ky = 3*mny
kz = 3*mnz
Mx = sparse(iix[1:kx], jjx[1:kx], sx[1:kx], (Nr+2)*(Ntheta+2)*(Nz+2), (Nr+2)*(Ntheta+2)*(Nz+2))
My = sparse(iiy[1:ky], jjy[1:ky], sy[1:ky], (Nr+2)*(Ntheta+2)*(Nz+2), (Nr+2)*(Ntheta+2)*(Nz+2))
Mz = sparse(iiz[1:kz], jjz[1:kz], sz[1:kz], (Nr+2)*(Ntheta+2)*(Nz+2), (Nr+2)*(Ntheta+2)*(Nz+2))
M = Mx + My + Mz;

(M, Mx, My, Mz)
end

function convectionUpwindTermCylindrical3D(u::FaceValue, u_upwind::FaceValue)
# u is a face variable
# extract data from the mesh structure
Nr = u.domain.dims[1]
Ntheta = u.domain.dims[2]
Nz = u.domain.dims[3]
G=reshape([1:(Nr+2)*(Ntheta+2)*(Nz+2);], Nr+2, Ntheta+2, Nz+2)
DRe = u.domain.cellsize.x[3:end]
DRw = u.domain.cellsize.x[1:end-2]
DRp = u.domain.cellsize.x[2:end-1]
DTHETA = zeros( 1, Ntheta+2)
DTHETA[:] .= u.domain.cellsize.y
DTHETAn = DTHETA[:,3:end]
DTHETAs = DTHETA[:,1:end-2]
DTHETAp = DTHETA[:,2:end-1]
DZ = zeros( 1, 1, Nz+2)
DZ[:] .= u.domain.cellsize.z
DZf = DZ[:,:,3:end]
DZb = DZ[:,:,1:end-2]
DZp = DZ[:,:,2:end-1]
rp = u.domain.cellcenters.x
rf = u.domain.facecenters.x

# define the vectors to stores the sparse matrix data
iix = zeros(Int64, 3*(Nr+2)*(Ntheta+2)*(Nz+2))
jjx = zeros(Int64, 3*(Nr+2)*(Ntheta+2)*(Nz+2))
sx = zeros(Float64, 3*(Nr+2)*(Ntheta+2)*(Nz+2))
iiy = zeros(Int64, 3*(Nr+2)*(Ntheta+2)*(Nz+2))
jjy = zeros(Int64, 3*(Nr+2)*(Ntheta+2)*(Nz+2))
sy = zeros(Float64, 3*(Nr+2)*(Ntheta+2)*(Nz+2))
iiz = zeros(Int64, 3*(Nr+2)*(Ntheta+2)*(Nz+2))
jjz = zeros(Int64, 3*(Nr+2)*(Ntheta+2)*(Nz+2))
sz = zeros(Float64, 3*(Nr+2)*(Ntheta+2)*(Nz+2))
mnx = Nr*Ntheta*Nz
mny = Nr*Ntheta*Nz
mnz = Nr*Ntheta*Nz

re = rf[2:Nr+1,:]
rw = rf[1:Nr,:]

# find the velocity direction for the upwind scheme
ue_min = u.xvalue[2:Nr+1,:,:]
ue_max = u.xvalue[2:Nr+1,:,:]
uw_min = u.xvalue[1:Nr,:,:]
uw_max = u.xvalue[1:Nr,:,:]
vn_min = u.yvalue[:,2:Ntheta+1,:]
vn_max = u.yvalue[:,2:Ntheta+1,:]
vs_min = u.yvalue[:,1:Ntheta,:]
vs_max = u.yvalue[:,1:Ntheta,:]
wf_min = u.zvalue[:,:,2:Nz+1]
wf_max = u.zvalue[:,:,2:Nz+1]
wb_min = u.zvalue[:,:,1:Nz]
wb_max = u.zvalue[:,:,1:Nz]

ue_min[u_upwind.xvalue[2:Nr+1,:,:].> 0.0] .= 0.0
ue_max[u_upwind.xvalue[2:Nr+1,:,:].< 0.0] .= 0.0
uw_min[u_upwind.xvalue[1:Nr,:,:].> 0.0] .= 0.0
uw_max[u_upwind.xvalue[1:Nr,:,:].< 0.0] .= 0.0
vn_min[u_upwind.yvalue[:,2:Ntheta+1,:].> 0.0] .= 0.0
vn_max[u_upwind.yvalue[:,2:Ntheta+1,:].< 0.0] .= 0.0
vs_min[u_upwind.yvalue[:,1:Ntheta,:].> 0.0] .= 0.0
vs_max[u_upwind.yvalue[:,1:Ntheta,:].< 0.0] .= 0.0
wf_min[u_upwind.zvalue[:,:,2:Nz+1].> 0.0] .= 0.0
wf_max[u_upwind.zvalue[:,:,2:Nz+1].< 0.0] .= 0.0
wb_min[u_upwind.zvalue[:,:,1:Nz].> 0.0] .= 0.0
wb_max[u_upwind.zvalue[:,:,1:Nz].< 0.0] .= 0.0

# calculate the coefficients for the internal cells
AE = re.*ue_min./(DRp.*rp)
AW = -rw.*uw_max./(DRp.*rp)
AN = vn_min./(DTHETAp.*rp)
AS = -vs_max./(DTHETAp.*rp)
AF = wf_min./DZp
AB = -wb_max./DZp
APx = (re.*ue_max.-rw.*uw_min)./(DRp.*rp)
APy = (vn_max.-vs_min)./(DTHETAp.*rp)
APz = (wf_max.-wb_min)./DZp

# Also correct for the boundary cells (not the ghost cells)
# Left boundary:
APx[1,:,:] .= APx[1,:,:].-rw[1,:,:].*uw_max[1,:,:]./(2.0*DRp[1]*rp[1,:,:])
AW[1,:,:] .= AW[1,:,:]/2.0
# Right boundary:
AE[end,:,:] .= AE[end,:,:]/2.0
APx[end,:,:] .= APx[end,:,:].+re[end,:,:].*ue_min[end,:,:]./(2.0*DRp[end]*rp[end,:,:])
# Bottom boundary:
APy[:,1,:] .= APy[:,1,:].-vs_max[:,1,:]./(2.0*DTHETAp[1]*rp[:,1,:])
AS[:,1,:] .= AS[:,1,:]/2.0
# Top boundary:
AN[:,end,:] .= AN[:,end,:]/2.0
APy[:,end,:] .= APy[:,end,:].+vn_min[:,end,:]./(2.0*DTHETAp[end]*rp[:,end,:])
# Back boundary:
APz[:,:,1] .= APz[:,:,1].-wb_max[:,:,1]/(2.0*DZp[1])
AB[:,:,1] .= AB[:,:,1]/2.0
# Front boundary:
AF[:,:,end] .= AF[:,:,end]/2.0
APz[:,:,end] .= APz[:,:,end].+wf_min[:,:,end]/(2.0*DZp[end])

AE = reshape(AE,mnx)
AW = reshape(AW,mnx)
AN = reshape(AN,mny)
AS = reshape(AS,mny)
AF = reshape(AF,mnz)
AB = reshape(AB,mnz)
APx = reshape(APx,mnx)
APy = reshape(APy,mny)
APz = reshape(APz,mnz)

# build the sparse matrix based on the numbering system
rowx_index = reshape(G[2:Nr+1,2:Ntheta+1,2:Nz+1],mnx)  # main diagonal x
iix[1:3*mnx] .= repeat(rowx_index,3)
rowy_index = reshape(G[2:Nr+1,2:Ntheta+1,2:Nz+1],mny)  # main diagonal y
iiy[1:3*mny] .= repeat(rowy_index,3)
rowz_index = reshape(G[2:Nr+1,2:Ntheta+1,2:Nz+1],mnz)  # main diagonal z
iiz[1:3*mnz] .= repeat(rowz_index,3)
jjx[1:3*mnx] .= [reshape(G[1:Nr,2:Ntheta+1,2:Nz+1],mnx); reshape(G[2:Nr+1,2:Ntheta+1,2:Nz+1],mnx); reshape(G[3:Nr+2,2:Ntheta+1,2:Nz+1],mnx)]
jjy[1:3*mny] .= [reshape(G[2:Nr+1,1:Ntheta,2:Nz+1],mny); reshape(G[2:Nr+1,2:Ntheta+1,2:Nz+1],mny); reshape(G[2:Nr+1,3:Ntheta+2,2:Nz+1],mny)]
jjz[1:3*mnz] .= [reshape(G[2:Nr+1,2:Ntheta+1,1:Nz],mnz); reshape(G[2:Nr+1,2:Ntheta+1,2:Nz+1],mnz); reshape(G[2:Nr+1,2:Ntheta+1,3:Nz+2],mnz)]
sx[1:3*mnx] .= [AW; APx; AE]
sy[1:3*mny] .= [AS; APy; AN]
sz[1:3*mnz] .= [AB; APz; AF]

# build the sparse matrix
kx = 3*mnx
ky = 3*mny
kz = 3*mnz
Mx = sparse(iix[1:kx], jjx[1:kx], sx[1:kx], (Nr+2)*(Ntheta+2)*(Nz+2), (Nr+2)*(Ntheta+2)*(Nz+2))
My = sparse(iiy[1:ky], jjy[1:ky], sy[1:ky], (Nr+2)*(Ntheta+2)*(Nz+2), (Nr+2)*(Ntheta+2)*(Nz+2))
Mz = sparse(iiz[1:kz], jjz[1:kz], sz[1:kz], (Nr+2)*(Ntheta+2)*(Nz+2), (Nr+2)*(Ntheta+2)*(Nz+2))
M = Mx + My + Mz;

(M, Mx, My, Mz)
end


# ============================= 3D Cylindrical Convection TVD Term ===============================
function convectionTvdTermCylindrical3D(u::FaceValue, phi::CellValue, FL::Function)
# u is a face variable
# a function to avoid division by zero
eps1 = 1.0e-20
fsign(phi_in) = (abs.(phi_in).>=eps1).*phi_in+eps1*(phi_in.==0.0)+eps1*(abs.(phi_in).<eps1).*sign.(phi_in)
# extract data from the mesh structure
Nr = u.domain.dims[1]
Ntheta = u.domain.dims[2]
Nz = u.domain.dims[3]
G=reshape([1:(Nr+2)*(Ntheta+2)*(Nz+2);], Nr+2, Ntheta+2, Nz+2)
DRe = u.domain.cellsize.x[3:end]
DRw = u.domain.cellsize.x[1:end-2]
DRp = u.domain.cellsize.x[2:end-1]
DTHETA = zeros( 1, Ntheta+2)
DTHETA[:] .= u.domain.cellsize.y
DTHETAn = DTHETA[:,3:end]
DTHETAs = DTHETA[:,1:end-2]
DTHETAp = DTHETA[:,2:end-1]
DZ = zeros( 1, 1, Nz+2)
DZ[:] .= u.domain.cellsize.z
DZf = DZ[:,:,3:end]
DZb = DZ[:,:,1:end-2]
DZp = DZ[:,:,2:end-1]
dr=0.5.*(u.domain.cellsize.x[1:end-1].+u.domain.cellsize.x[2:end])
dtheta = zeros( 1, Ntheta+1)
dtheta[:] .=0.5.*(u.domain.cellsize.y[1:end-1].+u.domain.cellsize.y[2:end])
dz = zeros( 1, 1, Nz+1)
dz[:] .=0.5.*(u.domain.cellsize.z[1:end-1].+u.domain.cellsize.z[2:end])
psiX_p = zeros(Nr+1,Ntheta,Nz)
psiX_m = zeros(Nr+1,Ntheta,Nz)
psiY_p = zeros(Nr,Ntheta+1,Nz)
psiY_m = zeros(Nr,Ntheta+1,Nz)
psiZ_p = zeros(Nr,Ntheta,Nz+1)
psiZ_m = zeros(Nr,Ntheta,Nz+1)
rp = u.domain.cellcenters.x
rf = u.domain.facecenters.x

# define the vectors to stores the sparse matrix data
iix = zeros(Int64, 3*(Nr+2)*(Ntheta+2)*(Nz+2))
jjx = zeros(Int64, 3*(Nr+2)*(Ntheta+2)*(Nz+2))
sx = zeros(Float64, 3*(Nr+2)*(Ntheta+2)*(Nz+2))
iiy = zeros(Int64, 3*(Nr+2)*(Ntheta+2)*(Nz+2))
jjy = zeros(Int64, 3*(Nr+2)*(Ntheta+2)*(Nz+2))
sy = zeros(Float64, 3*(Nr+2)*(Ntheta+2)*(Nz+2))
iiz = zeros(Int64, 3*(Nr+2)*(Ntheta+2)*(Nz+2))
jjz = zeros(Int64, 3*(Nr+2)*(Ntheta+2)*(Nz+2))
sz = zeros(Float64, 3*(Nr+2)*(Ntheta+2)*(Nz+2))
mnx = Nr*Ntheta*Nz
mny = Nr*Ntheta*Nz
mnz = Nr*Ntheta*Nz

# calculate the upstream to downstream gradient ratios for u>0 (+ ratio)
# x direction
dphiX_p = (phi.value[2:Nr+2, 2:Ntheta+1, 2:Nz+1].-phi.value[1:Nr+1, 2:Ntheta+1, 2:Nz+1])./dr
rX_p = dphiX_p[1:end-1,:,:]./fsign(dphiX_p[2:end,:,:])
psiX_p[2:Nr+1,:,:] .= 0.5.*FL.(rX_p).*(phi.value[3:Nr+2,2:Ntheta+1,2:Nz+1].-phi.value[2:Nr+1,2:Ntheta+1,2:Nz+1])
psiX_p[1,:,:] .= 0  # left boundary
# y direction
dphiY_p = (phi.value[2:Nr+1, 2:Ntheta+2, 2:Nz+1].-phi.value[2:Nr+1, 1:Ntheta+1, 2:Nz+1])./dtheta
rY_p = dphiY_p[:,1:end-1,:]./fsign(dphiY_p[:,2:end,:])
psiY_p[:,2:Ntheta+1,:] .= 0.5.*FL.(rY_p).*(phi.value[2:Nr+1,3:Ntheta+2,2:Nz+1].-phi.value[2:Nr+1, 2:Ntheta+1,2:Nz+1])
psiY_p[:,1,:] .= 0.0  # Bottom boundary
# z direction
dphiZ_p = (phi.value[2:Nr+1, 2:Ntheta+1, 2:Nz+2].-phi.value[2:Nr+1, 2:Ntheta+1, 1:Nz+1])./dz
rZ_p = dphiZ_p[:,:,1:end-1]./fsign(dphiZ_p[:,:,2:end])
psiZ_p[:,:,2:Nz+1] .= 0.5.*FL.(rZ_p).*(phi.value[2:Nr+1,2:Ntheta+1,3:Nz+2].-phi.value[2:Nr+1,2:Ntheta+1,2:Nz+1])
psiZ_p[:,:,1] .= 0.0  # Back boundary

# calculate the upstream to downstream gradient ratios for u<0 (- ratio)
# x direction
rX_m = dphiX_p[2:end,:,:]./fsign(dphiX_p[1:end-1,:,:])
psiX_m[1:Nr,:,:] .= 0.5.*FL.(rX_m).*(phi.value[1:Nr, 2:Ntheta+1, 2:Nz+1].-phi.value[2:Nr+1, 2:Ntheta+1, 2:Nz+1])
psiX_m[Nr+1,:,:] .= 0.0  # right boundary
# y direction
rY_m = dphiY_p[:,2:end,:]./fsign(dphiY_p[:,1:end-1,:])
psiY_m[:,1:Ntheta,:] .= 0.5.*FL.(rY_m).*(phi.value[2:Nr+1,1:Ntheta,2:Nz+1].-phi.value[2:Nr+1,2:Ntheta+1,2:Nz+1])
psiY_m[:,Ntheta+1,:] .= 0.0  # top boundary
# z direction
rZ_m = dphiZ_p[:,:,2:end]./fsign(dphiZ_p[:,:,1:end-1])
psiZ_m[:,:,1:Nz] .= 0.5.*FL.(rZ_m).*(phi.value[2:Nr+1,2:Ntheta+1,1:Nz].-phi.value[2:Nr+1,2:Ntheta+1,2:Nz+1])
psiZ_m[:,:,Nz+1] .= 0.0  # front boundary

re = rf[2:Nr+1]
rw = rf[1:Nr]

# find the velocity direction for the upwind scheme
ue_min = min.(u.xvalue[2:Nr+1,:,:],0)
ue_max = max.(u.xvalue[2:Nr+1,:,:],0)
uw_min = min.(u.xvalue[1:Nr,:,:],0)
uw_max = max.(u.xvalue[1:Nr,:,:],0)
vn_min = min.(u.yvalue[:,2:Ntheta+1,:],0)
vn_max = max.(u.yvalue[:,2:Ntheta+1,:],0)
vs_min = min.(u.yvalue[:,1:Ntheta,:],0)
vs_max = max.(u.yvalue[:,1:Ntheta,:],0)
wf_min = min.(u.zvalue[:,:,2:Nz+1],0)
wf_max = max.(u.zvalue[:,:,2:Nz+1],0)
wb_min = min.(u.zvalue[:,:,1:Nz],0)
wb_max = max.(u.zvalue[:,:,1:Nz],0)

# calculate the coefficients for the internal cells
AE = re.*ue_min./(DRp.*rp)
AW = -rw.*uw_max./(DRp.*rp)
AN = vn_min./(DTHETAp.*rp)
AS = -vs_max./(DTHETAp.*rp)
AF = wf_min./DZp
AB = -wb_max./DZp
APx = (re.*ue_max.-rw.*uw_min)./(DRp.*rp)
APy = (vn_max.-vs_min)./(DTHETAp.*rp)
APz = (wf_max.-wb_min)./DZp

# Also correct for the boundary cells (not the ghost cells)
# Left boundary:
APx[1,:,:] .= APx[1,:,:].-rw[1,:,:].*uw_max[1,:,:]./(2.0*DRp[1]*rp[1,:,:])
AW[1,:,:] .= AW[1,:,:]/2.0
# Right boundary:
AE[end,:,:] .= AE[end,:,:]/2.0
APx[end,:,:] .= APx[end,:,:].+re[end,:,:].*ue_min[end,:,:]./(2.0*DRp[end]*rp[end,:,:])
# Bottom boundary:
APy[:,1,:] .= APy[:,1,:].-vs_max[:,1,:]./(2.0*DTHETAp[1]*rp[:,1,:])
AS[:,1,:] .= AS[:,1,:]/2.0
# Top boundary:
AN[:,end,:] .= AN[:,end,:]/2.0
APy[:,end,:] .= APy[:,end,:].+vn_min[:,end,:]./(2.0*DTHETAp[end]*rp[:,end,:])
# Back boundary:
APz[:,:,1] .= APz[:,:,1].-wb_max[:,:,1]/(2.0*DZp[1])
AB[:,:,1] .= AB[:,:,1]/2.0
# Front boundary:
AF[:,:,end] .= AF[:,:,end]/2.0
APz[:,:,end] .= APz[:,:,end].+wf_min[:,:,end]/(2.0*DZp[end])

AE = reshape(AE,mnx)
AW = reshape(AW,mnx)
AN = reshape(AN,mny)
AS = reshape(AS,mny)
AF = reshape(AF,mnz)
AB = reshape(AB,mnz)
APx = reshape(APx,mnx)
APy = reshape(APy,mny)
APz = reshape(APz,mnz)

# build the sparse matrix based on the numbering system
rowx_index = reshape(G[2:Nr+1,2:Ntheta+1,2:Nz+1],mnx)  # main diagonal x
iix[1:3*mnx] .= repeat(rowx_index,3)
rowy_index = reshape(G[2:Nr+1,2:Ntheta+1,2:Nz+1],mny)  # main diagonal y
iiy[1:3*mny] .= repeat(rowy_index,3)
rowz_index = reshape(G[2:Nr+1,2:Ntheta+1,2:Nz+1],mnz)  # main diagonal z
iiz[1:3*mnz] .= repeat(rowz_index,3)
jjx[1:3*mnx] .= [reshape(G[1:Nr,2:Ntheta+1,2:Nz+1],mnx); reshape(G[2:Nr+1,2:Ntheta+1,2:Nz+1],mnx); reshape(G[3:Nr+2,2:Ntheta+1,2:Nz+1],mnx)]
jjy[1:3*mny] .= [reshape(G[2:Nr+1,1:Ntheta,2:Nz+1],mny); reshape(G[2:Nr+1,2:Ntheta+1,2:Nz+1],mny); reshape(G[2:Nr+1,3:Ntheta+2,2:Nz+1],mny)]
jjz[1:3*mnz] .= [reshape(G[2:Nr+1,2:Ntheta+1,1:Nz],mnz); reshape(G[2:Nr+1,2:Ntheta+1,2:Nz+1],mnz); reshape(G[2:Nr+1,2:Ntheta+1,3:Nz+2],mnz)]
sx[1:3*mnx] .= [AW; APx; AE]
sy[1:3*mny] .= [AS; APy; AN]
sz[1:3*mnz] .= [AB; APz; AF]

# calculate the TVD correction term
div_x = -(1.0./(DRp.*rp)).*(re.*(ue_max.*psiX_p[2:Nr+1,:,:].+ue_min.*psiX_m[2:Nr+1,:,:])-
              rw.*(uw_max.*psiX_p[1:Nr,:,:].+uw_min.*psiX_m[1:Nr,:,:]))
div_y = -(1.0./(DTHETAp.*rp)).*((vn_max.*psiY_p[:,2:Ntheta+1,:].+vn_min.*psiY_m[:,2:Ntheta+1,:])-
              (vs_max.*psiY_p[:,1:Ntheta,:].+vs_min.*psiY_m[:,1:Ntheta,:]))
div_z = -(1.0./DZp).*((wf_max.*psiZ_p[:,:,2:Nz+1].+wf_min.*psiZ_m[:,:,2:Nz+1])-
              (wb_max.*psiZ_p[:,:,1:Nz].+wb_min.*psiZ_m[:,:,1:Nz]))

# define the RHS Vector
RHS = zeros((Nr+2)*(Ntheta+2)*(Nz+2))
RHSx = zeros((Nr+2)*(Ntheta+2)*(Nz+2))
RHSy = zeros((Nr+2)*(Ntheta+2)*(Nz+2))
RHSz = zeros((Nr+2)*(Ntheta+2)*(Nz+2))

# assign the values of the RHS vector
row_index = rowx_index
RHS[row_index] .= reshape(div_x+div_y+div_z,Nr*Ntheta*Nz)
RHSx[rowx_index] .= reshape(div_x,Nr*Ntheta*Nz)
RHSy[rowy_index] .= reshape(div_y,Nr*Ntheta*Nz)
RHSz[rowz_index] .= reshape(div_z,Nr*Ntheta*Nz)

# build the sparse matrix
kx = 3*mnx
ky = 3*mny
kz = 3*mnz
Mx = sparse(iix[1:kx], jjx[1:kx], sx[1:kx], (Nr+2)*(Ntheta+2)*(Nz+2), (Nr+2)*(Ntheta+2)*(Nz+2))
My = sparse(iiy[1:ky], jjy[1:ky], sy[1:ky], (Nr+2)*(Ntheta+2)*(Nz+2), (Nr+2)*(Ntheta+2)*(Nz+2))
Mz = sparse(iiz[1:kz], jjz[1:kz], sz[1:kz], (Nr+2)*(Ntheta+2)*(Nz+2), (Nr+2)*(Ntheta+2)*(Nz+2))
M = Mx + My + Mz;

(M, RHS, Mx, My, Mz, RHSx, RHSy, RHSz)
end


# ============================= 3D Cylindrical Convection TVD RHS ===============================
function convectionTvdRHSCylindrical3D(u::FaceValue, phi::CellValue, FL::Function)
# u is a face variable
# a function to avoid division by zero
eps1 = 1.0e-20
fsign(phi_in) = (abs.(phi_in).>=eps1).*phi_in+eps1*(phi_in.==0.0)+eps1*(abs.(phi_in).<eps1).*sign.(phi_in)
# extract data from the mesh structure
Nr = u.domain.dims[1]
Ntheta = u.domain.dims[2]
Nz = u.domain.dims[3]
G=reshape([1:(Nr+2)*(Ntheta+2)*(Nz+2);], Nr+2, Ntheta+2, Nz+2)
DRe = u.domain.cellsize.x[3:end]
DRw = u.domain.cellsize.x[1:end-2]
DRp = u.domain.cellsize.x[2:end-1]
DTHETA = zeros( 1, Ntheta+2)
DTHETA[:] .= u.domain.cellsize.y
DTHETAn = DTHETA[:,3:end]
DTHETAs = DTHETA[:,1:end-2]
DTHETAp = DTHETA[:,2:end-1]
DZ = zeros( 1, 1, Nz+2)
DZ[:] .= u.domain.cellsize.z
DZf = DZ[:,:,3:end]
DZb = DZ[:,:,1:end-2]
DZp = DZ[:,:,2:end-1]
dr=0.5.*(u.domain.cellsize.x[1:end-1].+u.domain.cellsize.x[2:end])
dtheta = zeros( 1, Ntheta+1)
dtheta[:] .=0.5.*(u.domain.cellsize.y[1:end-1].+u.domain.cellsize.y[2:end])
dz = zeros( 1, 1, Nz+1)
dz[:] .=0.5.*(u.domain.cellsize.z[1:end-1].+u.domain.cellsize.z[2:end])
psiX_p = zeros(Nr+1,Ntheta,Nz)
psiX_m = zeros(Nr+1,Ntheta,Nz)
psiY_p = zeros(Nr,Ntheta+1,Nz)
psiY_m = zeros(Nr,Ntheta+1,Nz)
psiZ_p = zeros(Nr,Ntheta,Nz+1)
psiZ_m = zeros(Nr,Ntheta,Nz+1)
rp = u.domain.cellcenters.x
rf = u.domain.facecenters.x

# calculate the upstream to downstream gradient ratios for u>0 (+ ratio)
# x direction
dphiX_p = (phi.value[2:Nr+2, 2:Ntheta+1, 2:Nz+1].-phi.value[1:Nr+1, 2:Ntheta+1, 2:Nz+1])./dr
rX_p = dphiX_p[1:end-1,:,:]./fsign(dphiX_p[2:end,:,:])
psiX_p[2:Nr+1,:,:] .= 0.5.*FL.(rX_p).*(phi.value[3:Nr+2,2:Ntheta+1,2:Nz+1].-phi.value[2:Nr+1,2:Ntheta+1,2:Nz+1])
psiX_p[1,:,:] .= 0.0  # left boundary
# y direction
dphiY_p = (phi.value[2:Nr+1, 2:Ntheta+2, 2:Nz+1].-phi.value[2:Nr+1, 1:Ntheta+1, 2:Nz+1])./dtheta
rY_p = dphiY_p[:,1:end-1,:]./fsign(dphiY_p[:,2:end,:])
psiY_p[:,2:Ntheta+1,:] .= 0.5.*FL.(rY_p).*(phi.value[2:Nr+1,3:Ntheta+2,2:Nz+1].-phi.value[2:Nr+1, 2:Ntheta+1,2:Nz+1])
psiY_p[:,1,:] .= 0.0  # Bottom boundary
# z direction
dphiZ_p = (phi.value[2:Nr+1, 2:Ntheta+1, 2:Nz+2].-phi.value[2:Nr+1, 2:Ntheta+1, 1:Nz+1])./dz
rZ_p = dphiZ_p[:,:,1:end-1]./fsign(dphiZ_p[:,:,2:end])
psiZ_p[:,:,2:Nz+1] .= 0.5.*FL.(rZ_p).*(phi.value[2:Nr+1,2:Ntheta+1,3:Nz+2].-phi.value[2:Nr+1,2:Ntheta+1,2:Nz+1])
psiZ_p[:,:,1] .= 0.0  # Back boundary

# calculate the upstream to downstream gradient ratios for u<0 (- ratio)
# x direction
rX_m = dphiX_p[2:end,:,:]./fsign(dphiX_p[1:end-1,:,:])
psiX_m[1:Nr,:,:] .= 0.5.*FL.(rX_m).*(phi.value[1:Nr, 2:Ntheta+1, 2:Nz+1].-phi.value[2:Nr+1, 2:Ntheta+1, 2:Nz+1])
psiX_m[Nr+1,:,:] .= 0.0  # right boundary
# y direction
rY_m = dphiY_p[:,2:end,:]./fsign(dphiY_p[:,1:end-1,:])
psiY_m[:,1:Ntheta,:] .= 0.5.*FL.(rY_m).*(phi.value[2:Nr+1,1:Ntheta,2:Nz+1].-phi.value[2:Nr+1,2:Ntheta+1,2:Nz+1])
psiY_m[:,Ntheta+1,:] .= 0.0  # top boundary
# z direction
rZ_m = dphiZ_p[:,:,2:end]./fsign(dphiZ_p[:,:,1:end-1])
psiZ_m[:,:,1:Nz] .= 0.5.*FL.(rZ_m).*(phi.value[2:Nr+1,2:Ntheta+1,1:Nz].-phi.value[2:Nr+1,2:Ntheta+1,2:Nz+1])
psiZ_m[:,:,Nz+1] .= 0.0  # front boundary

re = rf[2:Nr+1]
rw = rf[1:Nr]

# find the velocity direction for the upwind scheme
ue_min = min.(u.xvalue[2:Nr+1,:,:],0)
ue_max = max.(u.xvalue[2:Nr+1,:,:],0)
uw_min = min.(u.xvalue[1:Nr,:,:],0)
uw_max = max.(u.xvalue[1:Nr,:,:],0)
vn_min = min.(u.yvalue[:,2:Ntheta+1,:],0)
vn_max = max.(u.yvalue[:,2:Ntheta+1,:],0)
vs_min = min.(u.yvalue[:,1:Ntheta,:],0)
vs_max = max.(u.yvalue[:,1:Ntheta,:],0)
wf_min = min.(u.zvalue[:,:,2:Nz+1],0)
wf_max = max.(u.zvalue[:,:,2:Nz+1],0)
wb_min = min.(u.zvalue[:,:,1:Nz],0)
wb_max = max.(u.zvalue[:,:,1:Nz],0)

# calculate the TVD correction term
div_x = -(1.0./(DRp.*rp)).*(re.*(ue_max.*psiX_p[2:Nr+1,:,:].+ue_min.*psiX_m[2:Nr+1,:,:])-
              rw.*(uw_max.*psiX_p[1:Nr,:,:].+uw_min.*psiX_m[1:Nr,:,:]))
div_y = -(1.0./(DTHETAp.*rp)).*((vn_max.*psiY_p[:,2:Ntheta+1,:].+vn_min.*psiY_m[:,2:Ntheta+1,:])-
              (vs_max.*psiY_p[:,1:Ntheta,:].+vs_min.*psiY_m[:,1:Ntheta,:]))
div_z = -(1.0./DZp).*((wf_max.*psiZ_p[:,:,2:Nz+1].+wf_min.*psiZ_m[:,:,2:Nz+1])-
              (wb_max.*psiZ_p[:,:,1:Nz].+wb_min.*psiZ_m[:,:,1:Nz]))

# define the RHS Vector
RHS = zeros((Nr+2)*(Ntheta+2)*(Nz+2))
RHSx = zeros((Nr+2)*(Ntheta+2)*(Nz+2))
RHSy = zeros((Nr+2)*(Ntheta+2)*(Nz+2))
RHSz = zeros((Nr+2)*(Ntheta+2)*(Nz+2))

# assign the values of the RHS vector
mnx = Nr*Ntheta*Nz
mny = Nr*Ntheta*Nz
mnz = Nr*Ntheta*Nz
rowx_index = reshape(G[2:Nr+1,2:Ntheta+1,2:Nz+1],mnx)  # main diagonal x
rowy_index = reshape(G[2:Nr+1,2:Ntheta+1,2:Nz+1],mny)  # main diagonal y
rowz_index = reshape(G[2:Nr+1,2:Ntheta+1,2:Nz+1],mnz)  # main diagonal z
row_index = rowx_index
RHS[row_index] .= reshape(div_x+div_y+div_z,Nr*Ntheta*Nz)
RHSx[rowx_index] .= reshape(div_x,Nr*Ntheta*Nz)
RHSy[rowy_index] .= reshape(div_y,Nr*Ntheta*Nz)
RHSz[rowz_index] .= reshape(div_z,Nr*Ntheta*Nz)

(RHS, RHSx, RHSy, RHSz)
end

function convectionTvdRHSCylindrical3D(u::FaceValue, phi::CellValue, FL::Function, u_upwind::FaceValue)
# u is a face variable
# a function to avoid division by zero
eps1 = 1.0e-20
fsign(phi_in) = (abs.(phi_in).>=eps1).*phi_in+eps1*(phi_in.==0.0)+eps1*(abs.(phi_in).<eps1).*sign.(phi_in)
# extract data from the mesh structure
Nr = u.domain.dims[1]
Ntheta = u.domain.dims[2]
Nz = u.domain.dims[3]
G=reshape([1:(Nr+2)*(Ntheta+2)*(Nz+2);], Nr+2, Ntheta+2, Nz+2)
DRe = u.domain.cellsize.x[3:end]
DRw = u.domain.cellsize.x[1:end-2]
DRp = u.domain.cellsize.x[2:end-1]
DTHETA = zeros( 1, Ntheta+2)
DTHETA[:] .= u.domain.cellsize.y
DTHETAn = DTHETA[:,3:end]
DTHETAs = DTHETA[:,1:end-2]
DTHETAp = DTHETA[:,2:end-1]
DZ = zeros( 1, 1, Nz+2)
DZ[:] .= u.domain.cellsize.z
DZf = DZ[:,:,3:end]
DZb = DZ[:,:,1:end-2]
DZp = DZ[:,:,2:end-1]
dr=0.5.*(u.domain.cellsize.x[1:end-1].+u.domain.cellsize.x[2:end])
dtheta = zeros( 1, Ntheta+1)
dtheta[:] .=0.5.*(u.domain.cellsize.y[1:end-1].+u.domain.cellsize.y[2:end])
dz = zeros( 1, 1, Nz+1)
dz[:] .=0.5.*(u.domain.cellsize.z[1:end-1].+u.domain.cellsize.z[2:end])
psiX_p = zeros(Nr+1,Ntheta,Nz)
psiX_m = zeros(Nr+1,Ntheta,Nz)
psiY_p = zeros(Nr,Ntheta+1,Nz)
psiY_m = zeros(Nr,Ntheta+1,Nz)
psiZ_p = zeros(Nr,Ntheta,Nz+1)
psiZ_m = zeros(Nr,Ntheta,Nz+1)
rp = u.domain.cellcenters.x
rf = u.domain.facecenters.x

# calculate the upstream to downstream gradient ratios for u>0 (+ ratio)
# x direction
dphiX_p = (phi.value[2:Nr+2, 2:Ntheta+1, 2:Nz+1].-phi.value[1:Nr+1, 2:Ntheta+1, 2:Nz+1])./dr
rX_p = dphiX_p[1:end-1,:,:]./fsign(dphiX_p[2:end,:,:])
psiX_p[2:Nr+1,:,:] .= 0.5.*FL.(rX_p).*(phi.value[3:Nr+2,2:Ntheta+1,2:Nz+1].-phi.value[2:Nr+1,2:Ntheta+1,2:Nz+1])
psiX_p[1,:,:] .= 0  # left boundary
# y direction
dphiY_p = (phi.value[2:Nr+1, 2:Ntheta+2, 2:Nz+1].-phi.value[2:Nr+1, 1:Ntheta+1, 2:Nz+1])./dtheta
rY_p = dphiY_p[:,1:end-1,:]./fsign(dphiY_p[:,2:end,:])
psiY_p[:,2:Ntheta+1,:] .= 0.5.*FL.(rY_p).*(phi.value[2:Nr+1,3:Ntheta+2,2:Nz+1].-phi.value[2:Nr+1, 2:Ntheta+1,2:Nz+1])
psiY_p[:,1,:] .= 0.0  # Bottom boundary
# z direction
dphiZ_p = (phi.value[2:Nr+1, 2:Ntheta+1, 2:Nz+2].-phi.value[2:Nr+1, 2:Ntheta+1, 1:Nz+1])./dz
rZ_p = dphiZ_p[:,:,1:end-1]./fsign(dphiZ_p[:,:,2:end])
psiZ_p[:,:,2:Nz+1] .= 0.5.*FL.(rZ_p).*(phi.value[2:Nr+1,2:Ntheta+1,3:Nz+2].-phi.value[2:Nr+1,2:Ntheta+1,2:Nz+1])
psiZ_p[:,:,1] .= 0.0  # Back boundary

# calculate the upstream to downstream gradient ratios for u<0 (- ratio)
# x direction
rX_m = dphiX_p[2:end,:,:]./fsign(dphiX_p[1:end-1,:,:])
psiX_m[1:Nr,:,:] .= 0.5.*FL.(rX_m).*(phi.value[1:Nr, 2:Ntheta+1, 2:Nz+1].-phi.value[2:Nr+1, 2:Ntheta+1, 2:Nz+1])
psiX_m[Nr+1,:,:] .= 0.0  # right boundary
# y direction
rY_m = dphiY_p[:,2:end,:]./fsign(dphiY_p[:,1:end-1,:])
psiY_m[:,1:Ntheta,:] .= 0.5.*FL.(rY_m).*(phi.value[2:Nr+1,1:Ntheta,2:Nz+1].-phi.value[2:Nr+1,2:Ntheta+1,2:Nz+1])
psiY_m[:,Ntheta+1,:] .= 0.0  # top boundary
# z direction
rZ_m = dphiZ_p[:,:,2:end]./fsign(dphiZ_p[:,:,1:end-1])
psiZ_m[:,:,1:Nz] .= 0.5.*FL.(rZ_m).*(phi.value[2:Nr+1,2:Ntheta+1,1:Nz].-phi.value[2:Nr+1,2:Ntheta+1,2:Nz+1])
psiZ_m[:,:,Nz+1] .= 0.0  # front boundary

re = rf[2:Nr+1]
rw = rf[1:Nr]

# find the velocity direction for the upwind scheme
ue_min = u.xvalue[2:Nr+1,:,:]
ue_max = u.xvalue[2:Nr+1,:,:]
uw_min = u.xvalue[1:Nr,:,:]
uw_max = u.xvalue[1:Nr,:,:]
vn_min = u.yvalue[:,2:Ntheta+1,:]
vn_max = u.yvalue[:,2:Ntheta+1,:]
vs_min = u.yvalue[:,1:Ntheta,:]
vs_max = u.yvalue[:,1:Ntheta,:]
wf_min = u.zvalue[:,:,2:Nz+1]
wf_max = u.zvalue[:,:,2:Nz+1]
wb_min = u.zvalue[:,:,1:Nz]
wb_max = u.zvalue[:,:,1:Nz]

ue_min[u_upwind.xvalue[2:Nr+1,:,:].> 0.0] .= 0.0
ue_max[u_upwind.xvalue[2:Nr+1,:,:].< 0.0] .= 0.0
uw_min[u_upwind.xvalue[1:Nr,:,:].> 0.0] .= 0.0
uw_max[u_upwind.xvalue[1:Nr,:,:].< 0.0] .= 0.0
vn_min[u_upwind.yvalue[:,2:Ntheta+1,:].> 0.0] .= 0.0
vn_max[u_upwind.yvalue[:,2:Ntheta+1,:].< 0.0] .= 0.0
vs_min[u_upwind.yvalue[:,1:Ntheta,:].> 0.0] .= 0.0
vs_max[u_upwind.yvalue[:,1:Ntheta,:].< 0.0] .= 0.0
wf_min[u_upwind.zvalue[:,:,2:Nz+1].> 0.0] .= 0.0
wf_max[u_upwind.zvalue[:,:,2:Nz+1].< 0.0] .= 0.0
wb_min[u_upwind.zvalue[:,:,1:Nz].> 0.0] .= 0.0
wb_max[u_upwind.zvalue[:,:,1:Nz].< 0.0] .= 0.0

# calculate the TVD correction term
div_x = -(1.0./(DRp.*rp)).*(re.*(ue_max.*psiX_p[2:Nr+1,:,:].+ue_min.*psiX_m[2:Nr+1,:,:])-
              rw.*(uw_max.*psiX_p[1:Nr,:,:].+uw_min.*psiX_m[1:Nr,:,:]))
div_y = -(1.0./(DTHETAp.*rp)).*((vn_max.*psiY_p[:,2:Ntheta+1,:].+vn_min.*psiY_m[:,2:Ntheta+1,:])-
              (vs_max.*psiY_p[:,1:Ntheta,:].+vs_min.*psiY_m[:,1:Ntheta,:]))
div_z = -(1.0./DZp).*((wf_max.*psiZ_p[:,:,2:Nz+1].+wf_min.*psiZ_m[:,:,2:Nz+1])-
              (wb_max.*psiZ_p[:,:,1:Nz].+wb_min.*psiZ_m[:,:,1:Nz]))

# define the RHS Vector
RHS = zeros((Nr+2)*(Ntheta+2)*(Nz+2))
RHSx = zeros((Nr+2)*(Ntheta+2)*(Nz+2))
RHSy = zeros((Nr+2)*(Ntheta+2)*(Nz+2))
RHSz = zeros((Nr+2)*(Ntheta+2)*(Nz+2))

# assign the values of the RHS vector
mnx = Nr*Ntheta*Nz
mny = Nr*Ntheta*Nz
mnz = Nr*Ntheta*Nz
rowx_index = reshape(G[2:Nr+1,2:Ntheta+1,2:Nz+1],mnx)  # main diagonal x
rowy_index = reshape(G[2:Nr+1,2:Ntheta+1,2:Nz+1],mny)  # main diagonal y
rowz_index = reshape(G[2:Nr+1,2:Ntheta+1,2:Nz+1],mnz)  # main diagonal z
row_index = rowx_index
RHS[row_index] .= reshape(div_x+div_y+div_z,Nr*Ntheta*Nz)
RHSx[rowx_index] .= reshape(div_x,Nr*Ntheta*Nz)
RHSy[rowy_index] .= reshape(div_y,Nr*Ntheta*Nz)
RHSz[rowz_index] .= reshape(div_z,Nr*Ntheta*Nz)

(RHS, RHSx, RHSy, RHSz)
end
