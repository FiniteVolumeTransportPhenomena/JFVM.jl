# ===============================
# Written by AAE
# TU Delft, Spring 2014
# simulkade.com
# ===============================

# ================================================================
# Changes:
#    2015-01-10 changed numberofcells to dims
# ================================================================

# ======================= Linear source term ========================
function linearSourceTerm(betta0::CellValue)
m = betta0.domain
cs = m.coordinatesystem
if is_1d(cs)
  Nx = m.dims[1]
  G = [1:Nx+2;]
  b = betta0.value[2:end-1]
  row_index = reshape(G[2:Nx+1],Nx)  # main diagonal (only internal cells)
  AP_diag = reshape(b,Nx)
  sparse(row_index, row_index, AP_diag, Nx+2, Nx+2)
elseif is_2d(cs)
  Nx = m.dims[1]
  Ny = m.dims[2]
  G=reshape([1:(Nx+2)*(Ny+2);], Nx+2, Ny+2)
  b = betta0.value[2:end-1,2:end-1]
  row_index = reshape(G[2:Nx+1,2:Ny+1],Nx*Ny)  # main diagonal (only internal cells)
  AP_diag = reshape(b,Nx*Ny)
  sparse(row_index, row_index, AP_diag, (Nx+2)*(Ny+2), (Nx+2)*(Ny+2))
elseif is_3d(cs)
  Nx = m.dims[1]
  Ny = m.dims[2]
  Nz = m.dims[3]
  G=reshape([1:(Nx+2)*(Ny+2)*(Nz+2);], Nx+2, Ny+2, Nz+2)
  b = betta0.value[2:end-1,2:end-1,2:end-1]
  row_index = reshape(G[2:Nx+1,2:Ny+1,2:Nz+1],Nx*Ny*Nz)  # main diagonal (only internal cells)
  AP_diag = reshape(b,Nx*Ny*Nz)
  sparse(row_index, row_index, AP_diag, (Nx+2)*(Ny+2)*(Nz+2), (Nx+2)*(Ny+2)*(Nz+2))
end
end


function linearSourceTerm(m::MeshStructure, betta0::Real)
linearSourceTerm(createCellVariable(m, betta0))
end


function linearSourceTerm(m::MeshStructure, betta0::AbstractArray{T}) where T<:Real
cs = m.coordinatesystem
if is_1d(cs)
  Nx = m.dims[1]
  G = [1:Nx+2;]
  row_index = reshape(G[2:Nx+1],Nx)  # main diagonal (only internal cells)
  AP_diag = reshape(betta0,Nx)
  sparse(row_index, row_index, AP_diag, Nx+2, Nx+2)
elseif is_2d(cs)
  Nx = m.dims[1]
  Ny = m.dims[2]
  G=reshape([1:(Nx+2)*(Ny+2);], Nx+2, Ny+2)
  row_index = reshape(G[2:Nx+1,2:Ny+1],Nx*Ny)  # main diagonal (only internal cells)
  AP_diag = reshape(betta0,Nx*Ny)
  sparse(row_index, row_index, AP_diag, (Nx+2)*(Ny+2), (Nx+2)*(Ny+2))
elseif is_3d(cs)
  Nx = m.dims[1]
  Ny = m.dims[2]
  Nz = m.dims[3]
  G=reshape([1:(Nx+2)*(Ny+2)*(Nz+2);], Nx+2, Ny+2, Nz+2)
  row_index = reshape(G[2:Nx+1,2:Ny+1,2:Nz+1],Nx*Ny*Nz)  # main diagonal (only internal cells)
  AP_diag = reshape(betta0,Nx*Ny*Nz)
  sparse(row_index, row_index, AP_diag, (Nx+2)*(Ny+2)*(Nz+2), (Nx+2)*(Ny+2)*(Nz+2))
end
end


# ================================== constant source term ================================
function constantSourceTerm(phi0::CellValue)
m = phi0.domain
cs = m.coordinatesystem
if is_1d(cs)
  Nx = m.dims[1]
  G = [1:Nx+2;]
  row_index = reshape(G[2:Nx+1],Nx)  # main diagonal (only internal cells)
  RHS = zeros(Nx+2)
  RHS[row_index] .= reshape(phi0.value[2:end-1],Nx)
elseif is_2d(cs)
  Nx = m.dims[1]
  Ny = m.dims[2]
  G=reshape([1:(Nx+2)*(Ny+2);], Nx+2, Ny+2)
  row_index = reshape(G[2:Nx+1,2:Ny+1],Nx*Ny)  # main diagonal (only internal cells)
  RHS = zeros((Nx+2)*(Ny+2))
  RHS[row_index] .= reshape(phi0.value[2:end-1,2:end-1],Nx*Ny)
elseif is_3d(cs)
  Nx = m.dims[1]
  Ny = m.dims[2]
  Nz = m.dims[3]
  G=reshape([1:(Nx+2)*(Ny+2)*(Nz+2);], Nx+2, Ny+2, Nz+2)
  row_index = reshape(G[2:Nx+1,2:Ny+1,2:Nz+1],Nx*Ny*Nz)  # main diagonal (only internal cells)
  RHS = zeros((Nx+2)*(Ny+2)*(Nz+2))
  RHS[row_index] .= reshape(phi0.value[2:end-1,2:end-1,2:end-1],Nx*Ny*Nz)
end
RHS
end


function constantSourceTerm(m::MeshStructure, phi0::Real)
constantSourceTerm(createCellVariable(m, phi0))
end


function constantSourceTerm(m::MeshStructure, phi0::AbstractArray{T}) where T<:Real
cs = m.coordinatesystem
if is_1d(cs)
  Nx = m.dims[1]
  G = [1:Nx+2;]
  row_index = reshape(G[2:Nx+1],Nx)  # main diagonal (only internal cells)
  RHS = zeros(Nx+2)
  RHS[row_index] .= reshape(phi0,Nx)
elseif is_2d(cs)
  Nx = m.dims[1]
  Ny = m.dims[2]
  G=reshape([1:(Nx+2)*(Ny+2);], Nx+2, Ny+2)
  row_index = reshape(G[2:Nx+1,2:Ny+1],Nx*Ny)  # main diagonal (only internal cells)
  RHS = zeros((Nx+2)*(Ny+2))
  RHS[row_index] .= reshape(phi0,Nx*Ny)
elseif is_3d(cs)
  Nx = m.dims[1]
  Ny = m.dims[2]
  Nz = m.dims[3]
  G=reshape([1:(Nx+2)*(Ny+2)*(Nz+2);], Nx+2, Ny+2, Nz+2)
  row_index = reshape(G[2:Nx+1,2:Ny+1,2:Nz+1],Nx*Ny*Nz)  # main diagonal (only internal cells)
  RHS = zeros((Nx+2)*(Ny+2)*(Nz+2))
  RHS[row_index] .= reshape(phi0,Nx*Ny*Nz)
end
RHS
end
