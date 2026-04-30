# ===============================
# Written by AAE
# TU Delft, Spring 2014
# simulkade.com
# ===============================

function transientTerm(phi_old::CellValue, dt::Real, alfa::Real)
transientTerm(phi_old, dt, alfa.*ones(Float64,tuple(phi_old.domain.dims...)))
end

function transientTerm(phi_old::CellValue, dt::Real)
transientTerm(phi_old, dt, 1.0)
end

function transientTerm(phi_old::CellValue, dt::Real, alfa::CellValue)
d = phi_old.domain.dimension
if d==1 || d==1.5
  transientTerm1D(phi_old, dt, alfa.value[2:end-1])
elseif d==2 || d==2.5 || d==2.8
  transientTerm2D(phi_old, dt, alfa.value[2:end-1,2:end-1])
elseif d==3 || d==3.2
  transientTerm3D(phi_old, dt, alfa.value[2:end-1,2:end-1,2:end-1])
end
end


function transientTerm(phi_old::CellValue, dt::Real, alfa::Array{T}) where T<:Real
d = phi_old.domain.dimension
if d==1 || d==1.5
  transientTerm1D(phi_old, dt, alfa)
elseif d==2 || d==2.5 || d==2.8
  transientTerm2D(phi_old, dt, alfa)
elseif d==3 || d==3.2
  transientTerm3D(phi_old, dt, alfa)
end

end

function transientTerm1D(phi_old::CellValue,
		    dt::Real, alfa::Array{T}) where T<:Real
# returns the matrix and RHS for a d(phi)/dt term

# extract data from the mesh structure
Nx = phi_old.domain.dims[1]
G = [1:Nx+2;]

# rearrange the matrix of k and build the sparse matrix for internal cells
row_index = reshape(G[2:Nx+1],Nx) # main diagonal (only internal cells)
AP_diag = reshape(alfa/dt,Nx)
M = sparse(row_index, row_index, AP_diag, Nx+2, Nx+2)

# define the RHS Vector
RHS = zeros(Nx+2)

# assign the values of the RHS vector
RHS[row_index] .= reshape(alfa.*phi_old.value[2:Nx+1]/dt,Nx)

(M, RHS)

end

function transientTerm2D(phi_old::CellValue,
		    dt::Real, alfa::Array{T}) where T<:Real
# returns the matrix and RHS for a d(phi)/dt term

# extract data from the mesh structure
Nx = phi_old.domain.dims[1]
Ny = phi_old.domain.dims[2]
G=reshape([1:(Nx+2)*(Ny+2);], Nx+2, Ny+2)

# rearrange the matrix of k and build the sparse matrix for internal cells
row_index = reshape(G[2:Nx+1,2:Ny+1],Nx*Ny) # main diagonal (only internal cells)
AP_diag = reshape(alfa/dt,Nx*Ny)
M = sparse(row_index, row_index, AP_diag, (Nx+2)*(Ny+2), (Nx+2)*(Ny+2))

# define the RHS Vector
RHS = zeros((Nx+2)*(Ny+2))

# assign the values of the RHS vector
RHS[row_index] .= reshape(alfa.*phi_old.value[2:Nx+1,2:Ny+1]/dt,Nx*Ny)

(M, RHS)

end

function transientTerm3D(phi_old::CellValue,
		    dt::Real, alfa::Array{T}) where T<:Real
# returns the matrix and RHS for a d(phi)/dt term

# extract data from the mesh structure
Nx = phi_old.domain.dims[1]
Ny = phi_old.domain.dims[2]
Nz = phi_old.domain.dims[3]
G=reshape([1:(Nx+2)*(Ny+2)*(Nz+2);], Nx+2, Ny+2, Nz+2)

# rearrange the matrix of k and build the sparse matrix for internal cells
row_index = reshape(G[2:Nx+1,2:Ny+1,2:Nz+1],Nx*Ny*Nz) # main diagonal (only internal cells)
AP_diag = reshape(alfa/dt,Nx*Ny*Nz)
M = sparse(row_index, row_index, AP_diag, (Nx+2)*(Ny+2)*(Nz+2), (Nx+2)*(Ny+2)*(Nz+2))

# define the RHS Vector
RHS = zeros((Nx+2)*(Ny+2)*(Nz+2))

# assign the values of the RHS vector
RHS[row_index] .= reshape(alfa.*phi_old.value[2:Nx+1,2:Ny+1,2:Nz+1]./dt,Nx*Ny*Nz)

(M, RHS)

end
