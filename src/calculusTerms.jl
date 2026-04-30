# ===============================
# Written by AAE
# TU Delft, Spring 2014
# simulkade.com
# ===============================

# ================================================================
# Changes:
#    2014-12-30 added 2D radial and 3D cylindrical grids
#    2014-12-31 many many bugs fixed
#    2015-01-10 extended to accept nonuniform grids
# ================================================================


function divergenceTerm(F::FaceValue)
  divergenceTerm(F, F.domain.coordinatesystem)

end

divergenceTerm(F::FaceValue, ::Cartesian1D) = divergenceTerm1D(F)
divergenceTerm(F::FaceValue, ::Cylindrical1D) = divergenceTermCylindrical1D(F)
divergenceTerm(F::FaceValue, ::Cartesian2D) = first(divergenceTerm2D(F))
divergenceTerm(F::FaceValue, ::Cylindrical2D) = first(divergenceTermCylindrical2D(F))
divergenceTerm(F::FaceValue, ::Radial2D) = first(divergenceTermRadial2D(F))
divergenceTerm(F::FaceValue, ::Cartesian3D) = first(divergenceTerm3D(F))
divergenceTerm(F::FaceValue, ::Cylindrical3D) = first(divergenceTermCylindrical3D(F))


# =============== Divergence 1D Term ============================
function divergenceTerm1D(F::FaceValue)
# This function calculates the divergence of a field
# using its face

# extract data from the mesh structure
Nx = F.domain.dims[1]
G = [1:Nx+2;]
DX = F.domain.cellsize.x[2:end-1]

# define the vector of cell index
row_index = reshape(G[2:Nx+1],Nx) # main diagonal

# compute the divergence
div_x = (F.xvalue[2:Nx+1].-F.xvalue[1:Nx])./DX

# define the RHS Vector
RHSdiv = zeros(Nx+2)

# assign the values of the RHS vector
RHSdiv[row_index] = reshape(div_x,Nx)

return RHSdiv

end


# =============== Divergence Cylindrical 1D Term ============================
function divergenceTermCylindrical1D(F::FaceValue)
# This function calculates the divergence of a field
# using its face

# extract data from the mesh structure
Nx = F.domain.dims[1]
G = [1:Nx+2;]
DX = F.domain.cellsize.x[2:end-1]
rp = F.domain.cellcenters.x
rf = F.domain.facecenters.x

# define the vector of cell index
row_index = reshape(G[2:Nx+1],Nx) # main diagonal

# reassign the east, west, north, and south flux vectors for the
# code readability
Fe = F.xvalue[2:Nx+1]
Fw = F.xvalue[1:Nx]
re = rf[2:Nx+1]
rw = rf[1:Nx]

# compute the divergence
div_x = (re.*Fe-rw.*Fw)./(DX.*rp)

# define the RHS Vector
RHSdiv = zeros(Nx+2)

# assign the values of the RHS vector
RHSdiv[row_index] = reshape(div_x,Nx)

return RHSdiv

end

# =============== Divergence 2D Term ============================
function divergenceTerm2D(F::FaceValue)
# This function calculates the divergence of a field
# using its face

# extract data from the mesh structure
Nx = F.domain.dims[1]
Ny = F.domain.dims[2]
G=reshape([1:(Nx+2)*(Ny+2);], Nx+2, Ny+2)
DX = F.domain.cellsize.x[2:end-1]
DY = zeros( 1, Ny)
DY[:] = F.domain.cellsize.y[2:end-1]

# define the vector of cell index
row_index = reshape(G[2:Nx+1,2:Ny+1],Nx*Ny) # main diagonal

# reassign the east, west, north, and south flux vectors for the
# code readability
Fe = F.xvalue[2:Nx+1,:]
Fw = F.xvalue[1:Nx,:]
Fn = F.yvalue[:,2:Ny+1]
Fs = F.yvalue[:,1:Ny]

# compute the divergence
div_x = (Fe - Fw)./DX
div_y = (Fn - Fs)./DY

# define the RHS Vector
RHSdiv = zeros((Nx+2)*(Ny+2))
RHSdivx = zeros((Nx+2)*(Ny+2))
RHSdivy = zeros((Nx+2)*(Ny+2))

# assign the values of the RHS vector
RHSdiv[row_index] = reshape(div_x+div_y,Nx*Ny)
RHSdivx[row_index] = reshape(div_x,Nx*Ny)
RHSdivy[row_index] = reshape(div_y,Nx*Ny)

(RHSdiv, RHSdivx, RHSdivy)

end



# =============== Divergence 2D Cylindrical Term ============================
function divergenceTermCylindrical2D(F::FaceValue)
# This function calculates the divergence of a field
# using its face

# extract data from the mesh structure
Nr = F.domain.dims[1]
Nz = F.domain.dims[2]
G=reshape([1:(Nr+2)*(Nz+2);], Nr+2, Nz+2)
dr = F.domain.cellsize.x[2:end-1]
dz= zeros( 1, Nz)
dz[:] = F.domain.cellsize.y[2:end-1]
rp = F.domain.cellcenters.x
rf = F.domain.facecenters.x

# define the vector of cell index
row_index = reshape(G[2:Nr+1,2:Nz+1],Nr*Nz) # main diagonal

# reassign the east, west, north, and south flux vectors for the
# code readability
Fe = F.xvalue[2:Nr+1,:]
Fw = F.xvalue[1:Nr,:]
Fn = F.yvalue[:,2:Nz+1]
Fs = F.yvalue[:,1:Nz]
re = rf[2:Nr+1]
rw = rf[1:Nr]

# compute the divergence
div_x = (re.*Fe - rw.*Fw)./(dr.*rp)
div_y = (Fn - Fs)./dz

# define the RHS Vector
RHSdiv = zeros((Nr+2)*(Nz+2))
RHSdivx = zeros((Nr+2)*(Nz+2))
RHSdivy = zeros((Nr+2)*(Nz+2))

# assign the values of the RHS vector
RHSdiv[row_index] = reshape(div_x+div_y,Nr*Nz)
RHSdivx[row_index] = reshape(div_x,Nr*Nz)
RHSdivy[row_index] = reshape(div_y,Nr*Nz)

(RHSdiv, RHSdivx, RHSdivy)

end



# =============== Divergence 2D Radial Term ============================
function divergenceTermRadial2D(F::FaceValue)
# This function calculates the divergence of a field
# using its face

# extract data from the mesh structure
Nr = F.domain.dims[1]
Ntheta = F.domain.dims[2]
G=reshape([1:(Nr+2)*(Ntheta+2);], Nr+2, Ntheta+2)
dr = F.domain.cellsize.x[2:end-1]
dtheta= zeros( 1, Ntheta)
dtheta[:]= F.domain.cellsize.y[2:end-1]
rp = F.domain.cellcenters.x
rf = F.domain.facecenters.x

# define the vector of cell index
row_index = reshape(G[2:Nr+1,2:Ntheta+1],Nr*Ntheta) # main diagonal

# reassign the east, west, north, and south flux vectors for the
# code readability
Fe = F.xvalue[2:Nr+1,:]
Fw = F.xvalue[1:Nr,:]
Fn = F.yvalue[:,2:Ntheta+1]
Fs = F.yvalue[:,1:Ntheta]
re = rf[2:Nr+1]
rw = rf[1:Nr]

# compute the divergence
div_x = (re.*Fe-rw.*Fw)./(dr.*rp)
div_y = (Fn-Fs)./(dtheta.*rp)

# define the RHS Vector
RHSdiv = zeros((Nr+2)*(Ntheta+2))
RHSdivx = zeros((Nr+2)*(Ntheta+2))
RHSdivy = zeros((Nr+2)*(Ntheta+2))

# assign the values of the RHS vector
RHSdiv[row_index] = reshape(div_x+div_y,Nr*Ntheta)
RHSdivx[row_index] = reshape(div_x,Nr*Ntheta)
RHSdivy[row_index] = reshape(div_y,Nr*Ntheta)

(RHSdiv, RHSdivx, RHSdivy)

end

# =============== Divergence 3D Term ============================
function divergenceTerm3D(F::FaceValue)
# This function calculates the divergence of a field
# using its face

# extract data from the mesh structure
Nx = F.domain.dims[1]
Ny = F.domain.dims[2]
Nz = F.domain.dims[3]
G=reshape([1:(Nx+2)*(Ny+2)*(Nz+2);], Nx+2, Ny+2, Nz+2)
dx = F.domain.cellsize.x[2:end-1]
dy = zeros( 1, Ny)
dy[:] = F.domain.cellsize.y[2:end-1]
dz = zeros( 1,1,Nz)
dz[:] = F.domain.cellsize.z[2:end-1]

# define the vector of cell index
row_index = reshape(G[2:Nx+1,2:Ny+1,2:Nz+1],Nx*Ny*Nz) # main diagonal

# reassign the east, west, north, and south flux vectors for the
# code readability
Fe = F.xvalue[2:Nx+1,:,:]
Fw = F.xvalue[1:Nx,:,:]
Fn = F.yvalue[:,2:Ny+1,:]
Fs = F.yvalue[:,1:Ny,:]
Ff = F.zvalue[:,:,2:Nz+1]
Fb = F.zvalue[:,:,1:Nz]

# compute the divergence
div_x = (Fe - Fw)./dx
div_y = (Fn - Fs)./dy
div_z = (Ff - Fb)./dz

# define the RHS Vector
RHSdiv = zeros((Nx+2)*(Ny+2)*(Nz+2))
RHSdivx = zeros((Nx+2)*(Ny+2)*(Nz+2))
RHSdivy = zeros((Nx+2)*(Ny+2)*(Nz+2))
RHSdivz = zeros((Nx+2)*(Ny+2)*(Nz+2))

# assign the values of the RHS vector
RHSdiv[row_index] = reshape(div_x+div_y+div_z,Nx*Ny*Nz)
RHSdivx[row_index] = reshape(div_x,Nx*Ny*Nz)
RHSdivy[row_index] = reshape(div_y,Nx*Ny*Nz)
RHSdivz[row_index] = reshape(div_z,Nx*Ny*Nz)

(RHSdiv, RHSdivx, RHSdivy, RHSdivz)

end



# =============== Divergence 3D Cylindrical Term ============================
function divergenceTermCylindrical3D(F::FaceValue)
# This function calculates the divergence of a field
# using its face

# extract data from the mesh structure
Nx = F.domain.dims[1]
Ny = F.domain.dims[2]
Nz = F.domain.dims[3]
G=reshape([1:(Nx+2)*(Ny+2)*(Nz+2);], Nx+2, Ny+2, Nz+2)
dx = F.domain.cellsize.x[2:end-1]
dy = zeros( 1, Ny)
dy[:] = F.domain.cellsize.y[2:end-1]
dz = zeros( 1,1,Nz)
dz[:] = F.domain.cellsize.z[2:end-1]
rp = F.domain.cellcenters.x

# define the vector of cell index
row_index = reshape(G[2:Nx+1,2:Ny+1,2:Nz+1],Nx*Ny*Nz) # main diagonal

# reassign the east, west, north, and south flux vectors for the
# code readability
Fe = F.xvalue[2:Nx+1,:,:]
Fw = F.xvalue[1:Nx,:,:]
Fn = F.yvalue[:,2:Ny+1,:]
Fs = F.yvalue[:,1:Ny,:]
Ff = F.zvalue[:,:,2:Nz+1]
Fb = F.zvalue[:,:,1:Nz]

# compute the divergence
div_x = (Fe - Fw)./dx
div_y = (Fn - Fs)./(dy.*rp)
div_z = (Ff - Fb)./dz

# define the RHS Vector
RHSdiv = zeros((Nx+2)*(Ny+2)*(Nz+2))
RHSdivx = zeros((Nx+2)*(Ny+2)*(Nz+2))
RHSdivy = zeros((Nx+2)*(Ny+2)*(Nz+2))
RHSdivz = zeros((Nx+2)*(Ny+2)*(Nz+2))

# assign the values of the RHS vector
RHSdiv[row_index] = reshape(div_x+div_y+div_z,Nx*Ny*Nz)
RHSdivx[row_index] = reshape(div_x,Nx*Ny*Nz)
RHSdivy[row_index] = reshape(div_y,Nx*Ny*Nz)
RHSdivz[row_index] = reshape(div_z,Nx*Ny*Nz)

(RHSdiv, RHSdivx, RHSdivy, RHSdivz)

end


# ===================== Gradient Term ============================
function gradientTerm(phi)
  # calculates the gradient of a variable
  # the output is a face variable
  cs = phi.domain.coordinatesystem
  if is_1d(cs)
    dx = 0.5*(phi.domain.cellsize.x[1:end-1].+phi.domain.cellsize.x[2:end])
    FaceValue(phi.domain,
      (phi.value[2:end].-phi.value[1:end-1])./dx,
      [1.0],
      [1.0])
  elseif (cs isa Cartesian2D) || (cs isa Cylindrical2D)
    dx = 0.5*(phi.domain.cellsize.x[1:end-1].+phi.domain.cellsize.x[2:end])
    Ny = phi.domain.dims[2]
    dy = zeros( 1, Ny+1)
    dy[:] = 0.5*(phi.domain.cellsize.y[1:end-1].+phi.domain.cellsize.y[2:end])
    FaceValue(phi.domain,
      (phi.value[2:end,2:end-1].-phi.value[1:end-1,2:end-1])./dx,
      (phi.value[2:end-1,2:end].-phi.value[2:end-1,1:end-1])./dy,
      [1.0])
  elseif cs isa Radial2D
    dx = 0.5*(phi.domain.cellsize.x[1:end-1].+phi.domain.cellsize.x[2:end])
    Ntheta = phi.domain.dims[2]
    dtheta = zeros( 1, Ntheta+1)
    dtheta[:] = 0.5*(phi.domain.cellsize.y[1:end-1].+phi.domain.cellsize.y[2:end])
    rp = phi.domain.cellcenters.x
    FaceValue(phi.domain,
      (phi.value[2:end,2:end-1].-phi.value[1:end-1,2:end-1])./dx,
      (phi.value[2:end-1,2:end].-phi.value[2:end-1,1:end-1])./(dtheta.*rp),
      [1.0])
  elseif cs isa Cartesian3D
    Ny = phi.domain.dims[2]
    Nz = phi.domain.dims[3]
    dx = 0.5*(phi.domain.cellsize.x[1:end-1].+phi.domain.cellsize.x[2:end])
    dy= zeros( 1, Ny+1)
    dy[:] = 0.5*(phi.domain.cellsize.y[1:end-1].+phi.domain.cellsize.y[2:end])
    dz= zeros( 1, 1, Nz+1)
    dz[:] = 0.5*(phi.domain.cellsize.z[1:end-1].+phi.domain.cellsize.z[2:end])
    FaceValue(phi.domain,
      (phi.value[2:end,2:end-1,2:end-1].-phi.value[1:end-1,2:end-1,2:end-1])./dx,
      (phi.value[2:end-1,2:end,2:end-1].-phi.value[2:end-1,1:end-1,2:end-1])./dy,
      (phi.value[2:end-1,2:end-1,2:end].-phi.value[2:end-1,2:end-1,1:end-1])./dz)
  elseif cs isa Cylindrical3D
    Ntheta = phi.domain.dims[2]
    Nz = phi.domain.dims[3]
    dx = 0.5*(phi.domain.cellsize.x[1:end-1].+phi.domain.cellsize.x[2:end])
    dy= zeros( 1, Ntheta+1)
    dy[:] = 0.5*(phi.domain.cellsize.y[1:end-1].+phi.domain.cellsize.y[2:end])
    dz= zeros( 1, 1, Nz+1)
    dz[:] = 0.5*(phi.domain.cellsize.z[1:end-1].+phi.domain.cellsize.z[2:end])
    rp = phi.domain.cellcenters.x
    FaceValue(phi.domain,
      (phi.value[2:end,2:end-1,2:end-1].-phi.value[1:end-1,2:end-1,2:end-1])./dx,
      (phi.value[2:end-1,2:end,2:end-1].-phi.value[2:end-1,1:end-1,2:end-1])./(dy.*rp),
      (phi.value[2:end-1,2:end-1,2:end].-phi.value[2:end-1,2:end-1,1:end-1])./dz)

  end
end


# ===================== Gradient Term ============================
function gradientCellTerm(phi)
  # calculates the gradient of a variable
  # the output is a face variable
  cs = phi.domain.coordinatesystem
  phi_face = linearMean(phi)
  if is_1d(cs)
    dx = phi.domain.cellsize.x[2:end-1]
    CellVector(phi.domain,
      (phi_face.xvalue[2:end] .- phi_face.xvalue[1:end-1])./dx,
      [1.0],
      [1.0])
  elseif (cs isa Cartesian2D) || (cs isa Cylindrical2D)
    dx = phi.domain.cellsize.x[2:end-1]
    Ny = phi.domain.dims[2]
    dy = zeros(1, Ny)
    dy[:] = phi.domain.cellsize.y[2:end-1]
    CellVector(phi.domain,
      (phi_face.xvalue[2:end, :] .- phi_face.xvalue[1:end-1, :])./dx,
      (phi_face.yvalue[:, 2:end] .- phi_face.yvalue[:, 1:end-1])./dy,
      [1.0])
  elseif cs isa Radial2D
    dx = phi.domain.cellsize.x[2:end-1]
    Ntheta = phi.domain.dims[2]
    dtheta = zeros(1, Ntheta)
    dtheta[:] = phi.domain.cellsize.y[2:end-1]
    rp = phi.domain.cellcenters.x
    CellVector(phi.domain,
      (phi_face.xvalue[2:end, :] .- phi_face.xvalue[1:end-1, :])./dx,
      (phi_face.yvalue[:, 2:end] .- phi_face.yvalue[:, 1:end-1])./(dtheta.*rp),
      [1.0])
  elseif cs isa Cartesian3D
    Ny = phi.domain.dims[2]
    Nz = phi.domain.dims[3]
    dx = phi.domain.cellsize.x[2:end-1]
    dy= zeros(1, Ny)
    dy[:] = phi.domain.cellsize.y[2:end-1]
    dz= zeros( 1, 1, Nz)
    dz[:] = phi.domain.cellsize.z[2:end-1]
    CellVector(phi.domain,
      (phi_face.xvalue[2:end, :, :] .- phi_face.xvalue[1:end-1, :, :])./dx,
      (phi_face.yvalue[:, 2:end, :] .- phi_face.yvalue[:, 1:end-1, :])./dy,
      (phi_face.zvalue[:, :, 2:end] .- phi_face.zvalue[:, :, 1:end-1])./dz)
  elseif cs isa Cylindrical3D
    Ntheta = phi.domain.dims[2]
    Nz = phi.domain.dims[3]
    dx = phi.domain.cellsize.x[2:end-1]
    dy= zeros(1, Ntheta)
    dy[:] = phi.domain.cellsize.y[2:end-1]
    dz= zeros(1, 1, Nz)
    dz[:] = phi.domain.cellsize.z[2:end-1]
    rp = phi.domain.cellcenters.x
    CellVector(phi.domain,
      (phi_face.xvalue[2:end, :, :] .- phi_face.xvalue[1:end-1, :, :])./dx,
      (phi_face.yvalue[:, 2:end, :] .- phi_face.yvalue[:, 1:end-1, :])./(dy.*rp),
      (phi_face.zvalue[:, :, 2:end] .- phi_face.zvalue[:, :, 1:end-1])./dz)
  end
end