# ===============================
# Written by AAE
# TU Delft, Spring 2014
# simulkade.com
# ===============================

# ================================================================
# Changes:
#    2014-12-30 added 2D radial and 3D cylindrical grids
#    2015-01-10 extended to accept nonuniform grids
# ================================================================

# ========================= DIFFUSION TERM ================================
function diffusionTerm(D::FaceValue)
	diffusionTerm(D, D.domain.coordinatesystem)
end

diffusionTerm(D::FaceValue, ::Cartesian1D) = diffusionTerm1D(D)
diffusionTerm(D::FaceValue, ::Cylindrical1D) = diffusionTermCylindrical1D(D)
diffusionTerm(D::FaceValue, ::Cartesian2D) = first(diffusionTerm2D(D))
diffusionTerm(D::FaceValue, ::Cylindrical2D) = first(diffusionTermCylindrical2D(D))
diffusionTerm(D::FaceValue, ::Radial2D) = first(diffusionTermRadial2D(D))
diffusionTerm(D::FaceValue, ::Cartesian3D) = first(diffusionTerm3D(D))
diffusionTerm(D::FaceValue, ::Cylindrical3D) = first(diffusionTermCylindrical3D(D))

# ======================== 1D CARTESIAN DIFFUSION =========================
function diffusionTerm1D(D::FaceValue)
# D is a face variable

# extract data from the mesh structure
Nx = D.domain.dims[1]
G = [1:Nx+2;]
DX = D.domain.cellsize.x
dx = 0.5*(DX[1:end-1].+DX[2:end])

# define the vectors to store the sparse matrix data
iix = zeros(Int64, 3*(Nx+2))
jjx = zeros(Int64, 3*(Nx+2))
sx = zeros(Float64, 3*(Nx+2))

# reassign the east, west for code readability
De = D.xvalue[2:Nx+1]./(dx[2:Nx+1].*DX[2:Nx+1])
Dw = D.xvalue[1:Nx]./(dx[1:Nx].*DX[2:Nx+1])

# calculate the coefficients for the internal cells
AE = reshape(De,Nx)
AW = reshape(Dw,Nx)
APx = -(AE+AW)

# build the sparse matrix based on the numbering system
rowx_index = reshape(G[2:Nx+1],Nx) # main diagonal x
iix[1:3*Nx] .= repeat(rowx_index,3)
jjx[1:3*Nx] .= [reshape(G[1:Nx],Nx); reshape(G[2:Nx+1],Nx); reshape(G[3:Nx+2],Nx)]
sx[1:3*Nx] .= [AW; APx; AE]

# build the sparse matrix
kx = 3*Nx
M = sparse(iix[1:kx], jjx[1:kx], sx[1:kx], Nx+2, Nx+2)
end

# ======================== 2D CARTESIAN DIFFUSION =========================
function diffusionTerm2D(D::FaceValue)
# D is a face variable
# extract data from the mesh structure
Nx = D.domain.dims[1]
Ny = D.domain.dims[2]
G=reshape([1:(Nx+2)*(Ny+2);], Nx+2, Ny+2)

DX = D.domain.cellsize.x
DY = zeros( 1, Ny+2)
DY[:] .= D.domain.cellsize.y
dx = 0.5*(DX[1:end-1].+DX[2:end])
dy = zeros( 1, Ny+1)
dy[:] .= 0.5*(DY[1:end-1].+DY[2:end])

# define the vectors to store the sparse matrix data
iix = zeros(Int64, 3*(Nx+2)*(Ny+2))
iiy = zeros(Int64, 3*(Nx+2)*(Ny+2))
jjx = zeros(Int64, 3*(Nx+2)*(Ny+2))
jjy = zeros(Int64, 3*(Nx+2)*(Ny+2))
sx = zeros(Float64, 3*(Nx+2)*(Ny+2))
sy = zeros(Float64, 3*(Nx+2)*(Ny+2))
mnx = Nx*Ny
mny = Nx*Ny

# reassign the east, west for code readability (use broadcasting in Julia)
De = D.xvalue[2:Nx+1,:]./(dx[2:Nx+1].*DX[2:Nx+1])
Dw = D.xvalue[1:Nx,:]./(dx[1:Nx].*DX[2:Nx+1])
Dn = D.yvalue[:,2:Ny+1]./(dy[:,2:Ny+1].*DY[:,2:Ny+1])
Ds = D.yvalue[:,1:Ny]./(dy[:,1:Ny].*DY[:,2:Ny+1])

# calculate the coefficients for the internal cells
AE = reshape(De,mnx)
AW = reshape(Dw,mnx)
AN = reshape(Dn,mny)
AS = reshape(Ds,mny)
APx = -(AE+AW)
APy = -(AN+AS)

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

# ======================== 3D CARTESIAN DIFFUSION =========================
function diffusionTerm3D(D::FaceValue)
# D is a face variable
# extract data from the mesh structure
Nx = D.domain.dims[1]
Ny = D.domain.dims[2]
Nz = D.domain.dims[3]
G=reshape([1:(Nx+2)*(Ny+2)*(Nz+2);], Nx+2, Ny+2, Nz+2)

# DX = repeat(D.domain.cellsize.x, outer=(1, Ny, Nz))
# DY = repeat(D.domain.cellsize.y.', outer=(Nx, 1, Nz))
# DZ = zeros(1,1,Nz+2);
# DZ(1,1,:) = D.domain.cellsize.z;
# DZ=repeat(DZ, Nx, Ny, 1);
# dx = 0.5*(DX(1:end-1,:,:)+DX(2:end,:,:));
# dy = 0.5*(DY(:,1:end-1,:)+DY(:,2:end,:));
# dz = 0.5*(DZ(:,:,1:end-1)+DZ(:,:,2:end));

DX = D.domain.cellsize.x
DY = zeros( 1, Ny+2)
DY[:] .= D.domain.cellsize.y
DZ = zeros( 1,1,Nz+2)
DZ[:] .= D.domain.cellsize.z
dx = 0.5*(DX[1:end-1].+DX[2:end])
dy = 0.5*(DY[:,1:end-1].+DY[:,2:end])
dz = 0.5*(DZ[:,:,1:end-1].+DZ[:,:,2:end])

# define the vectors to store the sparse matrix data
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
# code readability (use broadcasting)
De = D.xvalue[2:Nx+1,:,:]./(dx[2:Nx+1].*DX[2:Nx+1])
Dw = D.xvalue[1:Nx,:,:]./(dx[1:Nx].*DX[2:Nx+1])
Dn = D.yvalue[:,2:Ny+1,:]./(dy[:,2:Ny+1].*DY[:,2:Ny+1])
Ds = D.yvalue[:,1:Ny,:]./(dy[:,1:Ny].*DY[:,2:Ny+1])
Df = D.zvalue[:,:,2:Nz+1]./(dz[:,:,2:Nz+1].*DZ[:,:,2:Nz+1])
Db = D.zvalue[:,:,1:Nz]./(dz[:,:,1:Nz].*DZ[:,:,2:Nz+1])

# calculate the coefficients for the internal cells
AE = reshape(De,mnx)
AW = reshape(Dw,mnx)
AN = reshape(Dn,mny)
AS = reshape(Ds,mny)
AF = reshape(Df,mnz)
AB = reshape(Db,mnz)
APx = reshape(.-(De+Dw),mnx)
APy = reshape(.-(Dn+Ds),mny)
APz = reshape(.-(Df+Db),mnz)

# build the sparse matrix based on the numbering system
rowx_index = reshape(G[2:Nx+1,2:Ny+1,2:Nz+1],mnx)  # main diagonal x
iix[1:3*mnx] .= repeat(rowx_index,3)
rowy_index = reshape(G[2:Nx+1,2:Ny+1,2:Nz+1],mny)  # main diagonal y
iiy[1:3*mny] .= repeat(rowy_index,3);
rowz_index = reshape(G[2:Nx+1,2:Ny+1,2:Nz+1],mnz)  # main diagonal z
iiz[1:3*mnz] .= repeat(rowz_index,3)
jjx[1:3*mnx] .= [reshape(G[1:Nx,2:Ny+1,2:Nz+1],mnx);
		reshape(G[2:Nx+1,2:Ny+1,2:Nz+1],mnx);
		reshape(G[3:Nx+2,2:Ny+1,2:Nz+1],mnx)]
jjy[1:3*mny] .= [reshape(G[2:Nx+1,1:Ny,2:Nz+1],mny);
		reshape(G[2:Nx+1,2:Ny+1,2:Nz+1],mny);
		reshape(G[2:Nx+1,3:Ny+2,2:Nz+1],mny)]
jjz[1:3*mnz] .= [reshape(G[2:Nx+1,2:Ny+1,1:Nz],mnz);
		reshape(G[2:Nx+1,2:Ny+1,2:Nz+1],mnz);
		reshape(G[2:Nx+1,2:Ny+1,3:Nz+2],mnz)];
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



# ======================== 1D Cylindrical DIFFUSION =========================
function diffusionTermCylindrical1D(D::FaceValue)
# D is a face variable

# extract data from the mesh structure
Nx = D.domain.dims[1]
G = [1:Nx+2;]
DX = D.domain.cellsize.x
dx = 0.5.*(DX[1:end-1].+DX[2:end])
rp = D.domain.cellcenters.x
rf = D.domain.facecenters.x

# define the vectors to store the sparse matrix data
iix = zeros(Int64, 3*(Nx+2))
jjx = zeros(Int64, 3*(Nx+2))
sx = zeros(Float64, 3*(Nx+2))

# reassign the east, west for code readability
De = rf[2:Nx+1].*D.xvalue[2:Nx+1]./(rp.*dx[2:Nx+1].*DX[2:Nx+1])
Dw = rf[1:Nx].*D.xvalue[1:Nx]./(rp.*dx[1:Nx].*DX[2:Nx+1])

# calculate the coefficients for the internal cells
AE = reshape(De,Nx)
AW = reshape(Dw,Nx)
APx = reshape(-(De+Dw),Nx)

# build the sparse matrix based on the numbering system
rowx_index = reshape(G[2:Nx+1],Nx) # main diagonal x
iix[1:3*Nx] .= repeat(rowx_index,3)
jjx[1:3*Nx] .= [reshape(G[1:Nx],Nx); reshape(G[2:Nx+1],Nx); reshape(G[3:Nx+2],Nx)]
sx[1:3*Nx] .= [AW; APx; AE]

# build the sparse matrix
kx = 3*Nx
M = sparse(iix[1:kx], jjx[1:kx], sx[1:kx], Nx+2, Nx+2)
end



# ======================== 2D RADIAL DIFFUSION =========================
function diffusionTermRadial2D(D::FaceValue)
# D is a face variable
# extract data from the mesh structure

Nr = D.domain.dims[1]
Ntheta = D.domain.dims[2]
G=reshape([1:(Nr+2)*(Ntheta+2);], Nr+2, Ntheta+2)
DR = D.domain.cellsize.x
DTHETA = zeros( 1, Ntheta+2)
DTHETA[:] .= D.domain.cellsize.y
dr = 0.5.*(DR[1:end-1].+DR[2:end])
dtheta = zeros( 1, Ntheta+1)
dtheta[:] .= 0.5*(DTHETA[1:end-1].+DTHETA[2:end])
rp = D.domain.cellcenters.x
rf = D.domain.facecenters.x

# define the vectors to store the sparse matrix data
iix = zeros(Int64, 3*(Nr+2)*(Ntheta+2))
iiy = zeros(Int64, 3*(Nr+2)*(Ntheta+2))
jjx = zeros(Int64, 3*(Nr+2)*(Ntheta+2))
jjy = zeros(Int64, 3*(Nr+2)*(Ntheta+2))
sx = zeros(Float64, 3*(Nr+2)*(Ntheta+2))
sy = zeros(Float64, 3*(Nr+2)*(Ntheta+2))
mnx = Nr*Ntheta
mny = Nr*Ntheta

# reassign the east, west for code readability
De = rf[2:Nr+1].*D.xvalue[2:Nr+1,:]./(rp.*dr[2:Nr+1].*DR[2:Nr+1])
Dw = rf[1:Nr].*D.xvalue[1:Nr,:]./(rp.*dr[1:Nr].*DR[2:Nr+1])
Dn = D.yvalue[:,2:Ntheta+1]./(rp.*rp.*dtheta[:,2:Ntheta+1].*DTHETA[:,2:Ntheta+1])
Ds = D.yvalue[:,1:Ntheta]./(rp.*rp.*dtheta[:,1:Ntheta].*DTHETA[:,2:Ntheta+1])

# calculate the coefficients for the internal cells
AE = reshape(De,mnx)
AW = reshape(Dw,mnx)
AN = reshape(Dn,mny)
AS = reshape(Ds,mny)
APx = reshape(-(De+Dw),mnx)
APy = reshape(-(Dn+Ds),mny)

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



# ======================== 2D CYLINDRICAL DIFFUSION =========================
function diffusionTermCylindrical2D(D::FaceValue)
# D is a face variable
# extract data from the mesh structure
Nr = D.domain.dims[1]
Nz = D.domain.dims[2]
G=reshape([1:(Nr+2)*(Nz+2);], Nr+2, Nz+2)
DR = D.domain.cellsize.x
DZ = zeros( 1, Nz+2)
DZ[:] .= D.domain.cellsize.y
dr = 0.5*(DR[1:end-1].+DR[2:end])
dz = zeros( 1, Nz+1)
dz[:] .= 0.5*(DZ[1:end-1].+DZ[2:end])
rp = repeat(D.domain.cellcenters.x, 1, Nz)
rf = repeat(D.domain.facecenters.x, 1, Nz)

# define the vectors to store the sparse matrix data
iix = zeros(Int64, 3*(Nr+2)*(Nz+2))
iiy = zeros(Int64, 3*(Nr+2)*(Nz+2))
jjx = zeros(Int64, 3*(Nr+2)*(Nz+2))
jjy = zeros(Int64, 3*(Nr+2)*(Nz+2))
sx = zeros(Float64, 3*(Nr+2)*(Nz+2))
sy = zeros(Float64, 3*(Nr+2)*(Nz+2))
mnx = Nr*Nz
mny = Nr*Nz

# extract the diffusion coefficient data
Dx = D.xvalue
Dy = D.yvalue

# reassign the east, west for code readability
De = rf[2:Nr+1,:].*D.xvalue[2:Nr+1,:]./(rp.*dr[2:Nr+1].*DR[2:Nr+1])
Dw = rf[1:Nr,:].*D.xvalue[1:Nr,:]./(rp.*dr[1:Nr].*DR[2:Nr+1])
Dn = D.yvalue[:,2:Nz+1]./(dz[:,2:Nz+1].*DZ[:,2:Nz+1])
Ds = D.yvalue[:,1:Nz]./(dz[:,1:Nz].*DZ[:,2:Nz+1])

# calculate the coefficients for the internal cells
AE = reshape(De,mnx)
AW = reshape(Dw,mnx)
AN = reshape(Dn,mny)
AS = reshape(Ds,mny)
APx = reshape(.-(De+Dw),mnx)
APy = reshape(.-(Dn+Ds),mny)

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



# ======================== 3D CYLINDRICAL DIFFUSION =========================
function diffusionTermCylindrical3D(D::FaceValue)
# extract data from the mesh structure

Nr = D.domain.dims[1]
Ntheta = D.domain.dims[2]
Nz = D.domain.dims[3]
G=reshape([1:(Nr+2)*(Ntheta+2)*(Nz+2);], Nr+2, Ntheta+2, Nz+2)
DR = D.domain.cellsize.x
DTHETA = zeros( 1, Ntheta+2)
DTHETA[:] .= D.domain.cellsize.y
DZ = zeros( 1,1,Nz+2)
DZ[:] .= D.domain.cellsize.z
dr = 0.5.*(DR[1:end-1].+DR[2:end])
dtheta = 0.5.*(DTHETA[:,1:end-1].+DTHETA[:,2:end])
dz = 0.5.*(DZ[:,:,1:end-1].+DZ[:,:,2:end])
rp = D.domain.cellcenters.x # use broadcasting
rf = D.domain.facecenters.x

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

# extract the velocity data
# note: size(ux) = [1:m+1, 1:n] and size(uy) = [1:m, 1:n+1]
Dx = D.xvalue
Dy = D.yvalue
Dz = D.zvalue

# reassign the east, west, north, and south velocity vectors for the
# code readability
De = rf[2:Nr+1].*D.xvalue[2:Nr+1,:,:]./(rp.*dr[2:Nr+1].*DR[2:Nr+1])
Dw = rf[1:Nr].*D.xvalue[1:Nr,:,:]./(rp.*dr[1:Nr].*DR[2:Nr+1])
Dn = D.yvalue[:,2:Ntheta+1,:]./(rp.*rp.*dtheta[:,2:Ntheta+1].*DTHETA[:,2:Ntheta+1])
Ds = D.yvalue[:,1:Ntheta,:]./(rp.*rp.*dtheta[:,1:Ntheta].*DTHETA[:,2:Ntheta+1])
Df = D.zvalue[:,:,2:Nz+1]./(dz[:,:,2:Nz+1].*DZ[:,:,2:Nz+1])
Db = D.zvalue[:,:,1:Nz]./(dz[:,:,1:Nz].*DZ[:,:,2:Nz+1])

# calculate the coefficients for the internal cells
AE = reshape(De,mnx)
AW = reshape(Dw,mnx)
AN = reshape(Dn,mny)
AS = reshape(Ds,mny)
AF = reshape(Df,mnz)
AB = reshape(Db,mnz)
APx = reshape(.-(De+Dw),mnx)
APy = reshape(.-(Dn+Ds),mny)
APz = reshape(.-(Df+Db),mnz)

# build the sparse matrix based on the numbering system
rowx_index = reshape(G[2:Nr+1,2:Ntheta+1,2:Nz+1],mnx)  # main diagonal x
iix[1:3*mnx] .= repeat(rowx_index,3)
rowy_index = reshape(G[2:Nr+1,2:Ntheta+1,2:Nz+1],mny)  # main diagonal y
iiy[1:3*mny] .= repeat(rowy_index,3)
rowz_index = reshape(G[2:Nr+1,2:Ntheta+1,2:Nz+1],mnz)  # main diagonal z
iiz[1:3*mnz] .= repeat(rowz_index,3)
jjx[1:3*mnx] .= [reshape(G[1:Nr,2:Ntheta+1,2:Nz+1],mnx);
		reshape(G[2:Nr+1,2:Ntheta+1,2:Nz+1],mnx);
		reshape(G[3:Nr+2,2:Ntheta+1,2:Nz+1],mnx)]
jjy[1:3*mny] .= [reshape(G[2:Nr+1,1:Ntheta,2:Nz+1],mny);
		reshape(G[2:Nr+1,2:Ntheta+1,2:Nz+1],mny);
		reshape(G[2:Nr+1,3:Ntheta+2,2:Nz+1],mny)]
jjz[1:3*mnz] .= [reshape(G[2:Nr+1,2:Ntheta+1,1:Nz],mnz);
		reshape(G[2:Nr+1,2:Ntheta+1,2:Nz+1],mnz);
		reshape(G[2:Nr+1,2:Ntheta+1,3:Nz+2],mnz)]
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
