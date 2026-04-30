# ===============================
# Written by AAE
# TU Delft, Winter 2014
# simulkade.com
# ===============================


"""
fluxLimiter(flName::AbstractString)
This function returns a function handle to a flux limiter of user's
choice.
available flux limiters are: 'CHARM', 'HCUS', 'HQUICK', 'VanLeer',
'VanAlbada1', 'VanAlbada2', 'MinMod', 'SUPERBEE', 'Sweby', 'Osher',
'Koren', 'smart', 'MUSCL', 'QUICK', 'MC', and 'UMIST'.
Default limiter is 'SUPERBEE'. See:
<http://en.wikipedia.org/wiki/Flux_limiter>
"""
function fluxLimiter(flName::AbstractString)
# This function returns a function handle to a flux limiter of user's
# choice.
# available flux limiters are: 'CHARM', 'HCUS', 'HQUICK', 'VanLeer',
# 'VanAlbada1', 'VanAlbada2', 'MinMod', 'SUPERBEE', 'Sweby', 'Osher',
# 'Koren', 'smart', 'MUSCL', 'QUICK', 'MC', and 'UMIST'.
# Default limiter is 'SUPERBEE'. See:
# <http://en.wikipedia.org/wiki/Flux_limiter>

# find the flux limiter function

if flName=="CHARM"
  r->((r.>0.0).*r.*(3.0*r+1.0)./(((r+1.0).^2.0)+eps()*(r.==-1.0)))
elseif flName=="HCUS"
  r->(1.5*(r+abs.(r))./(r+2.0))
elseif flName=="HQUICK"
  r->(2.0*(r+abs.(r))./((r+3.0)+eps()*(r.==-3.0)))
elseif flName=="ospre"
  r->((1.5.*r.*(r+1.0))./(r.*(r+1.0)+1.0+eps()*((r.*(r+1.0)+1.0).==0.0)))
elseif flName=="VanLeer"
  r->((r.+abs.(r))./(1.0+abs.(r)))
elseif flName=="VanAlbada1"
  r->((r+r.*r)./(1.0+r.*r))
elseif flName=="VanAlbada2"
  r->(2.0*r./(1+r.*r))
elseif flName=="MinMod"
  r->((r>0.0).*min.(r,1.0))
elseif flName=="SUPERBEE"
  r->(max.(0.0, max.(min.(2.0*r,1.0), min.(r,2.0))))
elseif flName=="Osher"
  b=1.5
  r->(max.(0.0, min.(r,b)))
elseif flName=="Sweby"
  b=1.5
  r->(max.(0.0, max.(min.(b*r,1.0), min.(r,b))))
elseif flName=="smart"
  r->(max.(0.0, min.(4.0,min.(0.25+0.75*r, 2.0*r))))
elseif flName=="Koren"
  r->(max.(0.0, min.(2.0*r, min.((1.0+2.0*r)/3.0, 2.0))))
elseif flName=="MUSCL"
  r->(max.(0.0, min.(2.0*r, min.(0.5*(1+r), 2.0))))
elseif flName=="QUICK"
  r->(max.(0.0, min.(2.0, min.(2.0*r, (3.0+r)/4.0))))
elseif flName=="UMIST"
  r->(max.(0.0, min.(2.0, min.(2.0*r, min.((1.0+3.0*r)/4.0, (3.0+r)/4.0)))))
else
  println("The flux limiter of your choice is not available. The SUPERBEE flux limiter is used instead.")
  r->(max.(0.0, max.(min.(2.0*r,1.0), min.(r,2.0))))
end

end

"""
update!(phi_old::CellValue, phi_new::CellValue)
Copies the values of the new variable into the old cell variable:
phi_old.value .= phi_new.value
"""
function update!(phi_old::CellValue, phi_new::CellValue)
  phi_old.value .= phi_new.value
end

"""
update!(phi_old::FaceValue, phi_new::FaceValue)
Copies the values of the new variable into the old variable:
  phi_old.xvalue .= phi_new.xvalue
  phi_old.yvalue .= phi_new.yvalue
  phi_old.zvalue .= phi_new.zvalue
"""
function update!(phi_old::FaceValue, phi_new::FaceValue)
  phi_old.xvalue .= phi_new.xvalue
  phi_old.yvalue .= phi_new.yvalue
  phi_old.zvalue .= phi_new.zvalue
end

function faceEval(f::Function, x::FaceValue)
FaceValue(x.domain,
    f.(x.xvalue),
    f.(x.yvalue),
    f.(x.zvalue))
end

function faceEval(f::Function, x1::FaceValue, x2::FaceValue)
FaceValue(x1.domain,
    f.(x1.xvalue, x2.xvalue),
    f.(x1.yvalue, x2.yvalue),
    f.(x1.zvalue, x2.zvalue))
end

function faceEval(f::Function, x1::FaceValue, x2::FaceValue, x3::FaceValue)
FaceValue(x1.domain,
    f.(x1.xvalue, x2.xvalue, x3.xvalue),
    f.(x1.yvalue, x2.yvalue, x3.yvalue),
    f.(x1.zvalue, x2.zvalue, x3.zvalue))
end

function faceEval(f::Function, x1::FaceValue, x2::FaceValue, x3::FaceValue, x4::FaceValue)
FaceValue(x1.domain,
    f.(x1.xvalue, x2.xvalue, x3.xvalue, x4.xvalue),
    f.(x1.yvalue, x2.yvalue, x3.yvalue, x4.yvalue),
    f.(x1.zvalue, x2.zvalue, x3.zvalue, x4.zvalue))
end

function faceEval(f::Function, x1::FaceValue, x2::FaceValue, x3::FaceValue, x4::FaceValue, x5::FaceValue)
FaceValue(x1.domain,
    f.(x1.xvalue, x2.xvalue, x3.xvalue, x4.xvalue, x5.xvalue),
    f.(x1.yvalue, x2.yvalue, x3.yvalue, x4.yvalue, x5.yvalue),
    f.(x1.zvalue, x2.zvalue, x3.zvalue, x4.zvalue, x5.zvalue))
end

function faceEval(f::Function, x1::FaceValue, x2::FaceValue, x3::FaceValue, x4::FaceValue, x5::FaceValue, x6::FaceValue)
FaceValue(x1.domain,
    f.(x1.xvalue, x2.xvalue, x3.xvalue, x4.xvalue, x5.xvalue, x6.xvalue),
    f.(x1.yvalue, x2.yvalue, x3.yvalue, x4.yvalue, x5.yvalue, x6.yvalue),
    f.(x1.zvalue, x2.zvalue, x3.zvalue, x4.zvalue, x5.zvalue, x6.zvalue))
end

function cellEval(f::Function, x::CellValue)
CellValue(x.domain,
    f.(x.value))
end

function cellEval(f::Function, x1::CellValue, x2::CellValue)
CellValue(x1.domain,
    f.(x1.value, x2.value))
end

function cellEval(f::Function, x1::CellValue, x2::CellValue, x3::CellValue)
CellValue(x1.domain,
    f.(x1.value, x2.value, x3.value))
end

function cellEval(f::Function, x1::CellValue, x2::CellValue, x3::CellValue, x4::CellValue)
CellValue(x1.domain,
    f.(x1.value, x2.value, x3.value, x4.value))
end

function cellEval(f::Function, x1::CellValue, x2::CellValue, x3::CellValue, x4::CellValue, x5::CellValue)
CellValue(x1.domain,
    f.(x1.value, x2.value, x3.value, x4.value, x5.value))
end

function cellEval(f::Function, x1::CellValue, x2::CellValue, x3::CellValue, x4::CellValue, x5::CellValue, x6::CellValue)
CellValue(x1.domain,
    f.(x1.value, x2.value, x3.value, x4.value, x5.value, x6.value))
end

# ========================= Generate random perm field ======================
function permfieldlogrndg(Nx,k_avrg,V_dp,cl)
  # 1D random field generator:
  # hdf: Gaussian
  # acf: Gaussian
  Lx=1.0
  x = range(-Lx/2,Lx/2,length=Nx)
  s = -log(1-V_dp)
  mu = log(k_avrg)-s^2/2.0
  Z = s*randn(Nx)

  # Gaussian filter
  F = exp.(-x.^2/(cl^2/2.0))

  # correlation of surface using convolution (faltung), inverse
  # Fourier transform and normalizing prefactors
  f = sqrt(2.0/sqrt(pi))*sqrt(Lx/Nx/cl)*ifft(fft(Z).*fft(F))
  perm = exp.(mu .+ real(f))
end

function permfieldlogrnde(Nx,k_avrg,V_dp,cl)
  # 1D random field generator:
  # hdf: Gaussian
  # acf: Exponential
  Lx=1.0
  x = range(-Lx/2,Lx/2,length=Nx)
  s = -log(1-V_dp)
  mu = log(k_avrg)-s^2/2.0
  Z = s*randn(Nx)

  # Gaussian filter
  F = exp.(-abs.(x)/(cl/2.0))

  # correlation of surface using convolution (faltung), inverse
  # Fourier transform and normalizing prefactors
  f = sqrt(2.0)*sqrt(Lx/Nx/cl)*ifft(fft(Z).*fft(F))
  perm = exp.(mu .+ real.(f))
end


# ======================== 2D ==================================>
function permfieldlogrndg(Nx,Ny,k_avrg,V_dp,clx,cly)
# 2D random field generator:
# hdf: Gaussian
# acf: Gaussian
# The surface has a Gaussian height distribution function
# and Gaussian autocovariance functions
  Lx=1.0
  X = range(-Lx/2,Lx/2,length=Nx)
  Y = range(-Lx/2,Lx/2,length=Ny)'

  s = -log(1-V_dp)
  mu = log(k_avrg)-s^2/2
  Z = s*randn(Nx,Ny)
  # Gaussian filter
  F = exp.(-(X.^2/(clx^2/2.0).+Y.^2/(cly^2/2.0)))
  # correlated surface generation including convolution (faltning) and inverse
  # Fourier transform and normalizing prefactors
  f = 2.0/sqrt(pi)*Lx/sqrt(Nx*Ny)/sqrt(clx)/sqrt(cly)*ifft(fft(Z).*fft(F))
  perm = exp.(mu .+ real.(f))
end

function permfieldlogrnde(Nx,Ny,k_avrg,V_dp,clx,cly)
# 2D random field generator:
# hdf: Gaussian
# acf: Exponential
# The surface has a Gaussian height distribution function
# and Exponential autocovariance functions
  Lx=1.0
  X = range(-Lx/2,Lx/2,length=Nx)
  Y = range(-Lx/2,Lx/2,length=Ny)'

  s = -log(1-V_dp)
  mu = log(k_avrg)-s^2/2
  Z = s*randn(Nx,Ny)
  # Gaussian filter
  F = exp.(-(abs.(X)/(clx/2).+abs.(Y)/(cly/2)))
  # correlated surface generation including convolution (faltning) and inverse
  # Fourier transform and normalizing prefactors
  f = 2.0*Lx/sqrt(Nx*Ny)/sqrt(clx*cly)*ifft(fft(Z).*fft(F))
  perm = exp.(mu .+ real.(f))
end
# <========================= 2D ================================


# ======================== 3D ==================================>
function permfieldlogrndg(Nx,Ny,Nz,k_avrg,V_dp,clx,cly,clz)
# 2D random field generator:
# hdf: Gaussian
# acf: Gaussian
# The surface has a Gaussian height distribution function
# and Gaussian autocovariance functions
  Lx=1.0
  X = range(-Lx/2,Lx/2,length=Nx)
  Y = range(-Lx/2,Lx/2,length=Ny)'
  Z = zeros(1,1,Nz)
  Z[:]=range(-Lx/2,Lx/2,length=Nz)
  s = -log(1-V_dp)
  mu = log(k_avrg)-s^2/2
  z = s*randn(Nx,Ny,Nz)
  # Gaussian filter
  F = exp.(-(X.^2/(clx^2/2.0).+Y.^2/(cly^2/2.0).+Z.^2/(clz^2/2.0)))
  # correlated surface generation including convolution (faltning) and inverse
  # Fourier transform and normalizing prefactors
  f = 2.0/sqrt(pi)*Lx/(Nx*Ny*Nz)^(1/3)/(clx*cly*clz)^(1/3)*ifft(fft(z).*fft(F))
  perm = exp.(mu .+ real(f))
end

function permfieldlogrnde(Nx,Ny,Nz,k_avrg,V_dp,clx,cly,clz)
# 2D random field generator:
# hdf: Gaussian
# acf: Exponential
# The surface has a Gaussian height distribution function
# and Exponential autocovariance functions
  Lx=1.0
  X = range(-Lx/2,Lx/2,length=Nx)
  Y = range(-Lx/2,Lx/2,length=Ny)'
  Z = zeros(1,1,Nz)
  Z[:]=range(-Lx/2,Lx/2,length=Nz)
  s = -log(1-V_dp)
  mu = log(k_avrg)-s^2/2
  z = s*randn(Nx,Ny,Nz)
  # Gaussian filter
  F = exp.(-(abs.(X)/(clx/2.0).+abs.(Y)/(cly/2.0).+abs.(Z)/(clz/2.0)))
  # correlated surface generation including convolution (faltning) and inverse
  # Fourier transform and normalizing prefactors
  f = 2.0*Lx/(Nx*Ny*Nz)^(1/3)/(clx*cly*clz)^(1/3)*ifft(fft(z).*fft(F))
  perm = exp.(mu .+ real(f))
end
# <============================== 3D ==============================

"""
function cellvol = cellVolume(m::MeshStructure)
returns the volume of each cell in the form of a cell variable
"""
function cellVolume(m::MeshStructure)
  cs = m.coordinatesystem
  BC = createBC(m)
  if cs isa Cartesian1D
    c=m.cellsize.x[2:end-1]
  elseif cs isa Cylindrical1D
    c=2.0*pi()*m.cellsize.x[2:end-1].*m.cellcenters.x
  elseif cs isa Cartesian2D
    c=m.cellsize.x[2:end-1]*m.cellsize.y[2:end-1]'
  elseif cs isa Cylindrical2D # cylindrical
    c=2.0*pi*m.cellcenters.x.*m.cellsize.x[2:end-1]*m.cellsize.y[2:end-1]'
  elseif cs isa Radial2D # radial
    c=m.cellcenters.x.*m.cellsize.x[2:end-1]*m.cellsize.y[2:end-1]'
  elseif cs isa Cartesian3D
    z = zeros(1,1,m.dims[3])
    z[1,1,:] .= m.cellsize.z[2:end-1]
    c=m.cellsize.x[2:end-1]*m.cellsize.y[2:end-1]'.*z
  elseif cs isa Cylindrical3D
    z = zeros(1,1,m.dims[3])
    z[1,1,:] .= m.cellsize.z[2:end-1]
    c=m.cellcenters.x.*m.cellsize.x[2:end-1]*m.cellsize.y[2:end-1]'.*z
  end
  return createCellVariable(m, c, BC)
end

"""
this function reshapes a vetorized cell variable to its domain shape
matrix based on the mesh structure data; it is assumed that the phi
includes the ghost cell data as well.
"""
function reshapeCell(m, phi)
  reshape(Matrix(phi), tuple(m.dims+2...))
end

"""
this function reshapes a vetorized cell variable to its domain shape
matrix based on the mesh structure data; it is assumed that the phi
does NOT include the ghost cell data.
"""
function reshapeInternalCell(m, phi)
  reshape(Matrix(phi), tuple(m.dims...))
end

"""
returns the internal cells of a cell variable as an array of the same shape
"""
function internalCells(phi::CellValue)
  cs = phi.domain.coordinatesystem
  N = phi.domain.dims

  if is_1d(cs)
  	cellvar= phi.value[2:N[1]+1]
  elseif is_2d(cs)
  	cellvar= phi.value[2:N[1]+1, 2:N[2]+1]
  elseif is_3d(cs)
      cellvar= phi.value[2:N[1]+1, 2:N[2]+1, 2:N[3]+1]
  end
  return cellvar
end

"""
Integrate variable phi over the domain it is defined
"""
function domainInt(phi::CellValue)
  return sum(internalCells(phi).*internalCells(cellVolume(phi.domain)))
end
