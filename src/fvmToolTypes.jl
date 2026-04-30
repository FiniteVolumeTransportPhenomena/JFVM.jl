struct CellLocation{T<:Real}
  x::Array{T,1}
  y::Array{T,1}
  z::Array{T,1}
end

struct CellSize{T<:Real}
  x::Array{T,1}
  y::Array{T,1}
  z::Array{T,1}
end

struct FaceLocation{T<:Real}
  x::Array{T,1}
  y::Array{T,1}
  z::Array{T,1}
end

abstract type AbstractCoordinateSystem end

struct Cartesian1D <: AbstractCoordinateSystem end
struct Cylindrical1D <: AbstractCoordinateSystem end
struct Cartesian2D <: AbstractCoordinateSystem end
struct Cylindrical2D <: AbstractCoordinateSystem end
struct Radial2D <: AbstractCoordinateSystem end
struct Cartesian3D <: AbstractCoordinateSystem end
struct Cylindrical3D <: AbstractCoordinateSystem end

coordinatesystem_from_dimension(d::Real) =
  if d == 1
    Cartesian1D()
  elseif d == 1.5
    Cylindrical1D()
  elseif d == 2
    Cartesian2D()
  elseif d == 2.5
    Cylindrical2D()
  elseif d == 2.8
    Radial2D()
  elseif d == 3
    Cartesian3D()
  elseif d == 3.2
    Cylindrical3D()
  else
    error("JFVM: Unsupported mesh dimension code $(d).")
  end

dimension_from_coordinatesystem(::Cartesian1D) = 1.0
dimension_from_coordinatesystem(::Cylindrical1D) = 1.5
dimension_from_coordinatesystem(::Cartesian2D) = 2.0
dimension_from_coordinatesystem(::Cylindrical2D) = 2.5
dimension_from_coordinatesystem(::Radial2D) = 2.8
dimension_from_coordinatesystem(::Cartesian3D) = 3.0
dimension_from_coordinatesystem(::Cylindrical3D) = 3.2

is_1d(::Cartesian1D) = true
is_1d(::Cylindrical1D) = true
is_1d(::AbstractCoordinateSystem) = false

is_2d(::Cartesian2D) = true
is_2d(::Cylindrical2D) = true
is_2d(::Radial2D) = true
is_2d(::AbstractCoordinateSystem) = false

is_3d(::Cartesian3D) = true
is_3d(::Cylindrical3D) = true
is_3d(::AbstractCoordinateSystem) = false

struct MeshStructure{CS<:AbstractCoordinateSystem}
  coordinatesystem::CS
  dimension::Float64
  dims::Array{Int,1}
  cellsize::CellSize
  cellcenters::CellLocation
  facecenters::FaceLocation
  corner::Array{Int,1}
  edge::Array{Int,1}
end

function MeshStructure(dimension::Real,
                       dims::Array{Int,1},
                       cellsize::CellSize,
                       cellcenters::CellLocation,
                       facecenters::FaceLocation,
                       corner::Array{Int,1},
                       edge::Array{Int,1})
  cs = coordinatesystem_from_dimension(dimension)
  MeshStructure(cs, Float64(dimension), dims, cellsize, cellcenters, facecenters, corner, edge)
end

is_1d(m::MeshStructure) = is_1d(m.coordinatesystem)
is_2d(m::MeshStructure) = is_2d(m.coordinatesystem)
is_3d(m::MeshStructure) = is_3d(m.coordinatesystem)

struct CellValue # {T<:Real}
  domain::MeshStructure
  value::Union{Array{<:Real}, DenseArray{Bool}} #Union{Array{T}, BitArray{}}
end

struct CellVector # {T<:Real}
  domain::MeshStructure
  xvalue::Union{Array{<:Real}, DenseArray{Bool}} # Array{T}
  yvalue::Union{Array{<:Real}, DenseArray{Bool}} # Array{T}
  zvalue::Union{Array{<:Real}, DenseArray{Bool}} # Array{T}
end

struct FaceValue # {T<:Real}
  domain::MeshStructure
  xvalue::Union{Array{<:Real}, DenseArray{Bool}} # Array{T}
  yvalue::Union{Array{<:Real}, DenseArray{Bool}} # Array{T}
  zvalue::Union{Array{<:Real}, DenseArray{Bool}} # Array{T}
end

mutable struct BorderValue{T<:Real}
  a::Array{T}
  b::Array{T}
  c::Array{T}
  periodic::Bool
end

struct BoundaryCondition
  domain::MeshStructure
  left::BorderValue
  right::BorderValue
  bottom::BorderValue
  top::BorderValue
  back::BorderValue
  front::BorderValue
end
