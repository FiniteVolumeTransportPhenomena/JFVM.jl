__precompile__()

module JFVM

# global mumps_solver
# using PyPlot
try
  import MUMPS
  global mumps_solver = MUMPS
catch
  @info "MUMPS solver (optional) is not available."
end

using SparseArrays, FFTW
# using PyCall
# I prefer not to use the following command for the issues that it has on windows machines
# pygui_start(:wx)
# mayavis=0
# try
#   @pyimport mayavi.mlab as m
#   mayavis=m
# catch
#   warn("Mayavi is not installed or could not be imported.")
# end

import Base: +, -, *, /, ^, ==, >, >=, <, <=, broadcast, sin, cos, tan, cot, abs, exp, log, log10
export MeshStructure, BoundaryCondition, CellValue, FaceValue, CellVector,
  AbstractCoordinateSystem, Cartesian1D, Cylindrical1D, Cartesian2D,
  Cylindrical2D, Radial2D, Cartesian3D, Cylindrical3D,
  coordinatesystem_from_dimension, dimension_from_coordinatesystem,
  is_1d, is_2d, is_3d,
       arithmeticMean, geometricMean, harmonicMean, upwindMean, linearMean,
       tvdMean, createBC, boundaryConditionTerm, cellBoundary!, solvePDE,
       divergenceTerm, gradientTerm, convectionUpwindTerm, createCellVector,
       convectionTerm, convectionTvdTerm, diffusionTerm, createCellVariable,
       createFaceVariable, copyCell, fluxLimiter, createMesh1D,
       createMesh2D, createMesh3D, createMeshRadial2D, createMeshCylindrical2D,
       createMeshCylindrical3D, createMeshCylindrical1D, solveLinearPDE,
       linearSourceTerm, constantSourceTerm, transientTerm,
       solveMUMPSLinearPDE, faceEval, cellEval, permfieldlogrndg, permfieldlogrnde,
       JFVM_test, solveExplicitPDE, reshapeCell,
       cellVolume, reshapeInternalCell, internalCells, domainInt, convectionTvdRHS,
       linearMean!, update!, solveLinearPDE!
       
      #  visualizeCells, visualizeCellVectors, plot, imshow, xlabel, ylabel, figure, legend, pcolor, contour, colorbar,

include("fvmToolTypes.jl")
include("meshstructure.jl")
include("boundarycondition.jl")
include("domainVariables.jl")
include("diffusionterms.jl")
include("transientTerms.jl")
include("domainOperators.jl")
include("convectionTerms.jl")
include("averagingTerms.jl")
include("calculusTerms.jl")
include("sourceTerms.jl")
include("solveVisualizePDE.jl")
include("JFVMtools.jl")
include("jfvm_test.jl")

end # module
