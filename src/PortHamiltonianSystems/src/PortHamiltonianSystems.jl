module PortHamiltonianSystems

using LinearAlgebra, ControlSystemsBase, MatrixEquations, VectorizationTransformations
using JuMP, Hypatia

export PortHamiltonianStateSpace, phss # types/PortHamiltonianStateSpace.jl
export grampd, gram, prgrampd, prgram # gramians.jl
export kyp, kypare, kypmat, kypmin, kypmax # kyp.jl
export compose, decompose # convert.jl
export sym, skew, project_psd, issym, isskew, ispsd # utils.jl
export vech, unvech, unvec, duplication, sparsify! # utils.jl

import ControlSystemsBase: to_matrix, AbstractNumOrArray

include("types/PortHamiltonianStateSpace.jl")
include("analysis.jl")
include("convert.jl")
include("gramians.jl")
include("kyp.jl")
include("utils.jl")

end
