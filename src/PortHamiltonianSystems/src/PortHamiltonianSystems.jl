module PortHamiltonianSystems

using LinearAlgebra, ControlSystems, MatrixEquations
# using MATLAB
using QuadraticOutputSystems

export PortHamiltonianStateSpace, phss, hdss # types/PortHamiltonianStateSpace.jl
export grampd, gram, prgrampd, prgram # gramians.jl
export kyp, kypare, kypmat, kyp_min, kyp_max # kyp.jl
export compose, decompose, dedescriptorize # convert.jl
export sym, skew, project_psd # utils.jl

include("types/PortHamiltonianStateSpace.jl")
include("gramians.jl")
include("analysis.jl")
include("kyp.jl")
include("convert.jl")
include("utils.jl")

end
