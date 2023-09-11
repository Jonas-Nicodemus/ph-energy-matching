module EnergyMatching

using LinearAlgebra, SparseArrays, ControlSystems, MatrixEquations, JuMP, COSMO, Optim, LineSearches

# Require extra installation
# using Mosek, MosekTools, SeDuMi,

using PortHamiltonianSystems, QuadraticOutputSystems

export matchnrg, vech, unvech

include("sdp.jl")
include("bestricc.jl")
include("reshape.jl")
include("barrier.jl")

"""
    Σrem = matchnrg(Σ::PortHamiltonianStateSpace, Σr::PortHamiltonianStateSpace; solver=:BFGS, kwargs...)

Solves the energy matching problem for the solver `:BFGS` (default), `:COSMO` or `:BestRicc`.
* using `:BFGS` the problem is solved by the barrier method (see [`barrier`](@ref)).
* using `:COSMO` the problem is solved as equivalent SDP by the SDP solver [COSMO.jl](https://github.com/oxfordcontrol/COSMO.jl) (see [`sdp`](@ref)).
* using `:BestRicc` the problem is solved by returning the best solution to the algebraic Riccati equation (see [`bestricc`](@ref)).
"""
function matchnrg(Σ::PortHamiltonianStateSpace, Σr::PortHamiltonianStateSpace; solver=:BFGS, kwargs...)
    Σqo = hdss(Σ)
    Σr = ss(Σr)

    if solver == :COSMO
        Σrem = sdp(Σqo, Σr, COSMO.Optimizer; max_iter=100000, eps_abs=1e-7, eps_rel=1e-7, kwargs...)
    # elseif solver == :Mosek
    #     Σrem = sdp(Σqo, Σr, Mosek.Optimizer; INTPNT_CO_TOL_DFEAS=1e-7, kwargs...)
    # elseif solver == :SeDuMi
    #     Σrem = sdp(Σqo, Σr, SeDuMi.Optimizer; kwargs...)
    elseif solver == :BFGS
        Σrem = barrier(Σqo, Σr; kwargs...)
    elseif solver == :BestRicc
        Σrem = bestricc(Σqo, Σr; kwargs...)
    else
        error("Unknown solver")
    end

    return Σrem
end

end
