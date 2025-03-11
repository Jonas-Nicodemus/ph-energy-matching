module EnergyMatching

using LinearAlgebra, ControlSystemsBase, MatrixEquations
using JuMP, Hypatia, COSMO, Optim, LineSearches
# using MosekTools
using PortHamiltonianSystems, QuadraticOutputSystems

export minreal, matchnrg, hdss

include("sdp.jl")
include("bestricc.jl")
include("barrier.jl")
include("minreal.jl")

"""
    Σrem = matchnrg(Σ::PortHamiltonianStateSpace, Σr::PortHamiltonianStateSpace; solver=:BFGS, kwargs...)

Solves the energy matching problem.
"""
function matchnrg(Σ::PortHamiltonianStateSpace, Σr::PortHamiltonianStateSpace; solver=:BFGS, kwargs...)
    Σqo = hdss(Σ)

    if solver == :Hypatia
        Σrem = sdp(Σqo, ss(Σr), Hypatia.Optimizer; kwargs...)
    elseif solver == :COSMO
        Σrem = sdp(Σqo, ss(Σr), COSMO.Optimizer; max_iter=100000, eps_abs=1e-7, eps_rel=1e-7, kwargs...)
    # elseif solver == :Mosek
    #     Σrem = sdp(Σqo, ss(Σr), Mosek.Optimizer; INTPNT_CO_TOL_DFEAS=1e-7, kwargs...)
    elseif solver == :BFGS
        Σrem = barrier(Σqo, ss(Σr); Σr0 = Σr, kwargs...)
    elseif solver == :BR_BFGS
        Σrem = barrier(Σqo, ss(Σr); kwargs...)
    elseif solver == :BestRicc
        Σrem = bestricc(Σqo, ss(Σr); kwargs...)
    else
        error("Unknown solver $solver")
    end

    return Σrem
end

"""
    Σqo = hdss(J, R, Q, G, P, S, N)
    Σqo = hdss(Σph)

Converts a `PortHamiltonianStateSpace` to a `QuadraticOutputSystem` (Hamiltonian dynamic).
"""
function hdss(J, R, Q, G, P, S, N)
    A, B, _, _ = compose(J, R, Q, G, P, S, N)
    return qoss(A, B, 1/2*vec(Q)')
end

function hdss(Σph::PortHamiltonianStateSpace)
    return hdss(Σph.J, Σph.R, Σph.Q, Σph.G, Σph.P, Σph.S, Σph.N)
end

end
