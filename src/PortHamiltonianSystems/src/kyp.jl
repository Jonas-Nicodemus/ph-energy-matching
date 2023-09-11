function kyp(Σ; kwargs...)
    return prgram(Σ, :o; kwargs...)
end

"""
    Xmin = kyp_min(Σ; kwargs...)

Returns the minimal solution to the KYP inequality by solveing the Riccati equation for the stabilizing solution.
"""
function kyp_min(Σ; kwargs...)
    return prgram(Σ, :o; min=true, kwargs...)
end

"""
    Xmin = kyp_max(Σ; kwargs...)

Returns the maximal solution to the KYP inequality by solveing the Riccati equation for the anti-stabilizing solution.
"""
function kyp_max(Σ; kwargs...)
    return prgram(Σ, :o; min=false, kwargs...)
end

function kypare(Σ::PortHamiltonianStateSpace, X)
    return kypare(ss(Σ), X)
 end 

function kypare(Σ::StateSpace, X)
   return PortHamiltonianSystems.prare(Σ, :o, X)
end

"""
    W = kypmat(Σ, X) 

Returns the KYP matrix of the system `Σ` for the given matrix `X`.
"""
function kypmat(Σ::StateSpace, X)
    return [-Σ.A'*X-X*Σ.A  Σ.C'-X*Σ.B; Σ.C-Σ.B'*X Σ.D + Σ.D']
end

function kypmat(Σ::PortHamiltonianStateSpace, X)
    return kypmat(ss(Σ), X)
end