"""
    X = kyp(Σ; kwargs...)

Returns a solution to the KYP inequality by solving the corresponding linear matrix inequality.
"""
function kyp(Σ::StateSpace; kwargs...)
    model = Model(() -> Hypatia.Optimizer(verbose = false))
    for kwarg in keys(kwargs)
        set_optimizer_attribute(model, String(kwarg), kwargs[kwarg])
    end

    @variable(model, X[1:Σ.nx, 1:Σ.nx], PSD)
    @constraint(model, kypmat(Σ, X) in PSDCone())
    optimize!(model)

    return value.(X)
end

"""
    Xmin = kypmin(Σ; kwargs...)

Returns the minimal solution to the KYP inequality by solveing the Riccati equation for the stabilizing solution.
"""
function kypmin(Σ; kwargs...)
    return prgram(Σ, :o; min=true, kwargs...)
end

"""
    Xmin = kypmax(Σ; kwargs...)

Returns the maximal solution to the KYP inequality by solveing the Riccati equation for the anti-stabilizing solution.
"""
function kypmax(Σ; kwargs...)
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
