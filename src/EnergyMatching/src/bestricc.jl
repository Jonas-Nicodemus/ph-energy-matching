"""
    Σphr = bestricc(Σ::QuadraticOutputStateSpace, Σr::StateSpace; kwargs...)
    
Solves the energy matching problem, by returning the best solution of the algebraic Riccati equation.
"""
function bestricc(Σ::QuadraticOutputStateSpace, Σr::StateSpace; kwargs...)
    Xmin = kyp_min(Σr)
    Xmax = kyp_max(Σr)
    
    Ar = Σr.A
    Br = Σr.B
    Cr = Σr.C
    Dr = Σr.D

    Amin = Ar + Br * inv(-Dr - Dr') * (-Br' * Xmin + Cr)
    @assert all(real(eigvals(Amin)) .<= 0)
    
    BVs = projectors(Amin)
    Δ = Xmin - Xmax

    Σrmin = qoss(Ar, Br, 1/2*Xmin)
    Σrmax = qoss(Ar, Br, 1/2*Xmax)
    
    Vmin = norm(Σ - Σrmin)
    Vmax = norm(Σ - Σrmax)
    
    Vstar = Inf
    Xstar = zero(Xmin)
    if Vmin < Vmax
        Xstar .= Xmin
        Vstar = Vmin
    else
        Xstar .= Xmax
        Vstar = Vmax
    end

    for Z in combinations(BVs)
        P = Z * inv(Z' * Δ * Z) * Z' * Δ
        X = sym(Xmin * P + Xmax * (I - P))
        Σqor = qoss(Σr.A, Σr.B, 0.5*X)
        Vn =  norm(Σ - Σqor)
        if Vn < Vstar
            Xstar .= X
            Vstar = Vn
        end
    end
    return phss(Σr, Xstar)
end

function projectors(Ap)
    D, Z = eigen(Ap)
    BVs = Vector{Matrix{Float64}}(undef, 0)
    i = 1
    while i < length(D)
        if isreal(D[i])
            push!(BVs, real(Z[:, i][:, :]))
            i += 1
        else
            push!(
                BVs,
                hcat(
                    real(Z[:, i]) / norm(real(Z[:, i])),
                    imag(Z[:, i]) / norm(imag(Z[:, i])),
                ),
            )
            i += 2
        end
    end
    return BVs
end

function combinations(BVs)
    n = length(BVs)
    AC = Vector{Matrix{Float64}}(undef, 0)
    for i = 1:n
        for j = i:n
            Z, _, _ = svd(hcat(BVs[i:j]...))
            push!(AC, Z)
        end
    end
    return AC
end