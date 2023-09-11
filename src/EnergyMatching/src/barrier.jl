"""
    Σph = barrier(Σ, Σr; Σr0=nothing, kwargs...)

Solves the energy matching problem using the barrier method.
"""
function barrier(Σ::QuadraticOutputStateSpace, Σr::StateSpace; Σr0=nothing, kwargs...)
    if Σr0 === nothing
        Σr0 = EnergyMatching.bestricc(Σ, Σr)
    end
    
    x = vech(Σr0.Q)

    eps=1e-8
    W = kypmat(Σr, Σr0.Q)
    λ = minimum(eigvals(W))
    if λ < 0.0
        eps = -λ + 1e-8
    end
    @info "Using eps = $eps"
    
    fgo = objective(Σ, Σr)
    fgc = constraint(Σr; eps=eps)
        
    for α in exp10.(-3:-1:-15)
        res = Optim.optimize(Optim.only_fg!(combined(fgo, fgc, α)), x, Optim.BFGS(linesearch = LineSearches.BackTracking()),
            Optim.Options(
                g_tol = 1e-16,
                f_tol = 0.0,
                x_tol = 0.0,
                allow_f_increases = true,
                iterations = 50000,
                show_trace = true,
                show_every = 100,
            ))
        @info "α = $α, f = $(res.minimum),\n $(res))"
        x = res.minimizer
    end

    Qr = unvech(x)
    
    return phss(Σr, Qr)
end

"""
    fg = objective(Σ::QuadraticOutputStateSpace, Σr::StateSpace)

Returns the objective function and its gradient for the energy matching problem, in the standard format for [Optim.jl](https://github.com/JuliaNLSolvers/Optim.jl).
"""
function objective(Σ::QuadraticOutputStateSpace, Σr::StateSpace)
    # precomputations
    n = size(Σ.A, 1)
    r = size(Σr.A, 1)
    Pr = gram(Σr, :c)
    Y =  sylvc(Σ.A, Σr.A', -Σ.B * Σr.B')
    F = -(kron(I(r), Σ.A') + kron(Σr.A', I(n))) \ kron(I(r), 2 * Σ.M * Y)
    h2 = norm(Σ, 2)^2

    function fg(x)
        X = unvech(x)

        Z = reshape(F * reshape(X/4, :, 1), n, r)
        h2r = norm(qoss(Σr.A, Σr.B, 1/2 * X))^2
        
        # value
        f = h2 + h2r - 2 * tr(Σ.B' * Z * Σr.B)
        
        # gradient
        G = 1/2 * (Pr * X * Pr - 2 * Y' * Σ.M * Y)
        g = vech(2*G - diagm(diag(G)))

        return f, g
    end

    return fg
end

"""
    fg = constraint(Σ::QuadraticOutputStateSpace, Σr::StateSpace)

Returns the (barrier-)constraint function and its gradient for the energy matching problem, in the standard format for [Optim.jl](https://github.com/JuliaNLSolvers/Optim.jl).
"""
function constraint(Σr::StateSpace; eps=1e-8)
    
    function fg(x)
        X = unvech(x)
        W = sym(kypmat(Σr, X))
        W = W + eps * I(size(W, 1))
        
        if det(W) < 0.0 || any(eigvals(W) .< 0.0)
            f = Inf
            g = Inf * ones(length(x))
            return f, g
        else
            r, m = size(Σr.B)
            
            # value
            f = -log(det(W))

            # gradient
            G = -hcat(I(r), zeros(r, m)) * (W \ vcat(-Σr.A', -Σr.B')) -hcat(-Σr.A, -Σr.B) * (W \ vcat(I(r), zeros(m, r)))
            g = vech(2*G - diagm(diag(G)))
            
            return f, g
        end
    end

    return fg
end

"""
    fg = combined(fgo, fgc, α)

Returns the combined objective and constraint function and its gradient for given functions `fgo` and `fgc`, for a scalar weighting factor `α`.
The combined function is given by 

    f = fo + α * fc
    g = go + α * gc

where `fo` is the objective function `go` its gradient and `fc` is the constraint function with gradient `gc`.
For the energy matching problem `fgo` can be retrieved by [`objective`](@ref) and `fgc` by [`constraint`](@ref).
"""
function combined(fgo, fgc, α)
    function fg!(f, g, x)
        if g !== nothing
            f1, g1 = fgo(x)
            f2, g2 = fgc(x)
            g .= g1 + α * g2
            return f1+α*f2
        else
            g = zeros(length(x))
            f1, _ = fgo(x)
            f2, _ = fgc(x)
            return f1 + α*f2
        end
    end
end