"""
    Σph = barrier(Σ, Σr; Σr0=nothing, kwargs...)

Solves the energy matching problem using the barrier method.
"""
function barrier(Σ::QuadraticOutputStateSpace, Σr::StateSpace; Σr0=nothing, kwargs...)
    if Σr0 === nothing
        Σr0 = EnergyMatching.bestricc(Σ, Σr)
    end
    
    x = 1//2 * vech(Σr0.Q)

    eps=1e-8
    W = kypmat(Σr, Σr0.Q)
    λ = minimum(eigvals(W))
    if λ < 0.0
        eps = -λ + 1e-8
    end
    @info "Using eps = $eps"
    
    fo, go = objective(Σ, Σr)
    fc, gc = constraint(Σr; eps=eps)
        
    for α in exp10.(-3:-1:-15)
        
        # barrier method objective
        f = (x) -> fo(x) + α * fc(x)
        g!(g, x) = begin
            g .= go(x) + α * gc(x)
        end
        
        res = Optim.optimize(f, g!, x, Optim.BFGS(linesearch = LineSearches.BackTracking()),
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

    Qr = 2 * unvech(x)
    
    return phss(Σr, Qr)
end

"""
    fg = objective(Σ::QuadraticOutputStateSpace, Σr::StateSpace)

Returns the objective function and its gradient for the energy matching problem, in the standard format for [Optim.jl](https://github.com/JuliaNLSolvers/Optim.jl).
"""
function objective(Σ::QuadraticOutputStateSpace, Σr::StateSpace)
    r = Σr.nx
    
    # precomputations
    P = gram(Σ, :c)
    h2 = tr(P*unvec(Σ.M)*P*unvec(Σ.M))
    Pr = gram(Σr, :c)
    Y =  sylvc(Σ.A, Σr.A', -Σ.B * Σr.B')
    Dr = duplication(r)

    function f(x)
        X = unvech(x)
        
        # value
        f = h2 + tr(Pr * X * Pr * X) - 2*tr(Y' * unvec(Σ.M) * Y * X)
        # f = tr(Pr * X * Pr * X) - 2*tr(Y' * unvec(Σ.M) * Y * X)

        return f
    end

    function g(x)
        X = unvech(x)
        
        # gradient
        G = 2 * (Pr * X * Pr - Y' * unvec(Σ.M) * Y)
        g = Dr' * vec(G)
    end

    # function h(x)
    #     return 2 * Dr' * kron(Pr, Pr) * Dr
    # end

    return f, g
end


"""
    fg = constraint(Σ::QuadraticOutputStateSpace, Σr::StateSpace)

Returns the (barrier-)constraint function and its gradient for the energy matching problem, in the standard format for [Optim.jl](https://github.com/JuliaNLSolvers/Optim.jl).
"""
function constraint(Σr::StateSpace; eps=1e-8)
    r, m = Σr.nx, Σr.nu 
    Dr = duplication(r)

    function fg(x)
        X = unvech(x)
        W = sym(kypmat(Σr, 2*X))
        W = W + eps * I(size(W, 1))
        
        if det(W) < 0.0 || any(eigvals(W) .< 0.0)
            f = Inf
            g = Inf * ones(length(x))
            # H = Inf * ones(length(x), length(x))
            return f, g
        else
            # value
            f = -log(det(W))

            # gradient
            A = [Σr.A Σr.B]
            B = [I(r); zeros(m, r)]
            G = 2 * A * (W \ B) + 2 * B' * (W \ A') 
            g = Dr' * vec(G)

            return f, g
        end
    end

    f(x) = fg(x)[1]
    g(x) = fg(x)[2]

    return f, g
end
