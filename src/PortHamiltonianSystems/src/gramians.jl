import ControlSystemsBase: StateSpace, grampd, gram

"""
    L = grampd(Σph, opt; kwargs...)

Returns a Cholesky factor `L` of the Gramian of the system `Σph`. If `opt` is `:c` or `:o` 
it returns the controllability or observability Gramian, respectively by calling the 
`grampd` from [ControlSystems.jl](https://github.com/JuliaControl/ControlSystems.jl) package. 
If `opt` is `:pr_c` or `:pr_o` it returns a Cholesky factor `L` of the 
positive-real controllability or observability Gramian, respectively.
"""
function grampd(Σph::PortHamiltonianStateSpace, opt::Symbol; kwargs...)
    Σ = ss(Σph)
    if opt == :c || opt == :o
        return ControlSystemsBase.grampd(Σ, opt; kwargs...)
    elseif opt == :pr_c
        return prgrampd(Σ, :c; kwargs...)
    elseif opt == :pr_o
        return prgrampd(Σ, :o; kwargs...)
    else
        error("opt must be either :c, :o, :pr_c or :pr_o")
    end
end

"""
    X = gram(Σph, opt; kwargs...)

Returns the Gramian of the system `Σph`. If `opt` is `:c` or `:o` 
it returns the controllability or observability Gramian, respectively by calling the 
`gram` from [ControlSystems](https://github.com/JuliaControl/ControlSystems.jl) package. 
If `opt` is `:pr_c` or `:pr_o` it returns the 
positive-real controllability or observability Gramian, respectively (see [`prgram`](@ref)).
"""
function gram(Σph::PortHamiltonianStateSpace, opt::Symbol; kwargs...)
    Σ = ss(Σph)
    if opt == :c || opt == :o
        return ControlSystemsBase.gram(Σ, opt; kwargs...)
    elseif opt == :pr_c
        return prgram(Σ, :c; kwargs...)
    elseif opt == :pr_o
        return prgram(Σ, :o; kwargs...)
    else
        error("opt must be either :c, :o, :pr_c or :pr_o")
    end
end

"""
    X = prgram(Σ, opt; kwargs...)

Returns the positive-real Gramian of system `Σ`. If `opt` is `:c` or `:o` 
it returns the positive-real controllability or positive-real observability Gramian, respectively,
by solving the corresponding positive-real algebraic Riccati equation (see [`prare`](@ref)).
"""
function prgram(Σ::StateSpace, opt::Symbol; min=true)
    if opt == :o
        A = Array(Σ.A)
        B = Array(Σ.B)
        R = Array(-Σ.D - Σ.D')
        Q = zero(Σ.A)
        S = Array(-Σ.C')
    elseif opt == :c
        A = Array(Σ.A')
        B = Array(Σ.C')
        R = Array(-Σ.D - Σ.D')
        Q = zero(Σ.A)
        S = Array(-Σ.B)
    else
        error("opt must be either ':c' or ':o'")
    end
    
    X, _, _ = arec(A, B, R, Q, S, as=!min)
    
    return X
end

function prgram(Σph::PortHamiltonianStateSpace, opt::Symbol; kwargs...)
    Σ = ss(Σph)
    return prgram(Σ, opt; kwargs...)
end

"""
    L = prgrampd(Σ, opt; kwargs...)

Returns the Cholesky factor of the positive-real Gramian of system `Σ` (see [`prgram`](@ref)). If `opt` is `:c` or `:o` 
it returns the positive-real controllability or positive-real observability Gramian, respectively.

In the case that the solution of the positive-real algebraic Riccati equation is not positive definite (due to numerical errors), 
it is projected to the set of positive semi-definite matrices calling [`project_psd`](@ref).
"""
function prgrampd(Σ::StateSpace, opt::Symbol; kwargs...)
    X = prgram(Σ, opt; kwargs...)
    
    return cholesky(X; check=false).U
    # try
    #     return cholesky(X).U
    # catch # in case X is not positive definite
    #     @warn "Cholesky faild, since X is not positive definite, projecting to the nearest positive definite matrix"
    #     X = project_psd(X; eigtol=1e-8)
    #     return cholesky(X).U
    # end
end

function prgrampd(Σph::PortHamiltonianStateSpace, opt::Symbol; kwargs...)
    Σ = ss(Σph)
    return prgrampd(Σ, opt; kwargs...)
end

"""
    prare(Σ, opt, X; kwargs...)

Evaluates the positive-real algebraic Riccati equation for system `Σ` and candidat solution `X`.

If `opt` is `:o` the positive-real controllability algebraic Riccati equation is evaluated,

    A'X + XA + (C' - XB) inv(D + D') (C - B'X).

If `opt` is `:c` the positive-real observability algebraic Riccati equation is evaluated,

    AX + XA' + (B - XC') inv(D + D') (B' - CX).

"""
function prare(Σ::StateSpace, opt::Symbol, X)
    if opt == :o
        return Σ.A' * X + X * Σ.A + (Σ.C' - X * Σ.B) * inv(Σ.D + Σ.D') * (Σ.C - Σ.B' * X)
    elseif opt == :c
        return Σ.A * X + X * Σ.A' + (Σ.B - X * Σ.C') * inv(Σ.D + Σ.D') * (Σ.B' - Σ.C * X)
    end
end

function prare(Σph::PortHamiltonianStateSpace, opt::Symbol, X)
    Σ = ss(Σph)
    return prare(Σ, opt, X)
end
