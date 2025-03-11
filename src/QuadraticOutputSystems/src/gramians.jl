import ControlSystems: grampd, gram 

"""
    L = grampd(Σqo::QuadraticOutputStateSpace, opt::Symbol; kwargs...)

Returns a Cholesky factor `L` of the Gramian of system `Σqo`. If `opt` is `:c` 
it returns the controllability Gramian by calling the 
`ControlSystems.grampd` from [ControlSystems](https://github.com/JuliaControl/ControlSystems.jl) package.
If `opt` is `:o` it returns the Cholesky factor of the observability Gramian by solving the Lyapunov equation
    
    A'X + XA + MPM = 0,

where `P` is the controllability Gramian of `Σqo`.
"""
function grampd(Σqo::QuadraticOutputStateSpace, opt::Symbol; kwargs...)
    Σ = ss(Σqo)
    n = Σ.nx
    U = grampd(Σ, :c; kwargs...)
    if opt === :c 
        return U
    elseif opt === :o
        # return plyapc(Array(Σqo.A'), Σqo.M * U; kwargs...)
        return cholesky(lyapc(Σqo.A', Σqo.C'*Σqo.C + sum([reshape(Mi,n,n) * U*U' * reshape(Mi,n,n) for Mi ∈ eachrow(Σqo.M)]); kwargs...); check=false).U
    else
        error("opt must be either :c, :o")
    end
end

"""
    X = gram(Σqo::QuadraticOutputStateSpace, opt::Symbol; kwargs...)

Returns the Gramian of system `Σqo`. If `opt` is `:c` 
it returns the controllability Gramian by calling the 
`ControlSystems.gram` from [ControlSystems](https://github.com/JuliaControl/ControlSystems.jl) package.
If `opt` is `:o` it returns the observability Gramian by solving the Lyapunov equation
    
    A'X + XA + MPM = 0,

where `P` is the controllability Gramian of `Σqo`.
"""
function gram(Σqo::QuadraticOutputStateSpace, opt::Symbol; kwargs...)
    U = grampd(Σqo, opt; kwargs...)
    X = opt === :c ? U*U' : U'U
    Λ, U = eigen(hermitianpart(X))
    return hermitianpart(U * diagm(clamp.(Λ, 0, Inf)) * U')
end
