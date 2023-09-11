import ControlSystems: grampd, gram 

"""
    L = grampd(Σqo::QuadraticOutputStateSpace, opt::Symbol; kwargs...)

Returns a Cholesky factor `L` of the Gramian of system `Σqo`. If `opt` is `:c` 
it returns the controllability Gramian by calling the 
`ControlSystems.grampd` from [ControlSystems](https://github.com/JuliaControl/ControlSystems.jl) package.
If `opt` is `:o` it returns the Cholesky factor of the observability Gramian by solving the Lyapunov equation
    
    A'X + XA + MPM = 0,

where `P` is the controllability Gramian of `Σqo`.

**P. Benner, P. Goyal, and I. Pontes Duff.**
[Gramians, energy functionals, and balanced truncation for linear dynamical systems with quadratic outputs](https://doi.org/10.1109/TAC.2021.3086319).
IEEE Trans. Automat. Control, 67(2), 2022.
"""
function grampd(Σqo::QuadraticOutputStateSpace, opt::Symbol; kwargs...)
    Σ = ss(Σqo)
    U = grampd(Σ, :c; kwargs...)
    if opt === :c 
        return U
    elseif opt === :o
        return plyapc(Array(Σqo.A'), Σqo.M * U; kwargs...)
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

**P. Benner, P. Goyal, and I. Pontes Duff.**
[Gramians, energy functionals, and balanced truncation for linear dynamical systems with quadratic outputs](https://doi.org/10.1109/TAC.2021.3086319).
IEEE Trans. Automat. Control, 67(2), 2022.
"""
function gram(Σqo::QuadraticOutputStateSpace, opt::Symbol; kwargs...)
    if opt == :c || opt == :o
        U = grampd(Σqo, opt; kwargs...)
        return U * U'
    else
        error("opt must be either :c, :o")
    end
end