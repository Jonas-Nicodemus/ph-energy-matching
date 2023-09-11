"""
    norm(Σph, p=2; kwargs...)

Converts a `PortHamiltonianStateSpace` to a `ControlSystems.StateSpace` and calls `norm` on it.
For more details see [ControlSystems.jl](https://github.com/JuliaControl/ControlSystems.jl) package.
"""
function LinearAlgebra.norm(sys::PortHamiltonianStateSpace, p::Real=2; kwargs...)
    return LinearAlgebra.norm(ss(sys), p; kwargs...)
end