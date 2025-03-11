import LinearAlgebra: norm

"""
    norm(Σph, p=2; kwargs...)

Converts a `PortHamiltonianStateSpace` to a `StateSpace` and calls `norm` on it.
For more details see [ControlSystems.jl](https://github.com/JuliaControl/ControlSystems.jl) package.
"""
function norm(sys::PortHamiltonianStateSpace, p::Real=2; kwargs...)
    return norm(ss(sys), p; kwargs...)
end
