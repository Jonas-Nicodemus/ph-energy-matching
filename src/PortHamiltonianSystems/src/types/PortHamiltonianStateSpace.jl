import Base: +, -

"""
PortHamiltonianStateSpace{T}

An object representing a port-Hamiltonian state-space system.


    dx(t)/dt = (J-R)Qx(t) + (G-P)u(t)
    y(t)     = (G+P)'Qx(t) + (S-N)u(t)


See the function [`phss`](@ref) for a user facing constructor.

# Fields:
- `J::Matrix{T}`
- `R::Matrix{T}`
- `Q::Matrix{T}`
- `G::Matrix{T}`
- `P::Matrix{T}`
- `S::Matrix{T}`
- `N::Matrix{T}`
"""
struct PortHamiltonianStateSpace{T}
    J::Matrix{T}
    R::Matrix{T}
    Q::Matrix{T}
    G::Matrix{T}
    P::Matrix{T}
    S::Matrix{T}
    N::Matrix{T}
end

"""
    Σph = phss(J, R, Q, G, P, S, N)

Creates a port-Hamiltonian state-space model `Σph::PortHamiltonianStateSpace{T}`
with matrix element type `T`.
"""
phss(args...;kwargs...) = PortHamiltonianStateSpace(args...;kwargs...)

function -(sys1::PortHamiltonianStateSpace, sys2::PortHamiltonianStateSpace)
    return ss(sys1) - ss(sys2)
end

function -(sys1::ControlSystems.AbstractStateSpace, sys2::PortHamiltonianStateSpace)
    return sys1 - ss(sys)
end

function -(sys1::PortHamiltonianStateSpace, sys2::ControlSystems.AbstractStateSpace)
    return ss(sys1) - ss(sys2)
end

_string_mat_with_headers(X) = _string_mat_with_headers(Matrix(X))
function _string_mat_with_headers(X::Matrix)
    p = (io, m) -> Base.print_matrix(io, m)
    return replace(sprint(p, X), "\"" => " ")
end

Base.print(io::IO, sys::PortHamiltonianStateSpace) = show(io, sys)

function Base.show(io::IO, sys::PortHamiltonianStateSpace)
    println(io, typeof(sys))
    println(io, "J = \n", _string_mat_with_headers(sys.J))
    println(io, "R = \n", _string_mat_with_headers(sys.R))
    println(io, "Q = \n", _string_mat_with_headers(sys.Q))
    println(io, "G = \n", _string_mat_with_headers(sys.G))
    println(io, "P = \n", _string_mat_with_headers(sys.P))
    println(io, "S = \n", _string_mat_with_headers(sys.S))
    println(io, "N = \n", _string_mat_with_headers(sys.N))
end