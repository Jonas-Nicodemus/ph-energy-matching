import Base: +, -

"""
QuadraticOutputStateSpace{T}

An object representing a quadratic output state-space system.

    dx(t)/dt = Ax(t) + Bu(t)
    y(t)     = x(t)'Mx(t)

See the function [`qoss`](@ref) for a user facing constructor.

# Fields:
- `A::Matrix{T}`
- `B::Matrix{T}`
- `M::Matrix{T}`
"""
struct QuadraticOutputStateSpace{T}
    A::Matrix{T}
    B::Matrix{T}
    M::Matrix{T}
end

"""
    Σqo = qoss(A, B, M)

Creates a quadratic output state-space model `Σqo::QuadraticOutputStateSpace{T}`
with matrix element type `T`.
"""
qoss(args...;kwargs...) = QuadraticOutputStateSpace(args...;kwargs...)

function +(Σqo1::QuadraticOutputStateSpace{T1}, Σqo2::QuadraticOutputStateSpace{T2}) where {T1,T2}
    T = promote_type(T1, T2)

    n1 = size(Σqo1.A, 1)
    n2 = size(Σqo2.A, 1)
    A = [Σqo1.A zeros(T,n1,n2);
         zeros(T,n2,n1) Σqo2.A]
    B = [Σqo1.B ; Σqo2.B]
    M = [Σqo1.M zeros(T,n1,n2);
        zeros(T,n2,n1) Σqo2.M]
    return qoss(A, B, M)
end

function -(Σqo::QuadraticOutputStateSpace) 
    return qoss(Σqo.A, Σqo.B, -Σqo.M)
end

function -(Σqo1::QuadraticOutputStateSpace{T1}, Σqo2::QuadraticOutputStateSpace{T2}) where {T1,T2}
    return Σqo1 + (-Σqo2)
end

_string_mat_with_headers(X) = _string_mat_with_headers(Matrix(X))
function _string_mat_with_headers(X::Matrix)
    p = (io, m) -> Base.print_matrix(io, m)
    return replace(sprint(p, X), "\"" => " ")
end

Base.print(io::IO, sys::QuadraticOutputStateSpace) = show(io, sys)

function Base.show(io::IO, sys::QuadraticOutputStateSpace)
    println(io, typeof(sys)) 
    println(io, "A = \n", _string_mat_with_headers(sys.A))
    println(io, "B = \n", _string_mat_with_headers(sys.B))
    println(io, "M = \n", _string_mat_with_headers(sys.M))
end