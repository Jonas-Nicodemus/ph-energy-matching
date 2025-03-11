import Base: +, -
import ControlSystemsBase: to_matrix, AbstractNumOrArray

"""
QuadraticOutputStateSpace{T}

An object representing a quadratic output state space system.

    dx(t)/dt = Ax(t) + Bu(t)
    y(t)     = Cx(t) + M(x(t)⊗x(t))

See the function [`qoss`](@ref) for a user facing constructor.

# Fields:
- `A::Matrix{T}`
- `B::Matrix{T}`
- `C::Matrix{T}`
- `M::Matrix{T}`
"""
struct QuadraticOutputStateSpace{T}
    A::Matrix{T}
    B::Matrix{T}
    C::Matrix{T}
    M::Matrix{T}
end

function QuadraticOutputStateSpace(A::AbstractNumOrArray, B::AbstractNumOrArray, C::AbstractNumOrArray, M::AbstractNumOrArray)
    T = promote_type(eltype(A), eltype(B), eltype(C), eltype(M))
    if ndims(M) == 3
        M = reshape(M, size(M,1)*size(M,2), size(M,3))
    end

    return QuadraticOutputStateSpace{T}(to_matrix(T, A), to_matrix(T, B), to_matrix(T, C), to_matrix(T, M))
end

"""
    Σqo = qoss(A, B, C, M)
    Σqo = qoss(A, B, M)

Creates a quadratic output state-space model `Σqo::QuadraticOutputStateSpace{T}`
with matrix element type `T`.
"""
qoss(args...;kwargs...) = QuadraticOutputStateSpace(args...;kwargs...)

function qoss(A, B, M)
    n = size(A, 1)
    p = size(M, 1)
    return qoss(A, B, zeros(p, n), M)
end

function +(Σqo1::QuadraticOutputStateSpace{T1}, Σqo2::QuadraticOutputStateSpace{T2}) where {T1,T2}
    T = promote_type(T1, T2)

    n1 = size(Σqo1.A, 1)
    n2 = size(Σqo2.A, 1)
    A = [Σqo1.A zeros(T,n1,n2);
         zeros(T,n2,n1) Σqo2.A]
    B = [Σqo1.B ; Σqo2.B]
    C = [Σqo1.C Σqo2.C]
    
    p = size(Σqo1.M, 1)
    
    # not so nice
    M = zeros(p, (n1+n2)^2)
    for i in 1:p
        Me = [reshape(Σqo1.M[i,:], n1,n1) zeros(T,n1,n2);
              zeros(T,n2,n1) reshape(Σqo2.M[i,:], n2,n2)]
        M[i,:] = vec(Me)
    end

    return qoss(A, B, C, M)
end

function -(Σqo::QuadraticOutputStateSpace) 
    return qoss(Σqo.A, Σqo.B, -Σqo.C, -Σqo.M)
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
    println(io, "C = \n", _string_mat_with_headers(sys.C))
    println(io, "M = \n", _string_mat_with_headers(sys.M))
end
