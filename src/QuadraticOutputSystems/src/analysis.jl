import LinearAlgebra: norm

"""
    norm(Σqo, p=2; kwargs...)
    h2norm(Σqo; kwargs...)

Computes the H2 norm of the quadratic output state-space model `Σqo::QuadraticOutputStateSpace`.
"""
function norm(Σ::QuadraticOutputStateSpace, p::Real=2; kwargs...)
    if p == 2
        return h2norm(Σ; kwargs...)
    elseif p == Inf
        error("`Inf` norm not implemented yet.")
    else
        error("`p` must be either `2` or `Inf`")
    end
end

function h2norm(Σ::QuadraticOutputStateSpace, opt=:c; kwargs...)
    if !isstable(ss(Σ))
        return Inf
    end

    X = gram(Σ, opt; kwargs...)
    if opt == :o
        return sqrt(tr(Σ.B' * X * Σ.B))
    elseif opt == :c
        n = size(Σ.A, 1)
        return sqrt(tr(Σ.C * X * Σ.C') + tr(
            sum([X * reshape(Mi,n,n) * X * reshape(Mi,n,n) for Mi ∈ eachrow(Σ.M)])
            ))
    else
        error("Invalid option `opt`")
    end 
end

"""
    h2inner(Σ1, Σ2)

Computes the H2 inner product of the quadratic output state-space models 
`Σ1::QuadraticOutputStateSpace` and `Σ2::QuadraticOutputStateSpace`.
"""
function h2inner(Σ1::QuadraticOutputStateSpace, Σ2::QuadraticOutputStateSpace, opt=:c)
    Y = _h2inner_Y(Σ1, Σ2)
    if opt == :o
        Z = _h2inner_Z(Σ1, Σ2, Y)
        return tr(Σ1.B' * Z * Σ2.B)
    elseif opt == :c
        p, n = size(Σ1.C)
        r = size(Σ2.A, 1)
        return tr(Σ1.C*Y*Σ2.C') + tr(sum([Y' * reshape(Σ1.M[i,:],n,n) * Y * reshape(Σ2.M[i,:],r,r) for i ∈ 1:p]))
    else
        error("Invalid option `opt`")
    end
end

function _h2inner_Y(Σ1::QuadraticOutputStateSpace, Σ2::QuadraticOutputStateSpace)
    return sylvc(Σ1.A, Σ2.A', -Σ1.B * Σ2.B')
end

function _h2inner_Z(Σ1::QuadraticOutputStateSpace, Σ2::QuadraticOutputStateSpace)
    Y = _h2inner_Y(Σ1, Σ2)
    return _h2inner_Z(Σ1, Σ2, Y)
end

function _h2inner_Z(Σ1::QuadraticOutputStateSpace, Σ2::QuadraticOutputStateSpace, Y)
    p, n = size(Σ1.C)
    r = size(Σ2.A, 1)
    return sylvc(Σ1.A', Σ2.A, - sum([reshape(Σ1.M[i,:],n,n) * Y * reshape(Σ2.M[i,:],r,r) for i ∈ 1:p]) - Σ1.C' * Σ2.C)
end
