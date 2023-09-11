import LinearAlgebra: norm

"""
    norm(Σqo, p=2; kwargs...)
    h2norm(Σqo; kwargs...)

Computes the H2-norm of the quadratic output state-space model `Σqo::QuadraticOutputStateSpace{T}`.

**P. Benner, P. Goyal, and I. Pontes Duff.**
[Gramians, energy functionals, and balanced truncation for linear dynamical systems with quadratic outputs](https://doi.org/10.1109/TAC.2021.3086319).
IEEE Trans. Automat. Control, 67(2), 2022.
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

function h2norm(Σ::QuadraticOutputStateSpace; kwargs...)
    try
        O = gram(Σ, :o; kwargs...)
        return sqrt(tr(Σ.B' * O * Σ.B))
    catch
        return Inf
    end
end

"""
    h2inner(Σ1, Σ2)

Computes the H2-inner product of the quadratic output state-space models 
`Σ1::QuadraticOutputStateSpace` and `Σ2::QuadraticOutputStateSpace`.

**P. Benner, P. Goyal, and I. Pontes Duff.**
[Gramians, energy functionals, and balanced truncation for linear dynamical systems with quadratic outputs](https://doi.org/10.1109/TAC.2021.3086319).
IEEE Trans. Automat. Control, 67(2), 2022.
"""
function h2inner(Σ1::QuadraticOutputStateSpace, Σ2::QuadraticOutputStateSpace)
    Y = _h2inner_Y(Σ1, Σ2)
    Z = _h2inner_Z(Σ1, Σ2, Y)
    return tr(Σ1.B' * Z * Σ2.B)
end

function _h2inner_Y(Σ1::QuadraticOutputStateSpace, Σ2::QuadraticOutputStateSpace)
    return sylvc(Σ1.A, Σ2.A', -Σ1.B * Σ2.B')
end

function _h2inner_Z(Σ1::QuadraticOutputStateSpace, Σ2::QuadraticOutputStateSpace)
    Y = _h2inner_Y(Σ1, Σ2)
    return _h2inner_Z(Σ1, Σ2, Y)
end

function _h2inner_Z(Σ1::QuadraticOutputStateSpace, Σ2::QuadraticOutputStateSpace, Y)
    return sylvc(Σ1.A', Σ2.A, -Σ1.M * Y * Σ2.M)
end

function _h2inner_lin_sys(Σ1::QuadraticOutputStateSpace, Σ2::QuadraticOutputStateSpace, Y)
    n = size(Σ1.A, 1)
    r = size(Σ2.A, 1)
    return (kron(I(r), Σ1.A') + kron(Σ2.A', I(n))) \ kron(I(r), 2 * Σ1.M * Y)
end

function _h2inner_Z_kron(Σ1::QuadraticOutputStateSpace, Σ2::QuadraticOutputStateSpace, Y)
    n = size(Σ1.A, 1)
    r = size(Σ2.A, 1)
    F = _h2inner_lin_sys(Σ1, Σ2, Y) 
    return -1/4*reshape(F * reshape(2*Σ2.M, :, 1), n, r)
end

function _h2inner_kron(Σ1::QuadraticOutputStateSpace, Σ2::QuadraticOutputStateSpace)
    m = size(Σ1.B, 2)
    Y = _h2inner_Y(Σ1, Σ2)
    F = _h2inner_lin_sys(Σ1, Σ2, Y) 
    F2 = -1/4*kron(Σ1.B', Σ2.B') * F
    return reshape(F2 * reshape(2*Σ2.M, :, 1), m, m)
end