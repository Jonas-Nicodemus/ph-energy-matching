import ControlSystems: StateSpace, ss

"""
    Σ = ss(J, R, Q, G, P, S, N)
    Σ = ss(Σph)

Converts a `PortHamiltonianStateSpace` to a standard `ControlSystemsBase.StateSpace`.
"""
function ss(J, R, Q, G, P, S, N)
    A, B, C, D = compose(J, R, Q, G, P, S, N)
    return ss(A, B, C, D)
end

function ss(Σph::PortHamiltonianStateSpace)
    return ss(Σph.J, Σph.R, Σph.Q, Σph.G, Σph.P, Σph.S, Σph.N)
end

"""
    Σph = phss(Σ)
    Σph = phss(Σ, X)

Converts a `ControlSystemsBase.StateSpace` to a `PortHamiltonianStateSpace` by executing [`decompose`](@ref).
If `X` is not provided, the minimal solution of the KYP inequality is used.
"""
function phss(Σ::StateSpace)
    X = kyp(Σ; min=true)
    return phss(Σ, X) 
end

function phss(Σ::StateSpace, X)
    J, R, Q, G, P, S, N = decompose(Σ.A, Σ.B, Σ.C, Σ.D, X)
    return phss(J, R, Q, G, P, S, N)
end

"""
    Σqo = hdss(J, R, Q, G, P, S, N)
    Σqo = hdss(Σph)

Converts a `PortHamiltonianStateSpace` to a `QuadraticOutputSystem` (Hamiltonian dynamic).
"""
function hdss(J, R, Q, G, P, S, N)
    A, B, _, _ = compose(J, R, Q, G, P, S, N)
    return qoss(A, B, 1/2*Q)
end

function hdss(Σph::PortHamiltonianStateSpace)
    return hdss(Σph.J, Σph.R, Σph.Q, Σph.G, Σph.P, Σph.S, Σph.N)
end

"""
    A, B, C, D = compose(J, R, Q, G, P, S, N)

Composes the port-Hamiltonian matrices according to standard state-space matrices according to

    A = (J - R) * Q
    B = G - P
    C = (G + P)' * Q
    D = S - N.
"""
function compose(J, R, Q, G, P, S, N)
    A = (J - R) * Q
    B = G - P
    C = (G + P)' * Q
    D = S - N
    return A, B, C, D
end

"""
    J, R, Q, G, P, S, N = decompose(A, B, C, D, X)

Decomposes the standard state-space matrices to port-Hamiltonian matrices according to

    Q = X
    J = skew(A / X)
    R = -sym(A / X)
    G = 0.5 * (X \\ C' + B)
    P = 0.5 * (X \\ C' - B)
    S = sym(D)
    N = skew(D).
"""
function decompose(A, B, C, D, X)
    J = skew(A / X)
    R = -sym(A / X)
    Q = X
    G = 0.5 * (X \ C' + B)
    P = 0.5 * (X \ C' - B)
    S = sym(D)
    N = skew(D)
    return J, R, Q, G, P, S, N
end

