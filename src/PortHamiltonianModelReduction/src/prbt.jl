import ControlSystems: StateSpace

"""
    Σr = bt(Σ::StateSpace, r; Lx=grampd(Σ, :o), Ly=cholesky(gram(Σ, :c)).U)

Reduces the state dimension of the system `Σ` to `r` using standard square root balanced truncation.
The cholesky factors of the Gramians can be passed as optional arguments `Lx` and `Ly`.
The default values are the cholesky factors of the controllability and observability Gramians.
"""
function bt(Σ::StateSpace, r; Lx=grampd(Σ, :o), Ly=cholesky(gram(Σ, :c)).U)
    U, σ, Z = svd(Ly * Lx')
    U = U[:, 1:r]
    Z = Z[:, 1:r]
    σ = σ[1:r]
    V = Ly' * U * diagm(σ .^ (-1 / 2))
    W = Lx' * Z * diagm(σ .^ (-1 / 2))
    Ar = W' * Σ.A * V
    Br = W' * Σ.B
    Cr = Σ.C * V

    return ss(Ar, Br, Cr, Σ.D)
end

"""
    Σr = prbt(Σ::StateSpace, r)

Reduces the state dimension of the system `Σ` to `r` using positve real balanced truncation, which is passivity preserving.

**U. Desai and D.Pal.** 
[A transformation approach to stochastic model reduction](https://doi.org/10.1109/TAC.1984.1103438).
IEEE Trans. Automat. Control, 29(12):1097--1100, 1984.
"""
function prbt(Σ::StateSpace, r)
    Lx = prgrampd(Σ, :o)
    Ly = prgrampd(Σ, :c)
    return bt(Σ, r; Lx=Lx, Ly=Ly)
end

function prbt(Σph::PortHamiltonianStateSpace, r)
    Σr = prbt(ss(Σph), r)
    return phss(Σr) 
end