abstract type AbstractModel end

struct PortHamiltonianModel <: AbstractModel 
    name::String
    Σ::PortHamiltonianStateSpace
end

function dedescriptorize(E, J, R, Q, G, P, S, N)
    Q = sym(inv(E))
    return J, R, Q, G, P, S, N
end

function dedescriptorize(phsys::PHSystem)
    (; E, J, R, Q, G, P, S, N) = phsys
    J, R, Q, G, P, S, N = dedescriptorize(E, J, R, Q, G, P, S, N)
    return PHSystem(Matrix(1.0I, size(J)), J, R, Q, G, P, S, N)
end

function densify(phsys::PHSystem)
    E = Array(phsys.E)
    J = Array(phsys.J)
    R = Array(phsys.R)
    Q = Array(phsys.Q)
    G = Array(phsys.G)
    P = Array(phsys.P)
    S = Array(phsys.S)
    N = Array(phsys.N)
    return PHSystem(E, J, R, Q, G, P, S, N)
end

function addfeedthrough(phsys::PHSystem, ε)
    S = phsys.S + I * ε
    return PHSystem(phsys.E, phsys.J, phsys.R, phsys.Q, phsys.G, phsys.P, S, phsys.N)
end

function prepare(phsys::PHSystem; ε = 1e-6)
    phsys = densify(phsys)
    
    if !(phsys.E ≈ I)
        phsys = dedescriptorize(phsys)
    end
    
    phsys = addfeedthrough(phsys, ε)

    return phss(phsys.J, phsys.R, phsys.Q, phsys.G, phsys.P, phsys.S, phsys.N)
end

function msd(ε = 1e-6)
    return PortHamiltonianModel("msd", prepare(PHSystem(SingleMSDConfig()); ε=ε)) 
end
    
function poro(n=320, ε = 1e-6)
    return PortHamiltonianModel("poro", prepare(PHSystem(PoroElasticityConfig(n=n)); ε=ε))
end

function msd_Xmin(ε = 1e-6)
    fom = msd(ε)
    Xmin = kyp_min(fom.Σ)
    return PortHamiltonianModel("msd_Xmin", phss(ss(fom.Σ), Xmin))
end