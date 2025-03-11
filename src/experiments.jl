abstract type AbstractExperiment end

struct Experiment <: AbstractExperiment 
    name::String
    fom::PortHamiltonianStateSpace
end

function msd(ε = 1e-6)
    # E = I
    (; E, J, R, Q, G, P, S, N) = PHSystem(SingleMSDConfig())
    return Experiment("msd", phss(J, R, Q, G, P, S + I * ε, N))
end
    
function poro(n=980, ε = 1e-6)
    # Q = I
    (; E, J, R, Q, G, P, S, N) = PHSystem(PoroElasticityConfig(n=n))
    return Experiment("poro", phss(J, R, sym(E\I(n)), G, P, S + I * ε, N))
end

function msd_Xmin(ε = 1e-6)
    Σ = msd(ε).Σ
    Xmin = kypmin(Σ)
    return Experiment("msd_Xmin", phss(ss(Σ), Xmin))
end

function rcl(ε = 0)
    n = 5000
    N = Int(n/2)
    ri = 0.2 * ones(N + 1)
    ri[end] = 0.4
    Σph = rclnet(n)
    return Experiment("rcl", phss(Σph.J, Σph.R, Σph.Q, Σph.G, Σph.P, Σph.S + I * ε, Σph.N))
end

function rclnet(n, m=1; L=0.2, C=0.2, R=1.0)
    #  Structure preserving model reduction of port-Hamiltonian systems by moment matching at infinity
    #  Rostyslav V. Polyuga, Arjan van der Schaft
    if mod(n,2) != 0
        @error "$n must be even"
    end

    if length(L) != 1 && length(L) != Int(n/2)
        @error "L must be sclar or of length $(Int(n/2))"
    end
    if length(C) != 1 && length(C) != Int(n/2)
        @error "C must be sclar or of length $(Int(n/2))"
    end
    if length(R) != 1 && length(R) != Int(n/2) + 1
        @error "R must be sclar or of length $(Int(n/2) + 1)"
    end

    J = diagm(1 => ones(n-1)) + diagm(-1 => -1*ones(n-1))
        
    if length(R) == 1
        R = R * ones(Int(n/2) + 1)
    end

    diagR = zeros(n)
    index_diagR = 2:2:n
    diagR[index_diagR] = R[1:end-1]
    diagR[end] = diagR[end] + R[end]
    R = diagm(diagR)
    
    diagQ = zeros(n)
    index_diagQ_C = 1:2:n-1
    index_diagQ_L = 2:2:n
    diagQ[index_diagQ_C] .= 1 ./ C
    diagQ[index_diagQ_L] .= 1 ./ L
    Q = diagm(diagQ)

    # G
    if m == 1
        #  SISO
        G = zeros(n,1);
        G[1, 1] = 1;
    elseif m == 2
        #  MIMO
        G = zeros(n,2);
        G[1, 1] = 1;
        G[2, 2] = 1;
    else
        @error "Only m=1 and m=2 are supported, not m=$m"
    end
    
    return phss(J,R,Q,G)     
end 