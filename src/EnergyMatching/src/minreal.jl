function minreal(Σph::PortHamiltonianStateSpace)
    if !isposdef(Σph.Q)
        Σph = truncno(Σph)
    end

    return truncnc(Σph)
end

function truncno(Σph::PortHamiltonianStateSpace; ϵ=1e-12)
    F = eigen(Σph.Q)
    i = findall(x->abs(x)>ϵ, F.values)
    V = F.vectors[:,i]
    return phss(skew(V'*Σph.J*V), sym(V'*Σph.R*V), sym(V'*Σph.Q*V), V'*Σph.G, V'*Σph.P, Σph.S, Σph.N)  
end

function truncnc(Σph::PortHamiltonianStateSpace; ϵ=1e-12)
    L = cholesky(Σph.Q).L
    
    Σph1 = phss(L'*Σph.J*L, L'*Σph.R*L, I(size(L,1)), L'*Σph.G, L'*Σph.P, Σph.S, Σph.N)
    
    P = gram(Σph1, :c)
    F = eigen(P, sortby=-)
    i = findall(x->abs(x)>ϵ, F.values)
    r = length(i)
    V = F.vectors[:,i]
    
    return phss(V'*Σph1.J*V, V'*Σph1.R*V, I(r), V'*Σph1.G, V'*Σph1.P, Σph.S, Σph.N)
end
