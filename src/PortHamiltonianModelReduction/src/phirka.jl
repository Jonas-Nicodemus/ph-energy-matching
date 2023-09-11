"""
    Σr = phirka(Σph::PortHamiltonianStateSpace, r, num_runs; tol=1e-3, max_iter=50)

Calls [`phirka`](@ref) `num_runs` times and returns the best result.
"""
function phirka(Σph::PortHamiltonianStateSpace, r, num_runs; tol=1e-3, max_iter=50)
    h2 = Inf

    m = size(Σph.G, 2)
    Σphr = phss(zeros(r,r), zeros(r,r), zeros(r,r), zeros(r,m), zeros(r,m), zeros(m,m), zeros(m,m))

    for i in 1:num_runs
        Σphrᵢ = phirka(Σph, r; tol=tol, max_iter=max_iter)
        h2ᵢ = norm(Σph - Σphrᵢ)
        if h2ᵢ < h2
            h2 = h2ᵢ
            Σphr = Σphrᵢ
        end
    end

    return Σphr 
end


"""
    Σr = phirka(Σph::PortHamiltonianStateSpace, r; tol=1e-3, max_iter=50)

Reduces the state dimension of the port-Hamiltonian system `Σph` to `r` using the iterative rational Krylov algorithm.

**S. Gugercin, R. V. Polyuga, C. Beattie, and A. van der Schaft.**
[Structure-preserving tangential interpolation for model reduction of port-{Hamiltonian} systems](https://doi.org/10.1016/j.automatica.2012.05.052).
Automatica J. IFAC, 48(9):1963--1974, 2012.
"""
function phirka(Σph::PortHamiltonianStateSpace, r; tol=1e-3, max_iter=50)
    Σ = ss(Σph)
    A = Σ.A
    B = Σ.B
    Q = Σph.Q

    b = size(B,2)==1 ? ones(1,r) : randn(size(B,2),r)
    s = r==1 ? [(10+0*im)^0] : (10+0*im) .^ range(-1, stop=1, length=r)
    Ar, Br, V, W = interpolate(A, B, Q, s, b)

    s_prev = s
    s, b = interpolation_data(Ar, Br)
    i = 1

    while !converged(s, s_prev, tol) && max_iter > i
        Ar, Br, V, W = interpolate(A, B, Q, s, b)
        s_prev = s
        s, b = interpolation_data(Ar, Br)
        i += 1
    end
    
    if i == max_iter
        @warn "IRKA not converged after $max_iter iterations"
    else
        @debug "IRKA converged after $i iterations"
    end

    return project(Σph, V, W)
end

function interpolation_data(A, B)
    λ, U = eigen(A)
    s = -λ
    b = (U \ B)'
    return s, b
end

# function interpolation_data(A, B)
#     Λ, v = eigen(A')
#     b = transpose(v' * B)
#     return -Λ, b
# end

function interpolate(A, B, Q, s, b)
    V = construct_V(A, B, s, b)
    W = construct_W(V, Q)
    Ar = W' * A * V
    Br = W' * B
    return Ar, Br, V, W
end

function construct_V(A, B, s, b)
    n = size(A, 1)
    r = length(s)
    V = zeros(n, r)
    
    j = 1
    while j < r + 1
        x = (s[j] * I - A) \ (B * b[:,j])
        if abs(imag(s[j])) > 0
            V[:, j] = real(x)
            V[:, j+1] = imag(x)
            j += 2
        else
            V[:, j] = real(x)
            j += 1
        end
    end

    #  Orthonormalizing => worse results
    V = qr(V).Q[:,1:r]
    return V
end

function construct_W(V, Q)
    W = Q * V / (V' * Q * V)
    return W
end

function converged(s, s_prev, tol)
    sort_by = x -> (sign(real(x)), sign(imag(x)), real(x),imag(x))
    s = sort(s, by=sort_by)
    s_prev = sort(s_prev, by=sort_by)
    d = norm(s-s_prev)/norm(s_prev)
    return d < tol
end

function project(Σph, V, W)
    Jr = W' * Σph.J * W
    Rr = W' * Σph.R * W
    Qr = V' * Σph.Q * V
    Gr = W' * Σph.G
    Pr = W' * Σph.P
    
    return phss(Jr, Rr, Qr, Gr, Pr, Σph.S, Σph.N)
end