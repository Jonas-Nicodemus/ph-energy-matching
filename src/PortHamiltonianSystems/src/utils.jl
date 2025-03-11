"""
    lrcf(X, trunc_tol)

Computes an approximate low-rank Cholesky-like 
factorization of a symmetric positive semi-definite matrix ``X``
s.t. ``X = Z^T Z`` (up to a prescribed tolerance `trunc_tol`).
"""
function lrcf(X, trunc_tol)
    d,L = eigen(Symmetric(X))
    # remove negative eigenvalues (numerical errors)
    idx = findall(v -> v >= 0, d)
    dr, Lr = truncation(d[idx], L[:, idx], trunc_tol)
    return (Lr*diagm(sqrt.(dr)))'
end

"""
    truncation(d, L, trunc_tol) -> (dr, Lr)

Computes a rank revealing factorization for a given LDL-decomposition of
``S = L * \\mathrm{diag}(d) * L^T`` (up to a prescribed tolerance `trunc_tol`)
such that ``L_r * diag(d_r) * L_r^T \\approx S``.
"""
function truncation(d, L, trunc_tol)
    Q,R = qr(L)
    tmp = Symmetric(R*diagm(d)*R')
    d,U = eigen(tmp)
    p = sortperm(d, by=abs, rev=true)
    d = d[p]
    trunc_index = findlast(abs.(d/d[1]) .>= trunc_tol)
    return d[1:trunc_index], Q*U[:,p[1:trunc_index]]
end

"""
    Msym = sym(M)

Returns the symmetric part of a matrix `M`.
"""
function sym(M)
    return hermitianpart(M)
end

"""
    Mskew = skew(M)

Returns the skew-symmetric part of a matrix `M`.
"""   
function skew(M)
    return 0.5*(M-M')
end

"""
    Mpsd = project_psd(M; eigtol=1e-8)

Returns the nearest positive semi-definite matrix to `M` by setting negative eigenvalues to `eigtol`.
"""    
function project_psd(M; eigtol=1e-8)
    Λ, U = eigen(sym(M))
    return sym(U * diagm(clamp.(Λ, eigtol, Inf)) * U')
end

"""
    issym(M)

Returns `true` if `M` is symmetric, otherwise `false`.    
"""
function issym(M)
    return issymmetric(M)
end

"""
    isskew(M)

Returns `true` if `M` is skew-symmetric, otherwise `false`. 
"""    
function isskew(M)
    return isapprox(M, -M')
end

"""
    ispsd(M)

Returns `true` if `M` is positive semi-definite, otherwise `false`.
"""
function ispsd(M; ϵ=1e-8)
    return isposdef(M + ϵ*I) 
end

"""
    D = duplication(n)

Returns the duplication matrix of size `n` x `n`.
"""
function duplication(n)
    return duplication_matrix(n)
end

"""
    v = vech(M)

Returns the lower triangular part of `M` as a vector. Aka the half-vectorization of `M`.
For the inverse operation see [`unvech`](@ref).
"""
function vech(M)
    n = size(M, 1)
    @assert size(M, 2) == n
    return M[tril!(trues(n, n))]
end

"""
    M = unvech(v)

Returns the symmetric matrix `M` from the half-vectorized `v`, i.e., the inverse of [`vech`](@ref).
"""
function unvech(v)
    n = Int((sqrt(8 * length(v) + 1) - 1) / 2)
    D = duplication(n)
    return reshape(D * v, n, n)
end

"""
    M = unvec(v)

Returns the matrix `M` from the vectorized `v`, i.e., the inverse of `vec`.
"""
function unvec(v)
    n = Int(sqrt(length(v)))
    return reshape(v, n, n)
end

sparsify!(x; ε=1e-12) = x[abs.(x) .< ε] .= 0;
function sparsify(x; ε=1e-12)
    y = copy(x)
    sparsify!(y; ε=ε)
    return y
end

isinvertible(A::Matrix{Int}) = !isapprox(det(A), 0)
isinvertible(A::Matrix) = !isapprox(det(BigFloat.(A)), 0, atol = 1e-18)
