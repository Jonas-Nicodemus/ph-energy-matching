"""
    D = duplication(n)

Returns the duplication matrix of size `n` x `n` to transform the half-vectorization to a full vectorization.
"""
function duplication(n)
    i = collect(1:n*n)
    J = zeros(Int, n,n)
    z = 1
    for x = 1:n
      for y = x:n
        J[x,y] = J[y,x] = z  # duplication occurs here
        z += 1
      end
    end
    # pass 3 vectors {i-index, j-index, value} to sparse()
    return sparse(i, vec(J), 1)
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

function unvec(v)
    n = Int(sqrt(length(v)))
    return reshape(v, n, n)
end