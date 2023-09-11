"""
    Msym = sym(M)

Returns the symmetric part of a matrix `M`.
"""
function sym(M)
    return 0.5*(M+M')
end

"""
    Mskew = skew(M)

Returns the skew-symmetric part of a matrix `M`.
"""   
function skew(M)
    return 0.5*(M-M')
end

"""
    Mpsd = project_psd(M, eigtol)

Returns the nearest positive semi-definite matrix to `M` by setting negative eigenvalues to `eigtol`.
"""    
function project_psd(M, eigtol=1e-8)
    Λ, U = eigen(sym(M))
    return sym(U * diagm(clamp.(Λ, eigtol, Inf)) * U')
end

"""
    is_sym(M)

Returns `true` if `M` is symmetric, otherwise `false`.    
"""
function is_sym(M)
    return isapprox(M, M')
end

"""
    is_skew(M)

Returns `true` if `M` is skew-symmetric, otherwise `false`. 
"""    
function is_skew(M)
    return isapprox(M, -M')
end

"""
    is_psd(M)

Returns `true` if `M` is positive semi-definite, otherwise `false`.
"""
function is_psd(M)
    return isapprox(M, M') && all(eigvals(M) .>= 0)
end