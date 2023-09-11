"""
    Σphr = sdp(Σ::QuadraticOutputStateSpace, Σr::StateSpace; optimizer=COSMO.Optimizer, ε=1e-8, kwargs...)

Solves the energy matching problem using semidefinite programming.
"""
function sdp(Σ::QuadraticOutputStateSpace, Σr::StateSpace, optimizer=COSMO.Optimizer; ε=1e-8, kwargs...)
    n  = size(Σ.A, 1)
    r = size(Σr.A, 1)

    # Precompute
    Pr = gram(Σr, :c)
    Y =  sylvc(Σ.A, Σr.A', -Σ.B * Σr.B')
    h2 = norm(Σ, 2)^2

    model = JuMP.Model(optimizer)
    for kwarg in keys(kwargs)
        set_optimizer_attribute(model, String(kwarg), kwargs[kwarg])
    end

    @variable(model, Qr[1:r, 1:r], PSD)

    F = -(kron(I(r), Σ.A') + kron(Σr.A', I(n))) \ kron(I(r), 2 * Σ.M * Y)
    Z = reshape(F * reshape(Qr/4, :, 1), n, r)

    # Reduced quadratic output observability gramian for reduced system's H2 norm
    @variable(model, Or[1:r, 1:r], PSD)

    # Constraints
    W = Symmetric([-Σr.A'*Qr-Qr*Σr.A Σr.C'-Qr*Σr.B;Σr.C-Σr.B'*Qr Σr.D+Σr.D']) 
    O = Symmetric([-Σr.A' * Or - Or * Σr.A   -Qr/2; -Qr/2   inv(Pr)])
    @constraint(model, W - ε * I >= 0, PSDCone())
    @constraint(model, O - ε * I >= 0, PSDCone())

    # Objective and optimize   
    @objective(model, Min, h2 + tr(Σr.B' * Or * Σr.B) - 2 * tr(Σ.B' * Z * Σr.B))
    JuMP.optimize!(model)

    return phss(Σr, value.(Qr))
end