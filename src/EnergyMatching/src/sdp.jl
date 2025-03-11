"""
    Σphr = sdp(Σ::QuadraticOutputStateSpace, Σr::StateSpace; optimizer=COSMO.Optimizer, ε=1e-8, kwargs...)

Solves the energy matching problem using semidefinite programming.
"""
function sdp(Σ::QuadraticOutputStateSpace, Σr::StateSpace, optimizer=Hypatia.Optimizer; ε=1e-8, kwargs...)
    r = Σr.nx

    # Precompute
    Pr = gram(Σr, :c)
    Y =  sylvc(Σ.A, Σr.A', -Σ.B * Σr.B')
    h2 = h2norm(Σ)^2

    model = JuMP.Model(optimizer)
    for kwarg in keys(kwargs)
        set_optimizer_attribute(model, String(kwarg), kwargs[kwarg])
    end

    @variable(model, X[1:r, 1:r], PSD)
    Qr = 2 * X

    # Constraint
    W = [-Σr.A'*Qr-Qr*Σr.A Σr.C'-Qr*Σr.B;Σr.C-Σr.B'*Qr Σr.D+Σr.D']
    @constraint(model, W - ε * I >= 0, PSDCone())
    
    # Objective and optimize   
    @objective(model, Min, h2 + tr(Pr * X * Pr * X) - 2*tr(Y' * unvec(Σ.M) * Y * X))
    JuMP.optimize!(model)

    return phss(Σr, value.(Qr))
end
