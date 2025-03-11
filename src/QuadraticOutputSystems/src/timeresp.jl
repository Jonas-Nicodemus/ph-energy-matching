import ControlSystems: lsim, SimResult

"""
    result = lsim(Σqo, u, t; kwargs...)

Calculate the time response of the quadratic output state-space model `Σqo::QuadraticOutputStateSpace{T}` by
first treating it as standard state-space model and calling `ControlSystems.lsim` on it, and then
calculating the quadratic output `y_2(t)` as `y_2(t) = M(x(t)⊗x(t))` and add it to the linear output.
"""
function lsim(sys::QuadraticOutputStateSpace, u::Function, t::AbstractVector; kwargs...)
    res = lsim(ss(sys), u, t; kwargs...)
    n = size(sys.A, 1)
    y2 = zero(res.y)

    for i in 1:size(res.y, 1)
        y2[i,:] = sum(res.x .* (reshape(sys.M[i,:], n,n) * res.x), dims=1)
    end
    
    return SimResult(res.y .+ y2, res.t, res.x, res.u, sys)
end
