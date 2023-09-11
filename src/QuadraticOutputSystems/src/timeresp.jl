import ControlSystems: lsim, SimResult

"""
    result = lsim(Σqo, u, t; kwargs...)

Calculate the time response of the quadratic output state-space model `Σqo::QuadraticOutputStateSpace{T}` by
first treating it as standard state-space model and calling `ControlSystems.lsim` on it, and then
calculating the output `y(t)` as `y(t) = x(t)'Mx(t)`.
"""
function lsim(sys::QuadraticOutputStateSpace, u::Function, t::AbstractVector; kwargs...)
    res = lsim(ss(sys), u, t; kwargs...)
    y = sum(res.x .* (sys.M * res.x); dims=1)
    return SimResult(y, res.t, res.x, res.u, sys)
end