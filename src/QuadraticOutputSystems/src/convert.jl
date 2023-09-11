import ControlSystems: ss

function ss(Σqo::QuadraticOutputStateSpace)
    n, m = size(Σqo.B)
    return ss(Σqo.A, Σqo.B, zeros(0, n), zeros(0, m)) 
end