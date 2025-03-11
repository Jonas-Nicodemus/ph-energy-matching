import ControlSystems: ss

function ss(Σqo::QuadraticOutputStateSpace)
    return ss(Σqo.A, Σqo.B, Σqo.C, 0) 
end
