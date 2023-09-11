using Test

@testset "PortHamiltonianSystems.jl" begin
    include("test_convert.jl")
    include("test_gramians.jl")
    include("test_kyp.jl")
    include("test_utils.jl")
end
