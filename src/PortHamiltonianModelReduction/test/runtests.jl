using Test

@testset "PortHamiltonianModelReduction.jl" begin
    include("test_phirka.jl")
    include("test_prbt.jl")
end
