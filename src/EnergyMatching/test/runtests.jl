using EnergyMatching
using Test

@testset "EnergyMatching.jl" begin
    include("test_reshape.jl")
    include("test_bestricc.jl")
    include("test_barrier.jl")
    include("test_sdp.jl")
end
