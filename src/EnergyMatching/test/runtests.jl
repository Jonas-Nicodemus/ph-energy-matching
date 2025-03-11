using Test

@testset "EnergyMatching.jl" begin
    include("test_bestricc.jl")
    include("test_sdp.jl")
    include("test_barrier.jl")
    include("test_minreal.jl")
end
