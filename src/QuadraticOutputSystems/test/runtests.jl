using QuadraticOutputSystems
using Test

@testset "QuadraticOutputSystems.jl" begin
    include("test_gramians.jl")
    include("test_analysis.jl")
    include("test_timeresp.jl")
end
