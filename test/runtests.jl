using DrWatson, Test
@quickactivate "ph-energy-matching"

# Run test suite
println("Starting tests")
ti = time()

@testset "ph-energy-matching tests" begin
    include(srcdir("EnergyMatching", "test", "runtests.jl"))
    include(srcdir("PortHamiltonianModelReduction", "test", "runtests.jl"))
    include(srcdir("PortHamiltonianSystems", "test", "runtests.jl"))
    include(srcdir("QuadraticOutputSystems", "test", "runtests.jl"))
end

ti = time() - ti
println("\nTest took total time of:")
println(round(ti/60, digits = 3), " minutes")
