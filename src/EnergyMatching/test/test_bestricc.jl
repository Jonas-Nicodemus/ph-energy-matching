module Test_bestricc

using EnergyMatching
using ControlSystems
using QuadraticOutputSystems, PortHamiltonianSystems
using Test

@testset "test_bestricc.jl" begin
    A = [-2. 1.; -1. -1.]
    B = [6.; 0;;]
    C = B'
    D = [1;;]

    M = 1/2*[1. 0.; 0. 1.]
    
    Σ = qoss(A, B, M)
    Σr = ss(A[1:1,1:1], B[1:1,1:1], C[1:1,1:1], D)  
        
    @testset "bestricc" begin
        Σphr = EnergyMatching.bestricc(Σ, Σr)
        @test Σphr.Q ≈ kyp_min(Σr)
    end
end

end