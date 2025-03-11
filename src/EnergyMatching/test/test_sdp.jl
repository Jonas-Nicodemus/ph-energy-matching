module Test_sdp

using EnergyMatching
using ControlSystemsBase
using QuadraticOutputSystems, PortHamiltonianSystems
using Test

@testset "test_sdp.jl" begin
    A = [-2. 1.; -1. -1.]
    B = [6.; 0;;]
    C = B'
    D = [1;;]
    Q = [1. 0.; 0. 1.]
    
    Σ = phss(ss(A, B, C, D), Q)
    Σr = phss(ss(A[1:1,1:1], B[1:1,1:1], C[1:1,1:1], D), Q[1:1,1:1])  
        
    @testset "sdp" begin
        for solver in [:Hypatia, :COSMO]
            Σphr = matchnrg(Σ, Σr; solver=solver)
            @test Σphr.Q ≈ [160/169;;] atol = 1e-3
            @test norm(hdss(Σ) - hdss(Σphr))^2 ≈ 19 + 81/4 * (160/169)^2 - 2 * 3240/169 * (160/169) atol = 1e-3
        end
    end
end

end