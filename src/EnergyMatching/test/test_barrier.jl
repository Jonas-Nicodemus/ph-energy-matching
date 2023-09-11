module Test_barrier

using EnergyMatching
using LinearAlgebra, ControlSystems
using QuadraticOutputSystems, PortHamiltonianSystems
using Test

@testset "test_barrier.jl" begin
    A = [-2. 1.; -1. -1.]
    B = [6.; 0;;]
    C = B'
    D = [1;;]

    M = 1/2*[1. 0.; 0. 1.]
    
    Σ = qoss(A, B, M)

    @testset "objective" begin
        Σr = ss(A, B, C, D)
        fgo = EnergyMatching.objective(Σ, Σr)

        Q = kyp_min(Σr)
        x = vech(Q)

        @testset "f" begin
            fgo(x)[1] ≈ norm(Σ - qoss(Σr.A, Σr.B, 1/2*Q))^2
        end

        @testset "g" begin
            ε = 1e-8
            for (i,Δx) in enumerate(eachcol(I(length(x))))
                @test (fgo(x + ε * Δx)[1] - fgo(x)[1]) / ε ≈  fgo(x)[2][i] atol=1e-4 rtol=1e-4
            end
        end
    end

    @testset "barrier" begin
        Σr = ss(A[1:1,1:1], B[1:1,1:1], C[1:1,1:1], D) 
        Σphr = EnergyMatching.barrier(Σ, Σr)
        @test Σphr.Q ≈ [160/169;;]
        @test norm(Σ - hdss(Σphr))^2 ≈ 19 + 81/4 * (160/169)^2 - 2 * 3240/169 * (160/169) 

        Σr = ss(A, B, C, D)
        Σr0 = phss(Σr)
        Σphr = EnergyMatching.barrier(Σ, Σr; Σr0=Σr0)
        @test norm(Σphr.Q - 2*M) < 1e-6
    end
end

end