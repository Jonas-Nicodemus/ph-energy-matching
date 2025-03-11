module Test_barrier

using EnergyMatching
using LinearAlgebra, ControlSystemsBase
using QuadraticOutputSystems, PortHamiltonianSystems
using FiniteDifferences
using Test

@testset "test_barrier.jl" begin
    A = [-2. 1.; -1. -1.]
    B = [6.; 0;;]
    C = B'
    D = [1;;]

    M = 1//2*vec([1 0; 0 1])'
    
    Σ = qoss(A, B, M)
    Σr = ss(A[1:1,1:1], B[1:1,1:1], C[1:1,1:1], D)

    Xmin = kypmin(Σr)
    Xmax = kypmax(Σr)
    X0 = 1/2 * (Xmin + Xmax)
    M0 = 1/2 * vec(X0)'
    x = 1/2 * vech(X0) 

    @testset "objective" begin
        f, g = EnergyMatching.objective(Σ, Σr)

        @testset "f" begin
            @test f(x) ≈ norm(Σ - qoss(Σr.A, Σr.B, vec(M0)'))^2
        end

        @testset "g" begin
            @test norm(grad(central_fdm(5,1), f, x)[1] - g(x)) < 1e-10
        end

        # @testset "h" begin
        #     @test norm(jacobian(central_fdm(5, 1), g, x)[1] - h(x)) < 1e-10
        # end
    end

    @testset "constraint" begin
        f, g = EnergyMatching.constraint(Σr)

        @testset "g" begin
            @test norm(grad(central_fdm(5,1), f, x)[1] - g(x)) < 1e-10
        end

        # @testset "h" begin
        #     @test norm(jacobian(central_fdm(5, 1), g, x)[1] - H(x)) < 1e-12
        # end
    end

    @testset "combined" begin
        fo, go = EnergyMatching.objective(Σ, Σr)
        fc, gc = EnergyMatching.constraint(Σr)

        α = 1e-3
        f = (x) -> fo(x) + α * fc(x)
        g!(g, x) = begin
            g .= go(x) + α * gc(x)
        end

        @testset "g" begin
            g = zeros(length(x))
            g!(g, x)
            @test norm(grad(central_fdm(5,1), f, x)[1] - g) < 1e-10
        end
    end

    @testset "barrier" begin
        Σphr = EnergyMatching.barrier(Σ, Σr)
        @test Σphr.Q ≈ [160/169;;]
        @test norm(Σ - hdss(Σphr))^2 ≈ 19 + 81/4 * (160/169)^2 - 2 * 3240/169 * (160/169) 

        Σr = ss(A, B, C, D)
        Σr0 = phss(Σr)
        Σphr = EnergyMatching.barrier(Σ, Σr; Σr0=Σr0)
        @test norm(Σphr.Q - 2*unvec(M)) < 1e-6
    end
end

end