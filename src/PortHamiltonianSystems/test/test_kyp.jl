module Test_kyp

using Test

using PortHamiltonianSystems

using LinearAlgebra, ControlSystemsBase

@testset "test_kyp.jl" begin

    J = [0. 1.; -1. 0.]
    R = [2. 0.; 0. 1.]
    Q = [1. 0.; 0. 1.]
    G = [6.; 0.;;]
    P = zero(G)
    S = [1.;;]
    N = zero(S)

    Σph = phss(J, R, Q, G, P, S, N)
    Σ = ss(Σph)
        
    @testset "KYP are" begin
        for Σᵢ in [Σph, Σ]
            Xmin = kypmin(Σᵢ)
            @test norm(kypare(Σ, Xmin)) < 1e-8

            Xmax = kypmax(Σᵢ)
            @test norm(kypare(Σ, Xmax)) < 1e-8
        end
    end

    @testset "KYP lmi" begin
        X = kyp(Σ)

        @test ispsd(X) 
        @test ispsd(kypmat(Σ, X))
    end
end

end