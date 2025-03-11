using Test

using LinearAlgebra
using PortHamiltonianSystems

@testset "test_convert.jl" begin

    J = [0. 1.; -1. 0.]
    R = [2. 0.; 0. 1.]
    Q = [1. 0.; 0. 1.]
    G = [6.; 0.;;]
    P = zero(G)
    S = [1.;;]
    N = zero(S)

    Γ = [J G; G' N]
    W = [R P; P' S]

    @testset "phss" begin
        Σ1 = phss(J, R, Q, G, P, S, N)
        Σ2 = phss(J, R, Q, G)
        Σ3 = phss(Γ, W, Q)

        S = 1
        Σ4 = phss(J, R, Q, G, P, S, N)
        @test norm(Σ1 - Σ3) < 1e-12
    end
end