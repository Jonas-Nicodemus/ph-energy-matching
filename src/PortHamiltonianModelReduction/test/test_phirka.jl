module TestPHIRKA

using Test

using PortHamiltonianModelReduction

using LinearAlgebra, ControlSystems
using PortHamiltonianSystems

@testset "test_phirka.jl" begin
    J = [0. 1.; -1. 0.]
    R = [2. 0.; 0. 1.]
    Q = [1. 0.; 0. 1.]
    G = [6.; 0.;;]
    P = zero(G)
    S = [1.;;]
    N = zero(S)

    Σ = phss(J, R, Q, G, P, S, N)
    Σr = phirka(Σ, 1)
    @test norm(Σ - Σr) < 1e1

    Σr = phirka(Σ, 1, 3)
    @test norm(Σ - Σr) < 1e1
end
   
end