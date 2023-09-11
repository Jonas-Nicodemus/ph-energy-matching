module Test_kyp

using Test

using PortHamiltonianSystems

using ControlSystems

@testset "test_kyp.jl" begin

    J = [0. 1.; -1. 0.]
    R = [2. 0.; 0. 1.]
    Q = [1. 0.; 0. 1.]
    G = [6.; 0.;;]
    P = zero(G)
    S = [1.;;]
    N = zero(S)

    A = [-2. 1.; -1. -1.]
    B = [6.; 0;]
    C = B'
    D = [1.;;]

    Σph = phss(J, R, Q, G, P, S, N)
    Σ = ss(Σph)
        
    @testset "kyp" begin
        for min in [true, false]
            for Σᵢ in [Σph, Σ]
                X = kyp(Σᵢ; min=min)
                @test norm(kypare(Σ, X)) < 1e-8
            end
        end
    end

end

end