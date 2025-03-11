module Test_convert

using Test

using PortHamiltonianSystems

using ControlSystemsBase

@testset "test_convert.jl" begin

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

    @testset "compose" begin
        A, B, C, D = compose(J, R, Q, G, P, S, N)
        @test A == [-2. 1.; -1. -1.]
        @test B == [6.; 0.;;]
        @test C == [6. 0.;]
        @test D == [1.;;]
    end
    
    @testset "ss" begin
        Σ = ss(J, R, Q, G, P, S, N)
        @test Σ.A == [-2. 1.; -1. -1.]
        @test Σ.B == [6.; 0.;;]
        @test Σ.C == [6. 0.;]
        @test Σ.D == [1.;;]

        Σ = ss(Σph)
        @test Σ.A == [-2. 1.; -1. -1.]
        @test Σ.B == [6.; 0.;;]
        @test Σ.C == [6. 0.;]
        @test Σ.D == [1.;;]
    end

    @testset "decompose" begin
        X = [1. 0.; 0. 1.]
        J, R, Q, G, P, S, N = decompose(A, B, C, D, X)
        @test J == [0. 1.; -1. 0.]
        @test R == [2. 0.; 0. 1.]
        @test Q == [1. 0.; 0. 1.]
        @test G == [6.; 0.;;]
        @test P == [0.; 0.;;]
        @test S == [1.;;]
        @test N == [0.;;]
    end

    @testset "phss" begin
        Σ = ss(Σph)
        Σph_min = phss(Σ)
        @test norm(Σ - ss(Σph_min)) < 1e-8
    end
end

end
