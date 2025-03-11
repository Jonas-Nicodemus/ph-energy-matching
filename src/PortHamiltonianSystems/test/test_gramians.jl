module Test_gramians

using Test

using PortHamiltonianSystems

using ControlSystemsBase

@testset "test_gramians.jl" begin

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
        
    @testset "grampd" begin
        @testset ":c" begin
            opt = :c
            U = grampd(Σph, opt)
            X = U * U'
            @test norm(Σ.A * X + X * Σ.A' + Σ.B * Σ.B') < 1e-8
        end

        @testset ":o" begin
            opt = :o
            U = grampd(Σph, opt)
            X = U' * U
            @test norm(Σ.A' * X + X * Σ.A + Σ.C' * Σ.C) < 1e-8
        end

        @testset ":pr_c" begin
            opt = :pr_c
            U = grampd(Σph, opt)
            X = U' * U 
            @test norm(PortHamiltonianSystems.prare(Σph, :c, X)) < 1e-8
        end

        @testset ":pr_o" begin
            opt = :pr_o
            U = grampd(Σph, opt)
            X = U' * U 
            @test norm(PortHamiltonianSystems.prare(Σph, :o, X)) < 1e-8
        end
                    
    end

    @testset "gram" begin
        @testset ":c" begin
            opt = :c
            X = gram(Σph, opt)
            @test norm(Σ.A * X + X * Σ.A' + Σ.B * Σ.B') < 1e-8
        end

        @testset ":o" begin
            opt = :o
            X = gram(Σph, opt)
            @test norm(Σ.A' * X + X * Σ.A + Σ.C' * Σ.C) < 1e-8
        end

        @testset ":pr_c" begin
            opt = :pr_c
            X = gram(Σph, opt)
            @test norm(PortHamiltonianSystems.prare(Σph, :c, X)) < 1e-8
        end

        @testset ":pr_o" begin
            opt = :pr_o
            X = gram(Σph, opt) 
            @test norm(PortHamiltonianSystems.prare(Σph, :o, X)) < 1e-8
        end

    end

end

end