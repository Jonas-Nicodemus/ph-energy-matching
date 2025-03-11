module Test_prbt

using Test

using PortHamiltonianModelReduction

using LinearAlgebra, ControlSystemsBase
using PortHamiltonianSystems

import PortHamiltonianModelReduction: bt

@testset "test_prbt.jl" begin
    J = [0. 1.; -1. 0.]
    R = [2. 0.; 0. 1.]
    Q = [1. 0.; 0. 1.]
    G = [6.; 0.;;]
    P = zero(G)
    S = [1.;;]
    N = zero(S)

    Σph = phss(J, R, Q, G, P, S, N)
    Σ = ss(Σph)

    @testset "bt" begin
        @testset "r=2" begin
            r = 2
            Lx = grampd(Σ, :o)
            Ly = grampd(Σ, :c)'
            Σr = bt(Σ, r; Lx=Lx, Ly=Ly)
            P = gram(Σr, :c)
            O = gram(Σr, :o)
            @test norm(P - diagm(diag(P))) < 1e-10
            @test norm(O - diagm(diag(O))) < 1e-10
            @test norm(Σ - Σr) < 1e-10
        end

        @testset "r=1" begin
            r = 1
            Lx = grampd(Σ, :o)
            Ly = grampd(Σ, :c)'
            Σr = bt(Σ, r; Lx=Lx, Ly=Ly)
            P = gram(Σr, :c)
            O = gram(Σr, :o)
            @test norm(P - diagm(diag(P))) < 1e-10
            @test norm(O - diagm(diag(O))) < 1e-10
            @test norm(Σ - Σr) < 1e1
        end
    end

    @testset "prbt" begin
        r = 2
        Σr = prbt(Σ, r)
        P = prgram(Σr, :c)
        O = prgram(Σr, :o)

        @test norm(P - diagm(diag(P))) < 1e-10
        @test norm(O - diagm(diag(O))) < 1e-10
        @test norm(P - O) < 1e-10

        r = 1
        Σphr = prbt(Σph, r)
        Pr = prgram(Σphr, :c)
        Or = prgram(Σphr, :o)

        @test norm(Pr - diagm(diag(Pr))) < 1e-10
        @test norm(Or - diagm(diag(Or))) < 1e-10
        @test norm(Pr - Or) < 1e-10
    end

end

end
