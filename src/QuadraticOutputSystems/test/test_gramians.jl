module Test_gramians

using QuadraticOutputSystems
using ControlSystems
using Test

@testset "test_gramians.jl" begin
    A = [-2. 1.; -1. -1.]
    B = [6.; 0;;]
    M = [1/2 0.; 0. 1/2]
    
    Σqo = qoss(A, B, M)
        
    @testset "grampd" begin
        @testset ":c" begin
            opt = :c
            U = grampd(Σqo, opt)
            X = U * U'
            @test norm(Σqo.A * X + X * Σqo.A' + Σqo.B * Σqo.B') < 1e-8
        end

        @testset ":o" begin
            opt = :o
            U = grampd(Σqo, opt)
            X = U * U'
            
            UP = grampd(Σqo, :c)
            P = UP * UP'
            @test norm(Σqo.A' * X + X * Σqo.A + Σqo.M * P * Σqo.M) < 1e-8
        end       
    end

    @testset "gram" begin
        @testset ":c" begin
            opt = :c
            X = gram(Σqo, opt)
            @test norm(Σqo.A * X + X * Σqo.A' + Σqo.B * Σqo.B') < 1e-8
        end

        @testset ":o" begin
            opt = :o
            X = gram(Σqo, opt)
            
            P = gram(Σqo, :c)
            @test norm(Σqo.A' * X + X * Σqo.A + Σqo.M * P * Σqo.M) < 1e-8
        end
    end

end

end