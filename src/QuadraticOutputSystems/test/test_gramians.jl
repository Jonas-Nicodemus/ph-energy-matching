module Test_gramians

using QuadraticOutputSystems
using ControlSystems
using Test

@testset "test_gramians.jl" begin
    A = [-2. 1.; -1. -1.]
    B = [6.; 0;;]
    Q = [1. 0.; 0. 1.]
    M = vec(1/2*Q)'
    n = 2
    m = 1

    Σqo1 = qoss(A, B, M)
    Σqo2 = qoss(A, B, B', M)
    Σqo3 = qoss(A, B, [B'; zeros(1, n)], [zeros(m, n^2); M])
    
    for Σqo in [Σqo1, Σqo2, Σqo3]

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
                X = U' * U
                
                UP = grampd(Σqo, :c)
                P = UP * UP'
                @test norm(Σqo.A' * X + X * Σqo.A + Σqo.C'*Σqo.C + 1/4 * Q * P * Q) < 1e-8
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
                @test  norm(Σqo.A' * X + X * Σqo.A + Σqo.C'*Σqo.C + 1/4 * Q * P * Q) < 1e-8
            end
        end
    end
end

end