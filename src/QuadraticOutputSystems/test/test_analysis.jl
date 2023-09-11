module Test_analysis

using QuadraticOutputSystems
using ControlSystems
using Test

@testset "test_analysis.jl" begin
    A = [-2. 1.; -1. -1.]
    B = [6.; 0;;]
    M = 1/2*[1. 0.; 0. 1.]
    
    Σ = qoss(A, B, M)
    Σr = qoss(A[1:1,1:1], B[1:1,1:1], M[1:1,1:1])  
        
    @testset "norm" begin
        @testset "h2norm" begin
            @test h2norm(Σ)^2 ≈ 19
        end

        @testset "h2inner" begin
            @test h2inner(Σ, Σr) ≈ 19.17159763313609
        end
        
        @test norm(Σ - Σr)^2 ≈ norm(Σ)^2 + norm(Σr)^2 - 2*h2inner(Σ, Σr)
    end
end

end