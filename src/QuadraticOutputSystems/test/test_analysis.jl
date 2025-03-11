module Test_analysis

using QuadraticOutputSystems
using ControlSystems
using Test

@testset "test_analysis.jl" begin

    @testset "without C" begin
        A = [-2. 1.; -1. -1.]
        B = [6.; 0;;]
        Q = [1. 0.; 0. 1.]
        M = vec(1/2*Q)'

        Σ = qoss(A, B, M)
        Σr = qoss(A[1:1,1:1], B[1:1,1:1], vec(1/2*Q[1:1,1:1])')  
            
        @testset "norm" begin
            @testset "h2norm" begin
                for opt in (:c, :o)
                    @test h2norm(Σ, opt)^2 ≈ 19
                end
            end

            @testset "h2inner" begin
                for opt in (:c, :o)
                    @test h2inner(Σ, Σr, opt) ≈ 19.17159763313609
                end
            end
            
            @testset "h2error" begin
                for opt in (:c, :o)
                    @test h2norm(Σ - Σr, opt)^2 ≈ h2norm(Σ, opt)^2 + h2norm(Σr, opt)^2 - 2*h2inner(Σ, Σr, opt)
                end
            end
        end
    end

    @testset "with C" begin
        A = [-2. 1.; -1. -1.]
        B = [6.; 0;;]
        Q = [1. 0.; 0. 1.]
        M = vec(1/2*Q)'

        Σ = qoss(A, B, B', M)
        Σr = qoss(A[1:1,1:1], B[1:1,1:1], B'[1:1,1:1], vec(1/2*Q[1:1,1:1])')  
            
        @testset "norm" begin
            @testset "h2norm FOM" begin
                for opt in (:c, :o)
                    @test h2norm(Σ, opt)^2 ≈ 307
                end
            end

            @testset "h2norm ROM" begin
                for opt in (:c, :o)
                    @test h2norm(Σr, opt)^2 ≈ 344.25
                end
            end

            @testset "h2inner" begin
                for opt in (:c, :o)
                    @test h2inner(Σ, Σr, opt) ≈ 318.2485207100592
                end
            end
            
            @testset "h2error" begin
                for opt in (:c, :o)
                    @test h2norm(Σ - Σr, opt)^2 ≈ h2norm(Σ, opt)^2 + h2norm(Σr, opt)^2 - 2*h2inner(Σ, Σr, opt)
                end
            end
        end
    end
end

end