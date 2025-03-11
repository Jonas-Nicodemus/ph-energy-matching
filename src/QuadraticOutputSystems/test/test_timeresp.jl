module Test_timeresp

using QuadraticOutputSystems
using ControlSystems
using Test

@testset "test_timeresp.jl" begin
    A = [-2. 1.; -1. -1.]
    B = [6.; 0;;]
    C = B'
    M = 1/2*[1. 0.; 0. 1.]
    
    Σ = qoss(A, B, C, vec(M)')
        
    @testset "lsim" begin
        t = 0:0.1:10
        u(x,t) = [sin(t)]
        res = lsim(Σ, u, t) 
        @test isa(res.sys, QuadraticOutputStateSpace)
        @test size(res.y) == (1, length(t))
        @test size(res.x) == (2, length(t))
        @test size(res.u) == (1, length(t)) 
    end
end

end