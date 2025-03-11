module Test_minreal

using EnergyMatching
using LinearAlgebra, ControlSystemsBase
using PortHamiltonianSystems
using Test

@testset "minreal.jl" begin
    U = UpperTriangular([1.0 2.0 3.0; 0.0 4.0 5.0; 0.0 0.0 6.0])

    J = U' - U 
    R = U'*U

    G = [1; 0; 0]
    P = zero(G)
    S = [1.;;]
    N = zero(S)
    
    @testset "truncno" begin
        Q = U[1:2,:]'*U[1:2,:]
        Σph = phss(J, R, Q, G, P, S, N)
        Σphr = EnergyMatching.truncno(Σph)
        @test ss(Σphr).nx == 2
        @test norm(Σph - Σphr) < 1e-6
    end
        
    @testset "truncnc" begin
        Q = U'*U
        Σph = phss(J, R, Q, G, P, S, N)
        Σphr = EnergyMatching.truncnc(Σph; ϵ=1e-4)
        @test ss(Σphr).nx == 2
        @test norm(Σph - Σphr) < 1e-2
    end

    @testset "minreal" begin
        Q = U'*U
        Σph = phss(J, R, Q, G, P, S, N)
        Σphr = EnergyMatching.minreal(Σph)
        @test norm(Σph - Σphr) < 1e-2
    end
end

end