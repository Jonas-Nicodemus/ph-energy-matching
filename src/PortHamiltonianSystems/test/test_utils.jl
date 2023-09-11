module Test_utils

using Test

using PortHamiltonianSystems

using LinearAlgebra

@testset "test_utils.jl" begin

    @testset "sym" begin
        M = rand(10, 10)
        Msym = sym(M)
        @test norm(Msym - Msym') == 0
    end

    @testset "skew" begin
        M = rand(10, 10)
        Mskew = skew(M)
        @test norm(Mskew + Mskew') == 0
    end

    @testset "project_psd" begin
        eig_tol = 1e-12
        M = rand(10, 10)
        Mspsd = project_psd(M, eig_tol)
        @test norm(Mspsd - Mspsd') == 0
        @test all(eigvals(Mspsd) .>= -eig_tol)
    end

end

end