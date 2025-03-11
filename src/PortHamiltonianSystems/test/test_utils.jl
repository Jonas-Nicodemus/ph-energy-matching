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
        eigtol = 1e-12
        M = rand(10, 10)
        Mspsd = project_psd(M; eigtol=eigtol)
        @test norm(Mspsd - Mspsd') == 0
        @test all(eigvals(Mspsd) .>= -eigtol)
    end

    @testset "vech" begin
        S = [1 2; 2 3]
        s = vech(S)
        @test s == [1, 2, 3]

        S = [1 2 3; 2 4 5; 3 5 6]
        s = vech(S)
        @test s == [1, 2, 3, 4, 5, 6]
    end

    @testset "unvech" begin
        s = [1, 2, 3]
        S = unvech(s)
        @test S == [1 2; 2 3]

        s = [1, 2, 3, 4, 5, 6]
        S = unvech(s)
        @test S == [1 2 3; 2 4 5; 3 5 6]
    end

end

end