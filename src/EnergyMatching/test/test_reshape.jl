module Test_reshape

using EnergyMatching
using Test

@testset "reshape.jl" begin

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