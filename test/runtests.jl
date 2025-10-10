using OscarPuiseuxPolynomial
using Test

@testset "OscarPuiseuxPolynomial.jl" begin
    @testset "TrivialTests" begin
        @test 1+1==2
        @test 1+2==3
        # @test 1+1==3 # this one fails
    end
end
