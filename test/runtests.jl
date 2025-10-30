using OscarPuiseuxPolynomial
using Test
using Oscar

@testset "OscarPuiseuxPolynomial.jl" begin
    @testset "TrivialTests" begin
        @test 1+1==2
        @test 1+2==3
        # @test 1+1==3 # this one fails
    end

    @testset "Construction" begin
        K, (t1,t2,t3) = Oscar.polynomial_ring(QQ, ["t1","t2","t3"])
        K_p,(tp1,tp2,tp3) = puiseux_polynomial_ring(QQ, ["t1","t2","t3"])
        @test K_p.underlyingPolynomialRing == K
        @test K_p == OscarPuiseuxPolynomial.MPuiseuxPolyRing(K)

        h = 1+t1 + 2*t2+3*t1^4+t1*t2^4+t3^2
        g = OscarPuiseuxPolynomial.MPuiseuxPolyRingElem(K_p,h)
        @test h != g
        @test g.scale == 1
        @test g.shift == [0,0,0]

        g = OscarPuiseuxPolynomial.MPuiseuxPolyRingElem(K_p,h,[ZZ(1),ZZ(1),ZZ(1)],ZZ(3))
        @test g.scale==3
        @test g.shift==[ZZ(1),ZZ(1),ZZ(1)]

        h = t1^(2)
        g = puiseux_polynomial_ring_elem(K_p,h)
        @test g.scale == 1
        @test g.shift == [2,0,0]
        @test normalize!(g) == false

        h = t1^2*(1 + t1)
        g = puiseux_polynomial_ring_elem(K_p,h)
        @test g.scale == 1
        @test g.shift == [2,0,0]

        h = t1*t2^(2)*t3^(3)*(1+t1+t2+t3)
        g = puiseux_polynomial_ring_elem(K_p,h)
        @test g.scale == 1
        @test g.shift == [1,2,3]

        h = t1*t2^(2)*t3^(3)*(1+t1+t2+t3)
        g = puiseux_polynomial_ring_elem(K_p,h,skip_normalization=true)
        @test g.scale == 1
        @test g.shift == [0,0,0]
        @test normalize!(g) == true

        h = (1+t1+t2+t3)
        g = puiseux_polynomial_ring_elem(K_p,h,[ZZ(1),ZZ(2),ZZ(3)],ZZ(3),skip_normalization=true)
        @test g.scale == 3
        @test g.shift == [1,2,3]
        
        K, _ = polynomial_ring(QQ, ["t1","t2","t3"])
        Kt, _ = puiseux_polynomial_ring(QQ,["t1","t2","t3"])
        @test Kt.underlyingPolynomialRing == K
    end
    
    @testset "Getters" begin
        K, (t1,t2,t3) = polynomial_ring(QQ, ["t1","t2","t3"])
        Kp, (tp1,tp2,tp3) = puiseux_polynomial_ring(QQ,["t1","t2","t3"])

	    @test K == OscarPuiseuxPolynomial.underlying_polynomial_ring(Kp)
        @test QQ == OscarPuiseuxPolynomial.base_ring(Kp)
        @test QQ == OscarPuiseuxPolynomial.coefficient_ring(Kp)

        g = tp1^(1//2)+tp3^(1//3)
        @test OscarPuiseuxPolynomial.parent(g) == Kp
        @test OscarPuiseuxPolynomial.poly(g) == t1^3 + t3^2
        @test OscarPuiseuxPolynomial.scale(g) == 6
        @test OscarPuiseuxPolynomial.shift(g) == [0,0,0]


        
        g = tp1^(2//3)*tp1*tp2^(1//2)*tp3 + tp3^(3//7)*tp1*tp2^(1//2)*tp3 + tp2^(1//2)*tp1*tp2^(1//2)*tp3
        @test OscarPuiseuxPolynomial.parent(g) == Kp
        @test OscarPuiseuxPolynomial.poly(g) == t1^28 + t2^21 + t3^18
        @test OscarPuiseuxPolynomial.shift(g) == [42,21,42]
        @test OscarPuiseuxPolynomial.scale(g) == 2*3*7
    end

    @testset "Arithmetic" begin
        K, (u,v,w) = puiseux_polynomial_ring(QQ,["u","v","w"])
        g = u + v
        h = v + w
        @test 4*g == 4*u + 4*v
        @test h+g == u + 2*v + w
        @test h-g == w-u
        @test h*g == u*v + u*w + v^2 + v*w
        @test (g)^3 == u^3 + 3*u^2*v + 3*u*v^2 + v^3
        
        g = u^(1//2)*v^(2//3) + w^(1//4)
        h = u^(2//3)
        @test g*h == u^(7//6)*v^(2//3)+w^(1//4)*u^(2//3)

        @test (u^(1//2)+u^(-1//2))^2 == u+2+u^(-1)
   end

end


