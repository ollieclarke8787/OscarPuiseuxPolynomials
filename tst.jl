using Revise
using Oscar
using OscarPuiseuxPolynomial


Kt, t = puiseux_polynomial_ring(QQ, ["t"])
Ktx, x = polynomial_ring(Kt, ["x"])




length(x)
x^QQ(1//2)
S = base_ring(R)
t = first(gens(S))

f = FiniteSparsePuiseuxSeries.FiniteSparsePuiseuxSeriesRingElem(R, t^2, ZZ(0), ZZ(2))
f+one(R) # Danger!!

normalize!(f)
f
x
normalize!(x)
x

hash(f)
hash(x)
hash(one(R))
hash(zero(R))

one(R)==one(R)
iszero(zero(R))


one(R)+one(R)
one(R)+zero(R)
zero(R)+zero(R)
one(R)+x

one(R)-one(R)
x - one(R)
x - zero(R)

isone(one(R))
isone(zero(R))
isone(x)

x*x
(x^2-one(R))*(x^2+one(R))

f = x^3 + x^2 - x
normalize!(f)
f

f = FiniteSparsePuiseuxSeries.FiniteSparsePuiseuxSeriesRingElem(R, t^3+t+1, ZZ(0), ZZ(2))
