#################################################################################
#
#  Puiseux polynomials
#
#################################################################################


struct PuiseuxPolynomialRing <: Ring
    polynomialRing::MPolyRing

    function PuiseuxPolynomialRing(polynomialRing::MPolyRing)
        @assert ngens(polynomialRing) == 1 "underlying ring must be a multivariate polynomial ring with one variable"
        return new(polynomialRing)
    end
end

function Base.show(io::IO, R::PuiseuxPolynomialRing)
    print(io, "Puiseux polynomial ring over ", coefficient_field(R))
end

base_ring(R::PuiseuxPolynomialRing) = R.polynomialRing
coefficient_field(R::PuiseuxPolynomialRing) = base_ring(base_ring(R))


#################################################################################
#
# A Puiseux polynomial can be written of the form
#
#   t^{k//d} * (c_0 * t^(0//d) + ... + c_r * t^(r//d)).
#
# We store it as
# - a (sparse) poly: c_0 * t^0 + ... + c_r * t^r
# - a shift: k
# - a scale: d
#
# All outputs are normalized and all inputs (excepd for the `normalized!`
# function) are assumed to be normalized in the following sense, which makes
# the representation unique:
# - gcd(d, exponents of poly) = 1
# - c_0 != 0
#
###
mutable struct PuiseuxPolynomialRingElem <: RingElem
    parent::PuiseuxPolynomialRing
    poly::MPolyRingElem
    shift::ZZRingElem
    scale::ZZRingElem

    function PuiseuxPolynomialRingElem(Kt::PuiseuxPolynomialRing, f::MPolyRingElem, k::ZZRingElem = zero(ZZ), d::ZZRingElem = one(ZZ))
        @assert parent(f) == base_ring(Kt) "polynomial must be in the base ring"
        @assert d > 0 "scale must be positive"
        return new(Kt, f, k, d)
    end
end

parent(f::PuiseuxPolynomialRingElem) = f.parent
poly(f::PuiseuxPolynomialRingElem) = f.poly
shift(f::PuiseuxPolynomialRingElem) = f.shift
scale(f::PuiseuxPolynomialRingElem) = f.scale

iszero(f::PuiseuxPolynomialRingElem) = iszero(poly(f))
isone(f::PuiseuxPolynomialRingElem) = isone(poly(f)) && shift(f) == 0 && scale(f) == 1

function Base.show(io::IO, f::PuiseuxPolynomialRingElem)
    if iszero(f)
        print(io, "0")
        return
    elseif isone(f)
        print(io, "1")
        return
    end

    t = first(gens(base_ring(f)))
    terms = []
    for (c, e) in zip(coefficients(poly(f)), first.(exponents(poly(f))))
        push!(terms, "$(c)*$(t)^($( (e + shift(f)) // scale(f) ))")
    end
    print(io, join(terms, " + "))
end


function normalize!(f::PuiseuxPolynomialRingElem)
    if iszero(f)
        return false
    end

    baseRing = base_ring(f)

    # make sure scale is correct, i.e., gcd of exponents is 1
    polyf = poly(f)
    exponentsf = first.(exponents(polyf))
    gcdExponents = gcd(vcat(exponentsf, scale(f)))
    if gcdExponents > 1
        f.poly = sum([c*monomial(baseRing, [div(e, gcdExponents)]) for (c, e) in zip(coefficients(polyf), exponentsf)])
        f.shift = div(f.shift, gcdExponents)
        f.scale = div(f.scale, gcdExponents)
    end

    # make sure shift is correct, i.e., poly(f) has a constant term
    polyf = poly(f)
    exponentsf = first.(exponents(polyf))
    smallestExp = exponentsf[end]
    if smallestExp > 0
        f.poly = sum([c*monomial(baseRing, [e - smallestExp]) for (c, e) in zip(coefficients(polyf), exponentsf)])
        f.shift += smallestExp
    end
    return gcdExponents > 1 || smallestExp > 0
end

function ==(f::PuiseuxPolynomialRingElem, g::PuiseuxPolynomialRingElem)
    @assert parent(f) == parent(g) "elements must be in the same ring"

    # bring f and g in standard form
    normalize!(f)
    normalize!(g)

    return poly(f) == poly(g) && shift(f) == shift(g) && scale(f) == scale(g)
end

function hash(f::PuiseuxPolynomialRingElem, h::UInt)
    normalize!(f)
    return hash((parent(f), poly(f), shift(f), scale(f)), h)
end

gen(R::PuiseuxPolynomialRing) = PuiseuxPolynomialRingElem(R, first(gens(base_ring(R))))
zero(R::PuiseuxPolynomialRing) = PuiseuxPolynomialRingElem(R, zero(base_ring(R)))
one(R::PuiseuxPolynomialRing) = PuiseuxPolynomialRingElem(R, one(base_ring(R)))

coefficient_field(f::PuiseuxPolynomialRingElem) = coefficient_field(parent(f))
coefficients(f::PuiseuxPolynomialRingElem) = coefficients(poly(f))
function exponents(f::PuiseuxPolynomialRingElem)
    d = scale(f)
    k = shift(f)
    underlying_exponents = first.(exponents(poly(f)))
    return [(i + k) // d for i in underlying_exponents]
end

function puiseux_polynomial_ring(K::Field, variableName::String="t")
    base_ring, _ = Oscar.polynomial_ring(K, [variableName])
    Kt = PuiseuxPolynomialRing(base_ring)
    t = gen(Kt)
    return Kt, t
end


function rescale(f::PuiseuxPolynomialRingElem, newScale::ZZRingElem)
    # we assume f is normalized
    if newScale == scale(f)
        return f
    end
    @assert newScale > 0 "new scale must be positive"

    newScaleMultipleOfCurrentScale, scaleQuotient = divides(newScale, scale(f))
    @assert newScaleMultipleOfCurrentScale "new scale must be a multiple of the current scale"

    t = first(gens(base_ring(f)))
    newPoly = evaluate(poly(f), [t^scaleQuotient])
    newShift = shift(f) * scaleQuotient
    return PuiseuxPolynomialRingElem(parent(f), newPoly, newShift, newScale)
end


function +(f::PuiseuxPolynomialRingElem, g::PuiseuxPolynomialRingElem)
    @assert parent(f) == parent(g) "elements must be in the same ring"
    if iszero(f)
        return g
    elseif iszero(g)
        return f
    end

    newScale = lcm(scale(f), scale(g))

    frescaled = rescale(f, newScale)
    grescaled = rescale(g, newScale)

    newShift = min(shift(frescaled), shift(grescaled))
    t = first(gens(base_ring(f)))
    newPoly = t^(shift(frescaled)-newShift)*poly(frescaled) + t^(shift(grescaled)-newShift)*poly(grescaled)

    # normalize output
    fplusg = PuiseuxPolynomialRingElem(
        parent(f),
        newPoly,
        newShift,
        newScale)
    normalize!(fplusg)
    return fplusg
end

function -(f::PuiseuxPolynomialRingElem)
    return PuiseuxPolynomialRingElem(parent(f), -poly(f), shift(f), scale(f))
end

function -(f::PuiseuxPolynomialRingElem, g::PuiseuxPolynomialRingElem)
    return f + (-g)
end

function *(f::PuiseuxPolynomialRingElem, g::PuiseuxPolynomialRingElem)
    @assert parent(f) == parent(g) "elements must be in the same ring"
    if iszero(f) || iszero(g)
        return zero(parent(f))
    elseif isone(f)
        return g
    elseif isone(g)
        return f
    end

    newShift = shift(f) + shift(g)
    newScale = scale(f)*scale(g)
    newPoly = poly(f)*poly(g)

    # normalize output
    ftimesg = PuiseuxPolynomialRingElem(
        parent(f),
        newPoly,
        newShift,
        newScale)

    normalize!(ftimesg)
    return ftimesg
end

function length(f::PuiseuxPolynomialRingElem)
    return length(poly(f))
end

function ^(f::PuiseuxPolynomialRingElem, a::QQFieldElem)

    if denominator(a) == 1
        return f^numerator(a)
    end

    @assert length(f) == 1 "only single-term series can be exponentiated to rational powers"

    return PuiseuxPolynomialRingElem(
        parent(f),
        poly(f),
        shift(f)*numerator(a),
        scale(f)*denominator(a))
end

function ^(f::PuiseuxPolynomialRingElem, a::Rational{Int})
    return f^(QQ(a))
end
