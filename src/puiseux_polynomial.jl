#################################################################################
#
# Puiseux polynomials
# ===================
#
# !!! Currently only supports univariate Puiseux polynomials
#
# A Puiseux multivariate polynomial ring is just a wrapper around a normal
# multivariate polynomial ring for the normal poly below.
#
# A Puiseux polynomial can be written of the form
#
#   t^{k//d} * (c_0 * t^(0//d) + ... + c_r * t^(r//d)).
#
# We store it as
# - a normal poly: c_0 * t^0 + ... + c_r * t^r
# - a shift: k
# - a scale: d
#
# The representation above is normalized if
# - gcd(d, exponents of poly) = 1
# - c_0 != 0
# Every Puiseux polynomial has a unique normalized representation.
#
# With the exception of normalized! and rescale, all function inputs are assumed
# to be normalized and all function outputs will be normalized.
#
#################################################################################

struct PuiseuxPolynomialRing{T <: FieldElement} <: Ring
    underlyingPolynomialRing::MPolyRing

    function PuiseuxPolynomialRing(R::MPolyRing)
        @assert ngens(R) == 1 "underlying ring must be a multivariate polynomial ring with one variable"
        return new{elem_type(base_ring_type(R))}(R)
    end
end

mutable struct PuiseuxPolynomialRingElem{T <: FieldElement} <: RingElem
    parent::PuiseuxPolynomialRing{T}
    poly::MPolyRingElem
    shift::ZZRingElem
    scale::ZZRingElem

    function PuiseuxPolynomialRingElem(Kt::PuiseuxPolynomialRing, f::MPolyRingElem, k::ZZRingElem = zero(ZZ), d::ZZRingElem = one(ZZ))
        @assert parent(f) == underlying_polynomial_ring(Kt) "polynomial must be in the underlying polynomial ring"
        @assert d > 0 "scale must be positive"
        return new{elem_type(base_ring_type(parent(f)))}(Kt, f, k, d)
    end
end


#################################################################################
#
# Constructors
#
##################################################################################

function puiseux_polynomial_ring(K::Field, variableName::Vector{String}=["t"])
    @assert length(variableName)==1 "multivariate puiseux polynomials not supported"
    base_ring, _ = polynomial_ring(K, variableName)
    Kt = PuiseuxPolynomialRing(base_ring)
    t = gen(Kt)
    return Kt, t
end


#################################################################################
#
# Getters
#
#################################################################################

underlying_polynomial_ring(R::PuiseuxPolynomialRing) = R.underlyingPolynomialRing
base_ring(R::PuiseuxPolynomialRing) = base_ring(underlying_polynomial_ring(R))
coefficient_field(R::PuiseuxPolynomialRing) = base_ring(R)

parent(f::PuiseuxPolynomialRingElem) = f.parent
poly(f::PuiseuxPolynomialRingElem) = f.poly
shift(f::PuiseuxPolynomialRingElem) = f.shift
scale(f::PuiseuxPolynomialRingElem) = f.scale


#################################################################################
#
# Setters
#
#################################################################################

# WARNING: input may not be normalized
function normalize!(f::PuiseuxPolynomialRingElem)
    if iszero(f)
        return false
    end

    underlyingPolynomialRing = underlying_polynomial_ring(f)

    # make sure scale is correct, i.e., gcd of exponents is 1
    polyf = poly(f)
    exponentsf = first.(exponents(polyf))
    gcdExponents = gcd(vcat(exponentsf, scale(f)))
    if gcdExponents > 1
        f.poly = sum([c*monomial(underlyingPolynomialRing, [div(e, gcdExponents)]) for (c, e) in zip(coefficients(polyf), exponentsf)])
        f.shift = div(f.shift, gcdExponents)
        f.scale = div(f.scale, gcdExponents)
    end

    # make sure shift is correct, i.e., poly(f) has a constant term
    polyf = poly(f)
    exponentsf = first.(exponents(polyf))
    smallestExp = exponentsf[end]
    if smallestExp > 0
        f.poly = sum([c*monomial(underlyingPolynomialRing, [e - smallestExp]) for (c, e) in zip(coefficients(polyf), exponentsf)])
        f.shift += smallestExp
    end
    return gcdExponents > 1 || smallestExp > 0
end

# WARNING: output may not be normalized
function rescale(f::PuiseuxPolynomialRingElem, newScale::ZZRingElem)
    # we assume f is normalized
    if newScale == scale(f)
        return f
    end
    @assert newScale > 0 "new scale must be positive"

    newScaleMultipleOfCurrentScale, scaleQuotient = divides(newScale, scale(f))
    @assert newScaleMultipleOfCurrentScale "new scale must be a multiple of the current scale"

    t = first(gens(underlying_polynomial_ring(f)))
    newPoly = evaluate(poly(f), [t^scaleQuotient])
    newShift = shift(f) * scaleQuotient
    return PuiseuxPolynomialRingElem(parent(f), newPoly, newShift, newScale)
end


#################################################################################
#
# Converstions
#
#################################################################################

function (Kt::PuiseuxPolynomialRing)(i::Int)
    K = base_ring(Kt)
    return K(i)
end

#################################################################################
#
# Properties
#
#################################################################################

elem_type(::Type{PuiseuxPolynomialRing{T}}) where T <: FieldElement = PuiseuxPolynomialRingElem{T}
parent_type(::Type{PuiseuxPolynomialRingElem{T}}) where T <: FieldElement = PuiseuxPolynomialRing{T}
base_ring_type(::Type{PuiseuxPolynomialRing{T}}) where T <: FieldElement = parent_type(T)

gen(R::PuiseuxPolynomialRing) = PuiseuxPolynomialRingElem(R, first(gens(underlying_polynomial_ring(R))))
zero(R::PuiseuxPolynomialRing) = PuiseuxPolynomialRingElem(R, zero(underlying_polynomial_ring(R)))
one(R::PuiseuxPolynomialRing) = PuiseuxPolynomialRingElem(R, one(underlying_polynomial_ring(R)))

function hash(f::PuiseuxPolynomialRingElem, h::UInt)
    normalize!(f)
    return hash((parent(f), poly(f), shift(f), scale(f)), h)
end

iszero(f::PuiseuxPolynomialRingElem) = iszero(poly(f))
isone(f::PuiseuxPolynomialRingElem) = isone(poly(f)) && shift(f) == 0 && scale(f) == 1

function ==(f::PuiseuxPolynomialRingElem, g::PuiseuxPolynomialRingElem)
    @assert parent(f) == parent(g) "elements must be in the same ring"
    return poly(f) == poly(g) && shift(f) == shift(g) && scale(f) == scale(g)
end

coefficients(f::PuiseuxPolynomialRingElem) = coefficients(poly(f))
function exponents(f::PuiseuxPolynomialRingElem)
    d = scale(f)
    k = shift(f)
    underlying_exponents = first.(exponents(poly(f)))
    return [(i + k) // d for i in underlying_exponents]
end

function length(f::PuiseuxPolynomialRingElem)
    return length(poly(f))
end


#################################################################################
#
# Printing
#
#################################################################################

function Base.show(io::IO, R::PuiseuxPolynomialRing)
    print(io, "Puiseux polynomial ring over ", coefficient_field(R))
end

function Base.show(io::IO, f::PuiseuxPolynomialRingElem)
    if iszero(f)
        print(io, "0")
        return
    elseif isone(f)
        print(io, "1")
        return
    end

    t = first(gens(underlying_polynomial_ring(parent(f))))
    terms = []
    for (c, e) in zip(coefficients(poly(f)), first.(exponents(poly(f))))
        push!(terms, "$(c)*$(t)^($( (e + shift(f)) // scale(f) ))")
    end
    print(io, join(terms, " + "))
end


#################################################################################
#
# Operations
#
#################################################################################

function +(f::PuiseuxPolynomialRingElem, g::PuiseuxPolynomialRingElem)
    @assert parent(f) == parent(g) "elements must be in the same ring"
    if iszero(f)
        return g
    elseif iszero(g)
        return f
    end

    # rescale both to the lcm of their scales
    newScale = lcm(scale(f), scale(g))
    frescaled = rescale(f, newScale)
    grescaled = rescale(g, newScale)

    # add the polynomials, adjusting for shifts
    newShift = min(shift(frescaled), shift(grescaled))
    t = first(gens(underlying_polynomial_ring(f)))
    newPoly = t^(shift(frescaled)-newShift)*poly(frescaled) + t^(shift(grescaled)-newShift)*poly(grescaled)

    # normalize output, in case of cancellations
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

    # add shifts, multiply scales and polys
    newShift = shift(f)*scale(g) + shift(g)*scale(f)
    newScale = scale(f)*scale(g)
    newPoly = poly(f)*poly(g)

    # normalize output
    ftimesg = PuiseuxPolynomialRingElem(
        parent(f),
        newPoly,
        newShift,
        newScale)
    # todo: remove assertion after exhaustive testing
    @assert !(normalize!(ftimesg)) "product was not normalized"
    return ftimesg
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
