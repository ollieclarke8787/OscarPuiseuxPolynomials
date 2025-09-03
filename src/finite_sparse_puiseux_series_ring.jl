struct FiniteSparsePuiseuxSeriesRing <: Ring
    base_ring::MPolyRing

    function FiniteSparsePuiseuxSeriesRing(base_ring::MPolyRing)
        @assert ngens(base_ring) == 1 "underlying ring must be a multivariate polynomial ring with one variable"
        return new(base_ring)
    end
end


function Base.show(io::IO, R::FiniteSparsePuiseuxSeriesRing)
    print(io, "Finite sparse Puiseux series ring over ", coefficient_field(R))
end

base_ring(R::FiniteSparsePuiseuxSeriesRing) = R.base_ring
coefficient_field(R::FiniteSparsePuiseuxSeriesRing) = base_ring(base_ring(R))

# gen(R::FiniteSparsePuiseuxSeriesRing) = gens(base_ring(R))[1] ## this should be of type FiniteSparsePuiseuxSeriesRingElem


# given a puiseux series f = a_0 t^(k_0/d) + a_1 t^(k_1/d) + ... where k_0 < k_1 < ... integers
# we represent it as a polynomial: a_0 + a_1 t^(k_1 - k_0) + ...
# together with the shift k_0 and scale d
mutable struct FiniteSparsePuiseuxSeriesRingElem <: RingElem
    parent::FiniteSparsePuiseuxSeriesRing
    poly::MPolyRingElem
    shift::ZZRingElem
    scale::ZZRingElem

    function FiniteSparsePuiseuxSeriesRingElem(Kt::FiniteSparsePuiseuxSeriesRing, f::MPolyRingElem, k::ZZRingElem = zero(ZZ), d::ZZRingElem = one(ZZ))
        @assert parent(f) == base_ring(Kt) "polynomial must be in the underlying ring"
        @assert d > 0 "scale must be positive"
        return new(Kt, f, k, d)
    end
end

parent(f::FiniteSparsePuiseuxSeriesRingElem) = f.parent
poly(f::FiniteSparsePuiseuxSeriesRingElem) = f.poly
shift(f::FiniteSparsePuiseuxSeriesRingElem) = f.shift
scale(f::FiniteSparsePuiseuxSeriesRingElem) = f.scale

iszero(f::FiniteSparsePuiseuxSeriesRingElem) = iszero(poly(f))
isone(f::FiniteSparsePuiseuxSeriesRingElem) = isone(poly(f)) && shift(f) == 0 && scale(f) == 1

function Base.show(io::IO, f::FiniteSparsePuiseuxSeriesRingElem)
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


function normalize!(f::FiniteSparsePuiseuxSeriesRingElem)
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

function ==(f::FiniteSparsePuiseuxSeriesRingElem, g::FiniteSparsePuiseuxSeriesRingElem)
    @assert parent(f) == parent(g) "elements must be in the same ring"

    # bring f and g in standard form
    normalize!(f)
    normalize!(g)

    return poly(f) == poly(g) && shift(f) == shift(g) && scale(f) == scale(g)
end

function hash(f::FiniteSparsePuiseuxSeriesRingElem, h::UInt)
    normalize!(f)
    return hash((parent(f), poly(f), shift(f), scale(f)), h)
end

gen(R::FiniteSparsePuiseuxSeriesRing) = FiniteSparsePuiseuxSeriesRingElem(R, first(gens(base_ring(R))))
zero(R::FiniteSparsePuiseuxSeriesRing) = FiniteSparsePuiseuxSeriesRingElem(R, zero(base_ring(R)))
one(R::FiniteSparsePuiseuxSeriesRing) = FiniteSparsePuiseuxSeriesRingElem(R, one(base_ring(R)))

coefficient_field(f::FiniteSparsePuiseuxSeriesRingElem) = coefficient_field(parent(f))
coefficients(f::FiniteSparsePuiseuxSeriesRingElem) = coefficients(poly(f))
function exponents(f::FiniteSparsePuiseuxSeriesRingElem)
    d = scale(f)
    k = shift(f)
    underlying_exponents = first.(exponents(poly(f)))
    return [(i + k) // d for i in underlying_exponents]
end

function finite_sparse_puiseux_series_ring(K::Field, variableName::String="t")
    base_ring, _ = Oscar.polynomial_ring(K, [variableName])
    Kt = FiniteSparsePuiseuxSeriesRing(base_ring)
    t = gen(Kt)
    return Kt, t
end


function rescale(f::FiniteSparsePuiseuxSeriesRingElem, newScale::ZZRingElem)
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
    return FiniteSparsePuiseuxSeriesRingElem(parent(f), newPoly, newShift, newScale)
end


function +(f::FiniteSparsePuiseuxSeriesRingElem, g::FiniteSparsePuiseuxSeriesRingElem)
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
    fplusg = FiniteSparsePuiseuxSeriesRingElem(
        parent(f),
        newPoly,
        newShift,
        newScale)
    normalize!(fplusg)
    return fplusg
end

function -(f::FiniteSparsePuiseuxSeriesRingElem)
    return FiniteSparsePuiseuxSeriesRingElem(parent(f), -poly(f), shift(f), scale(f))
end

function -(f::FiniteSparsePuiseuxSeriesRingElem, g::FiniteSparsePuiseuxSeriesRingElem)
    return f + (-g)
end

function *(f::FiniteSparsePuiseuxSeriesRingElem, g::FiniteSparsePuiseuxSeriesRingElem)
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
    ftimesg = FiniteSparsePuiseuxSeriesRingElem(
        parent(f),
        newPoly,
        newShift,
        newScale)

    normalize!(ftimesg)
    return ftimesg
end

function length(f::FiniteSparsePuiseuxSeriesRingElem)
    return length(poly(f))
end

function ^(f::FiniteSparsePuiseuxSeriesRingElem, a::QQFieldElem)

    if denominator(a) == 1
        return f^numerator(a)
    end

    @assert length(f) == 1 "only single-term series can be exponentiated to rational powers"

    return FiniteSparsePuiseuxSeriesRingElem(
        parent(f),
        poly(f),
        shift(f)*numerator(a),
        scale(f)*denominator(a))
end

function ^(f::FiniteSparsePuiseuxSeriesRingElem, a::Rational{Int})
    return f^(QQ(a))
end
