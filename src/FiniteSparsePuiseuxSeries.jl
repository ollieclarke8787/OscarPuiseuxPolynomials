module FiniteSparsePuiseuxSeries

using Oscar

# functions we extend
import Base: parent, show, ==, +, -, *, //, hash, ^, length
import Oscar: zero, one, iszero, isone, coefficients, exponents, coefficient_field, gen, base_ring, normalize!

# source files
include("finite_sparse_puiseux_series_ring.jl")
include("exports.jl")


end
