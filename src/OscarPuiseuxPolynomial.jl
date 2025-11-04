module OscarPuiseuxPolynomial

using Oscar

# functions we extend
include("imports.jl")

# source files
include("puiseux_polynomial.jl")

# (non-generic) tropical features
include("semiring_map.jl")
include("hypersurface.jl")

# functions we define
include("exports.jl")


end
