#################################################################################
#
# tropical hypersurfaces from puiseux polynomials
#
#################################################################################

@doc raw"""
    tropical_hypersurface_over_fraction_field(f::AbstractAlgebra.Generic.MPoly{MPuiseuxPolyRingElem{K}}, nu::TropicalSemiringMap{MPuiseuxPolyRing{K}, MPuiseuxPolyRingElem{K}, MinOrMax}) where {K<:Ring, MinOrMax<:Union{typeof(min),typeof(max)}}

Return the tropical hypersurface of the tropical polynomial that is the image of
`f` under coefficient-wise `val`.  This is the tropical hypersurface of `f` over
the fraction field with respect to the valuation encoded by `nu`.  It is not the
tropical hypersurface of `f` according to Section 1.6 in Maclagan-Sturmfels.

# Examples
```jldoctest
julia> R,(x,y) = QQ[:x, :y];

julia> val = tropical_semiring_map(QQ,2)
Map into Min tropical semiring encoding the 2-adic valuation on Rational field

julia> f = x+y+2
x + y + 2

julia> tropical_hypersurface(f,val)
Min tropical hypersurface

```
"""
function tropical_hypersurface_over_fraction_field(
    f::AbstractAlgebra.Generic.MPoly{MPuiseuxPolyRingElem{K}},
    nu::TropicalSemiringMap{MPuiseuxPolyRing{K}, MPuiseuxPolyRingElem{K}, MinOrMax};
    weighted_polyhedral_complex_only::Bool=false
    ) where {K<:RingElem, MinOrMax<:Union{typeof(min),typeof(max)}}

    tropf = tropical_polynomial(f,nu)
    TropH = tropical_hypersurface(tropf,weighted_polyhedral_complex_only=weighted_polyhedral_complex_only)

    if !weighted_polyhedral_complex_only
        set_attribute!(TropH,:algebraic_polynomial,f)
        set_attribute!(TropH,:tropical_semiring_map,nu)
    end
    return TropH
end
