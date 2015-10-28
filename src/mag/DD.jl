@inline function isless(a::DD, b::DD)
    (a.hi < b.hi) || (a.hi==b.hi && a.lo<b.lo)
end
@inline function isless(a::DD, b::Float64)
    (a.hi < b) || (a.hi==b && a.lo<zero(Float64))
end
@inline function isless(a::Float64, b::DD)
    (a < b.hi) || (a==b.hi && b.lo<zero(Float64))
end
@inline function isequal(a::DD, b::DD)
    (a.hi == b.hi) && (a.lo == b.lo)
end
@inline function isequal(a::DD, b::Float64)
    (a.hi == b) && (a.lo == zero(Float64))
end
@inline function isequal(a::Float64, b::DD)
    (a == b.hi) && (b.lo == zero(Float64))
end

@inline (==)(a::DD,b::DD) = (a.hi == b.hi) && (a.lo == b.lo)
@inline (<)(a::DD,b::DD) = (a.hi < b.hi) || (a.hi==b.hi && a.lo<b.lo)
@inline (<=)(a::DD,b::DD) =  (a.hi < b.hi) || (a.hi==b.hi && a.lo<=b.lo)
@inline (>)(a::DD,b::DD) = (a.hi > b.hi) || (a.hi==b.hi && a.lo>b.lo)
@inline (>=)(a::DD,b::DD) = (a.hi > b.hi) || (a.hi==b.hi && a.lo>=b.lo)




function (floor)(a::DD)
    hi = floor(a.hi)
    lo = 0.0
    if (hi == a.hi)
        lo = floor(a.lo)
        hi,lo = eftSum2inOrder(hi,lo)
    end
    DD(hi,lo)
end

function (ceil)(a::DD)
    hi = ceil(a.hi)
    lo = 0.0
    if (hi == a.hi)
        lo = ceil(a.lo)
        hi,lo = eftSum2inOrder(hi,lo)
    end
    DD(hi,lo)
end

@inline function (trunc)(a::DD)
    a.hi >= zero(F) ? floor(a) : ceil(a)
end

function ldexp(a::DD,xp::Int)
    DD(ldexp(a.hi,xp),ldexp(a.lo,xp))
end

function frexp(a::DD)
    frhi, xphi = frexp(a.hi)
    frlo, xplo = frexp(a.lo)
    DD(frhi, ldexp(frlo,xplo-xphi)), xphi
end
