@inline function isless(a::TD, b::TD)
    (a.hi < b.hi) || (a.hi==b.hi && a.md<b.md) || (a.hi==b.hi && a.md==b.md && a.lo<b.lo)
end
@inline function isless(a::TD, b::Float64)
    (a.hi < b) || (a.hi==b && a.md<zero(Float64))
end
@inline function isless(a::Float64, b::TD)
    (a < b.hi) || (a==b.hi && b.lo>zero(Float64))
end
@inline function isequal(a::TD, b::TD)
    (a.hi == b.hi) && (a.md == b.md) && (a.lo == b.lo)
end
@inline function isequal(a::TD, b::Float64)
    (a.hi == b) && (a.md == zero(Float64))
end
@inline function isequal(a::Float64, b::TD)
    (a == b.hi) && (b.md == zero(Float64))
end

@inline (==)(a::TD,b::TD) = (a.hi == b.hi) && (a.md == b.md) && (a.lo == b.lo)
@inline (< )(a::TD,b::TD) = (a.hi < b.hi) || (a.hi==b.hi && a.md<b.md) || (a.hi==b.hi && a.md==b.md && a.lo<b.lo)
@inline (<=)(a::TD,b::TD) = (a.hi < b.hi) || (a.hi==b.hi && a.md<b.md) || (a.hi==b.hi && a.md==b.md && a.lo<=b.lo)
@inline (> )(a::TD,b::TD) = (a.hi > b.hi) || (a.hi==b.hi && a.md>b.md) || (a.hi==b.hi && a.md==b.md && a.lo>b.lo)
@inline (>=)(a::TD,b::TD) = (a.hi > b.hi) || (a.hi==b.hi && a.md>b.md) || (a.hi==b.hi && a.md==b.md && a.lo>=b.lo)

@inline (==)(a::TD,b::DD) = (==)(a,TD(b))
@inline (< )(a::TD,b::DD) = (< )(a,TD(b))
@inline (<=)(a::TD,b::DD) = (==)(a,TD(b))
@inline (> )(a::TD,b::DD) = (> )(a,TD(b))
@inline (>=)(a::TD,b::DD) = (>=)(a,TD(b))

@inline (==)(a::DD,b::TD) = (==)(TD(a),b)
@inline (< )(a::DD,b::TD) = (< )(TD(a),b)
@inline (<=)(a::DD,b::TD) = (==)(TD(a),b)
@inline (> )(a::DD,b::TD) = (> )(TD(a),b)
@inline (>=)(a::DD,b::TD) = (>=)(TD(a),b)


function (floor)(a::TD)
    hi = floor(a.hi)
    md = lo = 0.0
    if (hi == a.hi)
        md = floor(a.md)
        hi,md = eftSum2inOrder(hi,md)
        if md == a.md
            lo = floor(a.lo)
            md,lo = eftSum2inOrder(md,lo)
            hi,md = eftSum2inOrder(hi,md)
        end
    end
    TD(hi,md,lo)
end

function (ceil)(a::TD)
    hi = ceil(a.hi)
    md = lo = 0.0
    if (hi == a.hi)
        md = ceil(a.md)
        hi,md = eftSum2inOrder(hi,md)
        if md == a.md
            lo = ceil(a.lo)
            md,lo = eftSum2inOrder(md,lo)
            hi,md = eftSum2inOrder(hi,md)
        end
    end
end

@inline function (trunc)(a::TD)
    a.hi >= zero(F) ? floor(a) : ceil(a)
end

function ldexp(a::TD,xp::Int)
    TD(ldexp(a.hi,xp),ldexp(a.md,xp),ldexp(a.lo,xp))
end

function frexp(a::TD)
    frhi, xphi = frexp(a.hi)
    frmd, xpmd = frexp(a.md)
    frlo, xplo = frexp(a.lo)
    TD(frhi, ldexp(frmd,xpmd-xphi), ldexp(frlo,xplo-xphi)), xphi
end

