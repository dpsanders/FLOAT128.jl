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
