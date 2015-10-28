function (floor){F<:MachineFloat}(a::DD{F})
    hi = floor(a.hi)
    lo = 0.0
    if (hi == a.hi)
        lo = floor(a.lo)
        hi,lo = eftSum2inOrder(hi,lo)
    end
    DD(hi,lo)
end

function (ceil){F<:MachineFloat}(a::DD{F})
    hi = ceil(a.hi)
    lo = 0.0
    if (hi == a.hi)
        lo = ceil(a.lo)
        hi,lo = eftSum2inOrder(hi,lo)
    end
    DD(hi,lo)
end

@inline function (trunc){F<:MachineFloat}(a::DD{F})
    a.hi >= zero(F) ? floor(a) : ceil(a)
end

function ldexp{F<:MachineFloat}(a::DD{F},xp::Int)
    DD(ldexp(a.hi,xp),ldexp(a.lo,xp))
end

function frexp{F<:MachineFloat}(a::DD{F})
    frhi, xphi = frexp(a.hi)
    frlo, xplo = frexp(a.lo)
    DD(frhi, ldexp(frlo,xplo-xphi)), xphi
end
