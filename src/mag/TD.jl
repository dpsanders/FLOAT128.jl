function (floor){F<:MachineFloat}(a::TD{F})
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

function (ceil){F<:MachineFloat}(a::TD{F})
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

@inline function (trunc){F<:MachineFloat}(a::TD{F})
    a.hi >= zero(F) ? floor(a) : ceil(a)
end

function ldexp{F<:MachineFloat}(a::TD{F},xp::Int)
    TD(ldexp(a.hi,xp),ldexp(a.md,xp),ldexp(a.lo,xp))
end

function frexp{F<:MachineFloat}(a::TD{F})
    frhi, xphi = frexp(a.hi)
    frmd, xpmd = frexp(a.md)
    frlo, xplo = frexp(a.lo)
    TD(frhi, ldexp(frmd,xpmd-xphi), ldexp(frlo,xplo-xphi)), xphi
end

