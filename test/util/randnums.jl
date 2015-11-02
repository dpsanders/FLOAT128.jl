#=
random numbers for testing

bitgap is count of zero bits separating paired values
frval is the value used as frexp(val)[1]
expon is the value used as frexp(val)[2]
=#

if !isdefined(:DD)
   DD = FLOAT128.DD
end
if !isdefined(:TD)
   TD = FLOAT128.TD
end

using(Distributions)

baseBitgapDist = Geometric(0.435)
wideBitgapDist = Geometric(0.235)

baseBitgap() = rand(baseBitgapDist)
wideBitgap() = rand(wideBitgapDist)

function frval()
    s = rand()
    if (s < 0.5)
        s += 0.5
    end
    s
end

expon()  = rand(-5:5)

msval()  = ldexp(frval(),expon())
function lsval(frval::AbstractFloat, bitgap::Function)
    vexpon = frexp(frval)[2] - 54 - bitgap()
    ldexp(frval,vexpon)
end

function randhighlow(bitgap::Function=wideBitgap, expow2::Function=expon)
    fr = frval()
    xp = expow2()
    hi = ldexp(fr,xp)
    xp = xp - 54 - bitgap()
    lo = ldexp(frval(),xp)
    hi,lo
end

function randhighmedlow(bitgap::Function=wideBitgap, expow2::Function=expon)
    fr = frval()
    xp = expow2()
    hi = ldexp(fr,xp)
    xp = xp - 54 - bitgap()
    md = ldexp(frval(),xp)
    xp = xp - 54 - bitgap()
    lo = ldexp(frval(),xp)
    hi,md,lo
end

randhilo(bitgap::Function=wideBitgap, expow2::Function=expon) = DD(randhighlow(bitgap,expow2)...)
randhml(bitgap::Function=wideBitgap, expow2::Function=expon) = TD(randhighmedlow(bitgap,expow2)...)

rangedhilo(bg::Float64,en::Float64) = DD(rand(bg:en))
