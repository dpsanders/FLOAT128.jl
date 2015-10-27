#=
#    Internal Use Only
=#

isnan(a::TD) = isnan(a.hi)
isinf(a::TD) = isinf(a.hi)
isfinite(a::TD) = isfinite(a.hi)

sign(a::TD) = sign(a.hi)
signbit(a::TD) = signbit(a.hi)

(-)(a::TD) = TD(-a.hi,-a.md,-a.lo)
(abs)(a::TD) = (a.hi >= zero(Float64) ? a : -a)


@inline flipsign(a::TD,b::TD) = TD(flipsign(a.hi,b.hi), flipsign(a.md,b.hi), flipsign(a.lo,b.hi))
@inline flipsign(a::TD,b::DD) = TD(flipsign(a.hi,b.hi), flipsign(a.md,b.hi), flipsign(a.lo,b.hi))
flipsign(a::DD,b::TD) = DD(flipsign(a.hi,b.hi), flipsign(a.lo,b.hi))
flipsign(a::TD,b::Float64) = TD(flipsign(a.hi,b), flipsign(a.md,b.hi), flipsign(a.lo,b))
flipsign(a::TD,b::Integer) = TD(flipsign(a.hi,b), flipsign(a.md,b.hi), flipsign(a.lo,b))

@inline function copysign(a::TD,b::Float64)
    if (b < zero(Float64))
        a.hi < zero(Float64) ? a : -a
    else
        a.hi < zero(Float64) ? -a : a
end
copysign(a::TD,b::Integer) = copysign(a,convert(Float64,b))
copysign(a::TD,b::DD) = copysign(a,b.hi)
@inline copysign(a::TD,b::TD) = copysign(a, b.hi)
copysign(a::DD,b::TD) = copysign(a.hi,b.hi)


@inline function mulby2(a::TD)
    TD(a.hi*2.0, a.md*2.0, a.lo*2.0)
end

@inline function divby2(a::TD)
    TD(a.hi*0.5, a.md*0.5, a.lo*0.5)
end

@inline function mulbypow2(a::TD,p::Float64)
    TD(a.hi*p, a.md*p, a.lo*p)
end

@inline function divbypow2(a::TD,p::Float64)
    fr,xp = frexp(p)
    mulbypow2(a, ldexp(fr,-xp))
end


function (floorceil)(fn::Function, a::TD)
    hi = fn(a.hi)
    md = lo = zero(Float64)
    if (hi == a.hi)
        md = fn(a.md)
        hi,md = eftSum2inOrder(hi,md)
        if md == a.md
            lo = fn(a.lo)
            md,lo = eftSum2inOrder(md,lo)
            hi,md = eftSum2inOrder(hi,md)
        end
    end
    TD(hi,md,lo)
end

(floor)(a::TD) = floorceil(floor,a)
(ceil)(a::TD) = floorceil(ceil,a)


function round(a::TD)
    hi = round(a.hi)
    md = lo = zero(Float64)
    if (hi == a.hi)
        md = round(a.md)
        hi,md = eftSum2inOrder(hi,md)
        if md == a.md
            lo = round(a.lo)
            md,lo = eftSum2inOrder(md,lo)
            hi,md = eftSum2inOrder(hi,md)
        end
    end
    TD(hi,md,lo)
end
#(round)(a::TD) = floorceil(round,a)


@inline function (trunc)(a::TD)
    a.hi >= zero(Float64) ? floor(a) : ceil(a)
end

"""
stretch is the opposite of trunc()
it extends to the nearest integer away from zero
"""
@inline function (stretch)(a::TD)
    a.hi >= zero(Float64) ? ceil(a) : floor(a)
end

function fld{T<:TD}(a::T,b::T)
    if (b.hi == zero(Float64))
        throw(DomainError("denominator must be nonzero"))
    end
    floor( a/b )
end

function cld{T<:TD}(a::T,b::T)
    if (b.hi == zero(Float64))
        throw(DomainError("denominator must be nonzero"))
    end
    ceil( a/b )
end


function div{T<:TD}(a::T,b::T)
    if (b.hi == zero(Float64))
        throw(DomainError("denominator must be nonzero"))
    end
    temp = a/b
    trunc(temp)
end
div(a::TD,b::Float64) = div(a,TD(b))
div(a::Float64,b::TD) = div(TD(a),b)

function rem{T<:TD}(a::T,b::T)
    if (b.hi == zero(Float64))
        throw(DomainError("denominator must be nonzero"))
    end
    a - div(a,b)*b
end

%{T<:TD}(a::T,b::T) = a - trunc(a/b)
%(a::TD,b::Float64) = %(a,b)
%(a::Float64,b::TD) = %(a,b)

function divrem{T<:TD}(a::T,b::T)
    if (b.hi == zero(Float64))
        throw(DomainError("denominator must be nonzero"))
    end
    temp = a/b
    d = trunc(temp)
    r = a - d*b
    d,r
end


function mod{T<:TD}(a::T,b::T)
    if (b.hi == zero(Float64))
        throw(DomainError("modulus must be nonzero"))
    end
    if signbit(a.hi)==signbit(b.hi)
        rem(a,b)
    else
        d = floor(a/b)
        a - d*b
    end
end

mod(a::TD,b::Float64) = mod(a,convert(TD,b))
mod(a::Float64,b::TD) = mod(convert(TD,a),b)

function fldmod{T<:TD}(a::T,b::T)
    d = floor(a/b)
    a - d*b
    d,a
end
