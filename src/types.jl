typealias SignedInt Union{Int128,Int64,Int32,Int16}
typealias MachineFloat Union{Float64,Float32,Float16}

#=
The internal structure is 128 bits given as a magnitude-ordered pair of Float64 values.\\
The pair are *presumed* to be in cannonical form, where they do not overlap: hi⊕lo ≡ hi.\\
This representation is known as "double-double" (see refs.md)
=#

immutable DD <: Real
    hi::Float64
    lo::Float64
end

DD(a::Float64) = DD(a,zero(Float64))
DD{T<:Float64}(a::Tuple{T}) = DD(a[1])
DD{T<:Float64}(a::Tuple{T,T}) = DD(a[1],a[2])

convert(::Type{DD}, a::Float64) = DD(a)
convert(::Type{Float64}, a::DD) = a.hi
convert{T<:SignedInt}(::Type{DD}, a::T) = DD(convert(Float64,a))
convert{T<:SignedInt}(::Type{T}, a::DD) = convert(T,floor(a.hi))+convert(T,round(a.lo))

import Base:Tuple # needed??

convert(::Type{DD}, a::Tuple{Float64}) = DD(a[1])
convert(::Type{DD}, a::Tuple{Float64,Float64}) = DD(a[1],a[2])
convert(::Type{Tuple}, a::DD) = (a.hi,a.lo)
convert(::Type{Tuple{Float64,Float64}}, a::DD) = (a.hi,a.lo)

promote_rule(::Type{DD}, ::Type{Float64}) = DD
promote_rule(::Type{DD}, ::Type{Float32}) = DD
promote_rule(::Type{DD}, ::Type{Int64}) = DD
promote_rule(::Type{DD}, ::Type{Int32}) = DD


const epseps1 = eps(eps(1.0))

function clean(x::DD)
    if (abs(x.hi) <= epseps1)
        zero(DD)
    elseif (abs(x.lo) <= eps(eps(x.hi)))
        DD(x.hi,0.0)
    else
        x
    end
end

clean(v::Vector{DD}) = [clean(x) for x in v]
clean(m::Matrix{DD}) = reshape([clean(x) for x in m],size(m))
clean(a::Array{DD}) = reshape([clean(x) for x in a],size(a))


zero(::Type{DD}) = DD(zero(Float64),zero(Float64))
zero(a::DD)      = zero(DD)
one(::Type{DD})  = DD(one(Float64),zero(Float64))
one(a::DD)       = one(DD) 

#=
FL127(a::Float64) = DD(a, zero(Float64))

function FL127{T<:Float64}(a::T, b::T)
    hi = a + b
    lo = x - a
    lo = (a - (hi - lo)) + (b - lo)
    DD(hi,lo)
end

for T in (:AbstractFloat, Integer)
  @eval begin
    FL127(a::($T)) = Float127(convert(Float64,a))
    FL127(a::($T),b::($T)) = Float127(convert(Float64,a),convert(Float64,b))
  end
end
=#

typealias Float128 DD


#=
Some internal computations require higher precision to meet the targeted relative error.\\
The internal structure is 192 bits given as a magnitude-ordered triple of Float64 values.\\
The triple are *presumed* to be in cannonical form: hi⊕md ≡ hi,  md⊕lo ≡ md\\
This representation is known as "triple-double" (see refs.md)
=#

immutable TD <: Real
    hi::Float64
    md::Float64
    lo::Float64
end

TD(a::Float64) = TD(a,zero(Float64),zero(Float64))
TD{T<:Float64}(a::T,b::T) = TD(a,b,zero(Float64))

TD{T<:Float64}(a::Tuple{T}) = TD(a[1])
TD{T<:Float64}(a::Tuple{T,T}) = TD(a[1],a[2])
TD{T<:Float64}(a::Tuple{T,T}) = TD(a[1],a[2],a[3])


zero(::Type{TD}) = TD(zero(Float64),zero(Float64),zero(Float64))
zero(a::TD)      = zero(TD)
one(::Type{TD})  = TD(one(Float64),zero(Float64),zero(Float64))
one(a::TD)       = one(TD) 

convert(::Type{TD}, a::Float64) = TD(a)
convert(::Type{Float64}, a::TD) = a.hi
convert{T<:SignedInt}(::Type{TD}, a::T) = TD(convert(Float64,a))
convert{T<:SignedInt}(::Type{T}, a::TD) = convert(T,floor(a.hi))+convert(T,floor(a.md))+convert(T,round(a.lo))

convert(::Type{TD}, a::DD) = TD(a.hi, a.lo, zero(Float64))
convert(::Type{DD}, a::TD) = DD(a.hi, a.md)



convert(::Type{TD}, a::Tuple{Float64}) = TD(a[1])
convert(::Type{TD}, a::Tuple{Float64,Float64}) = TD(a[1],a[2])
convert(::Type{TD}, a::Tuple{Float64,Float64,Float64}) = TD(a[1],a[2],a[3])
convert(::Type{Tuple}, a::TD) = (a.hi,a.md,a.lo)
convert(::Type{Tuple{Float64,Float64,Float64}}, a::TD) = (a.hi,a.md,a.lo)

promote_rule(::Type{TD}, ::Type{Float64}) = TD
promote_rule(::Type{TD}, ::Type{Float32}) = TD
promote_rule(::Type{TD}, ::Type{Int64}) = TD
promote_rule(::Type{TD}, ::Type{Int32}) = TD


function clean(x::TD)
    ep = eps(eps(1.0))
    if (abs(x.hi) <= ep)
        zero(TD)
    elseif (abs(x.md) <= eps(eps(x.hi)))
        DD(x.hi,0.0,0.0)
    elseif (abs(x.lo) <= eps(eps(x.md)))
        DD(x.hi,x.md,0.0)
    else
        x
    end
end


clean(x::Float32) = abs(x) < eps(eps(one(Float32))) ? zero(Float32) : x
clean(x::Float64) = abs(x) < eps(eps(one(Float64))) ? zero(Float64) : x
clean(x::AbstractFloat) = abs(x) < eps(eps(one(typeof(x)))) ? zero(typeof(x)) : x
clean(x::Complex) = clean(x.re)+clean(x.im)*im
