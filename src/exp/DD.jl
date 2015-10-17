
function exp(a::DD)
   aa=TD(a.hi,a.lo,zero(Float64))
   aa=exp(aa)
   hi,lo=aa.hi,aa.md
   DD(hi,lo)
end



#   This is a natural logarithm (i.e., base e).
function log(a::DD)
  if (a == one(DD))
    return zero(DD)
  elseif (a.hi == zero(Float64))
    return DD(-Inf,-Inf)
  elseif (a.hi < zero(Float64))
    throw( ErrorException("log: Non-positive argument.") )
    return TD(NaN,NaN,NaN)
  end

  x = log(a.hi)
  xx = DD(x,zero(Float64))

  xx = xx + a * exp(-xx) - 1.0;
  #xx = xx + a * exp(-xx) - 1.0;
  y = TD(xx.hi,xx.lo,zero(Float64))
  y = y + a * exp(-y) - 1.0;
  DD(y.hi,y.md)
end



(^){T<:DD}(a::T,b::Int) = (b >= 0) ? npowr(a,b) : nroot(a,b)
(^){T<:DD}(a::T,b::T) = exp(b * log(a))
(^)(a::DD,b::Float64) = exp(DD(b) * log(a))
(^)(a::Float64,b::DD) = exp(b * log(DD(a)))
(^){T<:Signed}(a::DD,b::T) = exp(convert(DD,b) * log(a))
(^){T<:Signed}(a::T,b::DD) = exp(b * log(convert(DD,a)))
