const n_inv_fact = 16;
const inv_fact = DD[
DD(0.16666666666666666,9.25185853854297e-18),
DD(0.041666666666666664,2.3129646346357427e-18),
DD(0.008333333333333333,1.1564823173178714e-19),
DD(0.001388888888888889,-5.300543954373577e-20),
DD(0.0001984126984126984,1.7209558293420705e-22),
DD(2.48015873015873e-5,2.1511947866775882e-23),
DD(2.7557319223985893e-6,-1.858393274046472e-22),
DD(2.755731922398589e-7,2.3767714622250297e-23),
DD(2.505210838544172e-8,-1.448814070935912e-24),
DD(2.08767569878681e-9,-1.20734505911326e-25),
DD(1.6059043836821613e-10,1.2585294588752098e-26),
DD(1.1470745597729725e-11,2.0655512752830745e-28),
DD(7.647163731819816e-13,7.03872877733453e-30),
DD(4.779477332387385e-14,4.399205485834081e-31),
DD(2.8114572543455206e-15,1.6508842730861433e-31),
];


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
