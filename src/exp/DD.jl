
#const dd_n_inv_fact = 32;
const dd_inv_fact = DD[
  DD(1.0,0.0),
  DD(0.5,0.0),
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
  DD(1.5619206968586225e-16,1.1910679660273754e-32),
  DD(8.22063524662433e-18,2.2141894119604265e-34),
  DD(4.110317623312165e-19,1.4412973378659527e-36),
  DD(1.9572941063391263e-20,-1.3643503830087908e-36),
  DD(8.896791392450574e-22,-7.911402614872376e-38),
  DD(3.868170170630684e-23,-8.843177655482344e-40),
  DD(1.6117375710961184e-24,-3.6846573564509766e-41),
  DD(6.446950284384474e-26,-1.9330404233703465e-42),
  DD(2.4795962632247976e-27,-1.2953730964765229e-43),
  DD(9.183689863795546e-29,1.4303150396787322e-45),
  DD(3.279889237069838e-30,1.5117542744029879e-46),
  DD(1.1309962886447716e-31,1.0498015412959506e-47),
  DD(3.7699876288159054e-33,2.5870347832750324e-49),
  DD(1.216125041553518e-34,5.586290567888806e-51),
  DD(3.8003907548547434e-36,1.7457158024652518e-52)
];


const dd_n_exp_ints = 32;
const dd_exp_ints = DD[
    DD(2.718281828459045, 1.4456468917292502e-16),
    DD(7.38905609893065, -1.7971139497839148e-16),
    DD(20.085536923187668, -1.8275625525512858e-16),
    DD(54.598150033144236, 2.8741578015844115e-15),
    DD(148.4131591025766, 3.4863514900464198e-15),
    DD(403.4287934927351, 1.2359628024450387e-14),
    DD(1096.6331584284585, 9.869752640434095e-14),
    DD(2980.9579870417283, -2.7103295816873633e-14),
    DD(8103.083927575384, -2.1530877621067177e-13),
    DD(22026.465794806718, -1.3780134700517372e-12),
    DD(59874.14171519782, 1.7895764888916994e-12),
    DD(162754.79141900392, 5.30065881322063e-12),
    DD(442413.3920089205, 1.2118711752313224e-11),
    DD(1.2026042841647768e6, -1.5000525764327354e-11),
    DD(3.2690173724721107e6, -3.075806431120808e-11),
    DD(8.886110520507872e6, 5.321182483501564e-10),
    DD(2.41549527535753e7, -7.203995068362157e-10),
    DD(6.565996913733051e7, 1.4165536846555444e-9),
    DD(1.7848230096318725e8, 1.333018530234341e-8),
    DD(4.851651954097903e8, 4.880277289790406e-10),
    DD(1.3188157344832146e9, 8.043448618843281e-8),
    DD(3.584912846131592e9, -2.3519384005402157e-7),
    DD(9.744803446248903e9, -6.74501500127677e-7),
    DD(2.648912212984347e10, 7.670395527778119e-7),
    DD(7.200489933738588e10, -6.992440211033874e-6),
    DD(1.9572960942883878e11, -1.1364989227123904e-5),
    DD(5.3204824060179865e11, -2.8335783945658822e-5),
    DD(1.446257064291475e12, 7.602079742299693e-5),
    DD(3.931334297144042e12, 8.220112058084352e-5),
    DD(1.0686474581524463e13, -0.0007436345313492586),
    DD(2.9048849665247426e13, -0.0005501643178883202),
    DD(7.896296018268069e13, 0.007660978022635108),
];

function exp_taylor(a::DD)
  x = a
  x2 = x*x
  x3 = x*x2
  x4 = x2*x2
  x8 = x4*x4
  x12 = x4*x8
  x16 = x8*x8
  x20 = x8*x12
  x24 = x12*x12

  z = x + dd_inv_fact[2]*x2 + dd_inv_fact[3]*x3
  z2 = x4 * (dd_inv_fact[4] + x*dd_inv_fact[5] + x2*dd_inv_fact[6] + x3*dd_inv_fact[7])
  z3 = x8 * (dd_inv_fact[8] + x*dd_inv_fact[9] + x2*dd_inv_fact[10] + x3*dd_inv_fact[11])
  z4 = x12 * (dd_inv_fact[12] + x*dd_inv_fact[13] + x2*dd_inv_fact[14] + x3*dd_inv_fact[15])
  z5 = x16 * (dd_inv_fact[16] + x*dd_inv_fact[17] + x2*dd_inv_fact[18] + x3*dd_inv_fact[19])
  z6 = x20 * (dd_inv_fact[20] + x*dd_inv_fact[21] + x2*dd_inv_fact[22] + x3*dd_inv_fact[23])
  z7 = x24 * (dd_inv_fact[24] + x*dd_inv_fact[25] + x2*dd_inv_fact[26] + x3*dd_inv_fact[27])

  (((((z7+z6)+z5)+z4)+z3)+z2)+z + one(DD)
end


function exp01(a::DD)
  if a.hi >= 0.5
    sqr(exp_taylor(divby2(a)))
  else
    exp_taylor(a)
  end
end

function expGT0(a::DD)
    if a.hi < 1e-34
         return one(DD)
    elseif a.hi==1.0 && abs(a.lo) < 1e-34
         return dd_exp1
    end

    if a.hi > 1.0
      if a.hi <= 21.0
        ddidx = DD(floor(a.hi))
        idx = floor(Int,a.hi)
        expint = dd_exp_ints[idx]
        dlt = a - ddidx
        expdlt = exp01(dlt)
        expint*expdlt
      else
        t = TD(a)
        DD( exp(t) )
      end
    else
      exp01(a)
    end
end

function exp(a::DD)
    sb = signbit(a.hi)
    if abs(a.hi) == Inf
         return sb ? zero(DD) : a
    end

    if sb
         recip(expGT0(abs(a)))
    else
         expGT0(a)
    end
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
