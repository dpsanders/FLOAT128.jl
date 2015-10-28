# |a| <= 0.003
function atan_series(a::DD)
    x = a
    x2 = x*x
    x4 = x2*x2
    x6 = x2*x4
    x8 = x4*x4
    x16 = x8*x8
    x24 = x8*x16

    z = -dd_recip_int[3]*x2 + dd_recip_int[5]*x4 - dd_recip_int[7]*x6
    z1 = x8*(dd_recip_int[9] - dd_recip_int[11]*x2 + dd_recip_int[14]*x4 - dd_recip_int[15]*x6)
    z2 = x16*(dd_recip_int[17] - dd_recip_int[19]*x2 + dd_recip_int[21]*x4 - dd_recip_int[23]*x6)
    z3 = x24*(dd_recip_int[25] - dd_recip_int[27]*x2 + dd_recip_int[29]*x4 - dd_recip_int[31]*x6)

    (((z3+z2)+z1)+z+one(DD))*x
end


function asin(a::DD)
    isneg, abs_a = signbit(a.hi), abs(a)
    if abs_a.hi > 1.0
        throw(ErrorException("asin: Argument out of domain."))
        return dd_NaN
     elseif abs_a == one(DD)
        return (a.hi >= zero(Float64) ? dd_pi_over_2 : -dd_pi_over_2)
    end

    if (abs_a.hi <= 2.46552e-32)
        b = abs_a
    else
        #a = DD(asin(TD(abs_x)))
        y = TD(abs_a)
        x = y*y
        x = one(TD)-x
        x = sqrt(x)
        x = one(TD) + x
        z = atan2(y,x)
        z = mulby2(z)
        b = DD(z)
    end
    isneg ? -b : b
end

function acos(a::DD)
    isneg, abs_a = signbit(a.hi), abs(a)
    if abs_a.hi > 1.0
      throw(ErrorException("acos: Argument out of domain."))
      return dd_NaN
    elseif abs_a == one(DD)
      return (a.hi >= zero(Float64) ? dd_pi : zero(DD))
    elseif abs_a == zero(DD)
      return dd_pi_over_2
    end

    dd_pi_over_2 - asin(a)
end

function atan2(y::DD,x::DD)
    DD(atan2(TD(y),TD(x)))
end

function atan(x::DD)
    isneg, abs_x = signbit(x), abs(x)
    if (abs_x.hi <= 4.19745e-11)
        t = abs_x
    elseif abs_x.hi <= 0.003
        t = atan_series(abs_x)
    else
        t = atan2(abs_x, one(DD))
    end    
    isneg ? -t : t
end

# ok for |x| >= 1.005
# !!FIXME for 1 <= |x| < 1.005
function acsc(x::DD)
    isneg, abs_x = signbit(x), abs(x)
    if abs_a.hi < 1.0
      throw(ErrorException("acsc: Argument out of domain."))
      return dd_NaN
    end
    ac = asin(1.0/(-abs_x))
    isneg ? ac : -ac
end    
       
asec(x::DD) = acos(1.0/x)
acot(x::DD) = flipsign(dd_pi_over_2,x.hi) - atan(1.0/x)
