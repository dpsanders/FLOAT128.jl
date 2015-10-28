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

# for 1 <= |x| < 1+1/320 = 1.003125
function acscNear1(x::DD)
    dz = abs(x) - 1.0 
    sdz = sqrt(dz)
    s2  = FLOAT128.dd_sqrt2
    p2  = FLOAT128.dd_pi_over_2
    a = acscNear1inner(dz)
    b = a * sdz
    b = b * s2
    c = b / 28745297418623385600.0
    c = -c
    d = c + p2
    d
end



function acscNear1inner(z::DD)
    28745297418623385600+(-11977207257759744000+(7725298681255034880+(-5678479512384307200+(4471178803124633600+(-3678333901524172800+(3120230290312396800+(-2707335706710835200+(2390122732757760000+(-2139078913747776000+(1935592642417262400+(-1767382875804546000+1626037491584404899*dz)*dz)*dz)*dz)*dz)*dz)*dz)*dz)*dz)*dz)*dz)*dz
end


# recheck for 1.003125 <= |x| < 1.005
function acsc(x::DD)
    isneg, abs_x = signbit(x), abs(x)
    if abs_x.hi < 1.0
      throw(ErrorException("acsc: Argument out of domain."))
      return dd_NaN
    end
    # if abs_x.hi <= 1.003125 # 1+1/320
    #   ac = acscNear1(abs_x)
    if abs_x.hi > 1.005   
       ac = asin(1.0/abs_x)
    else
       ac = acscNear1(abs_x)
    end
    isneg ? -ac : ac
end    
   

function asecNear1(z::DD)
    dz = z - one(DD)
    v = 1.0+(-DD(5)/DD(12)+(DD(43)/DD(160)+(-DD(177)/DD(896)+(DD(2867)/DD(18432)+(-DD(11531)/DD(90112)+(DD(92479)/DD(851968)+(-DD(74069)/DD(786432)+(DD(11857475)/DD(142606336)+(-DD(47442055)/DD(637534208)+(DD(126527543)/DD(1879048192)+(-DD(1518418695)/DD(24696061952)+(DD(24295375159)/DD(429496729600)-DD(97182800711)*dz*(DD(1.0)/DD(1855425871872.0)))*dz)*dz)*dz)*dz)*dz)*dz)*dz)*dz)*dz)*dz)*dz)*dz
    w = v * dd_sqrt2
    w = w *sqrt(dz)
    w
end

function asec(x::DD)
    isneg, abs_x = signbit(x), abs(x)
    if abs_x.hi < 1.0
      throw(ErrorException("acsc: Argument out of domain."))
      return dd_NaN
    end
    if abs_x.hi > 1.005   
       ac = acos(1.0/abs_x)
    else
       ac = asecNear1(abs_x)
    end
    isneg ? -ac : ac
end    

 
acot(x::DD) = flipsign(dd_pi_over_2,x.hi) - atan(1.0/x)
