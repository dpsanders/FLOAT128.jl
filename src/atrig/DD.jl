# x=0..0.015
function asinconfrac(x::DD)
    x2 = x*x
    #x/(1+x2/(-6+x2/(10/17+x2/(-4046/549+x2/(602802/1173833+x2/(-104890816822/12889271025+1207527848162057*x2*(1/573251714022650)))))))
    x/(DD(1.0,0.0)+x2/(-DD(6.0,0.0)+x2/(DD(0.5882352941176471, -1.959217102279688e-17)+x2/(-DD(7.36976320582878, -4.044528322860315e-17)+x2/(DD(0.51353301534375, -3.4384261563723567e-17)+x2/(DD(-8.137839340840458, 8.011797547325665e-16)+x2*(DD((2.1064530966484747, 1.1717219329569228e-16)))))))))
end
      


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
    if abs_a.hi > one(Float64)
        throw(ErrorException("asin: Argument out of domain."))
        return dd_NaN
     elseif abs_a.hi == one(Float64) && abs_a.lo == zero(Float64)
        return (a.hi >= zero(Float64) ? dd_pi_over_2 : -dd_pi_over_2)
    end

    if (abs_a.hi <= 0.015)
        b = asinconfrac(abs_a)
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
    elseif abs_a.hi == one(Float64) && abs_a.lo == zero(Float64)
      return (a.hi >= zero(Float64) ? dd_pi : zero(DD))
    elseif abs_a.hi == zero(Float64) && abs.lo == zero(Float64)
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


#=
   arctan(x) for x in k/64 .. (k+1)/64 (k in 0..63)
   these are minimax results from maple over
     k/64 - 1/4096^2 .. (k+1)/64 + 1/4096^2
   they are not the best fit to DD coeffs
   but they are close enough
=#

function atan01numer(x::DD)
    DD(3.629893017141765e-35, 6.6406679896317375e-52) +
    (DD(0.9996754545703865, -4.598817090918266e-17) +
    (DD(0.02555812681060226, -9.958411852257287e-20) +
    (DD(1.0302585588731008, -7.183792179616217e-17) +
    (DD(0.01988046865875826, 7.896880383011565e-19) +
    (DD(0.2000869562087623, -8.47084451019246e-18) +
    DD(0.0017314342719200748, -7.732543321053218e-20) * x) *x)*x)*x)*x)*x
end
function atan01denom(x::DD)
    DD(0.9929146630650226, 1.896947743258429e-18),
    (DD(0.12687948850312847, -5.839018909578105e-18),
    (DD(1.3614681807580766, -4.280666382785087e-17),
    (DD(0.1412487658121559, 6.491929668608193e-18),
    (DD(0.4577547489764742, -2.5746404954730867e-18),
    (DD(0.03037246208354068, 9.645477620436761e-19),
    DD(0.02213623421350611, -1.1027549731601916e-18) * x) *x)*x)*x)*x)*x
end
@inline atan01of64(x::DD) = atan01numer(x) / atan01denom(x)

function atan12numer(x::DD)
    DD(1.6218277873895332e-25, 1.084747206314063e-41) +
    (DD(0.9974155586768761, 4.514103898642169e-18) + 
    (DD(0.0764916212253432, 3.500531525408491e-18) +
    (DD(1.0303435035210617, -6.890589883420933e-17) +
    (DD(0.05955192232775267, 2.077561630841645e-18) +
    (DD(0.20089995784387976, 4.140614508164837e-18) +
    DD(0.0051959782853776715, -1.158432657223438e-19)*x) *x)*x)*x)*x)*x
end
function atan12denom(x::DD)
    DD(0.9929146630650226, 1.896947743258429e-18) +
    (DD(0.12687948850312847, -5.839018909578105e-18) +
    (DD(1.3614681807580766, -4.280666382785087e-17) +
    (DD(0.1412487658121559, 6.491929668608193e-18) +
    (DD(0.4577547489764742, -2.5746404954730867e-18) +
    (DD(0.03037246208354068, 9.645477620436761e-19) +
    DD(0.02213623421350611, -1.1027549731601916e-18)*x) *x)*x)*x)*x)*x
end
@inline atan12of64(x::DD) = atan12numer(x) / atan12denom(x)


function atan2to3numer(x::DD)
   DD(1.5746491607703442e-22, -5.35357483834068e-39) +
   (DD(0.9929146630650226, 1.8434410645424224e-18) +
   (DD(0.12687948850312847, 2.538312664799496e-18) +
   (DD(1.0304966264030682, 8.1565115612327e-17) +
   (DD(0.09895560297783178, -4.74221417590191e-18) +
   (DD(0.20251495466768876, 4.4403079300469404e-18) +
   DD(0.00866543793078847, -6.0102146731863325e-19) * x) *x)*x)*x)*x)*x
end
function atan2to3denom(x::DD)
   DD(0.9929146630650226, 1.896947743258429e-18),
   (DD(0.12687948850312847, -5.839018909578105e-18),
   (DD(1.3614681807580766, -4.280666382785087e-17),
   (DD(0.1412487658121559, 6.491929668608193e-18),
   (DD(0.4577547489764742, -2.5746404954730867e-18),
   (DD(0.03037246208354068, 9.645477620436761e-19),
   DD(0.02213623421350611, -1.1027549731601916e-18)*x) *x)*x)*x)*x)*x
end
@inline atan2to3of64(x::DD) = atan2to3numer(x) / atan2to3denom(x)


function atan3to4numer(x::DD)
    DD(1.3187186340364065e-20, -6.470728923559978e-37) +
    (DD(0.9862102481582604, 2.7745960767391508e-17) +
    (DD(0.17636675197972745, 2.0310811075114496e-18) +
    (DD(1.0306848913144344, -2.8475702993687415e-17) +
    (DD(0.1379152536552334, -9.407544146202368e-18) +
    (DD(0.204910126648868, -5.721917057878642e-18) +
    DD(0.012142724623342087, 8.095456557657355e-19) * x) *x)*x)*x)*x)*x
end
function atan3to4denom(x::DD)
    DD(0.9862102481582604, 3.091568007464899e-17) +
    (DD(0.1763667519797271, 1.1498601703642871e-17) +
    (DD(1.359421640700545, -6.950382430252179e-17) +
    (DD(0.19670417098071252, 5.430012207236689e-19) +
    (DD(0.4608086239537724, -3.5080079009821014e-18) +
    (DD(0.04243743032802491, 2.7291468284337927e-18) +
    DD(0.022605741241803316, 6.327981810377862e-19)*x) *x)*x)*x)*x)*x
end
@inline atan3to4of64(x::DD) = atan3to4numer(x) / atan3to4denom(x)



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
    isneg ? DD(td_pi-ac) : ac
end    

 
acot(x::DD) = flipsign(dd_pi_over_2,x.hi) - atan(1.0/x)
