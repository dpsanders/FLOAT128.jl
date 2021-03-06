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
   DD(0.9929146630650226, 1.896947743258429e-18) +
   (DD(0.12687948850312847, -5.839018909578105e-18) +
   (DD(1.3614681807580766, -4.280666382785087e-17) +
   (DD(0.1412487658121559, 6.491929668608193e-18) +
   (DD(0.4577547489764742, -2.5746404954730867e-18) +
   (DD(0.03037246208354068, 9.645477620436761e-19) +
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



function atan4to5numer(x::DD)
DD(3.4986867510923863e-19, -4.80043199494538e-37) +
    (DD(0.9773577630727862, 4.222541831019509e-17) +
    (DD(0.224611377881578, -1.1083500370099093e-17) +
    (DD(1.0308599673016396, 2.2802867037807203e-17) +
    (DD(0.17625799155210703, 1.3836582834436884e-17) +
    (DD(0.20805322134678086, 6.198168426557029e-18) +
    DD(0.015630212957401312, 1.880151741876803e-19)*x) *x)*x)*x)*x)*x
end
function atan4to5denom(x::DD)
    DD(0.9773577630727863, -3.6247738309496816e-18) +
    (DD(0.2246113778815724, -5.126755876468263e-18) +
    (DD(1.3566458883261958, -3.8674444538119656e-17) +
    (DD(0.2511284508354467, -2.5250079300540653e-17) +
    (DD(0.4647969651117322, 2.5993613189837385e-17) +
    (DD(0.05441741583246426, 1.5779355560437385e-18) +
    DD(0.023225755314791163, 8.965019269474856e-19)*x) *x)*x)*x)*x)*x
end
@inline atan4to5of64(x::DD) = atan4to5numer(x) / atan4to5denom(x)

function atan5to6numer(x::DD)
    DD(4.724665895731685e-18, 3.197057622132749e-35) +
    (DD(0.9664297227423516, 4.7174145325688815e-17) +
    (DD(0.27128902188603604, -7.050713460875651e-18) +
    (DD(1.0309596465202269, 1.3983478441412103e-17) +
    (DD(0.21381562993695888, 7.849093812371462e-20) +
    (DD(0.21190211439467507, -1.3818021586137387e-17) +
    DD(0.01912954892540732, 1.2767502839033025e-18)*x) *x)*x)*x)*x)*x
end
function atan5to6denom(x::DD)
    DD(0.9664297227423524, -1.0951441920954721e-17) +
    (DD(0.2712890218859855, 1.419338629146824e-17) +
    (DD(1.3531028874365107, 6.755197286060292e-17) +
    (DD(0.3042453038355929, 2.6438459401187718e-17) +
    (DD(0.46965046699213003, -2.035931360386614e-17) +
    (DD(0.06628682504089098, 5.214516757590902e-18) +
    DD(0.02399120987040229, -2.370266164091456e-20)*x) *x)*x)*x)*x)*x
end
@inline atan5to6of64(x::DD) = atan5to6numer(x) / atan5to6denom(x)

#=
function atan6to7numer(x::DD)
    DD(4.08258884792285e-17, 1.9801226085996085e-33) +
    (DD(0.9535145413823523, -5.236875025807137e-19) +
    (DD(0.3160973631452324, 1.2845795833188462e-17) +
    (DD(1.0309096536691162, -7.345625130452875e-17) +
    (DD(0.25042597049018966, 7.039911798136601e-18) +
    (DD(0.21640553613902666, 1.2514546576627358e-18) +
    DD(0.022641480101284518, -1.4107627307355407e-18)*x) *x)*x)*x)*x)*x
end
function atan6to7denom(x::DD)
    DD(0.9535145413823575, 3.63393809486893e-17) +
    (DD(0.31609736314492026, -7.114941540056584e-18) +
    (DD(1.3487478341412331, 1.244385774499147e-17) +
    (DD(0.35579175792476014, 2.5929154529266115e-17) +
    (DD(0.4752852442377697, 1.9402468748279948e-17) +
    (DD(0.07801919424133187, -2.3738530590672923e-19) +
    DD(0.024895861851288902, -5.072156320295304e-19)*x) *x)*x)*x)*x)*x
end
@inline atan6to7of64(x::DD) = atan6to7numer(x) / atan6to7denom(x)
=#

function atan6to7numer(x::DD)
    DD(-1.0415085111722597e-19, -5.395102486091432e-36) +
    (DD(0.9462153358439911, -2.831401300378913e-17) +
    (DD(0.36126832267059383, 7.66314484633861e-18) +
    (DD(1.2775709776413058, -4.5201332076330683e-17) +
    (DD(0.3793593264332555, -1.9726308973556664e-17) +
    (DD(0.4252136862797208, 1.5736094412320355e-17) +
    (DD(0.07602702590216374, 3.85069313300002e-18) +
    DD(0.020647468580854707, 1.5575316780473585e-18)*x) *x)*x)*x)*x)*x)*x
end
function atan6to7denom(x::DD)
    DD(0.9462153358439911, -4.3782190107136284e-17) +
    (DD(0.3612683226705949, 2.496531220033659e-17) +
    (DD(1.5929760895892568, -8.615747353570366e-17) +
    (DD(0.4997821006581529, 2.632573078802357e-18) +
    (DD(0.7669626489442688, 1.5683319304926968e-17) +
    (DD(0.17036739541120313, 1.3669564178288488e-17) +
    (DD(0.09288008014809244, -3.2854050411161694e-19) +
    DD(0.008442533659821153, -4.3713113100221436e-21)*x) *x)*x)*x)*x)*x)*x
end
@inline atan6to7of64(x::DD) = atan6to7numer(x) / atan6to7denom(x)


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
