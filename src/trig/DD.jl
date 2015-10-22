function sin_taylor(a::DD)
  x = a
  x2 = x*x
  x3 = x*x2
  x4 = x2*x2
  x5 = x2*x3
  x6 = x3*x3
  x7 = x3*x4
  x8 = x4*x4
  x9 = x4*x5
  x17 = x8*x9
  x25 = x17*x8

  z = TD(x) + x*(-dd_inv_fact[3]*x2 + dd_inv_fact[5]*x4 - dd_inv_fact[7]*x6)
  z2 = x9 * (dd_inv_fact[9] - x2*dd_inv_fact[11] + x4*dd_inv_fact[13] - x6*dd_inv_fact[15])
  z3 = x17 * (dd_inv_fact[17] - x2*dd_inv_fact[19] + x4*dd_inv_fact[21] - x6*dd_inv_fact[23])
  z4 = x25 * (dd_inv_fact[25] - x2*dd_inv_fact[27] + x4*dd_inv_fact[29] - x6*dd_inv_fact[31])

  DD(z + ((z4+z3)+z2))
end

function cos_taylor(a::DD)
  x = a
  x2 = x*x
  x4 = x2*x2
  x6 = x2*x4
  x8 = x4*x4
  x16 = x8*x8
  x24 = x8*x16

  z = (-dd_inv_fact[2]*x2 + dd_inv_fact[4]*x4 - dd_inv_fact[6]*x6)
  z2 = x8 * (dd_inv_fact[8] - x2*dd_inv_fact[10] + x4*dd_inv_fact[12] - x6*dd_inv_fact[14])
  z3 = x16 * (dd_inv_fact[16] - x2*dd_inv_fact[18] + x4*dd_inv_fact[20] - x6*dd_inv_fact[22])
  z4 = x24 * (dd_inv_fact[24] - x2*dd_inv_fact[26] + x4*dd_inv_fact[28] - x6*dd_inv_fact[30])

  ((z4+z3)+z2)+z + 1.0
end


function sincos_taylor(a::DD)
    s = c = zero(DD)
    if abs(a.hi) < 0.39269908169872414 # pi/8
        s = sin_taylor(a)
        c = sqrt(1.0-s*s)
    else
        c = cos_taylor(a)
        s = sqrt(1.0-c*c)
    end
    s,c
end

const dd_n_tan_coeff = 31;
const dd_tan_coeff = DD[
   DD(1.0, 0.0),
   DD(0.9, 0.0),
   DD(0.3333333333333333, 1.850371707708594e-17),
   DD(0.0, 0.0),
   DD(0.13333333333333333, 1.8503717077085942e-18),
   DD(0.0, 0.0),
   DD(0.05396825396825397, -2.5552752154071065e-18),
   DD(0.0, 0.0),
   DD(0.021869488536155203, -1.7377829530067485e-19),
   DD(0.0, 0.0),
   DD(0.008863235529902197, -7.63300580171831e-19),
   DD(0.0, 0.0),
   DD(0.003592128036572481, -1.253823608406629e-19),
   DD(0.0, 0.0),
   DD(0.0014558343870513183, -6.214492640136062e-20),
   DD(0.0, 0.0),
   DD(0.000590027440945586, 3.478690842383652e-20),
   DD(0.0, 0.0),
   DD(0.00023912911424355248, 3.564613898329782e-21),
   DD(0.0, 0.0),
   DD(9.691537956929451e-5, -6.2386628755632464e-21),
   DD(0.0, 0.0),
   DD(3.927832388331683e-5, 1.3737015743076767e-21),
   DD(0.0, 0.0),
   DD(1.5918905069328964e-5, 1.0427554807190543e-21),
   DD(0.0, 0.0),
   DD(6.451689215655431e-6, 1.1519922496640058e-22),
   DD(0.0, 0.0),
   DD(2.6147711512907546e-6, -9.313685621299801e-23),
   DD(0.0, 0.0),
   DD(1.0597268320104654e-6, 2.3670525505213632e-24), # 31
];

# for a < 9/64 0.140625
function tan_taylor_series(a::DD)
  x = a
  x2 = x*x
  x3 = x*x2
  x4 = x2*x2
  x5 = x2*x3
  x6 = x3*x3
  x7 = x3*x4
  x8 = x4*x4
  x9 = x4*x5
  x17 = x8*x9
  x25 = x17*x8

  z = TD(x) + x*(dd_tan_coeff[3]*x2 + dd_tan_coeff[5]*x4 + dd_tan_coeff[7]*x6)
  z2 = x9 * (dd_tan_coeff[9] + x2*dd_tan_coeff[11] + x4*dd_tan_coeff[13] + x6*dd_tan_coeff[15])
  z3 = x17 * (dd_tan_coeff[17] + x2*dd_tan_coeff[19] + x4*dd_tan_coeff[21] + x6*dd_tan_coeff[23])
  z4 = x25 * (dd_tan_coeff[25] + x2*dd_tan_coeff[27] + x4*dd_tan_coeff[29] + x6*dd_tan_coeff[31])

  DD(z + ((z4+z3)+z2))
end

# a <= 0.15
function tan_taylor(a::DD)
    if (a.hi <= 9.0/64.0)
        tan_taylor_series(a)
    else
        s,c = sincos_taylor(a)
        s/c
    end
end

function sinpio2(a::DD)
    if a.hi == zero(Float64)
        return zero(DD)
    end
    if (a >= dd_pi_over_4)
        b = dd_pi_over_2 - a
        if (b.hi == zero(Float64))
           one(DD)
        else
           cos_taylor(b)
        end
    else
        sin_taylor(a)
    end
end

function sin02pi(a::DD)
    @assert(a >= zero(DD))
    sgnbit = false
    if (a.hi == zero(Float64)) | (a == dd_twopi)
        return zero(DD)
    end
    b = a
    if a >= dd_pi
        sgnbit = !sgnbit
        a = dd_twopi - a
    end
    # a in [0, pi)
    if a >= dd_pi_over_2
        a = dd_pi - a
    end
    # a in [0,pi/2)
    @assert(a >= zero(DD))
    @assert(a <= dd_pi_over_2)
    s = sinpio2(a)
    sgnbit ? -s : s
end

# !!make it better and put it in arith/DDsupp.jl!!
function mod2pi(a::DD)
    mod(a,dd_twopi)
end

function sin(a::DD)
    sgnbit, aa = signbit(a),abs(a)

    if aa.hi == zero(Float64)
        return aa
    end
    if aa < dd_pi_over_4
        aa = sin_taylor(aa)
    elseif aa >= dd_twopi
        aa = sin02pi( mod2pi(aa) )
    else
        aa = sin02pi(aa)
    end
    sgnbit ? -aa : aa
end
