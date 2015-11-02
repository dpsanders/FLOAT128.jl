@inline function mod2piAsTD(a::DD)
    b = TD(a)
    b1 = b * td_1_over_twopi
    f1 = floor(b1)
    chk = floor(b.hi * td_1_over_twopi_part2.hi)
    if (chk > zero(Float64))
       throw(ErrorException("arg ($a) is too large"))
    end
    d1 = f1*td_twopi
    d2 = f1*td_twopi_part2
  # d3 = f1*td_twopi_part3

    m = b - d1
    m = m - d2
  # m = m - d3

    m
end

function mod2pi(a::DD)
    DD(mod2piAsTD(a))
end

@inline function modPiOver2AsTD(a::DD)
    b = TD(a)
    b1 = b * td_2_over_pi
    f1 = floor(b1)
    chk = floor(b.hi * td_2_over_pi_part2.hi)
    if (chk > zero(Float64))
       throw(ErrorException("arg ($a) is too large"))
    end
    d1 = f1*td_pi_over_2
    d2 = f1*td_pi_over_2_part2
    # d3 = f1*td_pi_over_2_part3

    m = b - d1
    m = m - d2
    # m = m - d3

    m
end

# a --> [0 .. p/2)
function modPiOver2(a::DD)
    DD(modPiOver2AsTD(a))
end


@inline function modPiOver4AsTD(a::DD)
    b = TD(a)
    b1 = b * td_4_over_pi
    f1 = floor(b1)
    chk = floor(b.hi * td_4_over_pi_part2.hi)
    if (chk > zero(Float64))
       throw(ErrorException("arg ($a) is too large"))
    end
    d1 = f1*td_pi_over_4
    d2 = f1*td_pi_over_4_part2
    # d3 = f1*td_pi_over_4_part3

    m = b - d1
    m = m - d2
    # m = m - d3

    m
end

# a --> [0 .. p/2)
function modPiOver4(a::DD)
    DD(modPiOver4AsTD(a))
end

# a --> [-pi/4 .. p/4)
function modSignedPiOver4(a::DD)
    m = modPiOver2AsTD(a)
    if (m.hi >= pi/4.0)
        m = m - td_pi_over_2
        m = m - td_pi_over_2_part2
        # m = m - td_pi_over_2_part3
    end

    DD(m)
end
#=
more accurate
=#

function modSignedPiOver4(a::DD)
    m = modPiOver2AsTD(a)
    if (m.hi >= pi/4.0)
        m = m - td_pi_over_2
        m = m - td_pi_over_2_part2
        # m = m - td_pi_over_2_part3
    end

    DD(m)
end


# -2pi <= radian < 2pi
function resolveQuadrantInCircle(radian::DD)
    isneg, aradian = signbit(radian), abs(radian)
    if isneg
        aradian = DD(td_twopi - aradian)
    end
    # octant numbering for trig calcs are rotated by -pi/4
    # (see, e.g. Jack Crenshaw, Math Toolkit for Realtime Programming, p106)
    octant = trunc(Int, div(radian, dd_pi_over_4).hi)
    quadrant = (octant + isodd(octant)) >> 1
    quadrant = mod(quadrant,4)
    radi = modSignedPiOver4(aradian)
    quadrant, radi
end


const dd_n_tan_coeff = 16;
const dd_tan_coeff = DD[
   DD(1.0, 0.0),
   DD(0.3333333333333333, 1.850371707708594e-17),
   DD(0.13333333333333333, 1.8503717077085942e-18),
   DD(0.05396825396825397, -2.5552752154071065e-18),
   DD(0.021869488536155203, -1.7377829530067485e-19),
   DD(0.008863235529902197, -7.63300580171831e-19),
   DD(0.003592128036572481, -1.253823608406629e-19),
   DD(0.0014558343870513183, -6.214492640136062e-20),
   DD(0.000590027440945586, 3.478690842383652e-20),
   DD(0.00023912911424355248, 3.564613898329782e-21),
   DD(9.691537956929451e-5, -6.2386628755632464e-21),
   DD(3.927832388331683e-5, 1.3737015743076767e-21),
   DD(1.5918905069328964e-5, 1.0427554807190543e-21),
   DD(6.451689215655431e-6, 1.1519922496640058e-22),
   DD(2.6147711512907546e-6, -9.313685621299801e-23),
   DD(1.0597268320104654e-6, 2.3670525505213632e-24), # 31
];


function sin_taylor_series(a::DD)
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

function cos_taylor_series(a::DD)
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



# for a < 9/64 0.140625
function tan_taylor_series(radian::DD)
  x = radian
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

  z = TD(x) + x*(dd_tan_coeff[2]*x2 + dd_tan_coeff[3]*x4 + dd_tan_coeff[4]*x6)
  z2 = x9 * (dd_tan_coeff[5] + x2*dd_tan_coeff[6] + x4*dd_tan_coeff[7] + x6*dd_tan_coeff[8])
  z3 = x17 * (dd_tan_coeff[9] + x2*dd_tan_coeff[10] + x4*dd_tan_coeff[11] + x6*dd_tan_coeff[12])
  z4 = x25 * (dd_tan_coeff[13] + x2*dd_tan_coeff[14] + x4*dd_tan_coeff[15] + x6*dd_tan_coeff[16])

  DD(z + ((z4+z3)+z2))
end



function sin_taylor(radian::DD)
    hi = radian.hi
    if ((hi==zero(Float64)) | (hi==NaN) | (hi==pi/4.0) | (hi==pi/6.0))
        if hi == zero(Float64)
            return zero(DD)
        elseif hi == NaN
            return radian
        elseif hi == dd_pi_over_4.hi && radian.lo == dd_pi_over_4.lo
            return dd_1_over_sqrt2
        elseif hi == dd_pi_over_6.hi && radian.lo == dd_pi_over_6.lo
            return dd_half
        end
    end
    sin_taylor_series(radian)
end

function cos_taylor(radian::DD)
    hi = radian.hi
    if ((hi==zero(Float64)) | (hi==NaN) | (hi==pi/4.0) | (hi==pi/6.0))
        if hi == zero(Float64)
            return one(DD)
        elseif hi == NaN
            return radian
        elseif hi == dd_pi_over_4.hi && radian.lo == dd_pi_over_4.lo
            return dd_1_over_sqrt2
        elseif hi == dd_pi_over_6.hi && radian.lo == dd_pi_over_6.lo
            return dd_sqrt3_over_2
        end
    end
    cos_taylor_series(radian)
end


function sincos_taylor(radian::DD)
    s = c = zero(DD)
    if abs(radian.hi) < 0.39269908169872414 # pi/8
        s = sin_taylor(radian)
        ts = TD(s)
        c = DD(sqrt(one(TD)-ts*ts))
    else
        c = cos_taylor(radian)
        tc = TD(c)
        s = DD(sqrt(one(TD)-tc*tc))
    end
    s,c
end


# a <= 0.15
function tan_taylor(radian::DD)
    hi = radian.hi
    if ((hi==zero(Float64)) | (hi==NaN) | (hi==pi/4.0) | (hi==pi/6.0))
        if hi == zero(Float64)
            return zero(DD)
        elseif hi == NaN
            return radian
        elseif hi == dd_pi_over_4.hi && radian.lo == dd_pi_over_4.lo
            return one(DD)
        elseif hi == dd_pi_over_6.hi && radian.lo == dd_pi_over_6.lo
            return dd_sqrt3_over_3
        end
    end
    isneg, aradian = signbit(radian), abs(radian)
    if (1.0e-15 <= aradian.hi <= 9.0/64.0)
        t = tan_taylor_series(aradian)
    else
        s,c = sincos_taylor(aradian)
        t = s/c
    end
    isneg ? -t : t
end




#=
   for radian in [0,8*pi)
      relerr precise sin(radian) ~1.2e-32 (106 bits)
   for radian in [0,4096^3 *pi)
      relerr precise sin(radian) ~2e-32 (105 bits)
   for radian > 4096^3 *pi
      relerr 1.5e-29 (95 bits)
=#

function sinInCircle(radian::DD)
  quadrant, radi = resolveQuadrantInCircle(radian)
  isneg, aradi = signbit(radi), abs(radi)
  if (quadrant == 0)
    isneg ? -sin_taylor(aradi) : sin_taylor(aradi)
  elseif (quadrant == 1)
    cos_taylor(aradi)
  elseif (quadrant == 2)
    isneg ? sin_taylor(aradi) : -sin_taylor(aradi)
  else #if (quadrant == 3)
    -cos_taylor(aradi) # isneg ? -cos_taylor(aradi) : -cos_taylor(aradi)
  end
end

function sin(radian::DD)
  isneg, aradian = signbit(radian.hi), abs(radian)
  if aradian <= dd_pi_over_4
      s = sin_taylor(aradian)
   elseif aradian < dd_twopi
      s = sinInCircle(aradian)
   else
      r = mod2piAsTD(aradian)
      ddr = DD(r.hi,r.md)
      s = sinInCircle(ddr)*cos(r.lo)+sin(r.lo)*cos(ddr)
    end
    isneg ? -s : s
end

#=
  for radian in [0,8*pi)
      relerr precise sin(radian) ~1.2e-32 (106 bits)
   for radian in [0,4096^3 *pi)
      relerr precise sin(radian) ~2e-32 (105 bits)
   for radian > 4096^3 *pi
      relerr 1.5e-29 (95 bits)
=#
function cosInCircle(radian::DD)
  quadrant, radi = resolveQuadrantInCircle(radian)
  isneg, aradi = signbit(radi), abs(radi)
  if (quadrant == 0)
    cos_taylor(aradi)
  elseif (quadrant == 1)
    isneg ? sin_taylor(aradi) : -sin_taylor(aradi)
  elseif (quadrant == 2)
    -cos_taylor(aradi)
  else #if (quadrant == 3)
    isneg ? -sin_taylor(aradi) : sin_taylor(aradi)
  end
end

function cos(radian::DD)
   aradian = abs(radian)
   if aradian <= dd_pi_over_4
      cos_taylor(radian)
   elseif aradian < dd_twopi
      cosInCircle(aradian)
   else
      r = mod2piAsTD(aradian)
      ddr = DD(r.hi,r.md)
      cosInCircle(ddr)*cos(r.lo)-sinInCircle(ddr)*sin(r.lo)
    end
end



#=
function sincos(radian::DD)
  if abs(radian) <= dd_pi_over_4
      s = sin_taylor(radian)
      c = DD(sqrt(one(DD)-TD(s)*s))
   elseif abs(radian) < dd_twopi
      s = sinInCircle(radian)
      c = DD(sqrt(one(DD)-TD(s)*s))
      if (abs(radian) > dd_pi)
          c = -c
      end
   else
      r = mod2piAsTD(radian)
      ddr = DD(r.hi,r.md)
      s = sinInCircle(ddr)
      slo = sin(r.lo)
      c = cosInCircle(ddr)
      clo = cos(r.lo)
      s = s*clo + c*slo
      c = c*clo - s*slo
      if (abs(radian) > dd_pi)
          c = -c
      end
   end
   s,c
end
=#

function sincos(radian::DD)
   sin(radian),cos(radian)
end



# good for radian in 0..2pi(quadrants 0,1, 1/2 of q3)
#=
   for radian in [0,2*pi)
      relerr tan(radian) ~ 5e-3 (106 bits)
=#
function tanInCircle(radian::DD)
  quadrant, radi = resolveQuadrantInCircle(radian)
  isneg, aradi = signbit(radi), abs(radi)
  t = tan_taylor(aradi)
  if (quadrant == 0)
    isneg ? -t : t
  elseif (quadrant == 1)
    isneg ? recip(t) : -recip(t)
  elseif (quadrant == 2)
    isneg ? -t : t
  else #if (quadrant == 3)
    isneg ? recip(t) : -recip(t)
  end
end

#=
   for radian in [-8pi,8pi)
      relerr tan(radian) ~4.1e-32 (104 bits)
   for radian in [8pi,4096^3 *pi)
      relerr tan(radian) ~6.2e-32 (103 bits)
   for radian > 4096^3 *pi
      relerr 4.2e-29 (94 bits)
=#
#=
function tan(radian::DD)
   isneg, aradian = signbit(radian), abs(radian)
   if abs(radian.hi) < 9.0/64.0
      t = tan_taylor(aradian)
   elseif abs(radian < dd_twopi)
      t = tanInCircle(aradian)
   else
      r = mod2piAsTD(aradian)
      ddr = DD(r.hi,r.md)
      s = TD(sinInCircle(ddr))
      c = TD(cosInCircle(ddr))
      slo,clo = sin(r.lo),cos(r.lo)
      s = s*clo+c*slo; c=c*clo-s*slo;
      t = DD(s/c)
    end
    isneg ? -t : t
end
=#

function tan(radian::DD)
  isneg, aradian = signbit(radian.hi), abs(radian)
  s = c = zero(DD)
  if aradian <= dd_pi_over_4
      s = sin_taylor(aradian)
      c = cos_taylor(aradian)
   elseif aradian < dd_twopi
      s = sinInCircle(aradian)
      c = cosInCircle(aradian)
   else
      r = mod2piAsTD(aradian)
      ddr = DD(r.hi,r.md)
      s = sinInCircle(ddr)
      slo = sin(r.lo)
      c = cosInCircle(ddr)
      clo = cos(r.lo)
      s = s*clo+c*slo
      c = c*clo-s*slo
    end
    isneg ? -s/c : s/c
end

csc(radian::DD) = recip( sin(radian) )
sec(radian::DD) = recip( cos(radian) )
cot(radian::DD) = recip( tan(radian) )

