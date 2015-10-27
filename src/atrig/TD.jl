
function sincos(x::TD)
  s,c=sincos(DD(x))
  slo,clo = TD(sin(x.lo)),TD(cos(x.lo))
  s*clo+c*slo, c*clo-s*slo
end

# use Newton's iteration to solve
#    sin(z) = y/r    (1)
#    cos(z) = x/r    (2)
#  where r = sqrt(x^2 + y^2).
# The iteration is given by
#
#    z' = z + (y - sin(z)) / cos(z)          (for equation 1)
#    z' = z - (x - cos(z)) / sin(z)          (for equation 2)
#
#  Here, x and y are normalized so that x^2 + y^2 = 1.
#  If |x| > |y|, then first iteration is used since the
#  denominator is larger.  Otherwise, the second is used.


function atan2{T<:TD}(y::T,x::T)
  if x.hi == zero(Float64)
    if y.hi == zero(Float64)
      # Both x and y are zero
      throw(ErrorException("atan2: Both arguments zero."))
      return td_NaN
    end
    return (y.hi > zero(Float64) ? td_pi_over_2 : -td_pi_over_2)
  elseif y.hi == zero(Float64)
    return (x.hi > zero(Float64) ? zero(TD) : td_pi)
  end

  if (x == y)
    return (y.hi > zero(Float64)  ? td_pi_over_4 : -td_3pi_over_4)
  elseif (x == -y)
    return (y.hi > zero(Float64) ? td_3pi_over_4 : -td_pi_over_4)
  end

  r = sqrt(sqr(x) + sqr(y))
  xx = x / r
  yy = y / r

  # Compute double precision approximation to atan.
  z = TD(atan2(y.hi,x.hi))

  if (abs(xx.hi) > abs(yy.hi))
    # Use Newton iteration 1.  z' = z + (y - sin(z)) / cos(z)
    sin_z,cos_z = sincos(z)
    z += (yy - sin_z) / cos_z
    sin_z,cos_z = sincos(z)
    z += (yy - sin_z) / cos_z
    sin_z,cos_z = sincos(z)
    z += (yy - sin_z) / cos_z
  else
    # Use Newton iteration 2.  z' = z - (x - cos(z)) / sin(z)
    sin_z,cos_z = sincos(z)
    z -= (xx - cos_z) / sin_z;
    sin_z,cos_z = sincos(z)
    z -= (xx - cos_z) / sin_z;
    sin_z,cos_z = sincos(z)
    z -= (xx - cos_z) / sin_z;
  end

  z
end

atan(x::TD) = atan2(x,one(TD))

#=
  asin(x) = 2*atan(x / (1+sqrt(1-x*x)))
=#
function asin(a::TD)
    isneg, abs_a = signbit(a.hi), abs(a)
    if abs_a.hi > 1.0
      throw(ErrorException("asin: Argument out of domain."))
      return td_NaN
    elseif abs_a == one(TD)
      return (a.hi >= zero(Float64) ? td_pi_over_2 : -td_pi_over_2)
    end

    y = a
    x = a*a
    x = 1.0-x
    x = sqrt(x)
    x = 1.0 + x
    z = atan2(y,x)
    z = mulby2(z)
    z = abs(z)
    isneg ? -z : z
end

#=
  acos(x) = pi/2 - 2*atan(x / (1+sqrt(1-x*x)))
=#
function acos(a::TD)
    abs_a = abs(a)
    if abs_a > one(TD)
        throw(ErrorException("acos: abs(x) > 1.0"))
        return td_NaN
    elseif abs_a == one(TD)
        return (a.hi >= zero(Float64) ? zero(TD) : td_pi)
    end

#    td_pi_over_2 - asin(a)
    y = a
    x = a*a
    x = 1.0-x
    x = sqrt(x)
    x = 1.0 + x
    z = atan2(y,x)
    z = mulby2(z)
    td_pi_over_2 - z
end

function asinh(x::TD)
    if (x.hi >= zero(Float64))
       x2 = x*x
       z = sqrt(1.0+x2)
       z += x
       log(z)
    else
       z = asinh(-x)
       -z
    end
end

#=
function acosh(x::TD)
    xp1 = x*1.0
    xm1 = x-1.0
    xp1 = sqrt(0.5*xp1)
    xm1 = sqrt(0.5*xm1)
    2.0 * log(xp1+xm1)
end
=#
function acosh(x::TD)
    if x < one(TD)
      throw(ErrorException("acosh: x < 1.0"))
      return td_NaN
    elseif abs(x) < TD(1.0, 9.363352709384397e-97, 0.0)
      return zero(TD)
    end

    x2 = sqr(x)
    x2 = x2 - one(TD)
    x2 = sqrt(x2)
    a = -log(x - x2)
    b = log(x + x2)
    divby2(a+b)
end

function atanh(x::TD)
    if abs(x) > one(TD)
      throw(ErrorException("atanh: abs(x) > 1.0"))
      return td_NaN
    elseif abs(x) > TD(1.0,-9.363352709384397e-97, 0.0)
      return copysign(tn_Inf,x.hi)
    end

    onepx = 1.0+x
    onemx = 1.0-x
    onepx = log(onepx)
    onemx = log(onemx)
    0.5 * (onepx-onemx)
end

