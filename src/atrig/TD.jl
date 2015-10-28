
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

