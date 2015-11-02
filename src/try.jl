# Computes sin(a) and cos(a) using Taylor series.
#   Assumes |a| <= pi/2048

function sincos_taylor(a::DD)
  thresh = 0.5 * eps(a.hi)
  a=p=s=t=x=zero(DD)
  sin_a = zero(DD)
  cos_a = one(DD)

  if (a==zero(DD)) {
     return sin_a, cos_a
  }

  x = -sqr(a)
  s = a
  p = a
  i = 1
  
  p *= x
  t = p * dd_inv_fact[i]
  s += t
  i += 2
  
  while (i < dd_n_inv_fact && t.hi > thresh)
    p *= x
    t = p * inv_fact[i]
    s += t
    i += 2
  end
  
  sin_a = s;
  cos_a = sqrt(1.0 - sqr(s));
  sin_a, cos_a
end


function sincos(void sincos(a::DD)

  sin_a = zero(DD)
  cos_a = one(DD)
  if (a==zero(DD))
      return sin_a, cos_a
  end

  # approximately reduce by 2*pi
  z = trunc(a / dd_twopi)
  t = a - dd_twopi * z

  # approximately reduce by pi/2 and then by pi/1024.
  q = (t.hi / dd_pi_over_2.hi) + 0.5
  t -= dd_pi_over_2 * q
  j = trunc(Int,q)
  q = floor((t.hi / dd_pi_over_1024.hi) + 0.5)
  t -= dd_pi_over_1024 * q
  k = trunc(Int, q)
  abs_k = abs(k)

  if (j < -2 || j > 2)
    ErrorException("sincos: Cannot reduce modulo pi/2.")
    cos_a = sin_a = dd_NaN;
    return sin_a, cos_a
  end

  if (abs_k > 256)
    ErrorException("sincos: Cannot reduce modulo pi/1024.")
    cos_a = sin_a = dd_NaN;
    return sin_a, cos_a
  end

  sin_t, cos_t = sincos_taylor(t)

  if (k == 0)
    #if (j == 0)
    #  sin_a = sin_t
    #  cos_a = cos_t
    if (j == 1)
      sin_t, cos_t = cos_t, -sin_t
    elseif (j == -1)
      sin_t, cos_t = -cos_t, sin_t
    elseif (j != 0)
      sin_t,cos_t = -sin_t, -cos_t
    end
    return sin_t, cos_t
  end

  u = cos_table[abs_k]
  v = sin_table[abs_k]

  if (j == 0) {
    if (k > 0) {
      sin_a = u * sin_t + v * cos_t;
      cos_a = u * cos_t - v * sin_t;
    } else {
      sin_a = u * sin_t - v * cos_t;
      cos_a = u * cos_t + v * sin_t;
    }
  } else if (j == 1) {
    if (k > 0) {
      cos_a = - u * sin_t - v * cos_t;
      sin_a = u * cos_t - v * sin_t;
    } else {
      cos_a = v * cos_t - u * sin_t;
      sin_a = u * cos_t + v * sin_t;
    }
  } else if (j == -1) {
    if (k > 0) {
      cos_a = u * sin_t + v * cos_t;
      sin_a =  v * sin_t - u * cos_t;
    } else {
      cos_a = u * sin_t - v * cos_t;
      sin_a = - u * cos_t - v * sin_t;
    }
  } else {
    if (k > 0) {
      sin_a = - u * sin_t - v * cos_t;
      cos_a = v * sin_t - u * cos_t;
    } else {
      sin_a = v * cos_t - u * sin_t;
      cos_a = - u * cos_t - v * sin_t;
    }
  }
}
