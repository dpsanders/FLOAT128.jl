const dd_tanh_coeff = DD[
  TD(1.0, 0.0),
  TD(-0.3333333333333333, -1.850371707708594e-17),
  TD(0.13333333333333333, 1.8503717077085942e-18),
  TD(-0.05396825396825397, 2.5552752154071065e-18),
  TD(0.021869488536155203, -1.7377829530067485e-19),
  TD(-0.008863235529902197, 7.63300580171831e-19),
  TD(0.003592128036572481, -1.253823608406629e-19),
  TD(-0.0014558343870513183, 6.214492640136062e-20),
  TD(0.000590027440945586, 3.478690842383652e-20),
  TD(-0.00023912911424355248, -3.564613898329782e-21),
  TD(9.691537956929451e-5, -6.2386628755632464e-21),
  TD(-3.927832388331683e-5, -1.3737015743076767e-21),
  TD(1.5918905069328964e-5, 1.0427554807190543e-21),
  TD(-6.451689215655431e-6, -1.1519922496640058e-22),
  TD(2.6147711512907546e-6, -9.313685621299801e-23),
  TD(-1.0597268320104654e-6, -2.3670525505213632e-24),
  TD(4.294911078273806e-7, 1.1643520863702653e-23)
];

function tanh_taylor_series(radian::DD)
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

  z = TD(x) + x*(dd_tanh_coeff[2]*x2 + dd_tanh_coeff[3]*x4 + dd_tanh_coeff[4]*x6)
  z2 = x9 * (dd_tanh_coeff[5] + x2*dd_tanh_coeff[6] + x4*dd_tanh_coeff[7] + x6*dd_tanh_coeff[8])
  z3 = x17 * (dd_tanh_coeff[9] + x2*dd_tanh_coeff[10] + x4*dd_tanh_coeff[11] + x6*dd_tanh_coeff[12])
  z4 = x25 * (dd_tanh_coeff[13] + x2*dd_tanh_coeff[14] + x4*dd_tanh_coeff[15] + x6*dd_tanh_coeff[16])

  DD(z + ((z4+z3)+z2))
end


@inline function sinhAsTD(a::DD)
  isneg, abs_a = signbit(a), abs(a)
  x = TD(abs_a)
  s = sinh(x)
  isneg ? -s : s
end

sinh(a::DD) = DD(sinhAsTD(a))

@inline function coshAsTD(a::DD)
  x=TD(abs(a))
  cosh(x)
end

cosh(a::DD) = DD(coshAsTD(a))

@inline function tanhAsTD(a::DD)
  isneg, abs_a = signbit(a), abs(a)
  if abs_a.hi < 1.0e-15
      t = TD((exp(abs_a)-exp(-abs_a))/(exp(abs_a)+exp(-abs_a)))
  else 
      x = TD(abs_a)
      t = tanh(x)
  end    
  isneg ? -t : t
end

tanh(a::DD) = DD(tanhAsTD(a))

#tanh(a::DD) = sinh(a)/cosh(a) # DD(tanhAsTD(a))

@inline cschAsTD(a::DD) = recip(sinhAsTD(a))
csch(a::DD) = DD(cschAsTD(a))

@inline sechAsTD(a::DD) = recip(coshAsTD(a))
sech(a::DD) = DD(sechAsTD(a))

@inline cothAsTD(a::DD) = recip(cothAsTD(a))
coth(a::DD) = DD(cothAsTD(a))
