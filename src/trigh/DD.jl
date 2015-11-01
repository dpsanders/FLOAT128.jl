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

function sinh_taylor_series(a::DD)
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

  z = TD(x) + x*(dd_inv_fact[3]*x2 + dd_inv_fact[5]*x4 + dd_inv_fact[7]*x6)
  z2 = x9 * (dd_inv_fact[9] + x2*dd_inv_fact[11] + x4*dd_inv_fact[13] + x6*dd_inv_fact[15])
  z3 = x17 * (dd_inv_fact[17] + x2*dd_inv_fact[19] + x4*dd_inv_fact[21] + x6*dd_inv_fact[23])
  z4 = x25 * (dd_inv_fact[25] + x2*dd_inv_fact[27] + x4*dd_inv_fact[29] + x6*dd_inv_fact[31])

  DD(z + ((z4+z3)+z2))
end

function cosh_taylor_series(a::DD)
  x = a
  x2 = x*x
  x4 = x2*x2
  x6 = x2*x4
  x8 = x4*x4
  x16 = x8*x8
  x24 = x8*x16

  z = (dd_inv_fact[2]*x2 + dd_inv_fact[4]*x4 + dd_inv_fact[6]*x6)
  z2 = x8 * (dd_inv_fact[8] + x2*dd_inv_fact[10] + x4*dd_inv_fact[12] + x6*dd_inv_fact[14])
  z3 = x16 * (dd_inv_fact[16] + x2*dd_inv_fact[18] + x4*dd_inv_fact[20] + x6*dd_inv_fact[22])
  z4 = x24 * (dd_inv_fact[24] + x2*dd_inv_fact[26] + x4*dd_inv_fact[28] + x6*dd_inv_fact[30])

  ((z4+z3)+z2)+z + 1.0
end



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

  z = x + x*(dd_tanh_coeff[2]*x2 + dd_tanh_coeff[3]*x4 + dd_tanh_coeff[4]*x6)
  z2 = x9 * (dd_tanh_coeff[5] + x2*dd_tanh_coeff[6] + x4*dd_tanh_coeff[7] + x6*dd_tanh_coeff[8])
  z3 = x17 * (dd_tanh_coeff[9] + x2*dd_tanh_coeff[10] + x4*dd_tanh_coeff[11] + x6*dd_tanh_coeff[12])

  (z + (z3+z2))
end

#=
   tanh z = z / (1 + z^2/(3 + z^2/5 + ...))
=#   

function sinh(x::DD)
  isneg, abs_a = signbit(x), abs(x)
  if (abs_a.hi <= 0.1)
     s = sinh_taylor_series(abs_a)
  else
     epx = exp(abs_a)
     emx = 1.0/epx
     epx = epx - emx
     s = divby2(epx)
  end   
  isneg ? -s : s
end

function cosh(x::DD)
  abs_a = abs(x)
  if abs_a.hi <= 0.1
      c = cosh_taylor_series(abs_a)
  else    
      epx = exp(abs_a)
      emx = 1.0/epx
      epx = epx + emx
      c = divby2(epx)
  end
  c
end

function tanh(x::DD)
  isneg, abs_a = signbit(x), abs(x)
  if abs_a.hi < 0.0325
      t = tanh_taylor_series(abs_a)
  elseif abs_a.hi > 0.05
      ea = exp(abs_a)
      inv_ea = exp(-abs_a) # do not use 1/ea here
      t = (ea - inv_ea) / (ea + inv_ea)
  else    
      s = sinh(abs_a)
      c = s*s
      c = 1.0 - c
      c = sqrt(c)
      t = s/c
  end    
  isneg ? -t : t
end


#tanh(a::DD) = sinh(a)/cosh(a) # DD(tanhAsTD(a))

csch(a::DD) = recip(sinh(a))
sech(a::DD) = recip(cosh(a))
coth(a::DD) = recip(tanh(a))
