
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
  x = TD(abs_a)
  t = tanh(x)
  isneg ? -t : t
end
tanh(a::DD) = DD(tanhAsTD(a))

tanh(a::DD) = sinh(a)/cosh(a) # DD(tanhAsTD(a))

@inline cschAsTD(a::DD) = recip(sinhAsTD(a))
csch(a::DD) = DD(cschAsTD(a))

@inline sechAsTD(a::DD) = recip(coshAsTD(a))
sech(a::DD) = DD(sechAsTD(a))

@inline cothAsTD(a::DD) = recip(cothAsTD(a))
coth(a::DD) = DD(cothAsTD(a))
