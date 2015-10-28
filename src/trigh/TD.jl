
function sinh(x::TD)
  epx = exp(x)
  emx = 1.0/epx
  epx = epx - emx
  divby2(epx)
end

function cosh(x::TD)
  epx = exp(x)
  emx = 1.0/epx
  epx = epx + emx
  divby2(epx)
end

function tanh(x::TD)
  epx = exp(x)
  emx = 1.0/epx
  n = (epx - emx)
  d = (epx + emx)
  n/d
end
