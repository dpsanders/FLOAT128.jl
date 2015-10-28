
function sinh(x::TD)
  epx = exp(x)
  emx = exp(-x)
  epx = epx - emx
  divby2(epx)
end

function cosh(x::TD)
  epx = exp(x)
  emx = exp(-x)
  epx = epx + emx
  divby2(epx)
end

function tanh(x::TD)
  epx = exp(x)
  emx = exp(-x)
  n = (epx - emx)
  d = (epx + emx)
  n/d
end
