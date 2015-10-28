
function sinh(x::TD)
  divby2(exp(x) - exp(-x))
end

function cosh(x::TD)
  divby2(exp(x) + exp(-x))
end

function tanh(x::TD)
  epx = exp(x)
  emx = exp(-x)
  (epx - emx) / (epx + emx)
end
