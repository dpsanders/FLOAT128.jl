
function sinh(x::DD)
  0.5 * (exp(x) - exp(-x))
end

function cosh(x::DD)
  0.5 * (exp(x) + exp(-x))
end

function tanh(x::TD)
  epx = exp(x)
  emx = exp(-x)
  (epx - emx) / (epx + emx)
end
