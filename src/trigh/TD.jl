
function sinh(x::TD)
  0.5 * (exp(x) - exp(-x))
end

function cosh(x::TD)
  0.5 * (exp(x) + exp(-x))
end

function tanh(x::TD)
  epx = exp(x)
  emx = exp(-x)
  (epx - emx) / (epx + emx)
end
