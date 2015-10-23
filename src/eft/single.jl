@inline function eftSquare(a::Float64)
    x = a * a
    y = fma(a, a, -x)
    x,y
end

function eftCube(a::Float64)
    p = a*a; e = fma(a, a, -p)
    x = p*a; p = fma(p, a, -x)
    y = e*a
    x,y
end

# !!sassafrass!!
# presently 'y' must be negated to get the right result
# I do not know why.
#
@inline function eftRecip(a::Float64)
     x = one(Float64)/a
     y = fma(x,a,-1.0)/a
     x,-y
end
