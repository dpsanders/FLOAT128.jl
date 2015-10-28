
@inline function eftSum2{T<:FloatingPoint}(a::T, b::T)
  x = a + b
  t = x - a
  y = (a - (x - t)) + (b - t)
  x,y
end


@inline function eftSum2inOrder{T<:FloatingPoint}(a::T, b::T)
  x = a + b
  y = b - (x - a)
  x,y
end


@inline function eftDiff2{T<:FloatingPoint}(a::T, b::T)
  x = a - b
  t = x - a
  y = (a - (x - t)) - (b + t)
  x,y
end


@inline function eftDiff2inOrder{T<:FloatingPoint}(a::T, b::T)
  x = a - b
  y = (a - x) - b
  x,y
end

@inline function eftProd2{T<:FloatingPoint}(a::T, b::T)
    x = a * b
    y = fma(a, b, -x)
    x,y
end

# !!sassafrass!!
# presently 'y' must be negated to get the right result
# I do not know why.

@inline function eftDiv2{T<:FloatingPoint}(a::T,b::T)
     x = a/b
     y = fma(x,b,-a)/b
     x,-y
end

#=
"""
 (a,b) ↦ (x,y)\\
x⊕y ≖ a⊕b and x⊕y≖x
"""->eftSum2

"""
 (a,b) ↦ (x,y)\\
x⊕y ≖ a⊕b and x⊕y≖x\\
 *presumes* |a| ≥ |b|
"""->eftSum2InOrder

"""
 (a,b) ↦ (x,y)\\
x⊕y ≖ a⊝b and x⊕y≖x
"""->eftDiff2

"""
 (a,b) ↦ (x,y)\\
x⊕y ≖ a⊝b and x⊕y≖x\\
*presumes* |a| ≥ |b|
"""->eftDiff2inOrder
=#
