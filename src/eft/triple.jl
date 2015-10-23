function eftSum3{T<:Float64}(a::T,b::T,c::T)
    s,t = eftSum2(b, c)
    x,u = eftSum2(a, s)
    y,z = eftSum2(u, t)
    x,y = eftSum2inOrder(x, y)
    x,y,z
end

function eftSum3inOrder{T<:Float64}(a::T,b::T,c::T)
    s,t = eftSum2inOrder(b, c)
    x,u = eftSum2inOrder(a, s)
    y,z = eftSum2inOrder(u, t)
    x,y = eftSum2inOrder(x, y)
    x,y,z
end

function eftSum3as2{T<:Float64}(a::T,b::T,c::T)
    s,t = eftSum2(b, c)
    x,u = eftSum2(a, s)
    y = u+t
    x,y = eftSum2inOrder(x, y)
    x,y
end

function eftSum3inOrderAs2{T<:Float64}(a::T,b::T,c::T)
    s,t = eftSum2inOrder(b, c)
    x,u = eftSum2inOrder(a, s)
    y = u+t
    x,y = eftSum2inOrder(x, y)
    x,y
end

function eftProd3{T<:Float64}(a::T, b::T, c::T)
    p,e = eftProd2(a,b)
    x,p = eftProd2(p,c)
    y,z = eftProd2(e,c)
    x,y,z
end

function eftProd3as2{T<:Float64}(a::T, b::T, c::T)
    p,e = eftProd2(a,b)
    x,p = eftProd2(p,c)
    y = e*c
    x,y
end



"""
error free transformation for fused multiply add (fma)
 (a,b,c) ↦ (x,y,z)\\
x⊕y⊕z ≖ a⊗b⊕c and x==fma(a,b,c) and x⊕y≖x, y⊕z≖y

ref:
Some functions computable with a fused-mac
by Sylvie Boldo and Jean-Michel Muller
Proceedings of the 17th IEEE Symposium on Computer Arithmetic, 2005
"""
function eftFMA{T<:Float64}(a::T, b::T, c::T)
    x = fma(a,b,c)
    u1,u2 = eftProd2(a,b)
    a1,a2 = eftSum2(u2,c)
    b1,b2 = eftSum2(u1,a1)
    g = (b1-x)+b2
    y,z = eftSum2inOrder(g,a2)
    x,y,z
end

"""
error free transformation for fused multiply subtract (fms)
 (a,b,c) ↦ (x,y,z)\\
x⊕y⊕z ≖ a⊗b⊝c and x==fma(a,b,-c) and x⊕y≖x, y⊕z≖y
"""
function eftFMS{T<:Float64}(a::T, b::T, c::T)
    x = fma(a,b,c)
    u1,u2 = eftProd2(a,b)
    a1,a2 = eftDiff2(u2,c)
    b1,b2 = eftSum2(u1,a1)
    g = (b1-x)+b2
    y,z = eftSum2inOrder(g,a2)
    x,y,z
end
