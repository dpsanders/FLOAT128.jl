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
