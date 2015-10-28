function atanh(x::DD)
    isneg, abs_x = signbit(x.hi), abs(x)
    if abs_x.hi > one(DD)
        throw(ErrorException("atanh arg $(x) outside -1..+1"))
        return dd_NaN
    elseif abs_x >= (one(DD)-1.0e-96)
        return isneg ? -dd_Inf : dd_Inf
    end
    if abs_x.hi <= 0.00275
        z = atanh_series(abs_x)
    else
        nm = one(TD) + abs_x
        dn = one(TD) - abs_x
        q = nm/dn
        q = log(DD(q))
        z = divby2(q)
    end
    isneg ? -z : z
end

function asinh(x::DD)
    isneg, abs_x = signbit(x),abs(x)
    if isneg
        -asinh(abs_x)
    elseif abs_x.hi <= 0.045
        b = asinh_series(abs_x)
    else
        y = one(TD)+sqrt(TD(x)*x+one(TD))
        y = TD(x) / y
        y = atanh(DD(y))
        mulby2(y)
    end
end

function acosh(x::DD)
   isneg, abs_x = signbit(x.hi), abs(x)
    if isneg || abs_x.hi < one(Float64)
        throw(ErrorException("acosh arg $(x) outside +1.."))
        return dd_NaN
    elseif abs_x <= (one(DD)+1.0e-96)
        return zero(DD)
    end

    if abs_x.hi <= 1.01 #1.00275
        a = sqrt(TD(abs_x)*abs_x - one(TD))
        b = asinh(DD(a))
        b
    else
        a = sqrt(TD(abs_x)*abs_x - one(TD))
        b = a + abs_x
        b = log(DD(b))
        b
    end
end
