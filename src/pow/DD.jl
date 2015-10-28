
function (sqr){T<:DD}(a::T)
  t1,t2 = eftProd2(a.hi,a.hi)
  t3 = a.hi * a.lo
  t5 = t3 + t3
  t6 = t2 + t5
  t1,t6 = eftSum2inOrder(t1,t6)
  DD(t1,t6)
end

function (cub){T<:DD}(a::T)
  t1,t2 = eftProd2(a.hi,a.hi)
  t3 = a.hi * a.lo
  t5 = t3 + t3
  t6 = t2 + t5
  t1,t6 = eftSum2inOrder(t1,t6)
  a * DD(t1,t6)
end

function npow{T<:DD,N<:Integer}(a::T,n::N)
    if n == zero(N)
        if a.hi == zero(typeof(a.hi))
            throw(DomainError)
        else
            one(typeof(a))
        end
    end

    r = a
    s = one(typeof(a))
    m = abs(n)

    if m > one(N)
        # binary exponentiation
        while (m > 0)
            if m%2 == one(N)
                s = s*r
            end
            m >>= 1
            if m > zero(N)
                r = sqr(r)
            end
        end
    else
        s = r
    end

    if (n < zero(N))
        s = recip(s)
    end

    s
end

#=
function sqrt(a::DD)
    if a.hi <= zero(Float64)
       if a.hi == zero(Float64)
           return zero(DD)
       else
           throw(DomainError("sqrt expects a nonnegative base"))
       end
    end

    fl1 = one(Float64)
    r = one(Float64)/sqrt(a.hi)
    r += (r * (fl1 - (a * (r*r)))) * 0.5
    r += (r * (fl1 - (a * (r*r)))) * 0.5

    r*a
end
=#

function cubrt{T<:DD}(a::T)
    if a.hi <= zero(typeof(a.hi))
       if a.hi == zero(typeof(a.hi))
           return zero(typeof(a))
       else
           throw(DomainError("cubrt expects a nonnegative base"))
       end
    end

    fl1 = one(typeof(a.hi))
    r = a.hi^(-fl1/3.0)
    r += r * (fl1 - (a * (r*r))) / dd_three
    r += r * (fl1 - (a * (r*r))) / dd_three


    recip(r)
end

function nroot{T<:DD,N<:Integer}(a::T,p::N)
   if p <= one(N)
     if p == zero(N)
         return one(typeof(a))
     elseif p == one(N)
         return a
     end
   elseif a.hi <= zero(typeof(a.hi))
        if a.hi == zero(typeof(a.hi))
            return zero(typeof(a))
        else
            throw(DomainError("iroot expects a nonnegative base"))
        end
   end

   fl1 = one(typeof(a.hi))
   flp = convert(Float64,p)
   f2p = (T)(flp,zero(typeof(a.hi)))
   r = a.hi^(-fl1/flp)
   r += r * (fl1 - (a * (npowr(r,p)))) / f2p
   r += r * (fl1 - (a * (npowr(r,p)))) / f2p

   recip(r)
end
