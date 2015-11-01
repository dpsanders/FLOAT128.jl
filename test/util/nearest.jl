import Base:convert


function ieeedouble(xx::BigFloat)
    x=logabs=powermin=0.0
    mantissa=exponent=expmin=expmax=expmiddle=powermax=powermiddle=infmantissa=0.0
    sgn=0
    set_bigfloat_precision(320) ; #Digits = 100;
    x = xx;#$evalf(xx);
    if (x==0)
        sgn, exponent, mantissa = 1, -1022, 0
    else
      if (x < 0)
         sgn = -1
      else
        sgn = 1
      end
      x = abs(x);
      if x >=  2.0^(1023)*(2-2.0^(-53))
        mantissa = Inf
        exponent = 1023
      elseif x <= 2.0^(-1075)
        mantissa = 0; exponent = -1022
      elseif x <= 2.0^(-1022)
       exponent = -1022
      else
    # x is between 2^(-1022) and 2^(1024)
         powermin = 2.0^(-1022); expmin = -1022;
         powermax = 2^1024; expmax = 1024;
         while (expmax-expmin > 1)
            expmiddle = round(Int,(expmax+expmin)/2);
            powermiddle = 2.0^expmiddle;
            if x >= powermiddle
                powermin = powermiddle;
                expmin = expmiddle
            else
                powermax = powermiddle;
                expmax = expmiddle
            end
          end
          # now, expmax - expmin = 1 and powermin <= x < powermax,
          # powermin = 2^expmin and powermax = 2^expmax, so expmin is the exponent of x
          exponent = expmin;
         #end;
         infmantissa = x*2.0^(52-exponent);
         if infmantissa-round(Int,infmantissa) != 0.5
            mantissa = round(Int,infmantissa)
            else
              mantissa = floor(Int,infmantissa);
               if reinterpret(UInt64,mantissa) & one(UInt64) == one(UInt64)
                  mantissa = mantissa+1
               end
            end
         mantissa = mantissa*2.0^(-52);
      #end;
      end;
    end;
    sgn,exponent,mantissa;
end

function nearest(x::BigFloat)
    sgn,exponent,mantissa = ieeedouble(x)
    mantissa = convert(Float64,mantissa)
    sgn>=0 ? ldexp(mantissa,exponent) : -ldexp(mantissa,exponent)
end


function nnear(n::Int,x::BigFloat)
    bfprec = get_bigfloat_precision()
    set_bigfloat_precision(1280)
    z = zeros(Float64,n)
    z[1] = nearest(x)
    rm = x-BigFloat(z[1])
    for i in 2:n
       z[i] = nearest(rm)
       rm = rm - BigFloat(z[i])
    end
    set_bigfloat_precision(bfprec)
    z
 end

function nearest12(x::BigFloat)
    set_bigfloat_precision(4096)
    z = zeros(Float64,12)
    a = nnear(6,x)
    bfa = BigFloat(a[1])+BigFloat(a[2])+BigFloat(a[3])
    y = x-bfa
    b = nnear(6,y)
    bfb = BigFloat(b[1])+BigFloat(b[2])+BigFloat(b[3])
    y = x-(bfa+bfb)
    c = nnear(6,y)
    bfc = BigFloat(c[1])+BigFloat(c[2])+BigFloat(c[3])
    y = x-(bfa+bfb+bfc)
    d = nnear(6,y)

    z[1:3] = a[1:3]
    z[4:6] = b[1:3]
    z[7:9] = c[1:3]
    z[10:12] = d[1:3]
    z
end


function highlow(x::BigFloat)
    hi = nearest(x)
    lo = nearest(x-BigFloat(hi))
    hi,lo
end
function highlow(x::Real)
  if length(fieldnames(x)) == 0
    hilo(BigFloat(x))
  elseif length(fieldnames(x)) == 2
    hilo(BigFloat(x.hi)+BigFloat(x.lo))
  elseif length(fieldnames(x)) == 3
    hilo(BigFloat(x.hi)+BigFloat(x.md)+BigFloat(x.lo))
  else
    ErrorException("unsupported type")
  end
end

function highmedlow(x::BigFloat)
    hi = nearest(x)
    md = nearest(x-BigFloat(hi))
    lo = nearest(x-BigFloat(hi)-BigFloat(md))
    hi,md,lo
end
function highmedlow(x::Real)
  if length(fieldnames(x)) == 0
    hml(BigFloat(x))
  elseif length(fieldnames(x)) == 2
    hml(BigFloat(x.hi)+BigFloat(x.lo))
  elseif length(fieldnames(x)) == 3
    hml(BigFloat(x.hi)+BigFloat(x.md)+BigFloat(x.lo))
  else
    ErrorException("unsupported type")
  end
end

hilo(x::BigFloat) = DD(highlow(x))
hml(x::BigFloat) = TD(highmedlow(x)...)
hilo(x::AbstractFloat) = DD(highlow(x))
hml(x::AbstractFloat) = TD(highmedlow(x)...)


besthilo(fn::Function,a::DD) = hilo(fn(BigFloat(convert(Tuple,a))))
besthml(fn::Function,a::TD) = hml(fn(BigFloat(convert(Tuple,a))))
besthilo(fn::Function,a::TD) = besthilo(fn,DD(a.hi,a.md))
besthml(fn::Function,a::DD) = besthml(fn,TD(a.hi,a.lo))
besthilo(fn::Function,a::AbstractFloat) = besthilo(fn,DD(a))
besthml(fn::Function,a::AbstractFloat) = besthml(fn,TD(a))


convert(::Type{BigFloat}, a::Tuple{Float64,Float64}) = BigFloat(a[1])+BigFloat(a[2])
convert(::Type{BigFloat}, a::Tuple{Float64,Float64, Float64}) = BigFloat(a[1])+BigFloat(a[2])+BigFloat(a[3])

convert(::Type{BigFloat}, a::DD) = BigFloat(a.hi)+BigFloat(a.lo)
convert(::Type{BigFloat}, a::TD) = BigFloat(a.hi)+BigFloat(a.md)+BigFloat(a.lo)

convert(::Type{DD}, a::BigFloat) = hilo(a)
convert(::Type{TD}, a::BigFloat) = hml(a)
