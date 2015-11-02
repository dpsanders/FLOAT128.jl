#=
relerr testing
=#

set_bigfloat_precision(320);

@inline abserr(xIdeal::Real,xActual::Real) = abs(xIdeal-xActual)
relerr(xIdeal::Real,xActual::Real) = xIdeal==0 ? 0 : abserr(float(xIdeal),float(xActual))/abs(float(xIdeal))

function relerr1{T<:Float64}(xIdeal::T,xActual::T)
   a = abserr(xIdeal,xActual)
   x = abs(xActual)
   if x < eps(eps(1.0))
       zero(Float64)
   else
       a / xIdeal
   end
end

function relerr{T<:DD}(xIdeal::T,xActual::T)
    if (xIdeal.hi != xActual.hi)
        relerr1(xIdeal.hi, xActual.hi)
    else
        a = abs(xIdeal-xActual)
        x = abs(xActual)
        if (a.hi < 1e-90) || (x.hi < 1.e-90)
            zero(Float64)
        else
            (a / x).hi
        end
    end
end

relerr(xIdeal::TD,xActual::TD) = relerr(DD(xIdeal), DD(xActual))

#julia> r=[abs(tst3()) for i in 1:20000];sort(r)[end]

function runFnTestWide(fn::Function, n::Int=20_000, bitgap::Function=wideBitgap, expow2::Function=expon)
    src = zeros(DD,n)
    err = zeros(Float64,n)
    rerr = TD(0.0)
    serr = TD(0.0)
    for i in 1:n
        r = randhilo(bitgap,expow2)
        src[i] = r
        ideal = besthilo(fn,r)
        actual = fn(r)
        rerr1 = relerr(ideal,actual)
        if rerr1 > rerr
            rerr = rerr1
            serr = r
        end
        # err[i] = rerr
    end
    #src, err
    serr, rerr
end



function runFnTest(fn::Function, n::Int=20_000, lo::Float64=0.0, hi::Float64=1.0, bitgap::Function=baseBitgap)
    src = zeros(DD,n)
    err = zeros(Float64,n)
    rerr = DD(0.0)
    serr = DD(0.0)
    for i in 1:n
        r1 = rand(lo:hi)
        bg = bitgap()
        r2 = ldexp(rand(-0.99999:0.999999),-53-bg)
        r = DD(r1)+DD(r2)
        src[i] = r
        ideal = besthilo(fn,r)
        actual = fn(r)
        rerr1 = relerr(ideal,actual)
        #err[i] = rerr
        if rerr1 > rerr
            rerr = rerr1
            serr = r
        end

    end
    #src, err
    serr,rerr
end

function fnTestWide(fn::Function, n::Int=20_000, bitgap::Function=wideBitgap, expow2::Function=expon)
    #src, err = runFnTest(fn,n,bitgap,expow2)
    #maxerr = sort(err)[end]
    #srcmaxerr = src[ maxerr .== err ][1]
    srcmaxerr,maxerr=runFnTesWidet(fn,n,bitgap,expow2)
    frexp(maxerr)[2], maxerr, srcmaxerr
end

function fnTest(fn::Function, n::Int=20_000, lo::Float64=0.0, hi::Float64=1.0, bitgap::Function=baseBitgap)
    #src, err = runFnTest2(fn,n,lo,hi,bitgap)
    #maxerr = sort(err)[end]
    #srcmaxerr = src[ maxerr .== err ][1]
    #frexp(maxerr)[2], maxerr, srcmaxerr
    srcmaxerr,maxerr=runFnTest(fn,n,lo,hi,bitgap)
    frexp(maxerr)[2], maxerr, srcmaxerr
end




function runFnTestWideTD(fn::Function, n::Int=20_000, bitgap::Function=wideBitgap, expow2::Function=expon)
    src = zeros(TD,n)
    err = zeros(Float64,n)
    for i in 1:n
        r = randhml(bitgap,expow2)
        src[i] = r
        ideal = besthml(fn,r)
        actual = fn(r)
        rerr = relerr(ideal,actual)
        err[i] = rerr
    end
    src, err
end

function runFnTestTD(fn::Function, n::Int=20_000, bitgap::Function=wideBitgap, expow2::Function=expon)
    src = zeros(TD,n)
    err = zeros(Float64,n)
    for i in 1:n
        r = randhml(bitgap,expow2)
        src[i] = r
        ideal = besthml(fn,r)
        actual = fn(r)
        rerr = relerr(ideal,actual)
        err[i] = rerr
    end
    src, err
end

function fnTestWideTD(fn::Function, n::Int=20_000, bitgap::Function=wideBitgap, expow2::Function=expon)
    #src, err = runFnTest(fn,n,bitgap,expow2)
    #maxerr = sort(err)[end]
    #srcmaxerr = src[ maxerr .== err ][1]
    srcmaxerr,maxerr=runFnTesWideTD(fn,n,bitgap,expow2)
    frexp(maxerr)[2], maxerr, srcmaxerr
end

function fnTestTD(fn::Function, n::Int=20_000, lo::Float64=0.0, hi::Float64=1.0, bitgap::Function=baseBitgap)
    #src, err = runFnTest2(fn,n,lo,hi,bitgap)
    #maxerr = sort(err)[end]
    #srcmaxerr = src[ maxerr .== err ][1]
    #frexp(maxerr)[2], maxerr, srcmaxerr
    srcmaxerr,maxerr=runFnTestTD(fn,n,lo,hi,bitgap)
    frexp(maxerr)[2], maxerr, srcmaxerr
end

