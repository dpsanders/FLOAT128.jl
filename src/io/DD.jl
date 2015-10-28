const DD_typename = "Float128"

function show(io::IO, x::DD)
    hi=string(x.hi)
    lo=string(x.lo)
    print(io, string(DD_typename,"(",hi,", ",lo,")"))
end

function parse(::Type{DD}, str::AbstractString)
    if !(startswith(str,DD_typename))
        throw(ErrorException("$(DD_typename) not found in $(str)"))
    end
    s = str[10:(end-1)]
    hilo = split(s,',')
    if length(hilo)==1
        hilo = [hilo[1],"0.0"]
    end
    shi,slo = hilo[1],hilo[2]
    hi = parse(Float64,shi)
    lo = parse(Float64,slo)
    hi,lo = eftSum2(hi,lo)
    DD(hi,lo)
end

function hex(x::DD)
    hi,lo = x.hi, x.lo
    shi = @sprintf("%a",hi)
    slo = @sprintf("%a",lo)
    s = string(DD_typename,"(",shi,", ",slo,")")
    s
end
