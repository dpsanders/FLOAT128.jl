const TD_typename = "Float192"

function show(io::IO, x::TD)
    hi=string(x.hi)
    md=string(x.md)
    lo=string(x.lo)
    print(io, string(TD_typename,"(",hi,", ",md,", ",lo,")"))
end

function parse(::Type{TD}, str::AbstractString)
    if !(startswith(str,TD_typename))
        throw(ErrorException("$(TD_typename) not found in $(str)"))
    end
    s = str[10:(end-1)]
    himdlo = split(s,',')
    if length(hilo)==1
        himdlo = [himdlo[1],"0.0","0.0"]
    elseif length(himdlo)==2
        himdlo = [himdlo[1],himdlo[2],"0.0"]
    end

    shi,smd,slo = himdlo[1],himdlo[2],himdlo[3]
    hi = parse(Float64,shi)
    md = parse(Float64,smd)
    lo = parse(Float64,slo)
    hi,md,lo = eftSum3(hi,md,lo)
    TD(hi,md,lo)
end


function hex(x::TD)
    hi,md,lo = x.hi, x.md, x.lo
    shi = @sprintf("%a",hi)
    smd = @sprintf("%a",lo)
    slo = @sprintf("%a",lo)
    s = string(TD_typename,"(",shi,", ",smd,", ",slo,")")
    s
end
