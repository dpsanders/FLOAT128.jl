# FLOAT128.jl  
### *(work in progress -- incomplete)*
100+ valid significand bits for elementary functions with conventionally small values.

### exports 

type Float128

    convert from Float64|32|16, Int64|32|16, BigFloat
    
    ldexp, frexp, abs, sign, signbit, trunc, floor, ceil, round, 
    
    ==, !=, <, <=, >=, >, isless, isequal
    
    +, -, *, /, %, mod, rem, div, fld, cld, divrem, fldmod
    
    sqrt, exp, log, ^
    
    sin, cos, tan, csc, sec, cot, asin, acos, atan, acsc, asec, acot
    
    sinh, cosh, tanh, csch, sech, coth, asinh, acosh, atanh, acsch, asech, acoth 
    
    string, show, showcompact
    
dict Float128s{Symbol, Float128}

    :zero, :one, :two, :three, :four, :sqrt2, :log2, :log10, :exp1, :pi, :twopi
    
    reciprocals :half, :quarter, :_sqrt2, :_log2, :_log10, :_exp1, :_pi, :_twopi

clean(x::Float128)

    zeros the low order part when it is less than eps(eps(high order part))

### Use

