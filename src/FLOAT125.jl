VERSION >= v"0.4.0" && __precompile__()

module FLOAT125

# using

# import
import Base: convert, promote_rule,
             zero, one, isinf, isnan, isfinite,
             (-), abs, sign, signbit, copysign, flipsign,
             isequal, isless, (<),(<=),(==),(>=),(>),
             (+),(*),(/),(\), fma,
             round, floor, ceil, trunc, fld, cld,
             ldexp, frexp, modf,
             div, rem, divrem, mod, fldmod, (%),
             sqrt, hypot, (^),
             exp, log, log2, log10, expm1, log1p,
             sin, cos, tan, csc, sec, cot,
             sinh, cosh, tanh, csch, sech, coth,
             asin, acos, atan, acsc, asec, acot, atan2,
             asinh, acosh, atanh, acsch, asech, acoth,
             string, show, showcompact, parse, hex

# export

export Float125

# include
include("types.jl")

# error-free transformations by number of args in & vals out
include("eft/single.jl")
include("eft/double.jl")
include("eft/triple.jl")

# working types
include("const/DD.jl")
include("const/TD.jl")

# magnitude: comparison and subnumeric composition
include("mag/DD.jl")
include("mag/TD.jl")

# arithmetic and arithmetic support
include("arith/DD.jl")
include("arith/TD.jl")
include("arith/DDmanip.jl")
include("arith/TDmanip.jl")

# integer powers, reciprocal roots
# include("pow/DD.jl")
# include("pow/TD.jl")

# log, exp
include("exp/DD.jl")

# io
include("io/DD.jl")
include("io/TD.jl")


end # module
