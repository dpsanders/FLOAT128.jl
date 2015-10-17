VERSION >= v"0.4.0" && __precompile__()

module FLOAT125

# using

# import

# export

export Float125

# include

# error-free transformations by number of args in & vals out
include("eft/single.jl")
include("eft/double.jl")
include("eft/triple.jl")

# working types (double-double and triple-double)
include("types.jl")
include("const/DD.jl")
include("const/TD.jl")

# magnitude: comparison and subnumeric composition
include("mag/DD.jl")
include("mag/TD.jl")

# arithmetic and arithmetic support
include("arith/DDsupp.jl")
include("arith/TDsupp.jl")
include("arith/DD.jl")
include("arith/TD.jl")

# integer powers, reciprocal roots
include("pow/DD.jl")

# log, exp
include("exp/DD.jl")
include("exp/TD.jl")

# io
include("io/DD.jl")
include("io/TD.jl")


end # module
