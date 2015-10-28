#=
   The typealias *FloatingPoint* is used as a parameter for error free transformations.
   If you want multiple Floating Point types available, use a Union here.
   For most purposes, it should be either Float64 or, e.g. for CUDA Float32.
=#
typealias FloatingPoint Float64
