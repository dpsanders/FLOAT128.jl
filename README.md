# FLOAT128.jl  
### *(work in progress -- incomplete)*
100+ valid significand bits for elementary functions with conventionally small values.

exports type Float128

### Why?

typealias FP Float128

fp0=FP(0.0); fp1=FP(0.1); fp2=FP(0.2); fp3=FP(0.3);

##### this will work 

0.0  != 0.1+0.2+0.3 -0.3-0.2-0.1

fp0  == fp1+fp2+fp3 -fp3-fp2-fp1

##### this may matter

Starting with accurate numbers and after much computation, a valuation is obtained.
Some few or more of the least significant bits may reflect more noise than value.

While computing with a bunch of extra bits is not a panacea, if the computation
is well structured and the underlying code base avoids turning risky corners,
obtaining the valuation by keeping only the high order (Float64) constituents
is often better than presuming correctness of the lowest order octet(s).

