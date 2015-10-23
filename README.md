# FLOAT125.jl  
### *(work in progress -- incomplete)*
100+ valid significand bits for elementary functions with conventionally small values.

### Why?

typealias FP Float128

fp0=FP(0.0); fp1=FP(0.1); fp2=FP(0.2); fp3=FP(0.3);

##### more arithmetic just works

0.0  != 0.1+0.2+0.3 -0.3-0.2-0.1

fp0  == fp1+fp2+fp3 -fp3-fp2-fp1

