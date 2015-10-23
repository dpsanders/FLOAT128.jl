# FLOAT125.jl  
### *(work in progress -- incomplete)*
100+ valid significand bits for elementary functions with conventionally small values.

### Why?

typealias FP Float128

##### more arithmetic just works

fl0=0.0;     fl1=0.1;     fl2=0.2;     fl3=0.3;

fp0=FP(0.0); fp1=FP(0.1); fp2=FP(0.2); fp3=FP(0.3);

fl0  !=  fl1+fl2+fl3-fl3-fl2-fl1\\

fp0  ==  fp1+fp2+fp3-fp3-fp2-fp1

