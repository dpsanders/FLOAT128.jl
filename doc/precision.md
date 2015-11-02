# FLOAT128.jl  
### *(work in progress -- incomplete)*
100+ valid significand bits for elementary functions with conventionally small values.

  relative errors (*max found with 20,000 random double-double values each function*)
  relative time is fn(Float128)/fn(Float64), avg time over 20,000 evals (on 1 machine, so ymmv)
  
  relative big is fn(Float128)/fn(BigFloat_160) " 


| func | over | sig bits | rel err | rel time | rel big |
|------|------|-----------|---------|---------|---------|
| sqrt | 0..64G | 106 | 1.2e-32 | 1.35 | 1.06 |
|      |             |     |       |  | |
| exp  | -15..15   | 104 | 4.9e-32 |2.75  | 0.75 |
| exp  | -300..300   | 103 | 5.2e-32 |t | |
|      |             |     |       | |
| log  |    1..64G   | 105 | 1.9e-33 |8.0 | 1.6 |
|      |             |     |       | |
| sin, cos  | -2pi..2pi   | 106 | 1.2e-32 | 4.5 | 0.8 |
| tan  | -2pi..2pi   | 104 | 4.8e-32 | 10 | 3.0 |
| csc, sec, cot | -2pi..2pi | 104 | 4.8e-32 | t |
|      |             |     |       | |
| sin, cos  | -64Gpi..64Gpi   | 105 | 2.25-32 | t |
| tan  | -64Gpi..64Gpi   | 103 | 5.9e-32 | t |
| csc, sec, cot | -64Gpi..64Gpi | 103 | 8.8e-32 | t |
|      |             |     |       |
| asin, acos  | -1..1     | 106 | 8e-33 | t |
| atan, acot  | -1..1   | 106 | 1.2e-32 | t |
| acsc  | -64G..-1,1..64G | 106  | 2.5e-32  | t |
| asec  | -4096..-1,1..4096 | 104  | 2.7e-32  | t, !!improve!! |
|      |             |     |       |
| sinh, cosh, tanh  | -1..1   | 104 | 4.5e-32 | t |
| sinh, cosh, tanh  | -256..256   | 104 | 4.5e-32 | t |
| csch, sech, coth  | -1..1| 104  | 4.6e-32 | t |
| csch, sech, coth  | -256..256   | 103 | 6.7e-32 | t |
|      |            |     |       |
| asinh  | -1..1     | 105 | 2.2e-32 | t |
| acosh  |  1..16     | 105 | 2.0e-32 | t |
| atanh  | -1..1   | 109 | 1.5e-33 | t |
| acsch, asech, acoth  | | ?  |  |

** 64G = 4096^3 **

##### outliers
| function | arg | eval | true |
|----------|-----|------|------|
| tanh     |(5.551170634277014e-17, -6.162975822039155e-33) | (5.551170634277014e-17, -6.162975822039155e-33) | (5.551170634277013e-17, 6.162975822039155e-33) |
| asinh     |(1.1102341268554029e-16, -1.232595164407831e-32) | (1.1102341268554029e-16, -1.232595164407831e-32) | (1.1102341268554026e-16, 1.232595164407831e-32) |

  
