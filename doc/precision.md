# FLOAT128.jl  
### *(work in progress -- incomplete)*
100+ valid significand bits for elementary functions with conventionally small values.

  relative errors (*max found with 20,000 random double-double values each function*)
  
  ∝64 is time fn(Float128)/fn(Float64), @time (on 1 machine, so ymmv)
  
  ∝Big is time n(Float128)/fn(BigFloat_160)  similarly


| func | over | sig bits | rel err | ∝64 |  ∝Big |
|------|------|----------|---------|------|-------|
| sqrt | 0..64G | 106 | 1.2e-32 | 1.6| 0.6 |
|      |             |     |       |  | |
| exp  | -15..15   | 104 | 4.9e-32 |2.75  | 0.35 |
| exp  | -300..-15,15..300   | 103 | 5.2e-32 |2.35| 0.66 |
|      |             |     |       | |
| log  |    1..64G   | 105 | 1.9e-33 |4.8 | 0.6 |
|      |             |     |       | |
| sin, cos  | -2pi..2pi   | 106 | 1.2e-32 | 1.9 | 0.2 |
| tan  | -2pi..2pi   | 104 | 4.8e-32 | 2.3 | 0.4 |
| csc, sec, cot | -2pi..2pi | 104 | 4.8e-32 |  | |
|      |             |     |       | |
| asin, acos  | -1..1     | 106 | 8e-33 | 20 | 2.0 |
| atan, acot  | -1..1   | 106 | 1.2e-32 | 15 | 1.8 |
| acsc  | -64G..-1,1..64G | 106  | 2.5e-32  |  | |
| asec  | -4096..-1,1..4096 | 104  | 2.7e-32  |  | |
|      |             |     |       |
| sinh, cosh  | -1..1   | 104 | 4.5e-32 | 3.0 | 0.6 |
| tanh  | -1..1   | 104 | 4.5e-32 | 7.0 | 1.0 |
| csch, sech, coth  | -1..1| 104  | 4.6e-32 |  | |
|      |            |     |       |
| asinh  | -1..1     | 105 | 2.2e-32 | 15 | 3 |
| acosh  |  1..16     | 105 | 2.0e-32 | 20 | 3 |
| atanh  | -1..1   | 109 | 1.5e-33 | 9 | 2 |
| acsch, asech, acoth  | | ?  |  | |
|      |            |     |       |
|      |            |     |       |
| sin, cos  | -64Gpi..64Gpi   | 105 | 2.25-32 |  | |
| tan  | -64Gpi..64Gpi   | 103 | 5.9e-32 |  | |
| csc, sec, cot | -64Gpi..64Gpi | 103 | 8.8e-32 |  | |
|      |             |     |       |
| sinh, cosh, tanh  | -256..256   | 104 | 4.5e-32 |  | |
| csch, sech, coth  | -256..256   | 103 | 6.7e-32 |  | |

** 64G = 4096^3 **

##### outliers
| function | arg | eval | true |
|----------|-----|------|------|
| tanh     |(5.551170634277014e-17, -6.162975822039155e-33) | (5.551170634277014e-17, -6.162975822039155e-33) | (5.551170634277013e-17, 6.162975822039155e-33) |
| asinh     |(1.1102341268554029e-16, -1.232595164407831e-32) | (1.1102341268554029e-16, -1.232595164407831e-32) | (1.1102341268554026e-16, 1.232595164407831e-32) |

  
