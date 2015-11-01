# FLOAT128.jl  
### *(work in progress -- incomplete)*
100+ valid significand bits for elementary functions with conventionally small values.

  relative errors (*max found with 20,000 random double-double values each function*)


| func | over | sig bits | rel err | checked |
|------|------|-----------|---------|---------|
| sqrt | 0..64G | 106 | 1.2e-32 | t |
|      |             |     |       |  |
| exp  | -15..15   | 104 | 4.9e-32 |t  |
| exp  | -300..300   | 103 | 5.2e-32 |t |
|      |             |     |       | |
| log  |    1..64G   | 105 | 1.9e-33 |t  |
|      |             |     |       | |
| sin, cos  | -2pi..2pi   | 106 | 1.2e-32 | t |
| tan  | -2pi..2pi   | 104 | 4.8e-32 | t |
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
| sinh, cosh  | -1..1   | 104 | 4.5e-32 | t |
| tanh  | -1..1   |  | ? | |
| sinh, cosh  | -256..256   | 104 | 4.5e-32 | t |
| tanh  | -1024..1024   | ? | ? | |
| csch, sech, coth  | | ?  | ? |
|      |            |     |       |
|      |            |     |       |
| asinh  | -1..1     | ? | ? |
| acosh  |  1..16     | ? | ? |
| atanh  | -1..1   | ? | ? |
| acsch, asech, acoth  | | ?  |  |

** 64G = 4096^3 **
