# FLOAT125.jl  
### *(work in progress -- incomplete)*
100+ valid significand bits for elementary functions with conventionally small values.

  relative errors (*max found with 20,000 random double-double values each function*)


| func | over | good bits | rel err | 20,0000 |
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
| sinh, cosh, tanh  | -16..16   | 103 | 6e-32 |
| csch, sech, coth  | | accordingly  |  |
|      |            |     |       |
| asin, acos  | -1..1     | 106 | 2e-32 |
| atan  | -16..16   | 106 | 2e-32 |
| acsc, asec, acot  | | accordingly  |  |
|      |            |     |       |
| asinh  | -1024..1024     | 106 | 2e-32 |
| acosh  |  1..1024     | 106 | 2e-32 |
| atanh  | -1..1   | 106 | 2e-32 |
| acsch, asech, acoth  | | accordingly  |  |

** 64G = 4096^3 **
