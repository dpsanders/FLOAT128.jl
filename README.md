# FLOAT125.jl  
### *(work in progress -- incomplete)*
100+ valid significand bits for elementary functions with conventionally small values.

  relative errors (*max found with 24,000 random double-double values each function*)


| func | over | good bits | rel err |
|------|------|-----------|---------|
| sqrt | 1e-42..4e15 | 107 | 3e-33 |
|      |             |     |       |
| exp  | -512..512   | 107 | 3e-33 |
| log  |    1..512   | 107 | 3e-33 |
|      |             |     |       |
| exp  | -512..512   | 107 | 3e-33 |
| log  |    1..512   | 107 | 3e-33 |
|      |             |     |       |
| sin, cos  | -2pi..2pi   | 104 | 5e-32 |
| tan  | -2pi..2pi   | 103 | 6e-32 |
| csc, sec, cot  | | accordingly  |  |
|      |             |     |       |
| sin, cos  | -8pi..8pi   | 101 | 4e-31 |
| tan  | -8pi..8pi   | 100 | 8e-31 |
| csc, sec, cot  | | accordingly  |  |
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

