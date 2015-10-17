# FLOAT125.jl
100+ valid significand bits for elementary functions with conventionally small values.

  relative errors (max found in range)\\
  (24,000 rand valuess in given range)

                          good bits relerr

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
| sin  | -2pi..2pi   | 103 | 1e-31 |
| cos  | -2pi..2pi   | 103 | 1e-31 |
| tan  | -2pi..2pi   | 103 | 1e-31 |
| csc, sec, cot  | | accordingly  |  |
|      |             |     |       |
| sin  | -8pi..8pi   | 100 | 6e-31 |
| cos  | -8pi..8pi   | 100 | 6e-31 |
| tan  | -8pi..8pi   | 10- | 6e-31 |
| csc, sec, cot  | | accordingly  |  |
|      |             |     |       |
| sinh  | -16..16   | 103 | 6e-32 |
| cosh  | -16..16   | 103 | 6e-32 |
| tanh  | -16..16   | 103 | 6e-32 |
| csch, sech, coth  | | accordingly  |  |
|      |            |     |       |
| asin  | -1..1     | 106 | 2e-32 |
| acos  | -1..1     | 106 | 2e-32 |
| atan  | -16..16   | 106 | 2e-32 |
| acsc, asec, acot  | | accordingly  |  |
|      |            |     |       |
| asinh  | -1024..1024     | 106 | 2e-32 |
| acosh  |  1..1024     | 106 | 2e-32 |
| atanh  | -1..1   | 106 | 2e-32 |
| acsch, asech, acoth  | | accordingly  |  |

