# Cylindrical algebraic decomposition

## Full decomposition (seconds; cells / projection factors for Poly Mojo)

| Case | polymojo-lazard | polymojo-mccallum | polymojo-stream | qepcad-fullcad | PM/QEPCAD | Lazard cells | McCallum cells | factors (Lazard) | QEPCAD cells / factors | sample points checked |
|---|---|---|---|---|---|---|---|---|---|---|
| bank-a-real-implicitization-problem | 0.0615 | error | 0.12 | 0.008 | 7.69 | 919 | — | [1, 1, 2, 3, 2] | 1141 / [1, 1, 2, 3, 2] | 307 |
| bank-arnon-84 | 0.013 | 0.0127 | 0.0233 | 0.009 | 1.45 | 55 | 55 | [2, 1] | 66 / [2, 1] | 55 |
| bank-arnon-84-2 | 0.0164 | 0.0158 | 0.0281 | 0.005 | 3.28 | 61 | 61 | [5, 2] | 52 / [3, 2] | 61 |
| bank-ball-and-circular-cylinder | 0.0261 | 0.0253 | 0.0475 | 0.015 | 1.74 | 327 | 327 | [8, 5, 2] | 323 / [6, 4, 2] | 327 |
| bank-collins-and-johnson | 1.34 | 1.55 | 2.02 | 0.186 | 7.21 | 6223 | 6219 | [22, 9, 4] | 7025 / [22, 9, 4] | 389 |
| bank-collision-of-circle-and-square | 0.289 | 0.287 | 0.465 | 0.033 | 8.75 | 2283 | 2283 | [18, 6, 3] | 2083 / [12, 6, 3] | 381 |
| bank-concentric-circles | 0.00194 | 0.00189 | 0.00473 | 0.005 | 0.39 | 41 | 41 | [4, 2] | 50 / [4, 2] | 41 |
| bank-cyclic-3 | 0.0452 | 0.0428 | 0.0775 | 0.014 | 3.23 | 381 | 381 | [5, 5, 3] | 453 / [5, 5, 3] | 381 |
| bank-cyclic-4 | timeout | error | timeout | 182 | — | — | — | — | 104820 / [20, 10, 7, 4] | — |
| bank-davenport-and-heintz | 0.251 | 0.221 | 0.563 | 0.044 | 5.71 | 4949 | 4949 | [4, 6, 5, 3] | 5696 / [4, 6, 5, 3] | 381 |
| bank-edges-square-product | out-of-memory | error | timeout | error | — | — | — | — | — | — |
| bank-hong-90 | 0.0304 | error | 0.0593 | 0 | — | 639 | — | [1, 3, 3, 1, 1] | 41 / [0, 0, 1, 1, 1] | 320 |
| bank-mccallum-trivariate-random-polynomial | timeout | timeout | timeout | 7.13 | — | — | — | — | 1087 / [6, 2, 1] | — |
| bank-non-concentric-circles | 0.0029 | 0.00291 | 0.0077 | 0.005 | 0.58 | 59 | 59 | [6, 2] | 50 / [5, 2] | 59 |
| bank-off-center-ellipse | 7.2 | 8.8 | 12 | 0.611 | 11.78 | 4593 | 4997 | [14, 6, 2] | 5381 / [12, 6, 2] | 383 |
| bank-parametric-parabola | 0.00455 | error | 0.0183 | 0.01 | 0.45 | 115 | — | [1, 1, 2, 1] | 162 / [1, 1, 2, 1] | 115 |
| bank-piano-mover-s-problem-davenport | timeout | timeout | timeout | error | — | — | — | — | — | — |
| bank-positivity-of-quartic | 0.273 | error | 0.412 | 0.019 | 14.39 | 333 | — | [1, 3, 1, 1] | 280 / [1, 2, 1, 1] | 333 |
| bank-simplified-edges-square-product | 11.1 | 11 | 16.9 | 3.73 | 2.98 | 59067 | 61973 | [66, 16, 7] | 59641 / [63, 16, 7] | 400 |
| bank-simplified-putnum-example | 55.7 | 75.8 | 80 | 0.142 | 392.13 | 130829 | 164221 | [26, 13, 6, 2] | 12810 / [7, 8, 5, 2] | 399 |
| bank-simplified-yangxia | 0.988 | 1.82 | 1.95 | 0.678 | 1.46 | 6313 | 9793 | [1, 9, 7, 5] | 6824 / [1, 9, 7, 5] | 395 |
| bank-solotareff-3 | timeout | timeout | timeout | 9.41 | — | — | — | — | 266451 / [8, 6, 3, 5, 4] | — |
| bank-sphere-and-catastrophe | 0.252 | 0.181 | 0.298 | 0.018 | 13.98 | 753 | 753 | [7, 3, 2] | 637 / [5, 3, 2] | 377 |
| bank-term-rewrite | 0.101 | 0.103 | 0.201 | 0.019 | 5.33 | 1099 | 1099 | [7, 5, 2] | 1293 / [7, 5, 2] | 367 |
| bank-whitney-umbrella | 0.0516 | error | 0.124 | 0 | — | 907 | — | [1, 1, 2, 3, 2] | 1147 / [1, 1, 2, 3, 2] | 303 |
| bank-x-axis-ellipse-problem | 8.9 | error | 15.6 | 0.865 | 10.29 | 64625 | — | [5, 8, 5, 5, 2] | 82097 / [5, 8, 5, 5, 2] | 399 |

## Conjunctions (seconds; verdicts must agree)

| Case | polymojo-solve | z3 | qepcad-sat | PM/Z3 | verdicts | agree |
|---|---|---|---|---|---|---|
