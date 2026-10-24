# Cylindrical algebraic decomposition

## Full decomposition (seconds; cells / projection factors for Poly Mojo)

| Case | polymojo-lazard | polymojo-mccallum | polymojo-stream | qepcad-fullcad | PM/QEPCAD | Lazard cells | McCallum cells | factors (Lazard) | QEPCAD cells / factors | sample points checked |
|---|---|---|---|---|---|---|---|---|---|---|
| bank-a-real-implicitization-problem | 0.0599 | error | 0.126 | 0.008 | 7.49 | 919 | — | [1, 1, 2, 3, 2] | 1141 / [1, 1, 2, 3, 2] | 307 |
| bank-arnon-84 | 0.0123 | 0.0107 | 0.0259 | 0.009 | 1.37 | 55 | 55 | [2, 1] | 66 / [2, 1] | 55 |
| bank-arnon-84-2 | 0.0143 | 0.0131 | 0.0255 | 0.005 | 2.86 | 61 | 61 | [5, 2] | 52 / [3, 2] | 61 |
| bank-ball-and-circular-cylinder | 0.0244 | 0.0196 | 0.0616 | 0.015 | 1.63 | 327 | 327 | [8, 5, 2] | 323 / [6, 4, 2] | 327 |
| bank-collins-and-johnson | 1.24 | 1.02 | 1.96 | 0.186 | 6.66 | 6223 | 6219 | [22, 9, 4] | 7025 / [22, 9, 4] | 389 |
| bank-collision-of-circle-and-square | 0.263 | 0.238 | 0.527 | 0.033 | 7.96 | 2283 | 2283 | [18, 6, 3] | 2083 / [12, 6, 3] | 381 |
| bank-concentric-circles | 0.00166 | 0.00165 | 0.00454 | 0.005 | 0.33 | 41 | 41 | [4, 2] | 50 / [4, 2] | 41 |
| bank-cyclic-3 | 0.0446 | 0.0445 | 0.111 | 0.014 | 3.19 | 381 | 381 | [5, 5, 3] | 453 / [5, 5, 3] | 381 |
| bank-cyclic-4 | timeout | error | timeout | 182 | — | — | — | — | 104820 / [20, 10, 7, 4] | — |
| bank-davenport-and-heintz | 0.311 | 0.279 | 0.633 | 0.044 | 7.06 | 4949 | 4949 | [4, 6, 5, 3] | 5696 / [4, 6, 5, 3] | 381 |
| bank-edges-square-product | out-of-memory | error | timeout | error | — | — | — | — | — | — |
| bank-hong-90 | 0.0279 | error | 0.061 | 0 | — | 639 | — | [1, 3, 3, 1, 1] | 41 / [0, 0, 1, 1, 1] | 320 |
| bank-mccallum-trivariate-random-polynomial | 97.6 | timeout | 126 | 7.13 | 13.69 | 1439 | — | [8, 3, 1] | 1087 / [6, 2, 1] | 4 |
| bank-non-concentric-circles | 0.00284 | 0.00254 | 0.00708 | 0.005 | 0.57 | 59 | 59 | [6, 2] | 50 / [5, 2] | 59 |
| bank-off-center-ellipse | 5.2 | 6.03 | 6.48 | 0.611 | 8.52 | 4593 | 4997 | [14, 6, 2] | 5381 / [12, 6, 2] | 383 |
| bank-parametric-parabola | 0.00459 | error | 0.0112 | 0.01 | 0.46 | 115 | — | [1, 1, 2, 1] | 162 / [1, 1, 2, 1] | 115 |
| bank-piano-mover-s-problem-davenport | timeout | timeout | timeout | error | — | — | — | — | — | — |
| bank-positivity-of-quartic | 0.205 | error | 0.304 | 0.019 | 10.80 | 333 | — | [1, 3, 1, 1] | 280 / [1, 2, 1, 1] | 333 |
| bank-simplified-edges-square-product | 11.2 | 11.9 | 18.3 | 3.73 | 3.00 | 59067 | 61973 | [66, 16, 7] | 59641 / [63, 16, 7] | 400 |
| bank-simplified-putnum-example | 50.4 | 61.8 | 79.5 | 0.142 | 355.13 | 130829 | 164221 | [26, 13, 6, 2] | 12810 / [7, 8, 5, 2] | 399 |
| bank-simplified-yangxia | 1.07 | 2.12 | 2 | 0.678 | 1.57 | 6313 | 9793 | [1, 9, 7, 5] | 6824 / [1, 9, 7, 5] | 395 |
| bank-solotareff-3 | timeout | timeout | timeout | 9.41 | — | — | — | — | 266451 / [8, 6, 3, 5, 4] | — |
| bank-sphere-and-catastrophe | 0.221 | 0.198 | 0.308 | 0.018 | 12.26 | 753 | 753 | [7, 3, 2] | 637 / [5, 3, 2] | 377 |
| bank-term-rewrite | 0.118 | 0.0937 | 0.198 | 0.019 | 6.20 | 1099 | 1099 | [7, 5, 2] | 1293 / [7, 5, 2] | 367 |
| bank-whitney-umbrella | 0.0522 | error | 0.115 | 0 | — | 907 | — | [1, 1, 2, 3, 2] | 1147 / [1, 1, 2, 3, 2] | 303 |
| bank-x-axis-ellipse-problem | 8.86 | error | 15.9 | 0.865 | 10.24 | 64625 | — | [5, 8, 5, 5, 2] | 82097 / [5, 8, 5, 5, 2] | 399 |

## Conjunctions (seconds; verdicts must agree)

| Case | polymojo-solve | z3 | qepcad-sat | PM/Z3 | verdicts | agree |
|---|---|---|---|---|---|---|
