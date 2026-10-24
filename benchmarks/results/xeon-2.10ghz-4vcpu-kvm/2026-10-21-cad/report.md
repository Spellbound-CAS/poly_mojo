# Cylindrical algebraic decomposition

## Full decomposition (seconds; cells / projection factors for Poly Mojo)

| Case | polymojo-lazard | polymojo-mccallum | polymojo-stream | qepcad-fullcad | PM/QEPCAD | Lazard cells | McCallum cells | factors (Lazard) | QEPCAD cells / factors | sample points checked |
|---|---|---|---|---|---|---|---|---|---|---|
| bank-a-real-implicitization-problem | 0.131 | error | 0.29 | 0.018 | 7.27 | 919 | — | [1, 1, 2, 3, 2] | 1141 / [1, 1, 2, 3, 2] | 307 |
| bank-arnon-84 | 0.0465 | 0.0464 | 0.3 | 0.009 | 5.17 | 55 | 55 | [2, 1] | 66 / [2, 1] | 55 |
| bank-arnon-84-2 | 0.0434 | 0.034 | 0.0688 | 0 | — | 61 | 61 | [5, 2] | 52 / [3, 2] | 61 |
| bank-ball-and-circular-cylinder | 0.0806 | 0.0734 | 0.141 | 0 | — | 327 | 327 | [8, 5, 2] | 323 / [6, 4, 2] | 327 |
| bank-collins-and-johnson | 5.73 | 5.58 | 6.63 | 0.12 | 47.76 | 6223 | 6219 | [22, 9, 4] | 7025 / [22, 9, 4] | 389 |
| bank-collision-of-circle-and-square | 1.42 | 1.42 | 1.23 | 0.011 | 129.40 | 2283 | 2283 | [18, 6, 3] | 2083 / [12, 6, 3] | 381 |
| bank-concentric-circles | 0.00377 | 0.00375 | 0.00949 | 0.018 | 0.21 | 41 | 41 | [4, 2] | 50 / [4, 2] | 41 |
| bank-cyclic-3 | 0.207 | 0.213 | 0.169 | 0 | — | 381 | 381 | [5, 5, 3] | 453 / [5, 5, 3] | 381 |
| bank-cyclic-4 | timeout | error | timeout | 167 | — | — | — | — | 104820 / [20, 10, 7, 4] | — |
| bank-davenport-and-heintz | 0.625 | 0.611 | 1.29 | 0.029 | 21.54 | 4949 | 4949 | [4, 6, 5, 3] | 5696 / [4, 6, 5, 3] | 381 |
| bank-edges-square-product | out-of-memory | error | timeout | error | — | — | — | — | — | — |
| bank-hong-90 | 0.0782 | error | 0.197 | 0.013 | 6.02 | 639 | — | [1, 3, 3, 1, 1] | 41 / [0, 0, 1, 1, 1] | 320 |
| bank-mccallum-trivariate-random-polynomial | timeout | timeout | timeout | 6.2 | — | — | — | — | 1087 / [6, 2, 1] | — |
| bank-non-concentric-circles | 0.00882 | 0.00878 | 0.0174 | 0 | — | 59 | 59 | [6, 2] | 50 / [5, 2] | 59 |
| bank-off-center-ellipse | 126 | 167 | 250 | 0.429 | 292.57 | 4593 | 4997 | [14, 6, 2] | 5381 / [12, 6, 2] | 383 |
| bank-parametric-parabola | 0.0122 | error | 0.0328 | 0 | — | 115 | — | [1, 1, 2, 1] | 162 / [1, 1, 2, 1] | 115 |
| bank-piano-mover-s-problem-davenport | timeout | timeout | timeout | error | — | — | — | — | — | — |
| bank-positivity-of-quartic | 1.34 | error | 1 | 0 | — | 333 | — | [1, 3, 1, 1] | 280 / [1, 2, 1, 1] | 333 |
| bank-simplified-edges-square-product | timeout | timeout | 69.5 | 3.69 | — | — | — | — | 59641 / [63, 16, 7] | — |
| bank-simplified-putnum-example | timeout | timeout | timeout | 0.118 | — | — | — | — | 12810 / [7, 8, 5, 2] | — |
| bank-simplified-yangxia | 5.83 | 8.61 | timeout | 0.554 | 10.52 | 6313 | 9793 | [1, 9, 7, 5] | 6824 / [1, 9, 7, 5] | 395 |
| bank-solotareff-3 | timeout | timeout | timeout | 10.9 | — | — | — | — | 266451 / [8, 6, 3, 5, 4] | — |
| bank-sphere-and-catastrophe | 1.09 | 1.31 | 0.928 | 0 | — | 753 | 753 | [7, 3, 2] | 637 / [5, 3, 2] | 377 |
| bank-term-rewrite | 0.792 | 0.747 | 0.477 | 0.012 | 66.02 | 1099 | 1099 | [7, 5, 2] | 1293 / [7, 5, 2] | 367 |
| bank-whitney-umbrella | 0.0871 | error | 0.237 | 0.021 | 4.15 | 907 | — | [1, 1, 2, 3, 2] | 1147 / [1, 1, 2, 3, 2] | 303 |
| bank-x-axis-ellipse-problem | 39.1 | error | 49.5 | 0.695 | 56.28 | 64625 | — | [5, 8, 5, 5, 2] | 82097 / [5, 8, 5, 5, 2] | 399 |

## Conjunctions (seconds; verdicts must agree)

| Case | polymojo-solve | z3 | qepcad-sat | PM/Z3 | verdicts | agree |
|---|---|---|---|---|---|---|
| bank-a-real-implicitization-problem | 0.00124 | 0.000901 | 0.015 | 1.38 | sat/sat/sat | yes |
| bank-arnon-84 | 0.00466 | 0.00248 | 0.015 | 1.88 | sat/sat/sat | yes |
| bank-arnon-84-2 | 0.00567 | 0.00103 | 0.016 | 5.52 | sat/sat/sat | yes |
| bank-ball-and-circular-cylinder | 0.0212 | 0.00174 | 0.002 | 12.21 | sat/sat/sat | yes |
| bank-collins-and-johnson | 0.0234 | 0.00217 | 0 | 10.81 | sat/sat/sat | yes |
| bank-collision-of-circle-and-square | 0.041 | 0.00111 | 0 | 36.81 | unsat/unsat/unsat | yes |
| bank-cyclic-3 | 0.0657 | 0.00101 | 0 | 64.80 | unsat/unsat/unsat | yes |
| bank-cyclic-4 | 1.11 | 0.00112 | 0.003 | 994.53 | sat/sat/sat | yes |
| bank-edges-square-product | 2.71 | 0.00118 | 0.093 | 2288.85 | sat/sat/sat | yes |
| bank-hong-90 | 0.00116 | 0.000907 | 0.022 | 1.28 | sat/sat/sat | yes |
| bank-mccallum-trivariate-random-polynomial | 0.71 | 0.000999 | 0 | 710.35 | sat/sat/sat | yes |
| bank-parametric-parabola | 0.000486 | 0.000855 | 0 | 0.57 | sat/sat/sat | yes |
| bank-simplified-edges-square-product | 0.033 | 0.0013 | 0 | 25.30 | sat/sat/sat | yes |
| bank-simplified-putnum-example | 0.169 | 0.000934 | 0 | 180.93 | sat/sat/sat | yes |
| bank-simplified-yangxia | 0.551 | 0.0333 | 0.539 | 16.57 | unsat/unsat/unsat | yes |
| bank-solotareff-3 | 68.5 | 0.00122 | 0.08 | 56069.15 | unsat/unsat/unsat | yes |
| bank-sphere-and-catastrophe | 0.00689 | 0.00117 | 0 | 5.90 | sat/sat/sat | yes |
| bank-whitney-umbrella | 0.000996 | 0.000832 | 0.017 | 1.20 | sat/sat/sat | yes |
| hong-n2-sat | 0.00426 | 0.00097 | 0.017 | 4.39 | sat/sat/sat | yes |
| hong-n2-unsat | 0.00272 | 0.00133 | 0.013 | 2.05 | unsat/unsat/unsat | yes |
| hong-n3-sat | 0.0197 | 0.00141 | 0 | 13.98 | sat/sat/sat | yes |
| hong-n3-unsat | 0.0272 | 0.0027 | 0 | 10.08 | unsat/unsat/unsat | yes |
| hong-n4-sat | 0.131 | 0.00139 | 0.052 | 94.64 | sat/sat/sat | yes |
| hong-n4-unsat | 0.176 | 0.0184 | 0 | 9.56 | unsat/unsat/unsat | yes |
| motzkin-sat | 0.00248 | 0.00163 | 0 | 1.52 | sat/sat/sat | yes |
| motzkin-unsat | 0.00296 | 0.00293 | 0.014 | 1.01 | unsat/unsat/unsat | yes |
| sos-d4-sat | 8.51 | 0.00951 | 2.75 | 894.87 | sat/sat/sat | yes |
| sos-d4-unsat | 10.3 | 0.125 | 0.044 | 82.54 | unsat/unsat/unsat | yes |
| sos-d6-sat | timeout | 0.0198 | error | — | sat | yes |
| sos-d6-unsat | skipped | timeout | 21.2 | — | unsat | yes |
| sos-d8-sat | skipped | skipped | timeout | — |  | yes |
| sos-d8-unsat | skipped | skipped | skipped | — |  | yes |
| sphere-plane-n2-sat | 0.000719 | 0.00105 | 0 | 0.69 | sat/sat/sat | yes |
| sphere-plane-n2-unsat | 0.00144 | 0.000923 | 0 | 1.56 | unsat/unsat/unsat | yes |
| sphere-plane-n3-sat | 0.00189 | 0.000968 | 0 | 1.95 | sat/sat/sat | yes |
| sphere-plane-n3-unsat | 0.00888 | 0.00103 | 0 | 8.61 | unsat/unsat/unsat | yes |
| sphere-plane-n4-sat | 0.00364 | 0.00102 | 0 | 3.56 | sat/sat/sat | yes |
| sphere-plane-n4-unsat | 0.0432 | 0.00103 | 0 | 42.09 | unsat/unsat/unsat | yes |
| sphere-plane-n5-sat | 0.00813 | 0.00105 | 0.016 | 7.72 | sat/sat/sat | yes |
| sphere-plane-n5-unsat | 0.245 | 0.00106 | 0.013 | 231.51 | unsat/unsat/unsat | yes |
