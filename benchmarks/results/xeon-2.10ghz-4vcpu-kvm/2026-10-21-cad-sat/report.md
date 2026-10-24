# Cylindrical algebraic decomposition

## Full decomposition (seconds; cells / projection factors for Poly Mojo)

| Case | polymojo-lazard | polymojo-mccallum | polymojo-stream | qepcad-fullcad | PM/QEPCAD | Lazard cells | McCallum cells | factors (Lazard) | QEPCAD cells / factors | sample points checked |
|---|---|---|---|---|---|---|---|---|---|---|

## Conjunctions (seconds; verdicts must agree)

| Case | polymojo-solve | z3 | qepcad-sat | PM/Z3 | verdicts | agree |
|---|---|---|---|---|---|---|
| bank-a-real-implicitization-problem | 0.000223 | 0.000985 | — | 0.23 | sat/sat | yes |
| bank-arnon-84 | 0.000388 | 0.00271 | — | 0.14 | sat/sat | yes |
| bank-arnon-84-2 | 0.000299 | 0.00132 | — | 0.23 | sat/sat | yes |
| bank-ball-and-circular-cylinder | 0.00463 | 0.00229 | — | 2.02 | sat/sat | yes |
| bank-collins-and-johnson | 0.00327 | 0.00266 | — | 1.23 | sat/sat | yes |
| bank-collision-of-circle-and-square | 0.000158 | 0.00124 | — | 0.13 | unsat/unsat | yes |
| bank-cyclic-3 | 0.000821 | 0.00125 | — | 0.66 | unsat/unsat | yes |
| bank-cyclic-4 | 0.00154 | 0.00141 | — | 1.09 | sat/sat | yes |
| bank-edges-square-product | 0.00882 | 0.00147 | — | 6.00 | sat/sat | yes |
| bank-hong-90 | 0.000224 | 0.000991 | — | 0.23 | sat/sat | yes |
| bank-mccallum-trivariate-random-polynomial | 0.000284 | 0.00121 | — | 0.23 | sat/sat | yes |
| bank-parametric-parabola | 0.000122 | 0.000884 | — | 0.14 | sat/sat | yes |
| bank-simplified-edges-square-product | 0.00168 | 0.00167 | — | 1.00 | sat/sat | yes |
| bank-simplified-putnum-example | 0.00294 | 0.00116 | — | 2.54 | sat/sat | yes |
| bank-simplified-yangxia | 0.0547 | 0.0364 | — | 1.50 | unsat/unsat | yes |
| bank-solotareff-3 | 0.0365 | 0.00154 | — | 23.67 | unsat/unsat | yes |
| bank-sphere-and-catastrophe | 0.000197 | 0.00122 | — | 0.16 | sat/sat | yes |
| bank-whitney-umbrella | 0.000169 | 0.000985 | — | 0.17 | sat/sat | yes |
| hong-n2-sat | 0.000552 | 0.00126 | — | 0.44 | sat/sat | yes |
| hong-n2-unsat | 0.000873 | 0.00137 | — | 0.64 | unsat/unsat | yes |
| hong-n3-sat | 0.00137 | 0.00149 | — | 0.92 | sat/sat | yes |
| hong-n3-unsat | 0.0044 | 0.00327 | — | 1.34 | unsat/unsat | yes |
| hong-n4-sat | 0.00292 | 0.00165 | — | 1.78 | sat/sat | yes |
| hong-n4-unsat | 0.0121 | 0.0211 | — | 0.57 | unsat/unsat | yes |
| motzkin-sat | 0.000512 | 0.00202 | — | 0.25 | sat/sat | yes |
| motzkin-unsat | 0.00106 | 0.00329 | — | 0.32 | unsat/unsat | yes |
| sos-d4-sat | 0.000342 | 0.0121 | — | 0.03 | sat/sat | yes |
| sos-d4-unsat | 0.0689 | 0.181 | — | 0.38 | unsat/unsat | yes |
| sos-d6-sat | 0.000354 | 0.0221 | — | 0.02 | sat/sat | yes |
| sos-d6-unsat | timeout | timeout | — | — |  | yes |
| sos-d8-sat | skipped | skipped | — | — |  | yes |
| sos-d8-unsat | skipped | skipped | — | — |  | yes |
| sphere-plane-n2-sat | 0.000184 | 0.00109 | — | 0.17 | sat/sat | yes |
| sphere-plane-n2-unsat | 0.000392 | 0.00105 | — | 0.37 | unsat/unsat | yes |
| sphere-plane-n3-sat | 0.000179 | 0.00108 | — | 0.17 | sat/sat | yes |
| sphere-plane-n3-unsat | 0.00128 | 0.00117 | — | 1.09 | unsat/unsat | yes |
| sphere-plane-n4-sat | 0.000221 | 0.00112 | — | 0.20 | sat/sat | yes |
| sphere-plane-n4-unsat | 0.00239 | 0.00118 | — | 2.03 | unsat/unsat | yes |
| sphere-plane-n5-sat | 0.000295 | 0.00125 | — | 0.24 | sat/sat | yes |
| sphere-plane-n5-unsat | 0.00451 | 0.00125 | — | 3.61 | unsat/unsat | yes |
