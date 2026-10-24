# Cylindrical algebraic decomposition

## Full decomposition (seconds; cells / projection factors for Poly Mojo)

| Case | polymojo-lazard | polymojo-mccallum | polymojo-stream | qepcad-fullcad | PM/QEPCAD | Lazard cells | McCallum cells | factors (Lazard) | QEPCAD cells / factors | sample points checked |
|---|---|---|---|---|---|---|---|---|---|---|

## Conjunctions (seconds; verdicts must agree)

| Case | polymojo-solve | z3 | qepcad-sat | PM/Z3 | verdicts | agree |
|---|---|---|---|---|---|---|
| bank-a-real-implicitization-problem | 0.000132 | 0.000849 | 0.006 | 0.16 | sat/sat/sat | yes |
| bank-arnon-84 | 0.000294 | 0.00244 | 0 | 0.12 | sat/sat/sat | yes |
| bank-arnon-84-2 | 0.000152 | 0.00115 | 0.003 | 0.13 | sat/sat/sat | yes |
| bank-ball-and-circular-cylinder | 0.00256 | 0.00209 | 0.007 | 1.23 | sat/sat/sat | yes |
| bank-collins-and-johnson | 0.00203 | 0.00249 | 0 | 0.81 | sat/sat/sat | yes |
| bank-collision-of-circle-and-square | 0.00014 | 0.00111 | 0.009 | 0.13 | unsat/unsat/unsat | yes |
| bank-cyclic-3 | 0.000664 | 0.00112 | 0.02 | 0.59 | unsat/unsat/unsat | yes |
| bank-cyclic-4 | 0.00109 | 0.00124 | 0.017 | 0.87 | sat/sat/sat | yes |
| bank-edges-square-product | 0.00706 | 0.00132 | 0.151 | 5.33 | sat/sat/sat | yes |
| bank-hong-90 | 0.000129 | 0.000864 | 0 | 0.15 | sat/sat/sat | yes |
| bank-mccallum-trivariate-random-polynomial | 0.000163 | 0.00114 | 0.012 | 0.14 | sat/sat/sat | yes |
| bank-parametric-parabola | 7.52e-05 | 0.00079 | 0 | 0.10 | sat/sat/sat | yes |
| bank-simplified-edges-square-product | 0.00122 | 0.00149 | 0.013 | 0.82 | sat/sat/sat | yes |
| bank-simplified-putnum-example | 0.00244 | 0.00123 | 0.003 | 1.98 | sat/sat/sat | yes |
| bank-simplified-yangxia | 0.0294 | 0.0337 | 0.581 | 0.87 | unsat/unsat/unsat | yes |
| bank-solotareff-3 | 0.0183 | 0.00138 | 0.105 | 13.20 | unsat/unsat/unsat | yes |
| bank-sphere-and-catastrophe | 0.000113 | 0.00111 | 0.009 | 0.10 | sat/sat/sat | yes |
| bank-whitney-umbrella | 9.8e-05 | 0.00083 | 0 | 0.12 | sat/sat/sat | yes |
| hong-n2-sat | 0.000325 | 0.00112 | 0 | 0.29 | sat/sat/sat | yes |
| hong-n2-unsat | 0.000497 | 0.00124 | 0 | 0.40 | unsat/unsat/unsat | yes |
| hong-n3-sat | 0.000867 | 0.0013 | 0 | 0.67 | sat/sat/sat | yes |
| hong-n3-unsat | 0.0021 | 0.00284 | 0 | 0.74 | unsat/unsat/unsat | yes |
| hong-n4-sat | 0.00189 | 0.00145 | 0.048 | 1.31 | sat/sat/sat | yes |
| hong-n4-unsat | 0.00746 | 0.0204 | 0.012 | 0.37 | unsat/unsat/unsat | yes |
| motzkin-sat | 0.000378 | 0.00175 | 0 | 0.22 | sat/sat/sat | yes |
| motzkin-unsat | 0.000836 | 0.00282 | 0 | 0.30 | unsat/unsat/unsat | yes |
| sos-d4-sat | 0.000167 | 0.0102 | 3.13 | 0.02 | sat/sat/sat | yes |
| sos-d4-unsat | 0.0443 | 0.153 | 0.08 | 0.29 | unsat/unsat/unsat | yes |
| sos-d6-sat | 0.000179 | 0.0197 | error | 0.01 | sat/sat | yes |
| sos-d6-unsat | 33 | timeout | 22.8 | — | unsat/unsat | yes |
| sos-d8-sat | 0.000187 | skipped | timeout | — | sat | yes |
| sos-d8-unsat | 110 | skipped | skipped | — | unsat | yes |
| sphere-plane-n2-sat | 0.000139 | 0.00105 | 0.005 | 0.13 | sat/sat/sat | yes |
| sphere-plane-n2-unsat | 0.00033 | 0.000997 | 0 | 0.33 | unsat/unsat/unsat | yes |
| sphere-plane-n3-sat | 0.000117 | 0.00106 | 0.007 | 0.11 | sat/sat/sat | yes |
| sphere-plane-n3-unsat | 0.00106 | 0.00108 | 0 | 0.98 | unsat/unsat/unsat | yes |
| sphere-plane-n4-sat | 0.000141 | 0.00111 | 0 | 0.13 | sat/sat/sat | yes |
| sphere-plane-n4-unsat | 0.00215 | 0.00115 | 0 | 1.87 | unsat/unsat/unsat | yes |
| sphere-plane-n5-sat | 0.000168 | 0.00114 | 0.01 | 0.15 | sat/sat/sat | yes |
| sphere-plane-n5-unsat | 0.004 | 0.00124 | 0.008 | 3.22 | unsat/unsat/unsat | yes |
