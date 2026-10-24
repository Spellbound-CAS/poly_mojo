# Cylindrical algebraic decomposition

## Full decomposition (seconds; cells / projection factors for Poly Mojo)

| Case | polymojo-lazard | polymojo-mccallum | polymojo-stream | qepcad-fullcad | PM/QEPCAD | Lazard cells | McCallum cells | factors (Lazard) | QEPCAD cells / factors | sample points checked |
|---|---|---|---|---|---|---|---|---|---|---|

## Conjunctions (seconds; verdicts must agree)

| Case | polymojo-solve | z3 | qepcad-sat | PM/Z3 | verdicts | agree |
|---|---|---|---|---|---|---|
| bank-a-real-implicitization-problem | 0.000135 | 0.000849 | 0.006 | 0.16 | sat/sat/sat | yes |
| bank-arnon-84 | 0.000291 | 0.00244 | 0 | 0.12 | sat/sat/sat | yes |
| bank-arnon-84-2 | 0.000153 | 0.00115 | 0.003 | 0.13 | sat/sat/sat | yes |
| bank-ball-and-circular-cylinder | 0.00264 | 0.00209 | 0.007 | 1.26 | sat/sat/sat | yes |
| bank-collins-and-johnson | 0.00204 | 0.00249 | 0 | 0.82 | sat/sat/sat | yes |
| bank-collision-of-circle-and-square | 0.00014 | 0.00111 | 0.009 | 0.13 | unsat/unsat/unsat | yes |
| bank-cyclic-3 | 0.000742 | 0.00112 | 0.02 | 0.66 | unsat/unsat/unsat | yes |
| bank-cyclic-4 | 0.00113 | 0.00124 | 0.017 | 0.91 | sat/sat/sat | yes |
| bank-edges-square-product | 0.00731 | 0.00132 | 0.151 | 5.52 | sat/sat/sat | yes |
| bank-hong-90 | 0.00013 | 0.000864 | 0 | 0.15 | sat/sat/sat | yes |
| bank-mccallum-trivariate-random-polynomial | 0.000178 | 0.00114 | 0.012 | 0.16 | sat/sat/sat | yes |
| bank-parametric-parabola | 8.22e-05 | 0.00079 | 0 | 0.10 | sat/sat/sat | yes |
| bank-simplified-edges-square-product | 0.00122 | 0.00149 | 0.013 | 0.82 | sat/sat/sat | yes |
| bank-simplified-putnum-example | 0.00247 | 0.00123 | 0.003 | 2.00 | sat/sat/sat | yes |
| bank-simplified-yangxia | 0.0309 | 0.0337 | 0.581 | 0.92 | unsat/unsat/unsat | yes |
| bank-solotareff-3 | 0.0203 | 0.00138 | 0.105 | 14.69 | unsat/unsat/unsat | yes |
| bank-sphere-and-catastrophe | 0.000113 | 0.00111 | 0.009 | 0.10 | sat/sat/sat | yes |
| bank-whitney-umbrella | 9.7e-05 | 0.00083 | 0 | 0.12 | sat/sat/sat | yes |
| hong-n2-sat | 0.000326 | 0.00112 | 0 | 0.29 | sat/sat/sat | yes |
| hong-n2-unsat | 0.000556 | 0.00124 | 0 | 0.45 | unsat/unsat/unsat | yes |
| hong-n3-sat | 0.000899 | 0.0013 | 0 | 0.69 | sat/sat/sat | yes |
| hong-n3-unsat | 0.0022 | 0.00284 | 0 | 0.77 | unsat/unsat/unsat | yes |
| hong-n4-sat | 0.00213 | 0.00145 | 0.048 | 1.47 | sat/sat/sat | yes |
| hong-n4-unsat | 0.00763 | 0.0204 | 0.012 | 0.37 | unsat/unsat/unsat | yes |
| motzkin-sat | 0.000421 | 0.00175 | 0 | 0.24 | sat/sat/sat | yes |
| motzkin-unsat | 0.000872 | 0.00282 | 0 | 0.31 | unsat/unsat/unsat | yes |
| sos-d4-sat | 0.000183 | 0.0102 | 3.13 | 0.02 | sat/sat/sat | yes |
| sos-d4-unsat | 0.0481 | 0.153 | 0.08 | 0.31 | unsat/unsat/unsat | yes |
| sos-d6-sat | 0.000179 | 0.0197 | error | 0.01 | sat/sat | yes |
| sos-d6-unsat | 33.3 | timeout | 22.8 | — | unsat/unsat | yes |
| sos-d8-sat | 0.000186 | skipped | timeout | — | sat | yes |
| sos-d8-unsat | 121 | skipped | skipped | — | unsat | yes |
| sphere-plane-n2-sat | 0.000138 | 0.00105 | 0.005 | 0.13 | sat/sat/sat | yes |
| sphere-plane-n2-unsat | 0.000327 | 0.000997 | 0 | 0.33 | unsat/unsat/unsat | yes |
| sphere-plane-n3-sat | 0.000118 | 0.00106 | 0.007 | 0.11 | sat/sat/sat | yes |
| sphere-plane-n3-unsat | 0.00111 | 0.00108 | 0 | 1.03 | unsat/unsat/unsat | yes |
| sphere-plane-n4-sat | 0.000141 | 0.00111 | 0 | 0.13 | sat/sat/sat | yes |
| sphere-plane-n4-unsat | 0.00219 | 0.00115 | 0 | 1.90 | unsat/unsat/unsat | yes |
| sphere-plane-n5-sat | 0.00017 | 0.00114 | 0.01 | 0.15 | sat/sat/sat | yes |
| sphere-plane-n5-unsat | 0.00415 | 0.00124 | 0.008 | 3.34 | unsat/unsat/unsat | yes |
