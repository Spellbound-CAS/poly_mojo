# Cylindrical algebraic decomposition

## Full decomposition (seconds; cells / projection factors for Poly Mojo)

| Case | polymojo-lazard | polymojo-mccallum | polymojo-stream | qepcad-fullcad | PM/QEPCAD | Lazard cells | McCallum cells | factors (Lazard) | QEPCAD cells / factors | sample points checked |
|---|---|---|---|---|---|---|---|---|---|---|
| bank-a-real-implicitization-problem | 0.348 | error | 0.64 | — | — | 919 | — | [1, 1, 2, 3, 2] | — | 307 |
| bank-arnon-84 | timeout | timeout | timeout | — | — | — | — | — | — | — |
| bank-arnon-84-2 | timeout | timeout | timeout | — | — | — | — | — | — | — |
| bank-ball-and-circular-cylinder | 0.203 | 0.185 | 0.364 | — | — | 327 | 327 | [8, 5, 2] | — | 327 |
| bank-collins-and-johnson | timeout | timeout | timeout | — | — | — | — | — | — | — |
| bank-collision-of-circle-and-square | 4.21 | 3.64 | 6.31 | — | — | 2283 | 2283 | [18, 6, 3] | — | 381 |
| bank-concentric-circles | 0.00627 | 0.00618 | 0.0149 | — | — | 41 | 41 | [4, 2] | — | 41 |
| bank-cyclic-3 | 2.37 | 2.15 | 3.48 | — | — | 381 | 381 | [5, 5, 3] | — | 381 |
| bank-cyclic-4 | timeout | timeout | timeout | — | — | — | — | — | — | — |
| bank-davenport-and-heintz | 1.56 | 1.53 | 3.11 | — | — | 4949 | 4949 | [4, 6, 5, 3] | — | 381 |
| bank-mccallum-trivariate-random-polynomial | timeout | timeout | timeout | — | — | — | — | — | — | — |
| bank-non-concentric-circles | 0.0157 | 0.0164 | 0.0373 | — | — | 59 | 59 | [6, 2] | — | 59 |
| bank-off-center-ellipse | timeout | timeout | timeout | — | — | — | — | — | — | — |
| bank-parametric-parabola | 0.043 | error | 0.0537 | — | — | 115 | — | [1, 1, 2, 1] | — | 115 |
| bank-piano-mover-s-problem-davenport | timeout | timeout | timeout | — | — | — | — | — | — | — |
| bank-positivity-of-quartic | timeout | error | timeout | — | — | — | — | — | — | — |
| bank-simplified-edges-square-product | timeout | timeout | timeout | — | — | — | — | — | — | — |
| bank-simplified-putnum-example | timeout | timeout | timeout | — | — | — | — | — | — | — |
| bank-simplified-yangxia | timeout | timeout | timeout | — | — | — | — | — | — | — |
| bank-sphere-and-catastrophe | timeout | timeout | timeout | — | — | — | — | — | — | — |
| bank-term-rewrite | 2.33 | 2.37 | 4.52 | — | — | 1099 | 1099 | [7, 5, 2] | — | 367 |

## Conjunctions (seconds; verdicts must agree)

| Case | polymojo-solve | z3 | qepcad-sat | PM/Z3 | verdicts | agree |
|---|---|---|---|---|---|---|
