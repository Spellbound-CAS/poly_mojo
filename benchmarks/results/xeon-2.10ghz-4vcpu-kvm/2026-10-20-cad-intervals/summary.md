# CAD output intervals: precision floor (B) and canonical sectors (C)

Owner request: B, refine every irrational output coordinate to width at most
2^-K with K large enough that no interval is wider than nightly's on any case
nightly finishes; C, sector samples are the simplest rational of each gap.
Cost measured separately and together against the draft (31c157c) and
nightly (5a7e078).

- K = 830: nightly's narrowest output interval over the 17 full CAD cases it
  finishes is 2^-830 (collision-of-circle-and-square; term-rewrite 2^-693,
  sphere-n3-h2 2^-780).
- C: `_sector` returns the simplest rational of the true open gap between the
  neighbouring roots, found by refining until the inner and outer intervals
  agree, so it does not depend on isolating-interval depth.
- B: applied to output copies only (`cad.cells`, `cad.sample`); the
  certificate keeps the lifting intervals. An earlier variant that refined the
  stack samples themselves was 9x slower overall, because the 830-bit samples
  then fed every later lifting step.

Lazard decomposition, one CPU per variant, 300 s cap, 6 GB address limit.
Times are `cad.decompose` plus `cad.cells` (B's refinement happens in cells);
nightly times are `cad.decompose` only, 120 s cap, from its own run. Peak RSS
in MiB. cad.check is true on every finished case.

```
case                                nightly    draft        C        B      B+C | MiB nght drft    C    B  B+C
bank-a-real-implicitization-problem    0.227    0.230    0.292    3.029    3.055 |       53   53   52   53   53
bank-ball-and-circular-cylinder       0.143    0.220    0.123    1.033    1.055 |       29   28   28   29   29
bank-collision-of-circle-and-square    2.470    2.538    1.542    8.575    8.172 |      110  111  110  111  110
bank-concentric-circles               0.007    0.008    0.008    0.178    0.168 |       15   15   16   15   15
bank-cyclic-3                         0.545    0.354    0.222    2.402    1.747 |       28   28   28   28   28
bank-davenport-and-heintz             0.779    0.965    0.851    9.822   11.801 |      183  184  183  183  183
bank-hong-90                          0.087    0.146    0.102    0.091    0.095 |       51   51   51   51   51
bank-non-concentric-circles           0.011    0.012    0.013    0.263    0.259 |       16   16   16   16   16
bank-parametric-parabola              0.027    0.020    0.022    0.270    0.277 |       18   18   18   18   18
bank-term-rewrite                     1.882    2.108    0.992    7.822    7.124 |       47   48   47   48   47
bank-whitney-umbrella                 0.161    0.238    0.159    0.586    0.634 |       50   50   50   50   50
dense-n2-k1-d2                        0.007    0.010    0.008    0.147    0.102 |       15   15   15   15   15
dense-n2-k2-d2                        0.025    0.034    0.019    0.310    0.144 |       17   16   16   16   16
sphere-n2-h1                          0.024    0.035    0.018    0.285    0.229 |       16   16   15   16   16
sphere-n2-h2                          0.038    0.033    0.018    0.289    0.261 |       16   17   17   17   17
sphere-n3-h1                          0.742    0.785    0.367    3.821    3.621 |       31   32   31   31   30
sphere-n3-h2                          2.471    2.343    1.656   15.739   16.021 |      101  101   99  101   99
bank-simplified-yangxia                   -   20.552    6.915  183.893  161.381 |        -  433  428  434  428
bank-sphere-and-catastrophe               -    5.299    1.348   15.956   13.893 |        -   48   45   48   45
dense-n2-k1-d4                            -    1.470    0.927    7.193    7.038 |        -   35   35   35   35
dense-n2-k2-d3                            -    3.815    3.161   45.410   46.409 |        -   33   33   33   33
dense-n3-k1-d2                            -    1.882    1.136    5.702    5.327 |        -   27   27   27   27
dense-n2-k1-d3                            -    0.115    0.066    1.906    1.710 |        -   17   16   16   16
dense-n2-k3-d2                            -    0.465    0.289    9.274    8.575 |        -   28   28   28   28
dense-n2-k4-d2                            -    1.449    1.130   32.353   32.181 |        -   49   48   48   48
bank-arnon-84                             -    0.168    0.056    2.244    2.251 |        -   17   17   17   16
bank-arnon-84-2                           -    0.098    0.048    1.289    1.154 |        -   17   17   17   17
bank-collins-and-johnson                  -   11.125    7.259   27.352   22.197 |        -  291  278  291  278
bank-positivity-of-quartic                -    6.805    1.668   38.253   26.767 |        -   45   33   45   33
dense-n2-k3-d3                            -   19.009   18.659       TO       TO |        -   70   69    -    -
sphere-n4-h1                              -  107.883   21.400  180.215  123.459 |        -  411  376  411  376
Sum, 30 cases all variants finish: R 171.2 Cx2 51.8 Bx2 605.7 BC2 507.1
Sum, the 17 cases nightly finishes: nightly 9.6 R 10.1 Cx2 6.4 Bx2 54.7 BC2 54.8
```

Interval gate (`python -m harness.cad_gate NIGHTLY_ANSWERS BC_ANSWERS`): same
cells, dimensions and signs on all 17; no interval wider than nightly's on 16.
sphere-n3-h2 fails on 14 cells, where nightly's sector sample x0 = 2/3 made a
section coordinate rational (x1 = -1/3) and the canonical x0 = 1/2 makes it
irrational, refined to 2^-830 but no longer exact. Every other kind change
goes the other way (nightly algebraic, canonical exact) and passes.

## B with quadratic interval refinement

`polynomials.refine` keeps its Sturm validation but tightens by Abbott's
quadratic interval refinement (`_real_refine._tighten`) instead of bisecting
one bit per step, so B's 830 bits cost a few dozen steps. Same 31 cases,
lazard, decompose + cells, one CPU per run, merged with nightly 4db017f:

| Sum over | draft | C | B | B+C |
|---|---|---|---|---|
| 31 cases (all finish in every variant) | 190.2 s | 70.5 s | 196.1 s | 76.4 s |
| 17 cases nightly finishes (nightly 9.6 s) | 10.1 s | 6.4 s | 10.2 s | 6.6 s |

B now costs about 3 % (it was 3.5x with bisection: B+C 507 s); B+C is 2.5x
faster than the draft. Per case (draft / C / B / B+C, s):
simplified-yangxia 20.6 / 6.9 / 21.4 / 7.4, sphere-n4-h1 107.9 / 21.4 /
110.7 / 23.1, dense-n2-k3-d3 19.0 / 18.7 / 21.0 / 20.7, positivity-of-quartic
6.8 / 1.7 / 6.7 / 1.9, collision 2.5 / 1.5 / 2.5 / 1.5. The gate against
nightly is unchanged (16 of 17; sphere-n3-h2 as above).
