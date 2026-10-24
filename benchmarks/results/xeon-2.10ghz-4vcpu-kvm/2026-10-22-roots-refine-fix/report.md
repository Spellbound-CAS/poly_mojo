# Real roots

Host: Intel(R) Xeon(R) Processor @ 2.10GHz, 4 vCPUs, 16 GiB RAM, Ubuntu 24.04.5 LTS. Runs are pinned to one CPU, capped at 12G and 120 s.

Each answer is checked against an Arb reference: the number of distinct real roots, and every interval holding exactly one reference root (`wrong` otherwise). The `Poly Mojo / best` column divides Poly Mojo's time by the fastest other library on the case, each at its own achieved width.

Widths differ by tool: Poly Mojo and SymPy return intervals just under the requested width, while FLINT and msolve work at a fixed precision and return intervals many bits narrower. The isolate and refine tables therefore also time Poly Mojo at each one's achieved width (`Poly Mojo at FLINT width`, `Poly Mojo at msolve width`): 2^-b for the smallest power of two at or below that tool's widest interval on the case, never looser than requested. `matched / peer` divides that time by the peer's. The width tables give each tool's widest and median interval as log2 of width over 2^-bits (0 = exactly the requested width, -10 = 2^10 times narrower), and the longest endpoint in bits.

## Refinement of every root to width 2^-53 (Poly Mojo `roots.refine` of the isolating cover)

| Case | Poly Mojo | Poly Mojo (Sturm count) | FLINT | FLINT (Sturm count) | SymPy | msolve | Poly Mojo at FLINT width | Poly Mojo at msolve width | Poly Mojo / best | FLINT: matched / peer | msolve: matched / peer |
|---|---|---|---|---|---|---|---|---|---|---|---|
| chebyshev-100.refine53 | 0.0151 s | n/a | 0.0764 s | n/a | 7.54 s | 0.0505 s | 0.0226 s | 0.0192 s | 0.3 | 0.3 | 0.38 |
| chebyshev-200.refine53 | 0.0575 s | n/a | 0.862 s | n/a | 74.3 s | 0.213 s | 0.138 s | 0.0811 s | 0.27 | 0.16 | 0.38 |
| chebyshev-400.refine53 | 0.311 s | n/a | 11.4 s | n/a | timeout | 2.53 s | 0.758 s | 0.579 s | 0.12 | 0.067 | 0.23 |
| chebyshev-700.refine53 | 1.12 s | n/a | timeout | n/a | · | 30 s | n/a | 3.05 s | 0.037 |  | 0.1 |
| chebyshev-1000.refine53 | 3.48 s | n/a | · | n/a | · | timeout | n/a | n/a |  |  |  |
| mignotte-10.refine53 | 0.000527 s | n/a | 0.000783 s | n/a | 1.07 s | 0.00138 s | 0.000889 s | 0.00072 s | 0.67 | 1.1 | 0.52 |
| mignotte-20.refine53 | 0.00076 s | n/a | 0.00524 s | n/a | 0.00895 s | 0.00472 s | 0.00211 s | 0.00159 s | 0.16 | 0.4 | 0.34 |
| mignotte-50.refine53 | 0.00198 s | n/a | 0.149 s | n/a | 0.0572 s | 0.0356 s | 0.00933 s | 0.00519 s | 0.056 | 0.063 | 0.15 |
| mignotte-100.refine53 | 0.00355 s | n/a | 2.23 s | n/a | 0.163 s | 0.453 s | 0.0371 s | 0.00995 s | 0.022 | 0.017 | 0.022 |
| mignotte-200.refine53 | 0.0118 s | n/a | 47.6 s | n/a | 1.08 s | 14.1 s | 0.245 s | 0.0292 s | 0.011 | 0.0051 | 0.0021 |
| mignotte-400.refine53 | 0.0683 s | n/a | timeout | n/a | 13.9 s | timeout | n/a | n/a | 0.0049 |  |  |
| random1-100.refine53 | 0.000718 s | n/a | 0.0341 s | n/a | 0.296 s | 0.00564 s | 0.000997 s | 0.000941 s | 0.13 | 0.029 | 0.17 |
| random1-200.refine53 | 0.00082 s | n/a | 0.202 s | n/a | 0.776 s | 0.0265 s | 0.00157 s | 0.00117 s | 0.031 | 0.0078 | 0.044 |
| random1-500.refine53 | 0.00364 s | n/a | 1.74 s | n/a | 39.2 s | 0.141 s | 0.00906 s | 0.00923 s | 0.026 | 0.0052 | 0.066 |
| random1-1000.refine53 | 0.00781 s | n/a | 10.5 s | n/a | 50.4 s | 0.671 s | 0.0209 s | 0.0211 s | 0.012 | 0.002 | 0.031 |
| random1-2000.refine53 | 0.0434 s | n/a | 62.1 s | n/a | timeout | 12.4 s | 0.1 s | 0.114 s | 0.0035 | 0.0016 | 0.0092 |
| random1024-100.refine53 | 0.00412 s | n/a | 0.0374 s | n/a | 0.527 s | 0.128 s | 0.00485 s | 0.00476 s | 0.11 | 0.13 | 0.037 |
| random1024-200.refine53 | 0.00506 s | n/a | 0.201 s | n/a | 1.61 s | 0.203 s | 0.0077 s | 0.00653 s | 0.025 | 0.038 | 0.032 |
| random1024-500.refine53 | 0.0154 s | n/a | 1.73 s | n/a | 16.4 s | 0.603 s | 0.0265 s | 0.0241 s | 0.026 | 0.015 | 0.04 |
| random1024-1000.refine53 | 0.0308 s | n/a | 11.3 s | n/a | 85.9 s | 2.33 s | 0.0584 s | 0.0582 s | 0.013 | 0.0052 | 0.025 |
| random1024-2000.refine53 | 0.172 s | n/a | 105 s | n/a | timeout | 38.7 s | 0.535 s | 0.514 s | 0.0044 | 0.0051 | 0.013 |
| random64-100.refine53 | 2.72e-05 s | n/a | 0.0476 s | n/a | 0.00625 s | 0.0242 s | n/a | n/a | 0.0044 |  |  |
| random64-200.refine53 | 0.000965 s | n/a | 0.25 s | n/a | 0.782 s | 0.0314 s | 0.00208 s | 0.00151 s | 0.031 | 0.0083 | 0.048 |
| random64-500.refine53 | 0.00862 s | n/a | 1.9 s | n/a | 16.4 s | 0.176 s | 0.0192 s | 0.0194 s | 0.049 | 0.01 | 0.11 |
| random64-1000.refine53 | 0.0239 s | n/a | 12.5 s | n/a | timeout | 1.46 s | 0.0509 s | 0.0478 s | 0.016 | 0.0041 | 0.033 |
| random64-2000.refine53 | 0.0442 s | n/a | 47.2 s | n/a | · | 13.6 s | 0.105 s | 0.142 s | 0.0032 | 0.0022 | 0.01 |
| repeated-4.refine53 | 0.000233 s | n/a | error | n/a | 0.00336 s | 0.000112 s | n/a | 0.000319 s | 2.1 |  | 2.9 |
| repeated-8.refine53 | 0.000675 s | n/a | error | n/a | 0.0125 s | 0.0178 s | n/a | 0.000944 s | 0.054 |  | 0.053 |
| repeated-16.refine53 | 0.00148 s | n/a | error | n/a | 0.0531 s | 0.0229 s | n/a | 0.00203 s | 0.065 |  | 0.089 |
| repeated-32.refine53 | 0.00508 s | n/a | error | n/a | 0.274 s | 0.0432 s | n/a | 0.00701 s | 0.12 |  | 0.16 |
| repeated-64.refine53 | 0.0266 s | n/a | error | n/a | 1.93 s | 0.102 s | n/a | 0.0307 s | 0.26 |  | 0.3 |
| wilkinson-20.refine53 | 0.000117 s | n/a | 0.00301 s | n/a | 0.000873 s | 0.00267 s | 0.000126 s | n/a | 0.13 | 0.042 |  |
| wilkinson-40.refine53 | 0.000489 s | n/a | 0.0266 s | n/a | 0.00461 s | 0.0146 s | 0.000521 s | n/a | 0.11 | 0.02 |  |
| wilkinson-80.refine53 | 0.00322 s | n/a | 0.529 s | n/a | 0.0263 s | 0.0409 s | 0.00324 s | n/a | 0.12 | 0.0061 |  |
| wilkinson-160.refine53 | 0.0209 s | n/a | 7.36 s | n/a | 0.217 s | 0.119 s | 0.0204 s | n/a | 0.17 | 0.0028 |  |
| wilkinson-320.refine53 | 0.149 s | n/a | timeout | n/a | 2.05 s | 0.409 s | n/a | n/a | 0.36 |  |  |
| wilkinson-640.refine53 | 1.37 s | n/a | · | n/a | 20.8 s | 2.86 s | n/a | n/a | 0.48 |  |  |

Achieved widths (refine53): widest / median as log2 of width over 2^-53, longest endpoint in bits.

| Case | Poly Mojo | Poly Mojo (Sturm count) | FLINT | SymPy | msolve | Poly Mojo at FLINT width | Poly Mojo at msolve width |
|---|---|---|---|---|---|---|---|
| chebyshev-100.refine53 | -1.1 / -2.0, 125 b |  | -80.8 / -108.2, 562 b | -0.1 / -1.1, 61 b | -53.0 / -80.0, 281 b | -81.0 / -82.0, 287 b | -54.1 / -55.0, 231 b |
| chebyshev-200.refine53 | -1.0 / -2.0, 131 b |  | -210.2 / -269.4, 1075 b | -0.0 / -1.1, 67 b | -55.0 / -81.0, 297 b | -211.0 / -213.0, 551 b | -56.0 / -57.0, 241 b |
| chebyshev-400.refine53 | -1.0 / -2.0, 133 b |  | -465.4 / -587.5, 2100 b |  | -156.0 / -211.0, 567 b | -466.0 / -467.0, 1065 b | -157.0 / -158.0, 445 b |
| chebyshev-700.refine53 | -1.0 / -2.0, 133 b |  |  |  | -336.0 / -342.0, 1059 b |  | -337.0 / -338.0, 805 b |
| chebyshev-1000.refine53 | -1.0 / -2.0, 135 b |  |  |  |  |  |  |
| mignotte-10.refine53 | -0.1 / -0.9, 283 b |  | -170.3 / -170.3, 561 b | -1.6 / -1.8, 62 b | -68.0 / -76.0, 269 b | -172.1 / -172.1, 627 b | -69.1 / -69.1, 421 b |
| mignotte-20.refine53 | -1.4 / -1.6, 287 b |  | -392.1 / -392.1, 1070 b | -2.9 / -3.0, 126 b | -53.0 / -57.0, 274 b | -394.5 / -394.5, 1091 b | -54.5 / -54.5, 411 b |
| mignotte-50.refine53 | -1.1 / -1.4, 581 b |  | -804.7 / -805.0, 2096 b | -0.0 / -1.6, 324 b | -79.0 / -79.0, 339 b | -806.1 / -806.1, 2259 b | -80.0 / -80.1, 783 b |
| mignotte-100.refine53 | -1.1 / -1.4, 1157 b |  | -1661.3 / -1661.5, 4142 b | -1.1 / -3.7, 658 b | -81.0 / -83.0, 673 b | -1663.1 / -1663.3, 4771 b | -82.1 / -82.2, 1367 b |
| mignotte-200.refine53 | -1.1 / -1.2, 2327 b |  | -3376.7 / -3376.8, 8236 b | -0.5 / -1.6, 1322 b | -83.0 / -87.0, 1333 b | -3378.1 / -3378.2, 9439 b | -83.2 / -83.4, 2539 b |
| mignotte-400.refine53 | -1.0 / -1.4, 4661 b |  |  | -1.8 / -2.2, 2650 b |  |  |  |
| random1-100.refine53 | -1.9 / -2.0, 126 b |  | -65.2 / -65.8, 301 b | -0.7 / -0.7, 60 b | -66.0 / -83.0, 278 b | -67.9 / -68.0, 258 b | -67.9 / -68.0, 258 b |
| random1-200.refine53 | -1.4 / -1.4, 120 b |  | -192.7 / -192.7, 554 b | -0.1 / -0.1, 56 b | -83.0 / -83.0, 278 b | -194.2 / -194.2, 552 b | -84.2 / -84.2, 332 b |
| random1-500.refine53 | -1.1 / -1.7, 151 b |  | -447.5 / -447.7, 1066 b | -0.9 / -3.3, 60 b | -334.0 / -347.0, 809 b | -449.1 / -449.3, 1131 b | -335.1 / -335.3, 903 b |
| random1-1000.refine53 | -1.1 / -1.1, 127 b |  | -956.9 / -956.9, 2083 b | -3.1 / -3.1, 58 b | -723.0 / -723.0, 1686 b | -958.1 / -958.1, 2524 b | -724.1 / -724.1, 2056 b |
| random1-2000.refine53 | -1.1 / -1.4, 125 b |  | -954.1 / -955.9, 2081 b |  | -1013.0 / -1237.0, 3109 b | -956.1 / -956.4, 2528 b | -1014.1 / -1014.4, 2644 b |
| random1024-100.refine53 | -1.4 / -2.0, 114 b |  | -65.1 / -65.2, 300 b | -0.9 / -2.2, 58 b | -74.0 / -85.0, 288 b | -67.4 / -68.0, 246 b | -75.4 / -76.0, 262 b |
| random1024-200.refine53 | -1.1 / -1.2, 126 b |  | -191.4 / -191.6, 554 b | -1.2 / -1.8, 63 b | -78.0 / -82.0, 287 b | -193.1 / -193.2, 510 b | -79.1 / -79.2, 282 b |
| random1024-500.refine53 | -1.1 / -1.2, 125 b |  | -446.8 / -446.9, 1064 b | -2.2 / -2.3, 64 b | -243.0 / -343.0, 806 b | -448.1 / -448.2, 1072 b | -244.1 / -244.2, 664 b |
| random1024-1000.refine53 | -1.4 / -1.5, 307 b |  | -957.2 / -957.6, 2087 b | -0.3 / -0.5, 57 b | -723.0 / -729.0, 1581 b | -959.4 / -959.6, 2351 b | -724.4 / -724.6, 1881 b |
| random1024-2000.refine53 | -1.0 / -1.4, 158 b |  | -1976.5 / -1979.8, 4130 b |  | -1488.0 / -1496.0, 3350 b | -1978.0 / -1978.4, 4325 b | -1489.0 / -1489.4, 3347 b |
| random64-100.refine53 |  |  |  |  |  |  |  |
| random64-200.refine53 | -1.3 / -1.3, 122 b |  | -192.0 / -192.0, 549 b | -0.5 / -0.5, 57 b | -82.0 / -82.0, 282 b | -193.3 / -193.3, 557 b | -83.3 / -83.3, 337 b |
| random64-500.refine53 | -1.5 / -1.6, 124 b |  | -442.4 / -446.2, 1064 b | -0.1 / -0.6, 66 b | -335.0 / -340.0, 1046 b | -444.0 / -444.3, 1112 b | -336.0 / -336.3, 896 b |
| random64-1000.refine53 | -1.1 / -1.4, 304 b |  | -955.4 / -958.2, 2098 b |  | -716.0 / -723.0, 1584 b | -957.0 / -957.1, 2216 b | -717.0 / -717.1, 1736 b |
| random64-2000.refine53 | -1.1 / -1.8, 138 b |  | -952.8 / -954.4, 2081 b |  | -1489.0 / -1494.0, 3108 b | -954.1 / -954.4, 2134 b | -1490.4 / -1490.4, 4125 b |
| repeated-4.refine53 | -1.0 / -2.0, 123 b |  |  | -0.3 / -0.4, 56 b | -53.0 / -78.0, 266 b |  | -54.0 / -55.0, 229 b |
| repeated-8.refine53 | -1.1 / -1.5, 159 b |  |  | -0.3 / -0.4, 56 b | -53.0 / -74.0, 278 b |  | -54.1 / -54.5, 265 b |
| repeated-16.refine53 | -1.0 / -1.7, 127 b |  |  | -0.3 / -0.4, 57 b | -56.0 / -78.0, 278 b |  | -57.0 / -57.7, 239 b |
| repeated-32.refine53 | -1.0 / -1.8, 128 b |  |  | -0.1 / -1.3, 59 b | -61.0 / -79.0, 278 b |  | -62.0 / -62.8, 298 b |
| repeated-64.refine53 | -1.0 / -1.3, 131 b |  |  | -0.2 / -1.3, 62 b | -55.0 / -80.0, 283 b |  | -56.0 / -56.3, 284 b |
| wilkinson-20.refine53 | points, 6 b |  | -20.0 / -25.0, 290 b | points, 6 b | points, 6 b | points, 6 b |  |
| wilkinson-40.refine53 | points, 7 b |  | -96.3 / -107.3, 540 b | points, 7 b | points, 7 b | points, 7 b |  |
| wilkinson-80.refine53 | points, 8 b |  | -249.6 / -270.5, 1048 b | points, 8 b | points, 8 b | points, 8 b |  |
| wilkinson-160.refine53 | points, 9 b |  | -44.5 / -87.3, 1038 b | points, 9 b | points, 9 b | points, 9 b |  |
| wilkinson-320.refine53 | points, 10 b |  |  | points, 10 b | points, 10 b |  |  |
| wilkinson-640.refine53 | points, 11 b |  |  | points, 11 b | points, 11 b |  |  |

## Refinement of every root to width 2^-1000

| Case | Poly Mojo | Poly Mojo (Sturm count) | FLINT | FLINT (Sturm count) | SymPy | msolve | Poly Mojo at FLINT width | Poly Mojo at msolve width | Poly Mojo / best | FLINT: matched / peer | msolve: matched / peer |
|---|---|---|---|---|---|---|---|---|---|---|---|
| chebyshev-100.refine1000 | 0.0671 s | n/a | 0.118 s | n/a | timeout | 1.03 s | 0.11 s | 0.117 s | 0.57 | 0.93 | 0.11 |
| chebyshev-200.refine1000 | 0.278 s | n/a | 1.14 s | n/a | · | 8.77 s | 0.432 s | 0.464 s | 0.24 | 0.38 | 0.053 |
| chebyshev-400.refine1000 | 1.09 s | n/a | 11.5 s | n/a | · | 68.2 s | 1.56 s | 1.79 s | 0.094 | 0.14 | 0.026 |
| chebyshev-700.refine1000 | 3.82 s | n/a | timeout | n/a | · | timeout | n/a | n/a |  |  |  |
| chebyshev-1000.refine1000 | 9.2 s | n/a | · | n/a | · | · | n/a | n/a |  |  |  |
| mignotte-10.refine1000 | 0.00145 s | n/a | 0.00313 s | n/a | 1.2 s | 0.000415 s | 0.00208 s | 0.00253 s | 3.5 | 0.66 | 6.1 |
| mignotte-20.refine1000 | 0.00258 s | n/a | 0.011 s | n/a | timeout | 0.0165 s | 0.0033 s | 0.0033 s | 0.23 | 0.3 | 0.2 |
| mignotte-50.refine1000 | 0.00884 s | n/a | 0.143 s | n/a | · | 0.0494 s | 0.00994 s | 0.0107 s | 0.18 | 0.07 | 0.22 |
| mignotte-100.refine1000 | 0.0341 s | n/a | 2.24 s | n/a | · | 0.475 s | 0.0389 s | 0.0378 s | 0.072 | 0.017 | 0.08 |
| mignotte-200.refine1000 | 0.227 s | n/a | 46.9 s | n/a | · | 13.7 s | 0.261 s | 0.243 s | 0.017 | 0.0056 | 0.018 |
| mignotte-400.refine1000 | 1.79 s | n/a | timeout | n/a | · | timeout | n/a | n/a |  |  |  |
| random1-100.refine1000 | 0.0029 s | n/a | 0.213 s | n/a | 36.2 s | 0.0583 s | 0.00492 s | 0.00605 s | 0.05 | 0.023 | 0.1 |
| random1-200.refine1000 | 0.00307 s | n/a | 0.819 s | n/a | 116 s | 0.103 s | 0.00465 s | 0.00444 s | 0.03 | 0.0057 | 0.043 |
| random1-500.refine1000 | 0.0118 s | n/a | 4.94 s | n/a | timeout | 1.1 s | 0.0201 s | 0.0215 s | 0.011 | 0.0041 | 0.019 |
| random1-1000.refine1000 | 0.0202 s | n/a | 19.8 s | n/a | · | 3.81 s | 0.0324 s | 0.0374 s | 0.0053 | 0.0016 | 0.0098 |
| random1-2000.refine1000 | 0.108 s | n/a | 96.5 s | n/a | · | 24.6 s | 0.145 s | 0.175 s | 0.0044 | 0.0015 | 0.0071 |
| random1024-100.refine1000 | 0.00993 s | n/a | 0.257 s | n/a | 62.6 s | 0.199 s | 0.0126 s | 0.0137 s | 0.05 | 0.049 | 0.069 |
| random1024-200.refine1000 | 0.0109 s | n/a | 1.02 s | n/a | timeout | 0.372 s | 0.0157 s | 0.0167 s | 0.029 | 0.015 | 0.045 |
| random1024-500.refine1000 | 0.0303 s | n/a | 5.77 s | n/a | · | 1.49 s | 0.0426 s | 0.0462 s | 0.02 | 0.0074 | 0.031 |
| random1024-1000.refine1000 | 0.0626 s | n/a | 23.7 s | n/a | · | 7.75 s | 0.0884 s | 0.0975 s | 0.0081 | 0.0037 | 0.013 |
| random1024-2000.refine1000 | 0.348 s | n/a | 107 s | n/a | · | 101 s | 0.507 s | 0.642 s | 0.0034 | 0.0047 | 0.0064 |
| random64-100.refine1000 | 2.73e-05 s | n/a | 0.21 s | n/a | 0.0061 s | 0.0185 s | n/a | n/a | 0.0045 |  |  |
| random64-200.refine1000 | 0.00321 s | n/a | 0.88 s | n/a | 118 s | 0.115 s | 0.00533 s | 0.00535 s | 0.028 | 0.0061 | 0.046 |
| random64-500.refine1000 | 0.0235 s | n/a | 4.68 s | n/a | timeout | 1.57 s | 0.0355 s | 0.0336 s | 0.015 | 0.0076 | 0.021 |
| random64-1000.refine1000 | 0.0569 s | n/a | 20.1 s | n/a | · | 10.7 s | 0.0801 s | 0.1 s | 0.0053 | 0.004 | 0.0093 |
| random64-2000.refine1000 | 0.113 s | n/a | 84.4 s | n/a | · | 39.9 s | 0.15 s | 0.255 s | 0.0028 | 0.0018 | 0.0064 |
| repeated-4.refine1000 | 0.00123 s | n/a | error | n/a | 0.0794 s | 0.0202 s | n/a | 0.00291 s | 0.061 |  | 0.14 |
| repeated-8.refine1000 | 0.00333 s | n/a | error | n/a | 0.372 s | 0.0259 s | n/a | 0.00548 s | 0.13 |  | 0.21 |
| repeated-16.refine1000 | 0.00676 s | n/a | error | n/a | 1.84 s | 0.0338 s | n/a | 0.0108 s | 0.2 |  | 0.32 |
| repeated-32.refine1000 | 0.0201 s | n/a | error | n/a | 18.4 s | 0.11 s | n/a | 0.0354 s | 0.18 |  | 0.32 |
| repeated-64.refine1000 | 0.0753 s | n/a | error | n/a | timeout | 0.829 s | n/a | 0.131 s | 0.091 |  | 0.16 |
| wilkinson-20.refine1000 | 0.000129 s | n/a | 0.0119 s | n/a | 0.000873 s | 0.0113 s | 0.000128 s | n/a | 0.15 | 0.011 |  |
| wilkinson-40.refine1000 | 0.000536 s | n/a | 0.0545 s | n/a | 0.0055 s | 0.0236 s | 0.000563 s | n/a | 0.098 | 0.01 |  |
| wilkinson-80.refine1000 | 0.00343 s | n/a | 0.739 s | n/a | 0.0289 s | 0.0399 s | 0.00334 s | n/a | 0.12 | 0.0045 |  |
| wilkinson-160.refine1000 | 0.0237 s | n/a | 8.05 s | n/a | 0.217 s | 0.123 s | 0.0201 s | n/a | 0.19 | 0.0025 |  |
| wilkinson-320.refine1000 | 0.156 s | n/a | timeout | n/a | 2.05 s | 0.426 s | n/a | n/a | 0.37 |  |  |
| wilkinson-640.refine1000 | 1.29 s | n/a | · | n/a | 22 s | 3.32 s | n/a | n/a | 0.39 |  |  |

Achieved widths (refine1000): widest / median as log2 of width over 2^-1000, longest endpoint in bits.

| Case | Poly Mojo | Poly Mojo (Sturm count) | FLINT | SymPy | msolve | Poly Mojo at FLINT width | Poly Mojo at msolve width |
|---|---|---|---|---|---|---|---|
| chebyshev-100.refine1000 | -1.1 / -2.0, 2019 b |  | -923.8 / -950.6, 4150 b |  | -1051.0 / -1055.0, 4619 b | -925.1 / -926.0, 3867 b | -1052.1 / -1053.0, 4121 b |
| chebyshev-200.refine1000 | -1.0 / -2.0, 2025 b |  | -796.7 / -855.1, 4149 b |  | -1052.0 / -1054.0, 4620 b | -798.0 / -799.0, 3619 b | -1053.0 / -1054.0, 4129 b |
| chebyshev-400.refine1000 | -1.0 / -2.0, 2027 b |  | -542.4 / -664.5, 4148 b |  | -1053.0 / -1056.0, 4245 b | -544.0 / -545.0, 3113 b | -1054.0 / -1055.0, 4133 b |
| chebyshev-700.refine1000 | -1.0 / -2.0, 2027 b |  |  |  |  |  |  |
| chebyshev-1000.refine1000 | -1.0 / -2.0, 2029 b |  |  |  |  |  |  |
| mignotte-10.refine1000 | -1.1 / -1.1, 2179 b |  | -1014.8 / -1015.3, 4147 b | -0.2 / -1.6, 1010 b | -1049.0 / -1089.0, 4229 b | -1016.1 / -1016.1, 4209 b | -1050.1 / -1050.1, 4277 b |
| mignotte-20.refine1000 | -1.5 / -1.5, 2199 b |  | -982.3 / -982.8, 4144 b |  | -1056.0 / -1064.0, 4495 b | -984.5 / -984.5, 4165 b | -1057.5 / -1057.5, 4311 b |
| mignotte-50.refine1000 | -1.1 / -1.1, 2543 b |  | -881.4 / -881.9, 4142 b |  | -1084.0 / -1180.0, 4443 b | -883.1 / -883.1, 4307 b | -1085.1 / -1085.1, 4711 b |
| mignotte-100.refine1000 | -1.1 / -1.3, 3341 b |  | -714.3 / -714.5, 4142 b |  | -1054.0 / -1056.0, 4769 b | -716.1 / -716.3, 4771 b | -1055.1 / -1055.3, 5449 b |
| mignotte-200.refine1000 | -1.1 / -1.2, 4579 b |  | -2429.7 / -2429.8, 8236 b |  | -1056.0 / -1060.0, 4417 b | -2431.1 / -2431.2, 9439 b | -1056.2 / -1057.1, 6691 b |
| mignotte-400.refine1000 | -1.0 / -1.4, 7107 b |  |  |  |  |  |  |
| random1-100.refine1000 | -1.9 / -2.0, 2020 b |  | -1036.5 / -1037.3, 4138 b | -0.1 / -0.3, 1006 b | -1056.0 / -1057.0, 4144 b | -1038.9 / -1039.0, 4094 b | -1056.0 / -1057.9, 4132 b |
| random1-200.refine1000 | -1.2 / -1.2, 2060 b |  | -1036.6 / -1036.6, 4136 b | -0.1 / -0.1, 1003 b | -1056.0 / -1056.0, 4118 b | -1038.2 / -1038.2, 4134 b | -1057.2 / -1057.2, 4172 b |
| random1-500.refine1000 | -1.1 / -1.3, 2129 b |  | -1034.7 / -1035.5, 4131 b |  | -1051.0 / -1064.0, 4137 b | -1036.1 / -1036.3, 4199 b | -1052.1 / -1052.3, 4231 b |
| random1-1000.refine1000 | -1.1 / -1.1, 2504 b |  | -1034.0 / -1034.0, 4131 b |  | -1122.0 / -1122.0, 6161 b | -1035.1 / -1035.1, 4572 b | -1123.1 / -1123.1, 4748 b |
| random1-2000.refine1000 | -1.1 / -1.4, 2512 b |  | -1031.1 / -1033.2, 4131 b |  | -1066.0 / -1090.0, 6170 b | -1033.1 / -1033.4, 4576 b | -1067.1 / -1067.4, 4644 b |
| random1024-100.refine1000 | -1.1 / -2.0, 2259 b |  | -1036.8 / -1037.8, 4147 b | -0.6 / -1.2, 1007 b | -1047.0 / -1058.0, 4128 b | -1038.1 / -1039.0, 4333 b | -1048.1 / -1049.0, 4353 b |
| random1024-200.refine1000 | -1.1 / -1.2, 2020 b |  | -1035.7 / -1036.7, 4134 b |  | -1051.0 / -1055.0, 4368 b | -1037.1 / -1037.2, 4092 b | -1052.1 / -1052.2, 4122 b |
| random1024-500.refine1000 | -1.1 / -1.2, 2072 b |  | -1033.8 / -1034.5, 4134 b |  | -1060.0 / -1062.0, 4177 b | -1035.1 / -1035.2, 4140 b | -1061.1 / -1061.2, 4192 b |
| random1024-1000.refine1000 | -1.4 / -1.6, 2329 b |  | -1032.8 / -1032.9, 4133 b |  | -1062.0 / -1062.0, 6161 b | -1034.4 / -1034.6, 4395 b | -1063.4 / -1063.6, 4453 b |
| random1024-2000.refine1000 | -1.0 / -1.4, 2265 b |  | -1029.5 / -1032.8, 4130 b |  | -1089.0 / -2082.0, 6422 b | -1031.0 / -1031.4, 4325 b | -1090.0 / -1090.4, 4443 b |
| random64-100.refine1000 |  |  |  |  |  |  |  |
| random64-200.refine1000 | -1.3 / -1.3, 2067 b |  | -1034.7 / -1034.7, 4136 b | -0.7 / -0.7, 1004 b | -1055.0 / -1055.0, 4122 b | -1036.3 / -1036.3, 4137 b | -1056.3 / -1056.3, 4177 b |
| random64-500.refine1000 | -1.0 / -1.3, 2120 b |  | -1028.1 / -1031.9, 4135 b |  | -1052.0 / -1059.0, 6166 b | -1030.0 / -1030.3, 4178 b | -1053.0 / -1053.3, 4224 b |
| random64-1000.refine1000 | -1.0 / -1.1, 2198 b |  | -1032.3 / -1034.9, 4140 b |  | -1063.0 / -1071.0, 6161 b | -1034.0 / -1034.1, 5077 b | -1063.1 / -1064.1, 5137 b |
| random64-2000.refine1000 | -1.1 / -1.4, 2122 b |  | -1029.8 / -1031.3, 4131 b |  | -2078.0 / -2083.0, 6180 b | -1031.4 / -1031.4, 5101 b | -2079.4 / -2079.4, 7197 b |
| repeated-4.refine1000 | -1.0 / -2.0, 2017 b |  |  | -0.3 / -0.6, 1004 b | -1051.0 / -1053.0, 4118 b |  | -1052.0 / -1053.0, 4119 b |
| repeated-8.refine1000 | -1.1 / -1.5, 2053 b |  |  | -0.3 / -0.5, 1008 b | -1047.0 / -1052.0, 4135 b |  | -1048.0 / -1048.2, 4147 b |
| repeated-16.refine1000 | -1.0 / -1.7, 2271 b |  |  | -0.1 / -1.1, 1006 b | -1048.0 / -1052.0, 4226 b |  | -1049.0 / -1049.7, 4367 b |
| repeated-32.refine1000 | -1.0 / -1.8, 2070 b |  |  | -0.1 / -0.8, 1006 b | -1049.0 / -1052.0, 4619 b |  | -1050.0 / -1050.8, 4168 b |
| repeated-64.refine1000 | -1.0 / -1.3, 2068 b |  |  |  | -1050.0 / -1055.0, 4374 b |  | -1051.0 / -1051.3, 4168 b |
| wilkinson-20.refine1000 | points, 6 b |  | -991.6 / -996.5, 4128 b | points, 6 b | points, 6 b | points, 6 b |  |
| wilkinson-40.refine1000 | points, 7 b |  | -939.3 / -949.9, 4118 b | points, 7 b | points, 7 b | points, 7 b |  |
| wilkinson-80.refine1000 | points, 8 b |  | -836.0 / -856.9, 4116 b | points, 8 b | points, 8 b | points, 8 b |  |
| wilkinson-160.refine1000 | points, 9 b |  | -631.1 / -674.4, 4106 b | points, 9 b | points, 9 b | points, 9 b |  |
| wilkinson-320.refine1000 | points, 10 b |  |  | points, 10 b | points, 10 b |  |  |
| wilkinson-640.refine1000 | points, 11 b |  |  | points, 11 b | points, 11 b |  |  |

