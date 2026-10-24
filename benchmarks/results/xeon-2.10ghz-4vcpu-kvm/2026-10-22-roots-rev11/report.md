# Real roots

Host: Intel(R) Xeon(R) Processor @ 2.10GHz, 4 vCPUs, 16 GiB RAM, Ubuntu 24.04.5 LTS. Runs are pinned to one CPU, capped at 12G and 120 s.

Each answer is checked against an Arb reference: the number of distinct real roots, and every interval holding exactly one reference root (`wrong` otherwise). The `Poly Mojo / best` column divides Poly Mojo's time by the fastest other library on the case, each at its own achieved width.

Widths differ by tool: Poly Mojo and SymPy return intervals just under the requested width, while FLINT and msolve work at a fixed precision and return intervals many bits narrower. The isolate and refine tables therefore also time Poly Mojo at each one's achieved width (`Poly Mojo at FLINT width`, `Poly Mojo at msolve width`): 2^-b for the smallest power of two at or below that tool's widest interval on the case, never looser than requested. `matched / peer` divides that time by the peer's. The width tables give each tool's widest and median interval as log2 of width over 2^-bits (0 = exactly the requested width, -10 = 2^10 times narrower), and the longest endpoint in bits.

## Isolation (Poly Mojo `real_root_cover`, width 1/2)

| Case | Poly Mojo | Poly Mojo (Sturm count) | FLINT | FLINT (Sturm count) | SymPy | msolve | Poly Mojo at FLINT width | Poly Mojo at msolve width | Poly Mojo / best | FLINT: matched / peer | msolve: matched / peer |
|---|---|---|---|---|---|---|---|---|---|---|---|
| chebyshev-100.isolate | 0.00258 s | n/a | 0.0765 s | n/a | 0.383 s | 0.0121 s | 0.0263 s | 0.0177 s | 0.21 | 0.34 | 1.5 |
| chebyshev-200.isolate | 0.0337 s | n/a | 0.845 s | n/a | 3.32 s | 0.184 s | 0.175 s | 0.123 s | 0.18 | 0.21 | 0.67 |
| chebyshev-400.isolate | 0.19 s | n/a | 11.6 s | n/a | 40.1 s | 2.64 s | 1.1 s | 0.716 s | 0.072 | 0.095 | 0.27 |
| chebyshev-700.isolate | 1.82 s | n/a | timeout | n/a | timeout | 31.9 s | n/a | 4.76 s | 0.057 |  | 0.15 |
| chebyshev-1000.isolate | 5.37 s | n/a | · | n/a | · | timeout | n/a | n/a |  |  |  |
| mignotte-10.isolate | 5.99e-05 s | n/a | 0.000779 s | n/a | 0.000557 s | 0 s | 0.00162 s | 0.000167 s | 0.11 | 2.1 |  |
| mignotte-20.isolate | 7.94e-05 s | n/a | 0.00534 s | n/a | 0.000919 s | 0 s | 0.00226 s | 0.000282 s | 0.086 | 0.42 |  |
| mignotte-50.isolate | 0.000231 s | n/a | 0.138 s | n/a | 0.00353 s | 0.00881 s | 0.148 s | 0.00145 s | 0.065 | 1.1 | 0.16 |
| mignotte-100.isolate | 0.000684 s | n/a | 2.22 s | n/a | 0.0106 s | 0.462 s | 1.38 s | 0.00629 s | 0.065 | 0.62 | 0.014 |
| mignotte-200.isolate | 0.00309 s | n/a | 51.4 s | n/a | 0.0415 s | 14.8 s | 13.2 s | 0.0264 s | 0.075 | 0.26 | 0.0018 |
| mignotte-400.isolate | 0.0197 s | n/a | timeout | n/a | 0.179 s | timeout | n/a | n/a | 0.11 |  |  |
| random1-100.isolate | 0.0006 s | n/a | 0.027 s | n/a | 0.00791 s | 0.00302 s | 0.00144 s | 0.00136 s | 0.2 | 0.053 | 0.45 |
| random1-200.isolate | 0.00281 s | n/a | 0.177 s | n/a | 0.0305 s | 0.00758 s | 0.00368 s | 0.00398 s | 0.37 | 0.021 | 0.53 |
| random1-500.isolate | 0.0245 s | n/a | 1.93 s | n/a | 0.233 s | 0.112 s | 0.0349 s | 0.0345 s | 0.22 | 0.018 | 0.31 |
| random1-1000.isolate | 0.219 s | n/a | 8.3 s | n/a | 0.551 s | 0.633 s | 0.213 s | 0.235 s | 0.4 | 0.026 | 0.37 |
| random1-2000.isolate | 1.84 s | n/a | 62.4 s | n/a | 6.02 s | 11.9 s | 1.88 s | 2.03 s | 0.31 | 0.03 | 0.17 |
| random1024-100.isolate | 0.00209 s | n/a | 0.0309 s | n/a | 0.0191 s | 0.114 s | 0.00625 s | 0.00633 s | 0.11 | 0.2 | 0.056 |
| random1024-200.isolate | 0.00847 s | n/a | 0.155 s | n/a | 0.0671 s | 0.177 s | 0.0152 s | 0.0141 s | 0.13 | 0.098 | 0.08 |
| random1024-500.isolate | 0.0527 s | n/a | 1.66 s | n/a | 0.34 s | 0.692 s | 0.08 s | 0.0783 s | 0.16 | 0.048 | 0.11 |
| random1024-1000.isolate | 0.394 s | n/a | 6.95 s | n/a | 2.61 s | 2.56 s | 0.485 s | 0.475 s | 0.15 | 0.07 | 0.19 |
| random1024-2000.isolate | 5.44 s | n/a | 62 s | n/a | 37.8 s | 37.3 s | 5.63 s | 5.72 s | 0.15 | 0.091 | 0.15 |
| random64-100.isolate | 0.00171 s | n/a | 0.0396 s | n/a | 0.00616 s | 0.0148 s | n/a | n/a | 0.28 |  |  |
| random64-200.isolate | 0.00772 s | n/a | 0.189 s | n/a | 0.034 s | 0.0272 s | 0.00811 s | 0.00848 s | 0.28 | 0.043 | 0.31 |
| random64-500.isolate | 0.174 s | n/a | 2.15 s | n/a | 0.445 s | 0.143 s | 0.188 s | 0.22 s | 1.2 | 0.087 | 1.5 |
| random64-1000.isolate | 0.284 s | n/a | 8.75 s | n/a | 2.08 s | 1.66 s | 0.336 s | 0.331 s | 0.17 | 0.038 | 0.2 |
| random64-2000.isolate | 1.75 s | n/a | 48.3 s | n/a | 4.16 s | 13.6 s | 1.84 s | 1.82 s | 0.42 | 0.038 | 0.13 |
| repeated-4.isolate | 7.71e-05 s | n/a | error | n/a | wrong | 0.0039 s | n/a | 0.000166 s | 0.02 |  | 0.042 |
| repeated-8.isolate | 0.000262 s | n/a | error | n/a | wrong | 0 s | n/a | 0.000656 s |  |  |  |
| repeated-16.isolate | 0.00129 s | n/a | error | n/a | wrong | 0.0037 s | n/a | 0.00255 s | 0.35 |  | 0.69 |
| repeated-32.isolate | 0.0065 s | n/a | error | n/a | wrong | 0.022 s | n/a | 0.0117 s | 0.3 |  | 0.53 |
| repeated-64.isolate | 0.0455 s | n/a | error | n/a | wrong | 0.0838 s | n/a | 0.072 s | 0.54 |  | 0.86 |
| wilkinson-20.isolate | 0.000139 s | n/a | 0.00244 s | n/a | 0.000839 s | 0 s | 0.000295 s | n/a | 0.17 | 0.12 |  |
| wilkinson-40.isolate | 0.000535 s | n/a | 0.0225 s | n/a | 0.0041 s | 0 s | 0.00104 s | n/a | 0.13 | 0.046 |  |
| wilkinson-80.isolate | 0.00298 s | n/a | 0.565 s | n/a | 0.0271 s | 0.0144 s | 0.00657 s | n/a | 0.21 | 0.012 |  |
| wilkinson-160.isolate | 0.0272 s | n/a | 7.79 s | n/a | 0.192 s | 0.0819 s | 0.0434 s | n/a | 0.33 | 0.0056 |  |
| wilkinson-320.isolate | 0.196 s | n/a | timeout | n/a | 2.17 s | 0.462 s | n/a | n/a | 0.42 |  |  |
| wilkinson-640.isolate | 2.06 s | n/a | · | n/a | 23.5 s | 2.93 s | n/a | n/a | 0.7 |  |  |

Achieved widths (isolate): widest / median as log2 of width, longest endpoint in bits.

| Case | Poly Mojo | Poly Mojo (Sturm count) | FLINT | SymPy | msolve | Poly Mojo at FLINT width | Poly Mojo at msolve width |
|---|---|---|---|---|---|---|---|
| chebyshev-100.isolate | -7.0 / -8.0, 21 b |  | -133.8 / -161.2, 562 b | -4.5 / -6.3, 20 b | -57.0 / -68.0, 153 b | -134.0 / -135.0, 287 b | -58.1 / -59.0, 133 b |
| chebyshev-200.isolate | -8.0 / -9.0, 45 b |  | -263.2 / -322.4, 1075 b | -5.4 / -7.2, 24 b | -108.0 / -134.0, 297 b | -264.0 / -266.0, 551 b | -109.0 / -110.0, 241 b |
| chebyshev-400.isolate | -9.0 / -10.0, 49 b |  | -518.4 / -640.5, 2100 b | -6.1 / -8.2, 28 b | -209.0 / -264.0, 567 b | -519.0 / -520.0, 1065 b | -210.0 / -211.0, 445 b |
| chebyshev-700.isolate | -10.0 / -11.0, 51 b |  |  |  | -389.0 / -395.0, 1059 b |  | -390.0 / -391.0, 805 b |
| chebyshev-1000.isolate | -11.0 / -11.0, 53 b |  |  |  |  |  |  |
| mignotte-10.isolate | -2.0 / -2.1, 162 b |  | -223.3 / -223.3, 561 b | 0.0 / 0.0, 8 b | -9.0 / -9.0, 73 b | -225.1 / -225.3, 629 b | -10.1 / -10.4, 179 b |
| mignotte-20.isolate | -1.0 / -1.2, 149 b |  | -445.1 / -445.1, 1070 b | 0.0 / 0.0, 8 b | -14.0 / -16.0, 139 b | -447.5 / -447.5, 1091 b | -14.1 / -15.4, 185 b |
| mignotte-50.isolate | -1.0 / -2.4, 443 b |  | -857.7 / -858.0, 2096 b | 0.0 / 0.0, 8 b | -36.0 / -36.0, 339 b | -858.7 / -858.9, 4159 b | -37.1 / -37.2, 535 b |
| mignotte-100.isolate | -1.0 / -2.4, 1019 b |  | -1714.3 / -1714.5, 4142 b | 0.0 / 0.0, 8 b | -70.0 / -72.0, 673 b | -1715.3 / -1715.7, 8751 b | -71.4 / -71.6, 1199 b |
| mignotte-200.isolate | -1.0 / -2.4, 2189 b |  | -3429.7 / -3429.8, 8236 b | 0.0 / 0.0, 8 b | -136.0 / -140.0, 1333 b | -3430.1 / -3430.7, 18017 b | -136.2 / -137.1, 2515 b |
| mignotte-400.isolate | -1.0 / -2.4, 4523 b |  |  | 0.0 / 0.0, 8 b |  |  |  |
| random1-100.isolate | -2.0 / -3.0, 10 b |  | -54.6 / -55.0, 172 b | 0.0 / 0.0, 4 b | -58.0 / -71.0, 148 b | -56.9 / -57.0, 130 b | -59.9 / -60.0, 136 b |
| random1-200.isolate | -3.0 / -3.0, 8 b |  | -117.8 / -117.8, 300 b | 0.0 / 0.0, 3 b | -136.0 / -136.0, 278 b | -119.2 / -119.2, 296 b | -137.2 / -137.2, 332 b |
| random1-500.isolate | -1.4 / -3.0, 9 b |  | -500.5 / -500.7, 1066 b | 0.0 / 0.0, 6 b | -387.0 / -400.0, 809 b | -502.1 / -502.3, 1131 b | -388.1 / -388.3, 903 b |
| random1-1000.isolate | -1.4 / -1.4, 8 b |  | -499.9 / -499.9, 1062 b | 0.0 / 0.0, 3 b | -776.0 / -776.0, 1686 b | -501.1 / -501.1, 1021 b | -777.1 / -777.1, 2056 b |
| random1-2000.isolate | -3.0 / -3.0, 9 b |  | -1007.1 / -1008.9, 2081 b | 0.0 / 0.0, 5 b | -1066.0 / -1290.0, 3109 b | -1009.1 / -1009.4, 2528 b | -1067.1 / -1067.4, 2644 b |
| random1024-100.isolate | -3.0 / -4.0, 11 b |  | -54.1 / -54.4, 172 b | 0.0 / -1.6, 7 b | -63.0 / -74.0, 160 b | -56.4 / -57.0, 118 b | -64.4 / -65.0, 134 b |
| random1024-200.isolate | -3.0 / -3.0, 11 b |  | -116.6 / -116.6, 296 b | -1.0 / -1.0, 10 b | -131.0 / -135.0, 287 b | -118.1 / -118.2, 254 b | -132.1 / -132.2, 282 b |
| random1024-500.isolate | -1.4 / -3.0, 10 b |  | -499.8 / -499.9, 1064 b | 0.0 / 0.0, 5 b | -296.0 / -396.0, 806 b | -501.1 / -501.2, 1072 b | -297.1 / -297.2, 664 b |
| random1024-1000.isolate | -1.0 / -1.4, 189 b |  | -498.0 / -498.3, 1061 b | 0.0 / 0.0, 16 b | -776.0 / -782.0, 1581 b | -499.4 / -499.6, 1325 b | -777.4 / -777.6, 1881 b |
| random1024-2000.isolate | -1.0 / -5.0, 15 b |  | -1009.5 / -1009.9, 2084 b | 0.0 / -5.8, 14 b | -1541.0 / -1549.0, 3350 b | -1011.0 / -1011.4, 2285 b | -1542.0 / -1542.4, 3347 b |
| random64-100.isolate |  |  |  |  |  |  |  |
| random64-200.isolate | -1.0 / -1.0, 8 b |  | -116.6 / -116.6, 298 b | 0.0 / 0.0, 3 b | -135.0 / -135.0, 282 b | -118.3 / -118.3, 271 b | -136.3 / -136.3, 307 b |
| random64-500.isolate | -2.0 / -2.2, 16 b |  | -495.4 / -499.2, 1064 b | 0.0 / 0.0, 6 b | -388.0 / -393.0, 1046 b | -497.0 / -497.4, 1112 b | -389.0 / -389.4, 896 b |
| random64-1000.isolate | -1.0 / -4.0, 166 b |  | -498.8 / -499.9, 1063 b | 0.0 / -1.6, 14 b | -769.0 / -776.0, 1584 b | -500.0 / -500.1, 1196 b | -770.0 / -770.1, 1736 b |
| random64-2000.isolate | -2.0 / -2.0, 10 b |  | -1005.8 / -1007.4, 2081 b | 0.0 / 0.0, 4 b | -1542.0 / -1547.0, 3108 b | -1007.2 / -1007.4, 2134 b | -1543.2 / -1543.4, 3206 b |
| repeated-4.isolate | -1.0 / -3.0, 11 b |  |  |  | -8.0 / -9.0, 19 b |  | -9.0 / -9.0, 22 b |
| repeated-8.isolate | -2.0 / -2.0, 29 b |  |  |  | -13.0 / -15.0, 33 b |  | -13.7 / -14.3, 50 b |
| repeated-16.isolate | -2.0 / -4.0, 31 b |  |  |  | -24.0 / -27.0, 59 b |  | -25.0 / -25.7, 69 b |
| repeated-32.isolate | -2.0 / -6.0, 35 b |  |  |  | -41.0 / -52.0, 118 b |  | -42.0 / -42.8, 104 b |
| repeated-64.isolate | -2.0 / -7.0, 39 b |  |  |  | -75.0 / -101.0, 219 b |  | -76.0 / -76.3, 218 b |
| wilkinson-20.isolate | -2.0 / -2.0, 11 b |  | -8.9 / -14.1, 162 b | points, 6 b | points, 6 b | points, 6 b |  |
| wilkinson-40.isolate | -2.0 / -2.0, 12 b |  | -21.2 / -31.9, 280 b | points, 7 b | points, 7 b | points, 7 b |  |
| wilkinson-80.isolate | -2.0 / -2.0, 13 b |  | -46.3 / -66.8, 530 b | points, 8 b | points, 8 b | points, 8 b |  |
| wilkinson-160.isolate | -2.0 / -2.0, 14 b |  | -97.5 / -140.3, 1038 b | points, 9 b | points, 9 b | points, 9 b |  |
| wilkinson-320.isolate | -2.0 / -2.0, 15 b |  |  | points, 10 b | points, 10 b |  |  |
| wilkinson-640.isolate | -2.0 / -2.0, 16 b |  |  | points, 11 b | points, 11 b |  |  |

## Refinement of every root to width 2^-53 (Poly Mojo `roots.refine` of the isolating cover)

| Case | Poly Mojo | Poly Mojo (Sturm count) | FLINT | FLINT (Sturm count) | SymPy | msolve | Poly Mojo at FLINT width | Poly Mojo at msolve width | Poly Mojo / best | FLINT: matched / peer | msolve: matched / peer |
|---|---|---|---|---|---|---|---|---|---|---|---|
| chebyshev-100.refine53 | 0.0136 s | n/a | 0.0885 s | n/a | 8.09 s | 0.0301 s | 0.0236 s | 0.0188 s | 0.45 | 0.27 | 0.62 |
| chebyshev-200.refine53 | 0.0647 s | n/a | 0.887 s | n/a | 76.5 s | 0.206 s | 0.14 s | 0.0782 s | 0.31 | 0.16 | 0.38 |
| chebyshev-400.refine53 | 0.282 s | n/a | 12.2 s | n/a | timeout | 2.67 s | 0.749 s | 0.565 s | 0.11 | 0.061 | 0.21 |
| chebyshev-700.refine53 | 1.35 s | n/a | timeout | n/a | · | 30.9 s | n/a | 3.07 s | 0.044 |  | 0.099 |
| chebyshev-1000.refine53 | 3.18 s | n/a | · | n/a | · | timeout | n/a | n/a |  |  |  |
| mignotte-10.refine53 | 0.000576 s | n/a | 0.000783 s | n/a | 1.15 s | 0 s | 0.000956 s | 0.000742 s | 0.74 | 1.2 |  |
| mignotte-20.refine53 | 0.0008 s | n/a | 0.00549 s | n/a | 0.00951 s | 0 s | 0.00227 s | 0.00157 s | 0.15 | 0.41 |  |
| mignotte-50.refine53 | 0.00201 s | n/a | 0.195 s | n/a | 0.0621 s | 0.00909 s | 0.168 s | 0.00587 s | 0.22 | 0.86 | 0.65 |
| mignotte-100.refine53 | 0.00376 s | n/a | 2.41 s | n/a | 0.182 s | 0.492 s | 1.37 s | 0.007 s | 0.021 | 0.57 | 0.014 |
| mignotte-200.refine53 | 0.0127 s | n/a | 52.7 s | n/a | 1.17 s | 14.4 s | 14.4 s | 0.0222 s | 0.011 | 0.27 | 0.0015 |
| mignotte-400.refine53 | 0.0677 s | n/a | timeout | n/a | 17.6 s | timeout | n/a | n/a | 0.0038 |  |  |
| random1-100.refine53 | 0.000783 s | n/a | 0.0356 s | n/a | 0.295 s | 0 s | 0.00111 s | 0.00103 s | 0.022 | 0.031 |  |
| random1-200.refine53 | 0.000867 s | n/a | 0.226 s | n/a | 1.03 s | 0.00779 s | 0.00174 s | 0.00125 s | 0.11 | 0.0077 | 0.16 |
| random1-500.refine53 | 0.00377 s | n/a | 2.07 s | n/a | 42.1 s | 0.121 s | 0.00947 s | 0.00961 s | 0.031 | 0.0046 | 0.079 |
| random1-1000.refine53 | 0.00808 s | n/a | 11.8 s | n/a | 51.7 s | 0.573 s | 0.0213 s | 0.0212 s | 0.014 | 0.0018 | 0.037 |
| random1-2000.refine53 | 0.0417 s | n/a | 67.8 s | n/a | timeout | 11.5 s | 0.118 s | 0.112 s | 0.0036 | 0.0017 | 0.0097 |
| random1024-100.refine53 | 0.00419 s | n/a | 0.0396 s | n/a | 0.568 s | 0.135 s | 0.0055 s | 0.00522 s | 0.11 | 0.14 | 0.039 |
| random1024-200.refine53 | 0.00576 s | n/a | 0.215 s | n/a | 1.61 s | 0.192 s | 0.00766 s | 0.00653 s | 0.03 | 0.036 | 0.034 |
| random1024-500.refine53 | 0.0159 s | n/a | 1.73 s | n/a | 16.4 s | 0.585 s | 0.0264 s | 0.0252 s | 0.027 | 0.015 | 0.043 |
| random1024-1000.refine53 | 0.0361 s | n/a | 12 s | n/a | 85.2 s | 2.22 s | 0.0591 s | 0.0597 s | 0.016 | 0.0049 | 0.027 |
| random1024-2000.refine53 | 0.188 s | n/a | 108 s | n/a | timeout | 39 s | 0.541 s | 0.488 s | 0.0048 | 0.005 | 0.013 |
| random64-100.refine53 | 2.7e-05 s | n/a | 0.0473 s | n/a | 0.00633 s | 0.0115 s | n/a | n/a | 0.0043 |  |  |
| random64-200.refine53 | 0.00102 s | n/a | 0.254 s | n/a | 0.717 s | 0.0219 s | 0.00201 s | 0.00152 s | 0.047 | 0.0079 | 0.069 |
| random64-500.refine53 | 0.00822 s | n/a | 2.06 s | n/a | 16.4 s | 0.165 s | 0.0181 s | 0.0177 s | 0.05 | 0.0088 | 0.11 |
| random64-1000.refine53 | 0.0253 s | n/a | 12 s | n/a | timeout | 1.74 s | 0.0517 s | 0.0504 s | 0.015 | 0.0043 | 0.029 |
| random64-2000.refine53 | 0.0472 s | n/a | 50.6 s | n/a | · | 13.8 s | 0.0998 s | 0.14 s | 0.0034 | 0.002 | 0.01 |
| repeated-4.refine53 | 0.000249 s | n/a | error | n/a | 0.00318 s | 0 s | n/a | 0.000324 s | 0.078 |  |  |
| repeated-8.refine53 | 0.000705 s | n/a | error | n/a | 0.0126 s | 0.0063 s | n/a | 0.000965 s | 0.11 |  | 0.15 |
| repeated-16.refine53 | 0.0018 s | n/a | error | n/a | 0.0464 s | 0.0071 s | n/a | 0.0033 s | 0.25 |  | 0.47 |
| repeated-32.refine53 | 0.00539 s | n/a | error | n/a | 0.256 s | 0.0211 s | n/a | 0.00658 s | 0.26 |  | 0.31 |
| repeated-64.refine53 | 0.0251 s | n/a | error | n/a | 2.3 s | 0.0889 s | n/a | 0.0313 s | 0.28 |  | 0.35 |
| wilkinson-20.refine53 | 0.000171 s | n/a | 0.00298 s | n/a | 0.000962 s | 0 s | 0.000121 s | n/a | 0.18 | 0.041 |  |
| wilkinson-40.refine53 | 0.000521 s | n/a | 0.0357 s | n/a | 0.00457 s | 0 s | 0.000536 s | n/a | 0.11 | 0.015 |  |
| wilkinson-80.refine53 | 0.00348 s | n/a | 0.621 s | n/a | 0.0274 s | 0.00761 s | 0.00319 s | n/a | 0.46 | 0.0051 |  |
| wilkinson-160.refine53 | 0.0227 s | n/a | 7.88 s | n/a | 0.265 s | 0.0885 s | 0.0257 s | n/a | 0.26 | 0.0033 |  |
| wilkinson-320.refine53 | 0.183 s | n/a | timeout | n/a | 2.11 s | 0.408 s | n/a | n/a | 0.45 |  |  |
| wilkinson-640.refine53 | 1.33 s | n/a | · | n/a | 23.1 s | 2.93 s | n/a | n/a | 0.46 |  |  |

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
| mignotte-50.refine53 | -1.1 / -1.4, 581 b |  | -804.7 / -805.0, 2096 b | -0.0 / -1.6, 324 b | -79.0 / -79.0, 339 b | -805.7 / -805.9, 4159 b | -80.1 / -80.4, 759 b |
| mignotte-100.refine53 | -1.1 / -1.4, 1157 b |  | -1661.3 / -1661.5, 4142 b | -1.1 / -3.7, 658 b | -81.0 / -83.0, 673 b | -1662.3 / -1662.7, 8751 b | -82.1 / -82.4, 1339 b |
| mignotte-200.refine53 | -1.1 / -1.2, 2327 b |  | -3376.7 / -3376.8, 8236 b | -0.5 / -1.6, 1322 b | -83.0 / -87.0, 1333 b | -3377.1 / -3377.7, 18017 b | -83.2 / -84.1, 2515 b |
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
| chebyshev-100.refine1000 | 0.0716 s | n/a | 0.112 s | n/a | timeout | 1.12 s | 0.125 s | 0.136 s | 0.64 | 1.1 | 0.12 |
| chebyshev-200.refine1000 | 0.238 s | n/a | 1.28 s | n/a | · | 8.63 s | 0.43 s | 0.542 s | 0.19 | 0.33 | 0.063 |
| chebyshev-400.refine1000 | 1.12 s | n/a | 12.2 s | n/a | · | 71 s | 2.11 s | 1.72 s | 0.092 | 0.17 | 0.024 |
| chebyshev-700.refine1000 | 4.43 s | n/a | timeout | n/a | · | timeout | n/a | n/a |  |  |  |
| chebyshev-1000.refine1000 | 9.9 s | n/a | · | n/a | · | · | n/a | n/a |  |  |  |
| mignotte-10.refine1000 | 0.00151 s | n/a | 0.00357 s | n/a | 1.38 s | 0 s | 0.0022 s | 0.00267 s | 0.42 | 0.62 |  |
| mignotte-20.refine1000 | 0.00285 s | n/a | 0.0111 s | n/a | timeout | 0 s | 0.00386 s | 0.00361 s | 0.26 | 0.35 |  |
| mignotte-50.refine1000 | 0.207 s | n/a | 0.179 s | n/a | · | 0.0376 s | 0.929 s | 1.2 s | 5.5 | 5.2 | 32 |
| mignotte-100.refine1000 | 0.361 s | n/a | 2.41 s | n/a | · | 0.48 s | 1.24 s | 2.01 s | 0.75 | 0.52 | 4.2 |
| mignotte-200.refine1000 | 0.438 s | n/a | 54.6 s | n/a | · | 14.7 s | 14 s | 3.88 s | 0.03 | 0.26 | 0.26 |
| mignotte-400.refine1000 | 0.126 s | n/a | timeout | n/a | · | timeout | n/a | n/a |  |  |  |
| random1-100.refine1000 | 0.00295 s | n/a | 0.202 s | n/a | 36.3 s | 0.04 s | 0.00512 s | 0.00589 s | 0.074 | 0.025 | 0.15 |
| random1-200.refine1000 | 0.00282 s | n/a | 0.789 s | n/a | timeout | 0.0767 s | 0.00468 s | 0.00482 s | 0.037 | 0.0059 | 0.063 |
| random1-500.refine1000 | 0.0131 s | n/a | 5.25 s | n/a | · | 1.09 s | 0.0229 s | 0.0219 s | 0.012 | 0.0044 | 0.02 |
| random1-1000.refine1000 | 0.0211 s | n/a | 21.3 s | n/a | · | 3.66 s | 0.0321 s | 0.0356 s | 0.0058 | 0.0015 | 0.0097 |
| random1-2000.refine1000 | 0.101 s | n/a | 103 s | n/a | · | 27.6 s | 0.153 s | 0.179 s | 0.0037 | 0.0015 | 0.0065 |
| random1024-100.refine1000 | 0.00897 s | n/a | 0.274 s | n/a | 62.2 s | 0.175 s | 0.0126 s | 0.0139 s | 0.051 | 0.046 | 0.079 |
| random1024-200.refine1000 | 0.0111 s | n/a | 1.03 s | n/a | timeout | 0.35 s | 0.0152 s | 0.0166 s | 0.032 | 0.015 | 0.047 |
| random1024-500.refine1000 | 0.0307 s | n/a | 6.39 s | n/a | · | 1.71 s | 0.0446 s | 0.0474 s | 0.018 | 0.007 | 0.028 |
| random1024-1000.refine1000 | 0.0643 s | n/a | 24.7 s | n/a | · | 6.87 s | 0.0855 s | 0.1 s | 0.0094 | 0.0035 | 0.015 |
| random1024-2000.refine1000 | 0.371 s | n/a | 109 s | n/a | · | 99.5 s | 0.504 s | 0.657 s | 0.0037 | 0.0046 | 0.0066 |
| random64-100.refine1000 | 2.88e-05 s | n/a | 0.216 s | n/a | 0.00622 s | 0.00882 s | n/a | n/a | 0.0046 |  |  |
| random64-200.refine1000 | 0.00353 s | n/a | 0.816 s | n/a | timeout | 0.0928 s | 0.00549 s | 0.00565 s | 0.038 | 0.0067 | 0.061 |
| random64-500.refine1000 | 0.024 s | n/a | 5.32 s | n/a | · | 1.47 s | 0.0332 s | 0.0345 s | 0.016 | 0.0062 | 0.024 |
| random64-1000.refine1000 | 0.0503 s | n/a | 20.6 s | n/a | · | 11.5 s | 0.0888 s | 0.102 s | 0.0044 | 0.0043 | 0.0089 |
| random64-2000.refine1000 | 0.106 s | n/a | 87.7 s | n/a | · | 41.5 s | 0.162 s | 0.244 s | 0.0026 | 0.0018 | 0.0059 |
| repeated-4.refine1000 | 0.00133 s | n/a | error | n/a | 0.0808 s | 0.00364 s | n/a | 0.00289 s | 0.37 |  | 0.79 |
| repeated-8.refine1000 | 0.00313 s | n/a | error | n/a | 0.346 s | 0.00399 s | n/a | 0.00546 s | 0.79 |  | 1.4 |
| repeated-16.refine1000 | 0.00679 s | n/a | error | n/a | 1.94 s | 0.0283 s | n/a | 0.0113 s | 0.24 |  | 0.4 |
| repeated-32.refine1000 | 0.0215 s | n/a | error | n/a | 18.4 s | 0.121 s | n/a | 0.0336 s | 0.18 |  | 0.28 |
| repeated-64.refine1000 | 0.0906 s | n/a | error | n/a | timeout | 1.16 s | n/a | 0.122 s | 0.078 |  | 0.11 |
| wilkinson-20.refine1000 | 0.000126 s | n/a | 0.0107 s | n/a | 0.000945 s | 0 s | 0.000128 s | n/a | 0.13 | 0.012 |  |
| wilkinson-40.refine1000 | 0.000531 s | n/a | 0.0531 s | n/a | 0.00452 s | 0 s | 0.000545 s | n/a | 0.12 | 0.01 |  |
| wilkinson-80.refine1000 | 0.0034 s | n/a | 0.639 s | n/a | 0.0274 s | 0.00862 s | 0.00327 s | n/a | 0.39 | 0.0051 |  |
| wilkinson-160.refine1000 | 0.0206 s | n/a | 8.41 s | n/a | 0.292 s | 0.0772 s | 0.0217 s | n/a | 0.27 | 0.0026 |  |
| wilkinson-320.refine1000 | 0.164 s | n/a | timeout | n/a | 2.17 s | 0.461 s | n/a | n/a | 0.36 |  |  |
| wilkinson-640.refine1000 | 1.49 s | n/a | · | n/a | 23.6 s | 2.77 s | n/a | n/a | 0.54 |  |  |

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
| mignotte-50.refine1000 | -0.4 / -1.0, 4791 b |  | -881.4 / -881.9, 4142 b |  | -1084.0 / -1180.0, 4443 b | -882.3 / -882.6, 8667 b | -1084.2 / -1085.1, 9577 b |
| mignotte-100.refine1000 | -1.1 / -1.3, 5615 b |  | -714.3 / -714.5, 4142 b |  | -1054.0 / -1056.0, 4769 b | -715.3 / -715.7, 8751 b | -1054.5 / -1055.0, 10245 b |
| mignotte-200.refine1000 | -1.2 / -1.2, 6849 b |  | -2429.7 / -2429.8, 8236 b |  | -1056.0 / -1060.0, 4417 b | -2430.1 / -2430.7, 18017 b | -1056.2 / -1057.2, 11971 b |
| mignotte-400.refine1000 | -1.0 / -1.1, 6577 b |  |  |  |  |  |  |
| random1-100.refine1000 | -1.9 / -2.0, 2020 b |  | -1036.5 / -1037.3, 4138 b | -0.1 / -0.3, 1006 b | -1056.0 / -1057.0, 4144 b | -1038.9 / -1039.0, 4094 b | -1056.0 / -1057.9, 4132 b |
| random1-200.refine1000 | -1.2 / -1.2, 2060 b |  | -1036.6 / -1036.6, 4136 b |  | -1056.0 / -1056.0, 4118 b | -1038.2 / -1038.2, 4134 b | -1057.2 / -1057.2, 4172 b |
| random1-500.refine1000 | -1.1 / -1.3, 2129 b |  | -1034.7 / -1035.5, 4131 b |  | -1051.0 / -1064.0, 4137 b | -1036.1 / -1036.3, 4199 b | -1052.1 / -1052.3, 4231 b |
| random1-1000.refine1000 | -1.1 / -1.1, 2504 b |  | -1034.0 / -1034.0, 4131 b |  | -1122.0 / -1122.0, 6161 b | -1035.1 / -1035.1, 4572 b | -1123.1 / -1123.1, 4748 b |
| random1-2000.refine1000 | -1.1 / -1.4, 2512 b |  | -1031.1 / -1033.2, 4131 b |  | -1066.0 / -1090.0, 6170 b | -1033.1 / -1033.4, 4576 b | -1067.1 / -1067.4, 4644 b |
| random1024-100.refine1000 | -1.1 / -2.0, 2259 b |  | -1036.8 / -1037.8, 4147 b | -0.6 / -1.2, 1007 b | -1047.0 / -1058.0, 4128 b | -1038.1 / -1039.0, 4333 b | -1048.1 / -1049.0, 4353 b |
| random1024-200.refine1000 | -1.1 / -1.2, 2020 b |  | -1035.7 / -1036.7, 4134 b |  | -1051.0 / -1055.0, 4368 b | -1037.1 / -1037.2, 4092 b | -1052.1 / -1052.2, 4122 b |
| random1024-500.refine1000 | -1.1 / -1.2, 2072 b |  | -1033.8 / -1034.5, 4134 b |  | -1060.0 / -1062.0, 4177 b | -1035.1 / -1035.2, 4140 b | -1061.1 / -1061.2, 4192 b |
| random1024-1000.refine1000 | -1.4 / -1.6, 2329 b |  | -1032.8 / -1032.9, 4133 b |  | -1062.0 / -1062.0, 6161 b | -1034.4 / -1034.6, 4395 b | -1063.4 / -1063.6, 4453 b |
| random1024-2000.refine1000 | -1.0 / -1.4, 2265 b |  | -1029.5 / -1032.8, 4130 b |  | -1089.0 / -2082.0, 6422 b | -1031.0 / -1031.4, 4325 b | -1090.0 / -1090.4, 4443 b |
| random64-100.refine1000 |  |  |  |  |  |  |  |
| random64-200.refine1000 | -1.3 / -1.3, 2067 b |  | -1034.7 / -1034.7, 4136 b |  | -1055.0 / -1055.0, 4122 b | -1036.3 / -1036.3, 4137 b | -1056.3 / -1056.3, 4177 b |
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

## Exact real roots (Poly Mojo `real_roots`, SymPy `real_roots`)

| Case | Poly Mojo | Poly Mojo (Sturm count) | FLINT | FLINT (Sturm count) | SymPy | msolve | Poly Mojo / best |
|---|---|---|---|---|---|---|---|
| chebyshev-100.exact | 0.08 s | n/a | n/a | n/a | error | n/a |  |
| chebyshev-200.exact | 0.721 s | n/a | n/a | n/a | error | n/a |  |
| chebyshev-400.exact | 7.47 s | n/a | n/a | n/a | timeout | n/a |  |
| chebyshev-700.exact | 41.7 s | n/a | n/a | n/a | · | n/a |  |
| chebyshev-1000.exact | timeout | n/a | n/a | n/a | · | n/a |  |
| mignotte-10.exact | 0.000378 s | n/a | n/a | n/a | error | n/a |  |
| mignotte-20.exact | 0.000372 s | n/a | n/a | n/a | timeout | n/a |  |
| mignotte-50.exact | 0.00297 s | n/a | n/a | n/a | · | n/a |  |
| mignotte-100.exact | 0.035 s | n/a | n/a | n/a | · | n/a |  |
| mignotte-200.exact | 0.802 s | n/a | n/a | n/a | · | n/a |  |
| mignotte-400.exact | 20.6 s | n/a | n/a | n/a | · | n/a |  |
| random1-100.exact | 0.00723 s | n/a | n/a | n/a | error | n/a |  |
| random1-200.exact | 0.0298 s | n/a | n/a | n/a | error | n/a |  |
| random1-500.exact | 0.235 s | n/a | n/a | n/a | error | n/a |  |
| random1-1000.exact | 0.886 s | n/a | n/a | n/a | timeout | n/a |  |
| random1-2000.exact | 11.6 s | n/a | n/a | n/a | · | n/a |  |
| random1024-100.exact | 0.0127 s | n/a | n/a | n/a | error | n/a |  |
| random1024-200.exact | 0.118 s | n/a | n/a | n/a | error | n/a |  |
| random1024-500.exact | 0.408 s | n/a | n/a | n/a | timeout | n/a |  |
| random1024-1000.exact | 5.9 s | n/a | n/a | n/a | · | n/a |  |
| random1024-2000.exact | 27.6 s | n/a | n/a | n/a | · | n/a |  |
| random64-100.exact | 0.0112 s | n/a | n/a | n/a | 0.327 s | n/a | 0.034 |
| random64-200.exact | 0.0498 s | n/a | n/a | n/a | error | n/a |  |
| random64-500.exact | 0.6 s | n/a | n/a | n/a | error | n/a |  |
| random64-1000.exact | 5.12 s | n/a | n/a | n/a | timeout | n/a |  |
| random64-2000.exact | 9.77 s | n/a | n/a | n/a | · | n/a |  |
| repeated-4.exact | 0.000597 s | n/a | n/a | n/a | error | n/a |  |
| repeated-8.exact | 0.00181 s | n/a | n/a | n/a | error | n/a |  |
| repeated-16.exact | 0.0115 s | n/a | n/a | n/a | error | n/a |  |
| repeated-32.exact | 0.0639 s | n/a | n/a | n/a | error | n/a |  |
| repeated-64.exact | 0.576 s | n/a | n/a | n/a | error | n/a |  |
| wilkinson-20.exact | 0.0013 s | n/a | n/a | n/a | 0.0282 s | n/a | 0.046 |
| wilkinson-40.exact | 0.00547 s | n/a | n/a | n/a | 0.0748 s | n/a | 0.073 |
| wilkinson-80.exact | 0.0241 s | n/a | n/a | n/a | 0.291 s | n/a | 0.083 |
| wilkinson-160.exact | 0.105 s | n/a | n/a | n/a | 1.37 s | n/a | 0.077 |
| wilkinson-320.exact | 0.625 s | n/a | n/a | n/a | 9.58 s | n/a | 0.065 |
| wilkinson-640.exact | 3.73 s | n/a | n/a | n/a | 87.8 s | n/a | 0.042 |

## Distinct real roots (Poly Mojo `dense.count_roots`, Descartes, on a power-of-two root bound)

| Case | Poly Mojo | Poly Mojo (Sturm count) | FLINT | FLINT (Sturm count) | SymPy | msolve | Poly Mojo / best |
|---|---|---|---|---|---|---|---|
| chebyshev-100.count | 0.00269 s | 0.00214 s | 0.0307 s | 0.0309 s | 0.00512 s | n/a | 0.53 |
| chebyshev-200.count | 0.0213 s | 0.0109 s | 1.15 s | 1.25 s | 0.0214 s | n/a | 0.99 |
| chebyshev-400.count | 0.175 s | 0.0452 s | 30.5 s | 31.5 s | 0.0913 s | n/a | 1.9 |
| chebyshev-700.count | 1.94 s | 0.167 s | timeout | timeout | 0.397 s | n/a | 4.9 |
| chebyshev-1000.count | 5.59 s | 0.367 s | · | · | 0.889 s | n/a | 6.3 |
| mignotte-10.count | 6.02e-05 s | 2.41e-05 s | 1.76e-06 s | 1.09e-06 s | 0.000215 s | n/a | 55 |
| mignotte-20.count | 7.76e-05 s | 4.64e-05 s | 2.45e-06 s | 2.55e-06 s | 0.000334 s | n/a | 32 |
| mignotte-50.count | 0.00024 s | 0.000174 s | 8.82e-06 s | 8.62e-06 s | 0.000894 s | n/a | 28 |
| mignotte-100.count | 0.000729 s | 0.000593 s | 2.58e-05 s | 2.67e-05 s | 0.00199 s | n/a | 28 |
| mignotte-200.count | 0.00301 s | 0.0018 s | 8.73e-05 s | 8.62e-05 s | 0.0071 s | n/a | 35 |
| mignotte-400.count | 0.0173 s | 0.00748 s | 0.00033 s | 0.000363 s | 0.0241 s | n/a | 52 |
| random1-100.count | 0.000608 s | 0.168 s | 0.00156 s | 0.00148 s | 1.38 s | n/a | 0.41 |
| random1-200.count | 0.00261 s | 3.49 s | 0.0301 s | 0.0336 s | 66.5 s | n/a | 0.086 |
| random1-500.count | 0.0258 s | timeout | 0.954 s | 0.782 s | timeout | n/a | 0.033 |
| random1-1000.count | 0.207 s | · | 13.7 s | 12.9 s | · | n/a | 0.016 |
| random1-2000.count | 1.82 s | · | timeout | timeout | · | n/a |  |
| random1024-100.count | 0.00203 s | timeout | 6.14 s | 5.77 s | timeout | n/a | 0.00035 |
| random1024-200.count | 0.00842 s | · | 56.6 s | 54 s | · | n/a | 0.00016 |
| random1024-500.count | 0.0597 s | · | timeout | timeout | · | n/a |  |
| random1024-1000.count | 0.385 s | · | · | · | · | n/a |  |
| random1024-2000.count | 5.28 s | · | · | · | · | n/a |  |
| random64-100.count | 0.00165 s | 11.7 s | 0.0916 s | 0.104 s | 102 s | n/a | 0.018 |
| random64-200.count | 0.00736 s | timeout | 1.28 s | 1.25 s | timeout | n/a | 0.0059 |
| random64-500.count | 0.17 s | · | 32 s | 32 s | · | n/a | 0.0053 |
| random64-1000.count | 0.269 s | · | timeout | timeout | · | n/a |  |
| random64-2000.count | 1.85 s | · | · | · | · | n/a |  |
| repeated-4.count | 7.94e-05 s | 0.00111 s | 2.44e-05 s | 2.44e-05 s | 0.000535 s | n/a | 3.3 |
| repeated-8.count | 0.000282 s | 0.0108 s | 0.000365 s | 0.000417 s | 0.00365 s | n/a | 0.77 |
| repeated-16.count | 0.00136 s | 0.0996 s | 0.0137 s | 0.0158 s | 0.123 s | n/a | 0.099 |
| repeated-32.count | 0.00686 s | 5.47 s | 0.688 s | 0.546 s | 11.7 s | n/a | 0.013 |
| repeated-64.count | 0.0388 s | timeout | 19.9 s | 20.3 s | timeout | n/a | 0.002 |
| wilkinson-20.count | 0.00016 s | 0.000219 s | 3.47e-05 s | 4.1e-05 s | 0.000838 s | n/a | 4.6 |
| wilkinson-40.count | 0.000529 s | 0.00182 s | 0.00176 s | 0.00169 s | 0.00255 s | n/a | 0.31 |
| wilkinson-80.count | 0.00288 s | 0.0116 s | 0.0975 s | 0.104 s | 0.0112 s | n/a | 0.26 |
| wilkinson-160.count | 0.026 s | 0.0819 s | 4.15 s | 3.91 s | 0.0808 s | n/a | 0.32 |
| wilkinson-320.count | 0.196 s | 0.736 s | 118 s | 116 s | 0.652 s | n/a | 0.3 |
| wilkinson-640.count | 1.97 s | 5.07 s | timeout | timeout | 6.11 s | n/a | 0.32 |

