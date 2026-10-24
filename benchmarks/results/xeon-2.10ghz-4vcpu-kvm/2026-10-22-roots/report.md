# Real roots

Host: Intel(R) Xeon(R) Processor @ 2.10GHz, 4 vCPUs, 16 GiB RAM, Ubuntu 24.04.5 LTS. Runs are pinned to one CPU, capped at 12G and 120 s.

Each answer is checked against an Arb reference: the number of distinct real roots, and every interval holding exactly one reference root (`wrong` otherwise). The `Poly Mojo / best` column divides Poly Mojo's time by the fastest other library on the case, each at its own achieved width.

Widths differ by tool: Poly Mojo and SymPy return intervals just under the requested width, while FLINT and msolve work at a fixed precision and return intervals many bits narrower. The isolate and refine tables therefore also time Poly Mojo at each one's achieved width (`Poly Mojo at FLINT width`, `Poly Mojo at msolve width`): 2^-b for the smallest power of two at or below that tool's widest interval on the case, never looser than requested. `matched / peer` divides that time by the peer's. The width tables give each tool's widest and median interval as log2 of width over 2^-bits (0 = exactly the requested width, -10 = 2^10 times narrower), and the longest endpoint in bits.

## Isolation (Poly Mojo `real_root_cover`, width 1/2)

| Case | Poly Mojo | FLINT | FLINT (Sturm count) | SymPy | msolve | Poly Mojo at FLINT width | Poly Mojo at msolve width | Poly Mojo / best | FLINT: matched / peer | msolve: matched / peer |
|---|---|---|---|---|---|---|---|---|---|---|
| chebyshev-100.isolate | 0.00286 s | 0.0785 s | n/a | 0.367 s | 0.0303 s | 0.026 s | 0.0182 s | 0.094 | 0.33 | 0.6 |
| chebyshev-200.isolate | 0.023 s | 0.868 s | n/a | 3.21 s | 0.214 s | 0.186 s | 0.104 s | 0.11 | 0.21 | 0.49 |
| chebyshev-400.isolate | 0.189 s | 12 s | n/a | 35.9 s | 2.79 s | 1.02 s | 0.764 s | 0.068 | 0.085 | 0.27 |
| chebyshev-700.isolate | 1.86 s | timeout | n/a | timeout | 31.2 s | n/a | 4.64 s | 0.06 |  | 0.15 |
| chebyshev-1000.isolate | 5.33 s | · | n/a | · | timeout | n/a | n/a |  |  |  |
| mignotte-10.isolate | 6.12e-05 s | 0.000811 s | n/a | 0.000532 s | 0 s | 0.00106 s | 0.000154 s | 0.12 | 1.3 |  |
| mignotte-20.isolate | 6.14e-05 s | 0.00536 s | n/a | 0.000916 s | 0 s | 0.00225 s | 0.000266 s | 0.067 | 0.42 |  |
| mignotte-50.isolate | 0.000217 s | 0.14 s | n/a | 0.00351 s | 0.0216 s | 0.15 s | 0.00156 s | 0.062 | 1.1 | 0.072 |
| mignotte-100.isolate | 0.000731 s | 2.39 s | n/a | 0.0112 s | 0.478 s | 1.26 s | 0.00591 s | 0.065 | 0.53 | 0.012 |
| mignotte-200.isolate | 0.00304 s | 50.2 s | n/a | 0.0448 s | 14.5 s | 13.5 s | 0.0266 s | 0.068 | 0.27 | 0.0018 |
| mignotte-400.isolate | 0.0189 s | timeout | n/a | 0.195 s | timeout | n/a | n/a | 0.097 |  |  |
| random1-100.isolate | 0.000526 s | 0.0242 s | n/a | 0.00697 s | 0 s | 0.00121 s | 0.00121 s | 0.075 | 0.05 |  |
| random1-200.isolate | 0.00253 s | 0.159 s | n/a | 0.023 s | 0.000249 s | 0.00314 s | 0.0036 s | 10 | 0.02 | 15 |
| random1-500.isolate | 0.02 s | 1.62 s | n/a | 0.213 s | 0.0979 s | 0.0309 s | 0.0292 s | 0.2 | 0.019 | 0.3 |
| random1-1000.isolate | 0.185 s | 6.72 s | n/a | 0.497 s | 0.476 s | 0.192 s | 0.207 s | 0.39 | 0.029 | 0.43 |
| random1-2000.isolate | 1.6 s | 53 s | n/a | 4.86 s | 10.3 s | 1.79 s | 1.68 s | 0.33 | 0.034 | 0.16 |
| random1024-100.isolate | 0.00182 s | 0.0254 s | n/a | 0.0147 s | 0.0906 s | 0.00562 s | 0.00572 s | 0.12 | 0.22 | 0.063 |
| random1024-200.isolate | 0.0069 s | 0.153 s | n/a | 0.0604 s | 0.171 s | 0.0134 s | 0.0137 s | 0.11 | 0.087 | 0.08 |
| random1024-500.isolate | 0.0475 s | 1.57 s | n/a | 0.296 s | 0.5 s | 0.0721 s | 0.066 s | 0.16 | 0.046 | 0.13 |
| random1024-1000.isolate | 0.371 s | 5.65 s | n/a | 2 s | 2.04 s | 0.427 s | 0.434 s | 0.19 | 0.076 | 0.21 |
| random1024-2000.isolate | 5.16 s | 50 s | n/a | 31.1 s | 33.5 s | 5.49 s | 5.67 s | 0.17 | 0.11 | 0.17 |
| random64-100.isolate | 0.00149 s | 0.0381 s | n/a | 0.00611 s | 0 s | n/a | n/a | 0.24 |  |  |
| random64-200.isolate | 0.00669 s | 0.187 s | n/a | 0.0329 s | 0.00655 s | 0.00801 s | 0.00798 s | 1 | 0.043 | 1.2 |
| random64-500.isolate | 0.163 s | 1.89 s | n/a | 0.364 s | 0.152 s | 0.187 s | 0.189 s | 1.1 | 0.099 | 1.2 |
| random64-1000.isolate | 0.279 s | 7.77 s | n/a | 1.87 s | 1.39 s | 0.312 s | 0.316 s | 0.2 | 0.04 | 0.23 |
| random64-2000.isolate | 1.67 s | 42.8 s | n/a | 3.95 s | 12.2 s | 1.79 s | 1.89 s | 0.42 | 0.042 | 0.16 |
| repeated-4.isolate | 6.92e-05 s | error | n/a | wrong | 0 s | n/a | 0.000153 s |  |  |  |
| repeated-8.isolate | 0.000262 s | error | n/a | wrong | 0 s | n/a | 0.000633 s |  |  |  |
| repeated-16.isolate | 0.00128 s | error | n/a | wrong | 0 s | n/a | 0.00239 s |  |  |  |
| repeated-32.isolate | 0.00597 s | error | n/a | wrong | 0.009 s | n/a | 0.0111 s | 0.66 |  | 1.2 |
| repeated-64.isolate | 0.0416 s | error | n/a | wrong | 0.0694 s | n/a | 0.0674 s | 0.6 |  | 0.97 |
| wilkinson-20.isolate | 0.000137 s | 0.00272 s | n/a | 0.000916 s | 0 s | 0.000262 s | n/a | 0.15 | 0.096 |  |
| wilkinson-40.isolate | 0.000566 s | 0.0236 s | n/a | 0.00452 s | 0.00219 s | 0.00112 s | n/a | 0.26 | 0.047 |  |
| wilkinson-80.isolate | 0.00308 s | 0.547 s | n/a | 0.0282 s | 0.0252 s | 0.00641 s | n/a | 0.12 | 0.012 |  |
| wilkinson-160.isolate | 0.0235 s | 7.63 s | n/a | 0.217 s | 0.103 s | 0.0463 s | n/a | 0.23 | 0.0061 |  |
| wilkinson-320.isolate | 0.204 s | timeout | n/a | 2.07 s | 0.424 s | n/a | n/a | 0.48 |  |  |
| wilkinson-640.isolate | 1.98 s | · | n/a | 22.1 s | 2.96 s | n/a | n/a | 0.67 |  |  |

Achieved widths (isolate): widest / median as log2 of width, longest endpoint in bits.

| Case | Poly Mojo | FLINT | SymPy | msolve | Poly Mojo at FLINT width | Poly Mojo at msolve width |
|---|---|---|---|---|---|---|
| chebyshev-100.isolate | -7.0 / -8.0, 21 b | -133.8 / -161.2, 562 b | -4.5 / -6.3, 20 b | -57.0 / -68.0, 153 b | -134.0 / -135.0, 287 b | -58.1 / -59.0, 133 b |
| chebyshev-200.isolate | -8.0 / -9.0, 45 b | -263.2 / -322.4, 1075 b | -5.4 / -7.2, 24 b | -108.0 / -134.0, 297 b | -264.0 / -266.0, 551 b | -109.0 / -110.0, 241 b |
| chebyshev-400.isolate | -9.0 / -10.0, 49 b | -518.4 / -640.5, 2100 b | -6.1 / -8.2, 28 b | -209.0 / -264.0, 567 b | -519.0 / -520.0, 1065 b | -210.0 / -211.0, 445 b |
| chebyshev-700.isolate | -10.0 / -11.0, 51 b |  |  | -389.0 / -395.0, 1059 b |  | -390.0 / -391.0, 805 b |
| chebyshev-1000.isolate | -11.0 / -11.0, 53 b |  |  |  |  |  |
| mignotte-10.isolate | -2.0 / -2.1, 162 b | -223.3 / -223.3, 561 b | 0.0 / 0.0, 8 b | -9.0 / -9.0, 73 b | -225.1 / -225.3, 629 b | -10.1 / -10.4, 179 b |
| mignotte-20.isolate | -1.0 / -1.2, 149 b | -445.1 / -445.1, 1070 b | 0.0 / 0.0, 8 b | -14.0 / -16.0, 139 b | -447.5 / -447.5, 1091 b | -14.1 / -15.4, 185 b |
| mignotte-50.isolate | -1.0 / -2.4, 443 b | -857.7 / -858.0, 2096 b | 0.0 / 0.0, 8 b | -36.0 / -36.0, 339 b | -858.7 / -858.9, 4159 b | -37.1 / -37.2, 535 b |
| mignotte-100.isolate | -1.0 / -2.4, 1019 b | -1714.3 / -1714.5, 4142 b | 0.0 / 0.0, 8 b | -70.0 / -72.0, 673 b | -1715.3 / -1715.7, 8751 b | -71.4 / -71.6, 1199 b |
| mignotte-200.isolate | -1.0 / -2.4, 2189 b | -3429.7 / -3429.8, 8236 b | 0.0 / 0.0, 8 b | -136.0 / -140.0, 1333 b | -3430.1 / -3430.7, 18017 b | -136.2 / -137.1, 2515 b |
| mignotte-400.isolate | -1.0 / -2.4, 4523 b |  | 0.0 / 0.0, 8 b |  |  |  |
| random1-100.isolate | -2.0 / -3.0, 10 b | -54.6 / -55.0, 172 b | 0.0 / 0.0, 4 b | -58.0 / -71.0, 148 b | -56.9 / -57.0, 130 b | -59.9 / -60.0, 136 b |
| random1-200.isolate | -3.0 / -3.0, 8 b | -117.8 / -117.8, 300 b | 0.0 / 0.0, 3 b | -136.0 / -136.0, 278 b | -119.2 / -119.2, 296 b | -137.2 / -137.2, 332 b |
| random1-500.isolate | -1.4 / -3.0, 9 b | -500.5 / -500.7, 1066 b | 0.0 / 0.0, 6 b | -387.0 / -400.0, 809 b | -502.1 / -502.3, 1131 b | -388.1 / -388.3, 903 b |
| random1-1000.isolate | -1.4 / -1.4, 8 b | -499.9 / -499.9, 1062 b | 0.0 / 0.0, 3 b | -776.0 / -776.0, 1686 b | -501.1 / -501.1, 1021 b | -777.1 / -777.1, 2056 b |
| random1-2000.isolate | -3.0 / -3.0, 9 b | -1007.1 / -1008.9, 2081 b | 0.0 / 0.0, 5 b | -1066.0 / -1290.0, 3109 b | -1009.1 / -1009.4, 2528 b | -1067.1 / -1067.4, 2644 b |
| random1024-100.isolate | -3.0 / -4.0, 11 b | -54.1 / -54.4, 172 b | 0.0 / -1.6, 7 b | -63.0 / -74.0, 160 b | -56.4 / -57.0, 118 b | -64.4 / -65.0, 134 b |
| random1024-200.isolate | -3.0 / -3.0, 11 b | -116.6 / -116.6, 296 b | -1.0 / -1.0, 10 b | -131.0 / -135.0, 287 b | -118.1 / -118.2, 254 b | -132.1 / -132.2, 282 b |
| random1024-500.isolate | -1.4 / -3.0, 10 b | -499.8 / -499.9, 1064 b | 0.0 / 0.0, 5 b | -296.0 / -396.0, 806 b | -501.1 / -501.2, 1072 b | -297.1 / -297.2, 664 b |
| random1024-1000.isolate | -1.0 / -1.4, 189 b | -498.0 / -498.3, 1061 b | 0.0 / 0.0, 16 b | -776.0 / -782.0, 1581 b | -499.4 / -499.6, 1325 b | -777.4 / -777.6, 1881 b |
| random1024-2000.isolate | -1.0 / -5.0, 15 b | -1009.5 / -1009.9, 2084 b | 0.0 / -5.8, 14 b | -1541.0 / -1549.0, 3350 b | -1011.0 / -1011.4, 2285 b | -1542.0 / -1542.4, 3347 b |
| random64-100.isolate |  |  |  |  |  |  |
| random64-200.isolate | -1.0 / -1.0, 8 b | -116.6 / -116.6, 298 b | 0.0 / 0.0, 3 b | -135.0 / -135.0, 282 b | -118.3 / -118.3, 271 b | -136.3 / -136.3, 307 b |
| random64-500.isolate | -2.0 / -2.2, 16 b | -495.4 / -499.2, 1064 b | 0.0 / 0.0, 6 b | -388.0 / -393.0, 1046 b | -497.0 / -497.4, 1112 b | -389.0 / -389.4, 896 b |
| random64-1000.isolate | -1.0 / -4.0, 166 b | -498.8 / -499.9, 1063 b | 0.0 / -1.6, 14 b | -769.0 / -776.0, 1584 b | -500.0 / -500.1, 1196 b | -770.0 / -770.1, 1736 b |
| random64-2000.isolate | -2.0 / -2.0, 10 b | -1005.8 / -1007.4, 2081 b | 0.0 / 0.0, 4 b | -1542.0 / -1547.0, 3108 b | -1007.2 / -1007.4, 2134 b | -1543.2 / -1543.4, 3206 b |
| repeated-4.isolate | -1.0 / -3.0, 11 b |  |  | -8.0 / -9.0, 19 b |  | -9.0 / -9.0, 22 b |
| repeated-8.isolate | -2.0 / -2.0, 29 b |  |  | -13.0 / -15.0, 33 b |  | -13.7 / -14.3, 50 b |
| repeated-16.isolate | -2.0 / -4.0, 31 b |  |  | -24.0 / -27.0, 59 b |  | -25.0 / -25.7, 69 b |
| repeated-32.isolate | -2.0 / -6.0, 35 b |  |  | -41.0 / -52.0, 118 b |  | -42.0 / -42.8, 104 b |
| repeated-64.isolate | -2.0 / -7.0, 39 b |  |  | -75.0 / -101.0, 219 b |  | -76.0 / -76.3, 218 b |
| wilkinson-20.isolate | -2.0 / -2.0, 11 b | -8.9 / -14.1, 162 b | points, 6 b | points, 6 b | points, 6 b |  |
| wilkinson-40.isolate | -2.0 / -2.0, 12 b | -21.2 / -31.9, 280 b | points, 7 b | points, 7 b | points, 7 b |  |
| wilkinson-80.isolate | -2.0 / -2.0, 13 b | -46.3 / -66.8, 530 b | points, 8 b | points, 8 b | points, 8 b |  |
| wilkinson-160.isolate | -2.0 / -2.0, 14 b | -97.5 / -140.3, 1038 b | points, 9 b | points, 9 b | points, 9 b |  |
| wilkinson-320.isolate | -2.0 / -2.0, 15 b |  | points, 10 b | points, 10 b |  |  |
| wilkinson-640.isolate | -2.0 / -2.0, 16 b |  | points, 11 b | points, 11 b |  |  |

## Refinement of every root to width 2^-53 (Poly Mojo `roots.refine` of the isolating cover)

| Case | Poly Mojo | FLINT | FLINT (Sturm count) | SymPy | msolve | Poly Mojo at FLINT width | Poly Mojo at msolve width | Poly Mojo / best | FLINT: matched / peer | msolve: matched / peer |
|---|---|---|---|---|---|---|---|---|---|---|
| chebyshev-100.refine53 | 0.0155 s | 0.0779 s | n/a | 8.01 s | 0.0455 s | 0.0244 s | 0.0198 s | 0.34 | 0.31 | 0.44 |
| chebyshev-200.refine53 | 0.0641 s | 0.881 s | n/a | 75.1 s | 0.223 s | 0.158 s | 0.0829 s | 0.29 | 0.18 | 0.37 |
| chebyshev-400.refine53 | 0.313 s | 12.1 s | n/a | timeout | 2.71 s | 0.779 s | 0.539 s | 0.12 | 0.064 | 0.2 |
| chebyshev-700.refine53 | 1.28 s | timeout | n/a | · | 29.8 s | n/a | 2.62 s | 0.043 |  | 0.088 |
| chebyshev-1000.refine53 | 3.24 s | · | n/a | · | timeout | n/a | n/a |  |  |  |
| mignotte-10.refine53 | 0.000543 s | 0.000802 s | n/a | 1.16 s | 0 s | 0.000925 s | 0.000744 s | 0.68 | 1.2 |  |
| mignotte-20.refine53 | 0.000773 s | 0.00525 s | n/a | 0.00927 s | 0 s | 0.00216 s | 0.00158 s | 0.15 | 0.41 |  |
| mignotte-50.refine53 | 0.00197 s | 0.135 s | n/a | 0.0532 s | 0.0234 s | 0.16 s | 0.00365 s | 0.084 | 1.2 | 0.16 |
| mignotte-100.refine53 | 0.00379 s | 2.36 s | n/a | 0.179 s | 0.475 s | 1.29 s | 0.00714 s | 0.021 | 0.55 | 0.015 |
| mignotte-200.refine53 | 0.0128 s | 50.1 s | n/a | 1.13 s | 14.3 s | 13.3 s | 0.0235 s | 0.011 | 0.27 | 0.0016 |
| mignotte-400.refine53 | 0.0727 s | timeout | n/a | 14.4 s | timeout | n/a | n/a | 0.005 |  |  |
| random1-100.refine53 | 0.000614 s | 0.0317 s | n/a | 0.244 s | 0 s | 0.000844 s | 0.000811 s | 0.019 | 0.027 |  |
| random1-200.refine53 | 0.000725 s | 0.197 s | n/a | 0.667 s | 0.00767 s | 0.00126 s | 0.00108 s | 0.094 | 0.0064 | 0.14 |
| random1-500.refine53 | 0.00303 s | 1.49 s | n/a | 35 s | 0.0907 s | 0.00674 s | 0.00819 s | 0.033 | 0.0045 | 0.09 |
| random1-1000.refine53 | 0.00729 s | 9.51 s | n/a | 43.5 s | 0.549 s | 0.018 s | 0.0175 s | 0.013 | 0.0019 | 0.032 |
| random1-2000.refine53 | 0.0445 s | 53.6 s | n/a | timeout | 10.3 s | 0.0963 s | 0.116 s | 0.0043 | 0.0018 | 0.011 |
| random1024-100.refine53 | 0.00378 s | 0.0284 s | n/a | 0.453 s | 0.0999 s | 0.00452 s | 0.00443 s | 0.13 | 0.16 | 0.044 |
| random1024-200.refine53 | 0.00515 s | 0.16 s | n/a | 1.42 s | 0.172 s | 0.00694 s | 0.0061 s | 0.032 | 0.043 | 0.036 |
| random1024-500.refine53 | 0.015 s | 1.53 s | n/a | 13.4 s | 0.511 s | 0.025 s | 0.0246 s | 0.029 | 0.016 | 0.048 |
| random1024-1000.refine53 | 0.0316 s | 10.1 s | n/a | 76.6 s | 2.16 s | 0.0589 s | 0.0586 s | 0.015 | 0.0058 | 0.027 |
| random1024-2000.refine53 | 0.188 s | 94.6 s | n/a | timeout | 34.5 s | 0.5 s | 0.46 s | 0.0055 | 0.0053 | 0.013 |
| random64-100.refine53 | 2.38e-05 s | 0.0455 s | n/a | 0.00588 s | 0.000407 s | n/a | n/a | 0.058 |  |  |
| random64-200.refine53 | 0.000896 s | 0.22 s | n/a | 0.612 s | 0.0151 s | 0.00187 s | 0.0014 s | 0.06 | 0.0085 | 0.093 |
| random64-500.refine53 | 0.00678 s | 1.73 s | n/a | 14.2 s | 0.156 s | 0.024 s | 0.0184 s | 0.043 | 0.014 | 0.12 |
| random64-1000.refine53 | 0.0232 s | 11.2 s | n/a | timeout | 1.5 s | 0.0507 s | 0.0493 s | 0.015 | 0.0045 | 0.033 |
| random64-2000.refine53 | 0.0503 s | 43.3 s | n/a | · | 12.3 s | 0.107 s | 0.136 s | 0.0041 | 0.0025 | 0.011 |
| repeated-4.refine53 | 0.000211 s | error | n/a | 0.00291 s | 0 s | n/a | 0.000294 s | 0.073 |  |  |
| repeated-8.refine53 | 0.000565 s | error | n/a | 0.0116 s | 0 s | n/a | 0.000855 s | 0.049 |  |  |
| repeated-16.refine53 | 0.00134 s | error | n/a | 0.0393 s | 0 s | n/a | 0.00183 s | 0.034 |  |  |
| repeated-32.refine53 | 0.00482 s | error | n/a | 0.255 s | 0.016 s | n/a | 0.00613 s | 0.3 |  | 0.38 |
| repeated-64.refine53 | 0.0238 s | error | n/a | 1.76 s | 0.0897 s | n/a | 0.0282 s | 0.27 |  | 0.31 |
| wilkinson-20.refine53 | 0.000127 s | 0.00318 s | n/a | 0.00092 s | 0 s | 0.00013 s | n/a | 0.14 | 0.041 |  |
| wilkinson-40.refine53 | 0.000559 s | 0.0283 s | n/a | 0.00449 s | 0.00294 s | 0.000562 s | n/a | 0.19 | 0.02 |  |
| wilkinson-80.refine53 | 0.00349 s | 0.598 s | n/a | 0.0305 s | 0.0273 s | 0.00353 s | n/a | 0.13 | 0.0059 |  |
| wilkinson-160.refine53 | 0.0231 s | 7.66 s | n/a | 0.213 s | 0.0944 s | 0.0232 s | n/a | 0.24 | 0.003 |  |
| wilkinson-320.refine53 | 0.171 s | timeout | n/a | 2.21 s | 0.458 s | n/a | n/a | 0.37 |  |  |
| wilkinson-640.refine53 | 1.5 s | · | n/a | 22.9 s | 2.93 s | n/a | n/a | 0.51 |  |  |

Achieved widths (refine53): widest / median as log2 of width over 2^-53, longest endpoint in bits.

| Case | Poly Mojo | FLINT | SymPy | msolve | Poly Mojo at FLINT width | Poly Mojo at msolve width |
|---|---|---|---|---|---|---|
| chebyshev-100.refine53 | -1.1 / -2.0, 125 b | -80.8 / -108.2, 562 b | -0.1 / -1.1, 61 b | -53.0 / -80.0, 281 b | -81.0 / -82.0, 287 b | -54.1 / -55.0, 231 b |
| chebyshev-200.refine53 | -1.0 / -2.0, 131 b | -210.2 / -269.4, 1075 b | -0.0 / -1.1, 67 b | -55.0 / -81.0, 297 b | -211.0 / -213.0, 551 b | -56.0 / -57.0, 241 b |
| chebyshev-400.refine53 | -1.0 / -2.0, 133 b | -465.4 / -587.5, 2100 b |  | -156.0 / -211.0, 567 b | -466.0 / -467.0, 1065 b | -157.0 / -158.0, 445 b |
| chebyshev-700.refine53 | -1.0 / -2.0, 133 b |  |  | -336.0 / -342.0, 1059 b |  | -337.0 / -338.0, 805 b |
| chebyshev-1000.refine53 | -1.0 / -2.0, 135 b |  |  |  |  |  |
| mignotte-10.refine53 | -0.1 / -0.9, 283 b | -170.3 / -170.3, 561 b | -1.6 / -1.8, 62 b | -68.0 / -76.0, 269 b | -172.1 / -172.1, 627 b | -69.1 / -69.1, 421 b |
| mignotte-20.refine53 | -1.4 / -1.6, 287 b | -392.1 / -392.1, 1070 b | -2.9 / -3.0, 126 b | -53.0 / -57.0, 274 b | -394.5 / -394.5, 1091 b | -54.5 / -54.5, 411 b |
| mignotte-50.refine53 | -1.1 / -1.4, 581 b | -804.7 / -805.0, 2096 b | -0.0 / -1.6, 324 b | -79.0 / -79.0, 339 b | -805.7 / -805.9, 4159 b | -80.1 / -80.4, 759 b |
| mignotte-100.refine53 | -1.1 / -1.4, 1157 b | -1661.3 / -1661.5, 4142 b | -1.1 / -3.7, 658 b | -81.0 / -83.0, 673 b | -1662.3 / -1662.7, 8751 b | -82.1 / -82.4, 1339 b |
| mignotte-200.refine53 | -1.1 / -1.2, 2327 b | -3376.7 / -3376.8, 8236 b | -0.5 / -1.6, 1322 b | -83.0 / -87.0, 1333 b | -3377.1 / -3377.7, 18017 b | -83.2 / -84.1, 2515 b |
| mignotte-400.refine53 | -1.0 / -1.4, 4661 b |  | -1.8 / -2.2, 2650 b |  |  |  |
| random1-100.refine53 | -1.9 / -2.0, 126 b | -65.2 / -65.8, 301 b | -0.7 / -0.7, 60 b | -66.0 / -83.0, 278 b | -67.9 / -68.0, 258 b | -67.9 / -68.0, 258 b |
| random1-200.refine53 | -1.4 / -1.4, 120 b | -192.7 / -192.7, 554 b | -0.1 / -0.1, 56 b | -83.0 / -83.0, 278 b | -194.2 / -194.2, 552 b | -84.2 / -84.2, 332 b |
| random1-500.refine53 | -1.1 / -1.7, 151 b | -447.5 / -447.7, 1066 b | -0.9 / -3.3, 60 b | -334.0 / -347.0, 809 b | -449.1 / -449.3, 1131 b | -335.1 / -335.3, 903 b |
| random1-1000.refine53 | -1.1 / -1.1, 127 b | -956.9 / -956.9, 2083 b | -3.1 / -3.1, 58 b | -723.0 / -723.0, 1686 b | -958.1 / -958.1, 2524 b | -724.1 / -724.1, 2056 b |
| random1-2000.refine53 | -1.1 / -1.4, 125 b | -954.1 / -955.9, 2081 b |  | -1013.0 / -1237.0, 3109 b | -956.1 / -956.4, 2528 b | -1014.1 / -1014.4, 2644 b |
| random1024-100.refine53 | -1.4 / -2.0, 114 b | -65.1 / -65.2, 300 b | -0.9 / -2.2, 58 b | -74.0 / -85.0, 288 b | -67.4 / -68.0, 246 b | -75.4 / -76.0, 262 b |
| random1024-200.refine53 | -1.1 / -1.2, 126 b | -191.4 / -191.6, 554 b | -1.2 / -1.8, 63 b | -78.0 / -82.0, 287 b | -193.1 / -193.2, 510 b | -79.1 / -79.2, 282 b |
| random1024-500.refine53 | -1.1 / -1.2, 125 b | -446.8 / -446.9, 1064 b | -2.2 / -2.3, 64 b | -243.0 / -343.0, 806 b | -448.1 / -448.2, 1072 b | -244.1 / -244.2, 664 b |
| random1024-1000.refine53 | -1.4 / -1.5, 307 b | -957.2 / -957.6, 2087 b | -0.3 / -0.5, 57 b | -723.0 / -729.0, 1581 b | -959.4 / -959.6, 2351 b | -724.4 / -724.6, 1881 b |
| random1024-2000.refine53 | -1.0 / -1.4, 158 b | -1976.5 / -1979.8, 4130 b |  | -1488.0 / -1496.0, 3350 b | -1978.0 / -1978.4, 4325 b | -1489.0 / -1489.4, 3347 b |
| random64-100.refine53 |  |  |  |  |  |  |
| random64-200.refine53 | -1.3 / -1.3, 122 b | -192.0 / -192.0, 549 b | -0.5 / -0.5, 57 b | -82.0 / -82.0, 282 b | -193.3 / -193.3, 557 b | -83.3 / -83.3, 337 b |
| random64-500.refine53 | -1.5 / -1.6, 124 b | -442.4 / -446.2, 1064 b | -0.1 / -0.6, 66 b | -335.0 / -340.0, 1046 b | -444.0 / -444.3, 1112 b | -336.0 / -336.3, 896 b |
| random64-1000.refine53 | -1.1 / -1.4, 304 b | -955.4 / -958.2, 2098 b |  | -716.0 / -723.0, 1584 b | -957.0 / -957.1, 2216 b | -717.0 / -717.1, 1736 b |
| random64-2000.refine53 | -1.1 / -1.8, 138 b | -952.8 / -954.4, 2081 b |  | -1489.0 / -1494.0, 3108 b | -954.1 / -954.4, 2134 b | -1490.4 / -1490.4, 4125 b |
| repeated-4.refine53 | -1.0 / -2.0, 123 b |  | -0.3 / -0.4, 56 b | -53.0 / -78.0, 266 b |  | -54.0 / -55.0, 229 b |
| repeated-8.refine53 | -1.1 / -1.5, 159 b |  | -0.3 / -0.4, 56 b | -53.0 / -74.0, 278 b |  | -54.1 / -54.5, 265 b |
| repeated-16.refine53 | -1.0 / -1.7, 127 b |  | -0.3 / -0.4, 57 b | -56.0 / -78.0, 278 b |  | -57.0 / -57.7, 239 b |
| repeated-32.refine53 | -1.0 / -1.8, 128 b |  | -0.1 / -1.3, 59 b | -61.0 / -79.0, 278 b |  | -62.0 / -62.8, 298 b |
| repeated-64.refine53 | -1.0 / -1.3, 131 b |  | -0.2 / -1.3, 62 b | -55.0 / -80.0, 283 b |  | -56.0 / -56.3, 284 b |
| wilkinson-20.refine53 | points, 6 b | -20.0 / -25.0, 290 b | points, 6 b | points, 6 b | points, 6 b |  |
| wilkinson-40.refine53 | points, 7 b | -96.3 / -107.3, 540 b | points, 7 b | points, 7 b | points, 7 b |  |
| wilkinson-80.refine53 | points, 8 b | -249.6 / -270.5, 1048 b | points, 8 b | points, 8 b | points, 8 b |  |
| wilkinson-160.refine53 | points, 9 b | -44.5 / -87.3, 1038 b | points, 9 b | points, 9 b | points, 9 b |  |
| wilkinson-320.refine53 | points, 10 b |  | points, 10 b | points, 10 b |  |  |
| wilkinson-640.refine53 | points, 11 b |  | points, 11 b | points, 11 b |  |  |

## Refinement of every root to width 2^-1000

| Case | Poly Mojo | FLINT | FLINT (Sturm count) | SymPy | msolve | Poly Mojo at FLINT width | Poly Mojo at msolve width | Poly Mojo / best | FLINT: matched / peer | msolve: matched / peer |
|---|---|---|---|---|---|---|---|---|---|---|
| chebyshev-100.refine1000 | 0.0733 s | 0.13 s | n/a | timeout | 1.15 s | 0.123 s | 0.126 s | 0.56 | 0.94 | 0.11 |
| chebyshev-200.refine1000 | 0.281 s | 1.03 s | n/a | · | 9.38 s | 0.454 s | 0.474 s | 0.27 | 0.44 | 0.051 |
| chebyshev-400.refine1000 | 1.14 s | 12.6 s | n/a | · | 72.8 s | 1.89 s | 1.98 s | 0.09 | 0.15 | 0.027 |
| chebyshev-700.refine1000 | 4.46 s | timeout | n/a | · | timeout | n/a | n/a |  |  |  |
| chebyshev-1000.refine1000 | 9.82 s | · | n/a | · | · | n/a | n/a |  |  |  |
| mignotte-10.refine1000 | 0.00154 s | 0.00334 s | n/a | 1.3 s | 0 s | 0.00224 s | 0.00261 s | 0.46 | 0.67 |  |
| mignotte-20.refine1000 | 0.00276 s | 0.0118 s | n/a | timeout | 0 s | 0.00351 s | 0.00379 s | 0.23 | 0.3 |  |
| mignotte-50.refine1000 | 0.219 s | 0.165 s | n/a | · | 0.0353 s | 0.91 s | 1.16 s | 6.2 | 5.5 | 33 |
| mignotte-100.refine1000 | 0.363 s | 2.47 s | n/a | · | 0.529 s | 1.29 s | 1.97 s | 0.69 | 0.52 | 3.7 |
| mignotte-200.refine1000 | 0.492 s | 50.6 s | n/a | · | 14.7 s | 13.9 s | 3.81 s | 0.033 | 0.27 | 0.26 |
| mignotte-400.refine1000 | 0.131 s | timeout | n/a | · | timeout | n/a | n/a |  |  |  |
| random1-100.refine1000 | 0.00261 s | 0.188 s | n/a | 30.2 s | 0.0312 s | 0.00446 s | 0.00492 s | 0.084 | 0.024 | 0.16 |
| random1-200.refine1000 | 0.00264 s | 0.748 s | n/a | 101 s | 0.0802 s | 0.00438 s | 0.00437 s | 0.033 | 0.0059 | 0.055 |
| random1-500.refine1000 | 0.0114 s | 4.22 s | n/a | timeout | 1 s | 0.0191 s | 0.0195 s | 0.011 | 0.0045 | 0.019 |
| random1-1000.refine1000 | 0.0177 s | 17.2 s | n/a | · | 3.16 s | 0.0272 s | 0.0346 s | 0.0056 | 0.0016 | 0.011 |
| random1-2000.refine1000 | 0.0894 s | 85.5 s | n/a | · | 22.4 s | 0.141 s | 0.171 s | 0.004 | 0.0016 | 0.0076 |
| random1024-100.refine1000 | 0.00771 s | 0.234 s | n/a | 53.4 s | 0.147 s | 0.0124 s | 0.0113 s | 0.053 | 0.053 | 0.077 |
| random1024-200.refine1000 | 0.00997 s | 1 s | n/a | timeout | 0.331 s | 0.0147 s | 0.0142 s | 0.03 | 0.015 | 0.043 |
| random1024-500.refine1000 | 0.0302 s | 5.43 s | n/a | · | 1.33 s | 0.0421 s | 0.0452 s | 0.023 | 0.0078 | 0.034 |
| random1024-1000.refine1000 | 0.0589 s | 21.5 s | n/a | · | 6.14 s | 0.0713 s | 0.0939 s | 0.0096 | 0.0033 | 0.015 |
| random1024-2000.refine1000 | 0.335 s | 94.7 s | n/a | · | 85.1 s | 0.47 s | 0.621 s | 0.0039 | 0.005 | 0.0073 |
| random64-100.refine1000 | 2.38e-05 s | 0.209 s | n/a | 0.00581 s | 0 s | n/a | n/a | 0.0041 |  |  |
| random64-200.refine1000 | 0.00319 s | 0.801 s | n/a | 106 s | 0.0911 s | 0.00504 s | 0.00496 s | 0.035 | 0.0063 | 0.054 |
| random64-500.refine1000 | 0.024 s | 4.7 s | n/a | timeout | 1.37 s | 0.034 s | 0.0327 s | 0.018 | 0.0072 | 0.024 |
| random64-1000.refine1000 | 0.0518 s | 18.9 s | n/a | · | 9.14 s | 0.0829 s | 0.104 s | 0.0057 | 0.0044 | 0.011 |
| random64-2000.refine1000 | 0.108 s | 76.4 s | n/a | · | 35.9 s | 0.157 s | 0.245 s | 0.003 | 0.0021 | 0.0068 |
| repeated-4.refine1000 | 0.00116 s | error | n/a | 0.0731 s | 0 s | n/a | 0.00272 s | 0.016 |  |  |
| repeated-8.refine1000 | 0.00292 s | error | n/a | 0.348 s | 0.0028 s | n/a | 0.00497 s | 1 |  | 1.8 |
| repeated-16.refine1000 | 0.00628 s | error | n/a | 1.7 s | 0.0193 s | n/a | 0.0106 s | 0.33 |  | 0.55 |
| repeated-32.refine1000 | 0.0187 s | error | n/a | 15.9 s | 0.0821 s | n/a | 0.0312 s | 0.23 |  | 0.38 |
| repeated-64.refine1000 | 0.0714 s | error | n/a | timeout | 0.695 s | n/a | 0.116 s | 0.1 |  | 0.17 |
| wilkinson-20.refine1000 | 0.000135 s | 0.0119 s | n/a | 0.000985 s | 0 s | 0.000137 s | n/a | 0.14 | 0.011 |  |
| wilkinson-40.refine1000 | 0.00057 s | 0.0571 s | n/a | 0.00445 s | 0 s | 0.000577 s | n/a | 0.13 | 0.01 |  |
| wilkinson-80.refine1000 | 0.00351 s | 0.697 s | n/a | 0.0292 s | 0.0277 s | 0.00341 s | n/a | 0.13 | 0.0049 |  |
| wilkinson-160.refine1000 | 0.0225 s | 7.85 s | n/a | 0.214 s | 0.0939 s | 0.0228 s | n/a | 0.24 | 0.0029 |  |
| wilkinson-320.refine1000 | 0.182 s | timeout | n/a | 2.09 s | 0.425 s | n/a | n/a | 0.43 |  |  |
| wilkinson-640.refine1000 | 1.46 s | · | n/a | 22 s | 2.95 s | n/a | n/a | 0.49 |  |  |

Achieved widths (refine1000): widest / median as log2 of width over 2^-1000, longest endpoint in bits.

| Case | Poly Mojo | FLINT | SymPy | msolve | Poly Mojo at FLINT width | Poly Mojo at msolve width |
|---|---|---|---|---|---|---|
| chebyshev-100.refine1000 | -1.1 / -2.0, 2019 b | -923.8 / -950.6, 4150 b |  | -1051.0 / -1055.0, 4619 b | -925.1 / -926.0, 3867 b | -1052.1 / -1053.0, 4121 b |
| chebyshev-200.refine1000 | -1.0 / -2.0, 2025 b | -796.7 / -855.1, 4149 b |  | -1052.0 / -1054.0, 4620 b | -798.0 / -799.0, 3619 b | -1053.0 / -1054.0, 4129 b |
| chebyshev-400.refine1000 | -1.0 / -2.0, 2027 b | -542.4 / -664.5, 4148 b |  | -1053.0 / -1056.0, 4245 b | -544.0 / -545.0, 3113 b | -1054.0 / -1055.0, 4133 b |
| chebyshev-700.refine1000 | -1.0 / -2.0, 2027 b |  |  |  |  |  |
| chebyshev-1000.refine1000 | -1.0 / -2.0, 2029 b |  |  |  |  |  |
| mignotte-10.refine1000 | -1.1 / -1.1, 2179 b | -1014.8 / -1015.3, 4147 b | -0.2 / -1.6, 1010 b | -1049.0 / -1089.0, 4229 b | -1016.1 / -1016.1, 4209 b | -1050.1 / -1050.1, 4277 b |
| mignotte-20.refine1000 | -1.5 / -1.5, 2199 b | -982.3 / -982.8, 4144 b |  | -1056.0 / -1064.0, 4495 b | -984.5 / -984.5, 4165 b | -1057.5 / -1057.5, 4311 b |
| mignotte-50.refine1000 | -0.4 / -1.0, 4791 b | -881.4 / -881.9, 4142 b |  | -1084.0 / -1180.0, 4443 b | -882.3 / -882.6, 8667 b | -1084.2 / -1085.1, 9577 b |
| mignotte-100.refine1000 | -1.1 / -1.3, 5615 b | -714.3 / -714.5, 4142 b |  | -1054.0 / -1056.0, 4769 b | -715.3 / -715.7, 8751 b | -1054.5 / -1055.0, 10245 b |
| mignotte-200.refine1000 | -1.2 / -1.2, 6849 b | -2429.7 / -2429.8, 8236 b |  | -1056.0 / -1060.0, 4417 b | -2430.1 / -2430.7, 18017 b | -1056.2 / -1057.2, 11971 b |
| mignotte-400.refine1000 | -1.0 / -1.1, 6577 b |  |  |  |  |  |
| random1-100.refine1000 | -1.9 / -2.0, 2020 b | -1036.5 / -1037.3, 4138 b | -0.1 / -0.3, 1006 b | -1056.0 / -1057.0, 4144 b | -1038.9 / -1039.0, 4094 b | -1056.0 / -1057.9, 4132 b |
| random1-200.refine1000 | -1.2 / -1.2, 2060 b | -1036.6 / -1036.6, 4136 b | -0.1 / -0.1, 1003 b | -1056.0 / -1056.0, 4118 b | -1038.2 / -1038.2, 4134 b | -1057.2 / -1057.2, 4172 b |
| random1-500.refine1000 | -1.1 / -1.3, 2129 b | -1034.7 / -1035.5, 4131 b |  | -1051.0 / -1064.0, 4137 b | -1036.1 / -1036.3, 4199 b | -1052.1 / -1052.3, 4231 b |
| random1-1000.refine1000 | -1.1 / -1.1, 2504 b | -1034.0 / -1034.0, 4131 b |  | -1122.0 / -1122.0, 6161 b | -1035.1 / -1035.1, 4572 b | -1123.1 / -1123.1, 4748 b |
| random1-2000.refine1000 | -1.1 / -1.4, 2512 b | -1031.1 / -1033.2, 4131 b |  | -1066.0 / -1090.0, 6170 b | -1033.1 / -1033.4, 4576 b | -1067.1 / -1067.4, 4644 b |
| random1024-100.refine1000 | -1.1 / -2.0, 2259 b | -1036.8 / -1037.8, 4147 b | -0.6 / -1.2, 1007 b | -1047.0 / -1058.0, 4128 b | -1038.1 / -1039.0, 4333 b | -1048.1 / -1049.0, 4353 b |
| random1024-200.refine1000 | -1.1 / -1.2, 2020 b | -1035.7 / -1036.7, 4134 b |  | -1051.0 / -1055.0, 4368 b | -1037.1 / -1037.2, 4092 b | -1052.1 / -1052.2, 4122 b |
| random1024-500.refine1000 | -1.1 / -1.2, 2072 b | -1033.8 / -1034.5, 4134 b |  | -1060.0 / -1062.0, 4177 b | -1035.1 / -1035.2, 4140 b | -1061.1 / -1061.2, 4192 b |
| random1024-1000.refine1000 | -1.4 / -1.6, 2329 b | -1032.8 / -1032.9, 4133 b |  | -1062.0 / -1062.0, 6161 b | -1034.4 / -1034.6, 4395 b | -1063.4 / -1063.6, 4453 b |
| random1024-2000.refine1000 | -1.0 / -1.4, 2265 b | -1029.5 / -1032.8, 4130 b |  | -1089.0 / -2082.0, 6422 b | -1031.0 / -1031.4, 4325 b | -1090.0 / -1090.4, 4443 b |
| random64-100.refine1000 |  |  |  |  |  |  |
| random64-200.refine1000 | -1.3 / -1.3, 2067 b | -1034.7 / -1034.7, 4136 b | -0.7 / -0.7, 1004 b | -1055.0 / -1055.0, 4122 b | -1036.3 / -1036.3, 4137 b | -1056.3 / -1056.3, 4177 b |
| random64-500.refine1000 | -1.0 / -1.3, 2120 b | -1028.1 / -1031.9, 4135 b |  | -1052.0 / -1059.0, 6166 b | -1030.0 / -1030.3, 4178 b | -1053.0 / -1053.3, 4224 b |
| random64-1000.refine1000 | -1.0 / -1.1, 2198 b | -1032.3 / -1034.9, 4140 b |  | -1063.0 / -1071.0, 6161 b | -1034.0 / -1034.1, 5077 b | -1063.1 / -1064.1, 5137 b |
| random64-2000.refine1000 | -1.1 / -1.4, 2122 b | -1029.8 / -1031.3, 4131 b |  | -2078.0 / -2083.0, 6180 b | -1031.4 / -1031.4, 5101 b | -2079.4 / -2079.4, 7197 b |
| repeated-4.refine1000 | -1.0 / -2.0, 2017 b |  | -0.3 / -0.6, 1004 b | -1051.0 / -1053.0, 4118 b |  | -1052.0 / -1053.0, 4119 b |
| repeated-8.refine1000 | -1.1 / -1.5, 2053 b |  | -0.3 / -0.5, 1008 b | -1047.0 / -1052.0, 4135 b |  | -1048.0 / -1048.2, 4147 b |
| repeated-16.refine1000 | -1.0 / -1.7, 2271 b |  | -0.1 / -1.1, 1006 b | -1048.0 / -1052.0, 4226 b |  | -1049.0 / -1049.7, 4367 b |
| repeated-32.refine1000 | -1.0 / -1.8, 2070 b |  | -0.1 / -0.8, 1006 b | -1049.0 / -1052.0, 4619 b |  | -1050.0 / -1050.8, 4168 b |
| repeated-64.refine1000 | -1.0 / -1.3, 2068 b |  |  | -1050.0 / -1055.0, 4374 b |  | -1051.0 / -1051.3, 4168 b |
| wilkinson-20.refine1000 | points, 6 b | -991.6 / -996.5, 4128 b | points, 6 b | points, 6 b | points, 6 b |  |
| wilkinson-40.refine1000 | points, 7 b | -939.3 / -949.9, 4118 b | points, 7 b | points, 7 b | points, 7 b |  |
| wilkinson-80.refine1000 | points, 8 b | -836.0 / -856.9, 4116 b | points, 8 b | points, 8 b | points, 8 b |  |
| wilkinson-160.refine1000 | points, 9 b | -631.1 / -674.4, 4106 b | points, 9 b | points, 9 b | points, 9 b |  |
| wilkinson-320.refine1000 | points, 10 b |  | points, 10 b | points, 10 b |  |  |
| wilkinson-640.refine1000 | points, 11 b |  | points, 11 b | points, 11 b |  |  |

## Exact real roots (Poly Mojo `real_roots`, SymPy `real_roots`)

| Case | Poly Mojo | FLINT | FLINT (Sturm count) | SymPy | msolve | Poly Mojo / best |
|---|---|---|---|---|---|---|
| chebyshev-100.exact | 0.0794 s | n/a | n/a | error | n/a |  |
| chebyshev-200.exact | 0.707 s | n/a | n/a | error | n/a |  |
| chebyshev-400.exact | 6.74 s | n/a | n/a | error | n/a |  |
| chebyshev-700.exact | timeout | n/a | n/a | timeout | n/a |  |
| chebyshev-1000.exact | · | n/a | n/a | · | n/a |  |
| mignotte-10.exact | 0.000239 s | n/a | n/a | error | n/a |  |
| mignotte-20.exact | 0.000391 s | n/a | n/a | timeout | n/a |  |
| mignotte-50.exact | 0.00292 s | n/a | n/a | · | n/a |  |
| mignotte-100.exact | 0.0364 s | n/a | n/a | · | n/a |  |
| mignotte-200.exact | 0.802 s | n/a | n/a | · | n/a |  |
| mignotte-400.exact | 20.5 s | n/a | n/a | · | n/a |  |
| random1-100.exact | 0.00653 s | n/a | n/a | error | n/a |  |
| random1-200.exact | 0.0281 s | n/a | n/a | error | n/a |  |
| random1-500.exact | 0.226 s | n/a | n/a | error | n/a |  |
| random1-1000.exact | 0.733 s | n/a | n/a | timeout | n/a |  |
| random1-2000.exact | 10.1 s | n/a | n/a | · | n/a |  |
| random1024-100.exact | 0.0116 s | n/a | n/a | error | n/a |  |
| random1024-200.exact | 0.0824 s | n/a | n/a | error | n/a |  |
| random1024-500.exact | 0.347 s | n/a | n/a | timeout | n/a |  |
| random1024-1000.exact | 5.59 s | n/a | n/a | · | n/a |  |
| random1024-2000.exact | 26.2 s | n/a | n/a | · | n/a |  |
| random64-100.exact | 0.013 s | n/a | n/a | 0.34 s | n/a | 0.038 |
| random64-200.exact | 0.0485 s | n/a | n/a | error | n/a |  |
| random64-500.exact | 0.538 s | n/a | n/a | error | n/a |  |
| random64-1000.exact | 5.04 s | n/a | n/a | timeout | n/a |  |
| random64-2000.exact | 9.25 s | n/a | n/a | · | n/a |  |
| repeated-4.exact | 0.000517 s | n/a | n/a | error | n/a |  |
| repeated-8.exact | 0.00168 s | n/a | n/a | error | n/a |  |
| repeated-16.exact | 0.0111 s | n/a | n/a | error | n/a |  |
| repeated-32.exact | 0.0624 s | n/a | n/a | error | n/a |  |
| repeated-64.exact | 0.516 s | n/a | n/a | error | n/a |  |
| wilkinson-20.exact | 0.00128 s | n/a | n/a | 0.0272 s | n/a | 0.047 |
| wilkinson-40.exact | 0.0056 s | n/a | n/a | 0.0839 s | n/a | 0.067 |
| wilkinson-80.exact | 0.0243 s | n/a | n/a | 0.31 s | n/a | 0.078 |
| wilkinson-160.exact | 0.116 s | n/a | n/a | 1.43 s | n/a | 0.081 |
| wilkinson-320.exact | 0.605 s | n/a | n/a | 9.67 s | n/a | 0.063 |
| wilkinson-640.exact | 4.92 s | n/a | n/a | 85.6 s | n/a | 0.058 |

## Distinct real roots (Poly Mojo `sturm_sequence` + `count_roots` on a power-of-two root bound)

| Case | Poly Mojo | FLINT | FLINT (Sturm count) | SymPy | msolve | Poly Mojo / best |
|---|---|---|---|---|---|---|
| chebyshev-100.count | 0.0022 s | 0.033 s | 0.0324 s | 0.00555 s | n/a | 0.4 |
| chebyshev-200.count | 0.0102 s | 1.05 s | 1.15 s | 0.022 s | n/a | 0.46 |
| chebyshev-400.count | 0.046 s | 30.5 s | 30.5 s | 0.0989 s | n/a | 0.47 |
| chebyshev-700.count | 0.166 s | timeout | timeout | 0.363 s | n/a | 0.46 |
| chebyshev-1000.count | 0.375 s | · | · | 0.896 s | n/a | 0.42 |
| mignotte-10.count | 2.55e-05 s | 1.15e-06 s | 1.16e-06 s | 0.000225 s | n/a | 22 |
| mignotte-20.count | 4.61e-05 s | 2.66e-06 s | 2.65e-06 s | 0.000337 s | n/a | 17 |
| mignotte-50.count | 0.000162 s | 9.18e-06 s | 9.19e-06 s | 0.000862 s | n/a | 18 |
| mignotte-100.count | 0.000659 s | 2.81e-05 s | 2.8e-05 s | 0.00214 s | n/a | 24 |
| mignotte-200.count | 0.00192 s | 9.21e-05 s | 9.24e-05 s | 0.0069 s | n/a | 21 |
| mignotte-400.count | 0.00764 s | 0.000345 s | 0.000342 s | 0.0254 s | n/a | 22 |
| random1-100.count | 0.16 s | 0.00155 s | 0.00158 s | 1.47 s | n/a | 1e+02 |
| random1-200.count | 3.02 s | 0.0225 s | 0.0227 s | 62.8 s | n/a | 1.3e+02 |
| random1-500.count | timeout | 0.782 s | 0.8 s | timeout | n/a |  |
| random1-1000.count | · | 11.7 s | 11.5 s | · | n/a |  |
| random1-2000.count | · | timeout | timeout | · | n/a |  |
| random1024-100.count | timeout | 5.05 s | 4.91 s | timeout | n/a |  |
| random1024-200.count | · | 46.7 s | 46.5 s | · | n/a |  |
| random1024-500.count | · | timeout | timeout | · | n/a |  |
| random1024-1000.count | · | · | · | · | n/a |  |
| random1024-2000.count | · | · | · | · | n/a |  |
| random64-100.count | 10.4 s | 0.0901 s | 0.0892 s | 91.8 s | n/a | 1.2e+02 |
| random64-200.count | timeout | 1.12 s | 1.12 s | timeout | n/a |  |
| random64-500.count | · | 29 s | 28.1 s | · | n/a |  |
| random64-1000.count | · | timeout | timeout | · | n/a |  |
| random64-2000.count | · | · | · | · | n/a |  |
| repeated-4.count | 0.000926 s | 2.18e-05 s | 2.21e-05 s | 0.000463 s | n/a | 42 |
| repeated-8.count | 0.0106 s | 0.00034 s | 0.000348 s | 0.00332 s | n/a | 31 |
| repeated-16.count | 0.0904 s | 0.0131 s | 0.0131 s | 0.108 s | n/a | 6.9 |
| repeated-32.count | 4.83 s | 0.511 s | 0.548 s | 10.5 s | n/a | 9.5 |
| repeated-64.count | timeout | 17.7 s | 18.3 s | timeout | n/a |  |
| wilkinson-20.count | 0.000202 s | 3.66e-05 s | 3.66e-05 s | 0.00086 s | n/a | 5.5 |
| wilkinson-40.count | 0.00179 s | 0.00188 s | 0.00185 s | 0.00264 s | n/a | 0.97 |
| wilkinson-80.count | 0.0124 s | 0.0879 s | 0.0874 s | 0.0115 s | n/a | 1.1 |
| wilkinson-160.count | 0.0878 s | 3.73 s | 3.81 s | 0.0726 s | n/a | 1.2 |
| wilkinson-320.count | 0.662 s | 115 s | 114 s | 0.632 s | n/a | 1 |
| wilkinson-640.count | 5.05 s | timeout | timeout | 6.31 s | n/a | 0.8 |

