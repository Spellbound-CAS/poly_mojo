# Arithmetic results

Host: Intel(R) Xeon(R) Processor @ 2.80GHz, 4 vCPUs, 16 GiB RAM. Runs are pinned to one CPU, capped at 3G and 600 s.

Versions: Poly Mojo e88c776db002; FLINT python-flint 0.9.0, FLINT 3.6.0; SymPy 1.14.0 (gmpy ground types); Singular 4.4.1.

## Summary

Geometric mean of competitor time ÷ Poly Mojo time over cases both solved with a confirmed answer and where the slower took at least 10 ms (above 1 means Poly Mojo was faster).

| Ladder | FLINT | SymPy | Singular |
|---|---|---|---|
| dense-div.GF251 | 1.54 (*n* = 1) | 4.98e+03 (*n* = 3) | 2 (*n* = 1) |
| dense-div.GF4611686018427387847 | 1.94 (*n* = 2) | 1.76e+03 (*n* = 3) | — |
| dense-div.ZZ-1024b | 2.56 (*n* = 3) | 104 (*n* = 2) | 38.1 (*n* = 2) |
| dense-div.ZZ-64b | 4.02 (*n* = 2) | 797 (*n* = 2) | 255 (*n* = 3) |
| dense-div.ZZ-8b | 2.04 (*n* = 1) | 2.32e+03 (*n* = 2) | 428 (*n* = 2) |
| dense-mul.GF251 | 1.84 (*n* = 1) | 8.36e+03 (*n* = 3) | 3.85 (*n* = 1) |
| dense-mul.GF4611686018427387847 | 2.45 (*n* = 1) | 4.73e+03 (*n* = 3) | — |
| dense-mul.ZZ-1024b | 1.6 (*n* = 2) | 230 (*n* = 2) | 1.57 (*n* = 2) |
| dense-mul.ZZ-64b | 1.67 (*n* = 1) | 2.31e+03 (*n* = 2) | 2.87 (*n* = 2) |
| dense-mul.ZZ-8b | 1.57 (*n* = 1) | 7.15e+03 (*n* = 2) | 2.87 (*n* = 1) |
| dense-sqr.GF251 | — | — | 5.52 (*n* = 1) |
| dense-sqr.GF4611686018427387847 | 2.17 (*n* = 1) | — | — |
| dense-sqr.ZZ-1024b | 1.35 (*n* = 2) | — | 2.22 (*n* = 2) |
| dense-sqr.ZZ-64b | 1.36 (*n* = 1) | — | 3.11 (*n* = 2) |
| dense-sqr.ZZ-8b | 1.15 (*n* = 1) | — | 3.25 (*n* = 1) |

## Answer agreement

Every answer that two or more libraries computed agrees.

## Per-ladder results

Each cell gives the median time per call and the incremental peak memory; `·` marks an instance skipped after a smaller one on the same ladder hit a limit; `n/a` marks an unsupported field.

### dense-div.GF251

| n | Poly Mojo | FLINT | SymPy | Singular |
|---|---|---|---|---|
| 100 | < 1 ms<br>< 1 MiB | < 1 ms<br>< 1 MiB | 15.5 ms<br>< 1 MiB | < 1 ms<br>< 1 MiB |
| 1000 | < 1 ms<br>< 1 MiB | < 1 ms<br>< 1 MiB | 1.51 s<br>< 1 MiB | < 1 ms<br>< 1 MiB |
| 10000 | 2.76 ms<br>2.01 MiB | 3.21 ms<br>1.07 MiB | 255 s<br>5.77 MiB | 5.09 ms<br>1.8 MiB |
| 100000 | 35.9 ms<br>9.46 MiB | 55.3 ms<br>6.01 MiB | timeout | 71.8 ms<br>12 MiB |

### dense-div.GF4611686018427387847

| n | Poly Mojo | FLINT | SymPy | Singular |
|---|---|---|---|---|
| 100 | < 1 ms<br>< 1 MiB | < 1 ms<br>< 1 MiB | 14.7 ms<br>< 1 MiB | n/a |
| 1000 | < 1 ms<br>< 1 MiB | < 1 ms<br>< 1 MiB | 1.72 s<br>< 1 MiB | n/a |
| 10000 | 10.9 ms<br>3.54 MiB | 19 ms<br>1.86 MiB | 165 s<br>2.92 MiB | n/a |
| 100000 | 118 ms<br>13.9 MiB | 253 ms<br>26 MiB | timeout | n/a |

### dense-div.ZZ-1024b

| n | Poly Mojo | FLINT | SymPy | Singular |
|---|---|---|---|---|
| 100 | 1.95 ms<br>< 1 MiB | 1.17 ms<br>< 1 MiB | 5.4 ms<br>< 1 MiB | 2.38 ms<br>< 1 MiB |
| 1000 | 15.7 ms<br>< 1 MiB | 24.9 ms<br>2.73 MiB | 647 ms<br>< 1 MiB | 228 ms<br>1.12 MiB |
| 10000 | 241 ms<br>3.15 MiB | 636 ms<br>41.7 MiB | 63.3 s<br>7.43 MiB | 24.2 s<br>6.47 MiB |
| 100000 | 2.72 s<br>25.3 MiB | 11 s<br>384 MiB | timeout | timeout |

### dense-div.ZZ-64b

| n | Poly Mojo | FLINT | SymPy | Singular |
|---|---|---|---|---|
| 100 | < 1 ms<br>< 1 MiB | < 1 ms<br>< 1 MiB | 3.22 ms<br>< 1 MiB | < 1 ms<br>< 1 MiB |
| 1000 | 1.11 ms<br>< 1 MiB | 3.01 ms<br>< 1 MiB | 350 ms<br>< 1 MiB | 41.5 ms<br>< 1 MiB |
| 10000 | 19.6 ms<br>2.57 MiB | 68.9 ms<br>4.28 MiB | 39.5 s<br>3.25 MiB | 4.15 s<br>2.6 MiB |
| 100000 | 233 ms<br>11 MiB | 1.07 s<br>44.6 MiB | timeout | 490 s<br>22.1 MiB |

### dense-div.ZZ-8b

| n | Poly Mojo | FLINT | SymPy | Singular |
|---|---|---|---|---|
| 100 | < 1 ms<br>< 1 MiB | < 1 ms<br>< 1 MiB | 3.23 ms<br>< 1 MiB | < 1 ms<br>< 1 MiB |
| 1000 | < 1 ms<br>< 1 MiB | < 1 ms<br>< 1 MiB | 360 ms<br>< 1 MiB | 9.86 ms<br>< 1 MiB |
| 10000 | 6.23 ms<br>3.26 MiB | 8.86 ms<br>< 1 MiB | 34.7 s<br>3.13 MiB | 864 ms<br>1.46 MiB |
| 100000 | 80.7 ms<br>8.79 MiB | 165 ms<br>4.56 MiB | timeout | 107 s<br>14.4 MiB |

### dense-mul.GF251

| n | Poly Mojo | FLINT | SymPy | Singular |
|---|---|---|---|---|
| 100 | < 1 ms<br>< 1 MiB | < 1 ms<br>< 1 MiB | 13.4 ms<br>< 1 MiB | < 1 ms<br>1.11 MiB |
| 1000 | < 1 ms<br>< 1 MiB | < 1 ms<br>< 1 MiB | 1.63 s<br>< 1 MiB | < 1 ms<br>1.46 MiB |
| 10000 | < 1 ms<br>< 1 MiB | < 1 ms<br>1.02 MiB | 146 s<br>6.64 MiB | 2.63 ms<br>2.71 MiB |
| 100000 | 7.93 ms<br>10.9 MiB | 14.6 ms<br>5.46 MiB | timeout | 30.5 ms<br>16.4 MiB |

### dense-mul.GF4611686018427387847

| n | Poly Mojo | FLINT | SymPy | Singular |
|---|---|---|---|---|
| 100 | < 1 ms<br>< 1 MiB | < 1 ms<br>< 1 MiB | 13.2 ms<br>< 1 MiB | n/a |
| 1000 | < 1 ms<br>< 1 MiB | < 1 ms<br>< 1 MiB | 1.49 s<br>< 1 MiB | n/a |
| 10000 | 2.06 ms<br>1.82 MiB | 4.69 ms<br>1.7 MiB | 157 s<br>6.63 MiB | n/a |
| 100000 | 28.4 ms<br>10.8 MiB | 69.6 ms<br>13.7 MiB | timeout | n/a |

### dense-mul.ZZ-1024b

| n | Poly Mojo | FLINT | SymPy | Singular |
|---|---|---|---|---|
| 100 | < 1 ms<br>< 1 MiB | < 1 ms<br>< 1 MiB | 4.71 ms<br>< 1 MiB | < 1 ms<br>1.59 MiB |
| 1000 | 6.1 ms<br>3.98 MiB | 4.6 ms<br>2.33 MiB | 501 ms<br>1.1 MiB | 5.82 ms<br>3.61 MiB |
| 10000 | 86.5 ms<br>44.4 MiB | 145 ms<br>34.4 MiB | 55.9 s<br>11.3 MiB | 152 ms<br>39.3 MiB |
| 100000 | 1.35 s<br>450 MiB | 2.05 s<br>357 MiB | timeout | 1.88 s<br>396 MiB |

### dense-mul.ZZ-64b

| n | Poly Mojo | FLINT | SymPy | Singular |
|---|---|---|---|---|
| 100 | < 1 ms<br>< 1 MiB | < 1 ms<br>< 1 MiB | 2.5 ms<br>< 1 MiB | < 1 ms<br>1.35 MiB |
| 1000 | < 1 ms<br>< 1 MiB | < 1 ms<br>< 1 MiB | 287 ms<br>< 1 MiB | < 1 ms<br>1.66 MiB |
| 10000 | 4.58 ms<br>3.46 MiB | 7.49 ms<br>2.07 MiB | 31.9 s<br>5.76 MiB | 14.5 ms<br>5.34 MiB |
| 100000 | 57.6 ms<br>32.8 MiB | 96.4 ms<br>23.8 MiB | timeout | 149 ms<br>41.9 MiB |

### dense-mul.ZZ-8b

| n | Poly Mojo | FLINT | SymPy | Singular |
|---|---|---|---|---|
| 100 | < 1 ms<br>< 1 MiB | < 1 ms<br>< 1 MiB | 2.41 ms<br>< 1 MiB | < 1 ms<br>1.19 MiB |
| 1000 | < 1 ms<br>< 1 MiB | < 1 ms<br>< 1 MiB | 282 ms<br>< 1 MiB | < 1 ms<br>1.53 MiB |
| 10000 | 1.61 ms<br>1.51 MiB | 1.67 ms<br>1.11 MiB | 29.7 s<br>5.26 MiB | 3.4 ms<br>2.7 MiB |
| 100000 | 16.4 ms<br>12.9 MiB | 25.8 ms<br>6.46 MiB | timeout | 47.1 ms<br>14.7 MiB |

### dense-sqr.GF251

| n | Poly Mojo | FLINT | SymPy | Singular |
|---|---|---|---|---|
| 100 | < 1 ms<br>< 1 MiB | < 1 ms<br>< 1 MiB |  | < 1 ms<br>< 1 MiB |
| 1000 | < 1 ms<br>< 1 MiB | < 1 ms<br>< 1 MiB |  | < 1 ms<br>1.43 MiB |
| 10000 | < 1 ms<br>< 1 MiB | < 1 ms<br>< 1 MiB |  | 2.47 ms<br>3.01 MiB |
| 100000 | 6.08 ms<br>6.7 MiB | 9.96 ms<br>5.21 MiB |  | 33.6 ms<br>19.4 MiB |

### dense-sqr.GF4611686018427387847

| n | Poly Mojo | FLINT | SymPy | Singular |
|---|---|---|---|---|
| 100 | < 1 ms<br>< 1 MiB | < 1 ms<br>< 1 MiB |  | n/a |
| 1000 | < 1 ms<br>< 1 MiB | < 1 ms<br>< 1 MiB |  | n/a |
| 10000 | 1.52 ms<br>1.66 MiB | 3.21 ms<br>1.57 MiB |  | n/a |
| 100000 | 21.4 ms<br>5.28 MiB | 46.4 ms<br>11.6 MiB |  | n/a |

### dense-sqr.ZZ-1024b

| n | Poly Mojo | FLINT | SymPy | Singular |
|---|---|---|---|---|
| 100 | < 1 ms<br>< 1 MiB | < 1 ms<br>< 1 MiB |  | < 1 ms<br>1.61 MiB |
| 1000 | 5.07 ms<br>3.42 MiB | 2.94 ms<br>1.5 MiB |  | 5.85 ms<br>3.91 MiB |
| 10000 | 68.7 ms<br>39.2 MiB | 92.5 ms<br>19.7 MiB |  | 166 ms<br>41.2 MiB |
| 100000 | 937 ms<br>411 MiB | 1.27 s<br>203 MiB |  | 1.92 s<br>416 MiB |

### dense-sqr.ZZ-64b

| n | Poly Mojo | FLINT | SymPy | Singular |
|---|---|---|---|---|
| 100 | < 1 ms<br>< 1 MiB | < 1 ms<br>< 1 MiB |  | < 1 ms<br>1.4 MiB |
| 1000 | < 1 ms<br>< 1 MiB | < 1 ms<br>< 1 MiB |  | < 1 ms<br>1.77 MiB |
| 10000 | 4.28 ms<br>3.17 MiB | 5.24 ms<br>1.68 MiB |  | 11.7 ms<br>6.14 MiB |
| 100000 | 46.9 ms<br>27.8 MiB | 63.7 ms<br>16.5 MiB |  | 165 ms<br>49.6 MiB |

### dense-sqr.ZZ-8b

| n | Poly Mojo | FLINT | SymPy | Singular |
|---|---|---|---|---|
| 100 | < 1 ms<br>< 1 MiB | < 1 ms<br>< 1 MiB |  | < 1 ms<br>1.27 MiB |
| 1000 | < 1 ms<br>< 1 MiB | < 1 ms<br>< 1 MiB |  | < 1 ms<br>1.62 MiB |
| 10000 | 1.6 ms<br>1.65 MiB | < 1 ms<br>1.05 MiB |  | 3.51 ms<br>3.03 MiB |
| 100000 | 15.2 ms<br>7.14 MiB | 17.4 ms<br>3.95 MiB |  | 49.4 ms<br>17.8 MiB |

