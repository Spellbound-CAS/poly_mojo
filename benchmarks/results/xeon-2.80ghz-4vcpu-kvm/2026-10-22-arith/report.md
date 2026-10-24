# Arithmetic results

Host: Intel(R) Xeon(R) Processor @ 2.80GHz, 4 vCPUs, 16 GiB RAM. Runs are pinned to one CPU, capped at 3G and 600 s.

Versions: Poly Mojo d3e021c61ab6; FLINT python-flint 0.9.0, FLINT 3.6.0; SymPy 1.14.0 (gmpy ground types); Singular 4.4.1.

## Summary

Geometric mean of competitor time ÷ Poly Mojo time over cases both solved with a confirmed answer and where the slower took at least 10 ms (above 1 means Poly Mojo was faster).

| Ladder | FLINT | SymPy | Singular |
|---|---|---|---|
| dense-mul.GF4611686018427387847 | — | — | — |
| dense-mul.ZZ-1024b | 0.914 (*n* = 2) | 230 (*n* = 2) | 1.16 (*n* = 2) |
| dense-mul.ZZ-64b | 1.84 (*n* = 1) | 2.49e+03 (*n* = 2) | 2.67 (*n* = 2) |
| dense-mul.ZZ-8b | 1.49 (*n* = 1) | 7.58e+03 (*n* = 2) | 2.89 (*n* = 1) |
| dense-sqr.GF251 | 1.83 (*n* = 1) | 9.38e+03 (*n* = 3) | 5.41 (*n* = 1) |
| dense-sqr.GF4611686018427387847 | 2.05 (*n* = 1) | 5.9e+03 (*n* = 3) | — |
| dense-sqr.ZZ-1024b | 0.794 (*n* = 2) | 333 (*n* = 2) | 1.64 (*n* = 2) |
| dense-sqr.ZZ-64b | 1.43 (*n* = 1) | 2.77e+03 (*n* = 2) | 3.09 (*n* = 2) |
| dense-sqr.ZZ-8b | 1.22 (*n* = 1) | 8.87e+03 (*n* = 2) | 3.73 (*n* = 1) |

## Answer agreement

Every answer that two or more libraries computed agrees.

## Per-ladder results

Each cell gives the median time per call and the incremental peak memory; `·` marks an instance skipped after a smaller one on the same ladder hit a limit; `n/a` marks an unsupported field.

### dense-mul.GF4611686018427387847

| n | Poly Mojo | FLINT | SymPy | Singular |
|---|---|---|---|---|
| 100000 |  |  | timeout |  |

### dense-mul.ZZ-1024b

| n | Poly Mojo | FLINT | SymPy | Singular |
|---|---|---|---|---|
| 100 | < 1 ms<br>< 1 MiB | < 1 ms<br>< 1 MiB | 4.67 ms<br>< 1 MiB | < 1 ms<br>1.59 MiB |
| 1000 | 6.49 ms<br>4.15 MiB | 4.64 ms<br>2.24 MiB | 540 ms<br>1.09 MiB | 5.54 ms<br>3.66 MiB |
| 10000 | 95 ms<br>44.7 MiB | 153 ms<br>34.4 MiB | 60.5 s<br>11.3 MiB | 152 ms<br>39.3 MiB |
| 100000 | 3.75 s<br>450 MiB | 1.94 s<br>357 MiB | timeout | 3.14 s<br>396 MiB |

### dense-mul.ZZ-64b

| n | Poly Mojo | FLINT | SymPy | Singular |
|---|---|---|---|---|
| 100 | < 1 ms<br>< 1 MiB | < 1 ms<br>< 1 MiB | 2.64 ms<br>< 1 MiB | < 1 ms<br>1.36 MiB |
| 1000 | < 1 ms<br>< 1 MiB | < 1 ms<br>< 1 MiB | 307 ms<br>< 1 MiB | < 1 ms<br>1.67 MiB |
| 10000 | 4.85 ms<br>3.43 MiB | 7.64 ms<br>1.95 MiB | 35.3 s<br>5.75 MiB | 12.2 ms<br>5.31 MiB |
| 100000 | 55.4 ms<br>32.7 MiB | 102 ms<br>23.9 MiB | timeout | 157 ms<br>42 MiB |

### dense-mul.ZZ-8b

| n | Poly Mojo | FLINT | SymPy | Singular |
|---|---|---|---|---|
| 100 | < 1 ms<br>< 1 MiB | < 1 ms<br>< 1 MiB | 2.39 ms<br>< 1 MiB | < 1 ms<br>1.2 MiB |
| 1000 | < 1 ms<br>< 1 MiB | < 1 ms<br>< 1 MiB | 311 ms<br>< 1 MiB | < 1 ms<br>1.54 MiB |
| 10000 | 1.67 ms<br>1.55 MiB | 1.62 ms<br>1.14 MiB | 31.4 s<br>5.38 MiB | 3.73 ms<br>2.76 MiB |
| 100000 | 16.2 ms<br>12.8 MiB | 24 ms<br>6.56 MiB | timeout | 46.8 ms<br>14.7 MiB |

### dense-sqr.GF251

| n | Poly Mojo | FLINT | SymPy | Singular |
|---|---|---|---|---|
| 100 | < 1 ms<br>< 1 MiB | < 1 ms<br>< 1 MiB | 12.9 ms<br>< 1 MiB | < 1 ms<br>< 1 MiB |
| 1000 | < 1 ms<br>< 1 MiB | < 1 ms<br>< 1 MiB | 1.5 s<br>< 1 MiB | < 1 ms<br>1.36 MiB |
| 10000 | < 1 ms<br>< 1 MiB | < 1 ms<br>< 1 MiB | 167 s<br>6.62 MiB | 2.75 ms<br>3.03 MiB |
| 100000 | 7.01 ms<br>6.64 MiB | 12.8 ms<br>5.01 MiB | timeout | 37.9 ms<br>19.5 MiB |

### dense-sqr.GF4611686018427387847

| n | Poly Mojo | FLINT | SymPy | Singular |
|---|---|---|---|---|
| 100 | < 1 ms<br>< 1 MiB | < 1 ms<br>< 1 MiB | 14 ms<br>< 1 MiB | n/a |
| 1000 | < 1 ms<br>< 1 MiB | < 1 ms<br>< 1 MiB | 1.76 s<br>< 1 MiB | n/a |
| 10000 | 1.61 ms<br>1.05 MiB | 3.35 ms<br>1.56 MiB | 170 s<br>6.62 MiB | n/a |
| 100000 | 25 ms<br>5.29 MiB | 51.2 ms<br>11.5 MiB | timeout | n/a |

### dense-sqr.ZZ-1024b

| n | Poly Mojo | FLINT | SymPy | Singular |
|---|---|---|---|---|
| 100 | < 1 ms<br>1.04 MiB | < 1 ms<br>< 1 MiB | 4.95 ms<br>< 1 MiB | < 1 ms<br>1.67 MiB |
| 1000 | 5.25 ms<br>3.54 MiB | 3 ms<br>1.54 MiB | 576 ms<br>1.1 MiB | 6.52 ms<br>3.91 MiB |
| 10000 | 68.4 ms<br>39.3 MiB | 94.5 ms<br>19.7 MiB | 69 s<br>11.3 MiB | 158 ms<br>41.2 MiB |
| 100000 | 3.07 s<br>411 MiB | 1.4 s<br>203 MiB | timeout | 3.55 s<br>416 MiB |

### dense-sqr.ZZ-64b

| n | Poly Mojo | FLINT | SymPy | Singular |
|---|---|---|---|---|
| 100 | < 1 ms<br>< 1 MiB | < 1 ms<br>< 1 MiB | 2.49 ms<br>< 1 MiB | < 1 ms<br>1.41 MiB |
| 1000 | < 1 ms<br>< 1 MiB | < 1 ms<br>< 1 MiB | 310 ms<br>< 1 MiB | < 1 ms<br>1.76 MiB |
| 10000 | 4.21 ms<br>3.33 MiB | 5.17 ms<br>1.72 MiB | 37.8 s<br>5.75 MiB | 12 ms<br>6.13 MiB |
| 100000 | 46 ms<br>27.7 MiB | 65.6 ms<br>16.4 MiB | timeout | 154 ms<br>49.7 MiB |

### dense-sqr.ZZ-8b

| n | Poly Mojo | FLINT | SymPy | Singular |
|---|---|---|---|---|
| 100 | < 1 ms<br>< 1 MiB | < 1 ms<br>< 1 MiB | 2.61 ms<br>< 1 MiB | < 1 ms<br>1.21 MiB |
| 1000 | < 1 ms<br>< 1 MiB | < 1 ms<br>< 1 MiB | 343 ms<br>< 1 MiB | < 1 ms<br>1.62 MiB |
| 10000 | 1.7 ms<br>1.14 MiB | 1.08 ms<br>1.05 MiB | 38.9 s<br>5.39 MiB | 3.49 ms<br>3.02 MiB |
| 100000 | 14.3 ms<br>6.98 MiB | 17.4 ms<br>3.94 MiB | timeout | 53.1 ms<br>17.7 MiB |

