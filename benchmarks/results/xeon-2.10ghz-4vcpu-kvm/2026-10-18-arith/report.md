# Arithmetic results

Host: Intel(R) Xeon(R) Processor @ 2.10GHz, 4 vCPUs, 16 GiB RAM. Runs are pinned to one CPU, capped at 3G and 600 s.

Versions: Poly Mojo fdba30fb4669; FLINT python-flint 0.9.0, FLINT 3.6.0; SymPy 1.14.0 (gmpy ground types); Singular 4.4.1.

## Summary

Geometric mean of competitor time ÷ Poly Mojo time over cases both solved with a confirmed answer and where the slower took at least 10 ms (above 1 means Poly Mojo was faster).

| Ladder | FLINT | SymPy | Singular |
|---|---|---|---|
| dense-div.GF251 | 1.85 (*n* = 1) | 3.45e+03 (*n* = 3) | 2.53 (*n* = 1) |
| dense-div.GF4611686018427387847 | 2.38 (*n* = 2) | 1.56e+03 (*n* = 3) | — |
| dense-div.ZZ-1024b | 1.96 (*n* = 3) | 71.4 (*n* = 2) | 29.7 (*n* = 2) |
| dense-div.ZZ-64b | 2.78 (*n* = 2) | 403 (*n* = 2) | 170 (*n* = 3) |
| dense-div.ZZ-8b | 1.13 (*n* = 2) | 1.17e+03 (*n* = 2) | 226 (*n* = 2) |
| dense-mul.GF251 | 2.51 (*n* = 1) | 7.51e+03 (*n* = 3) | 4.49 (*n* = 1) |
| dense-mul.GF4611686018427387847 | 2.93 (*n* = 1) | 4.32e+03 (*n* = 3) | — |
| dense-mul.ZZ-1024b | 1.74 (*n* = 2) | 202 (*n* = 2) | 1.61 (*n* = 2) |
| dense-mul.ZZ-64b | 1.46 (*n* = 1) | 1.59e+03 (*n* = 2) | 2.31 (*n* = 2) |
| dense-mul.ZZ-8b | 1.76 (*n* = 1) | 7.38e+03 (*n* = 2) | 3.45 (*n* = 1) |
| dense-sqr.GF251 | 2.98 (*n* = 1) | 7.92e+03 (*n* = 3) | 6.62 (*n* = 1) |
| dense-sqr.GF4611686018427387847 | 3.07 (*n* = 1) | 5.46e+03 (*n* = 3) | — |
| dense-sqr.ZZ-1024b | 1.37 (*n* = 2) | 270 (*n* = 2) | 2.61 (*n* = 2) |

## Answer agreement

Every answer that two or more libraries computed agrees.

## Per-ladder results

Each cell gives the median time per call and the incremental peak memory; `·` marks an instance skipped after a smaller one on the same ladder hit a limit; `n/a` marks an unsupported field.

### dense-div.GF251

| n | Poly Mojo | FLINT | SymPy | Singular |
|---|---|---|---|---|
| 100 | < 1 ms<br>< 1 MiB | < 1 ms<br>< 1 MiB | 12.1 ms<br>< 1 MiB | < 1 ms<br>< 1 MiB |
| 1000 | < 1 ms<br>< 1 MiB | < 1 ms<br>< 1 MiB | 1.39 s<br>< 1 MiB | < 1 ms<br>1 MiB |
| 10000 | 3.32 ms<br>2.2 MiB | 2.93 ms<br>1.05 MiB | 153 s<br>5.82 MiB | 5.12 ms<br>2.02 MiB |
| 100000 | 28.3 ms<br>11.2 MiB | 52.5 ms<br>7.3 MiB | timeout | 71.5 ms<br>12.2 MiB |

### dense-div.GF4611686018427387847

| n | Poly Mojo | FLINT | SymPy | Singular |
|---|---|---|---|---|
| 100 | < 1 ms<br>< 1 MiB | < 1 ms<br>< 1 MiB | 11 ms<br>< 1 MiB | n/a |
| 1000 | < 1 ms<br>< 1 MiB | < 1 ms<br>< 1 MiB | 1.21 s<br>< 1 MiB | n/a |
| 10000 | 7.76 ms<br>3.36 MiB | 17.2 ms<br>1.96 MiB | 122 s<br>3.06 MiB | n/a |
| 100000 | 97.8 ms<br>12.6 MiB | 250 ms<br>25.9 MiB | timeout | n/a |

### dense-div.ZZ-1024b

| n | Poly Mojo | FLINT | SymPy | Singular |
|---|---|---|---|---|
| 100 | 1.96 ms<br>< 1 MiB | 1.59 ms<br>< 1 MiB | 4.49 ms<br>< 1 MiB | 2.21 ms<br>< 1 MiB |
| 1000 | 18.4 ms<br>2.96 MiB | 23.3 ms<br>2.61 MiB | 543 ms<br>< 1 MiB | 206 ms<br>1.12 MiB |
| 10000 | 302 ms<br>51.9 MiB | 598 ms<br>41.6 MiB | 52.1 s<br>7.41 MiB | 23.8 s<br>6.45 MiB |
| 100000 | 3.44 s<br>425 MiB | 10.3 s<br>384 MiB | timeout | timeout |

### dense-div.ZZ-64b

| n | Poly Mojo | FLINT | SymPy | Singular |
|---|---|---|---|---|
| 100 | < 1 ms<br>< 1 MiB | < 1 ms<br>< 1 MiB | 2.39 ms<br>< 1 MiB | < 1 ms<br>< 1 MiB |
| 1000 | 1.75 ms<br>< 1 MiB | 2.88 ms<br>< 1 MiB | 269 ms<br>< 1 MiB | 37.8 ms<br>< 1 MiB |
| 10000 | 25.8 ms<br>8.29 MiB | 62.2 ms<br>4.27 MiB | 27.3 s<br>2.99 MiB | 4.05 s<br>2.67 MiB |
| 100000 | 275 ms<br>54.8 MiB | 885 ms<br>44.6 MiB | timeout | 396 s<br>22.1 MiB |

### dense-div.ZZ-8b

| n | Poly Mojo | FLINT | SymPy | Singular |
|---|---|---|---|---|
| 100 | < 1 ms<br>< 1 MiB | < 1 ms<br>< 1 MiB | 2.6 ms<br>< 1 MiB | < 1 ms<br>< 1 MiB |
| 1000 | < 1 ms<br>< 1 MiB | < 1 ms<br>< 1 MiB | 295 ms<br>< 1 MiB | 6.91 ms<br>< 1 MiB |
| 10000 | 10.5 ms<br>5.84 MiB | 9.45 ms<br>< 1 MiB | 29.6 s<br>3.05 MiB | 731 ms<br>1.37 MiB |
| 100000 | 119 ms<br>33.1 MiB | 166 ms<br>4.14 MiB | timeout | 87 s<br>14.4 MiB |

### dense-mul.GF251

| n | Poly Mojo | FLINT | SymPy | Singular |
|---|---|---|---|---|
| 100 | < 1 ms<br>< 1 MiB | < 1 ms<br>< 1 MiB | 13 ms<br>< 1 MiB | < 1 ms<br>1.19 MiB |
| 1000 | < 1 ms<br>< 1 MiB | < 1 ms<br>< 1 MiB | 1.24 s<br>< 1 MiB | < 1 ms<br>1.36 MiB |
| 10000 | < 1 ms<br>1.37 MiB | < 1 ms<br>< 1 MiB | 128 s<br>6.62 MiB | 2.11 ms<br>2.66 MiB |
| 100000 | 6.06 ms<br>10.8 MiB | 15.2 ms<br>5.07 MiB | timeout | 27.2 ms<br>16.3 MiB |

### dense-mul.GF4611686018427387847

| n | Poly Mojo | FLINT | SymPy | Singular |
|---|---|---|---|---|
| 100 | < 1 ms<br>< 1 MiB | < 1 ms<br>< 1 MiB | 11 ms<br>< 1 MiB | n/a |
| 1000 | < 1 ms<br>< 1 MiB | < 1 ms<br>< 1 MiB | 1.21 s<br>< 1 MiB | n/a |
| 10000 | 1.72 ms<br>1.15 MiB | 4.61 ms<br>1.59 MiB | 127 s<br>6.51 MiB | n/a |
| 100000 | 21.6 ms<br>10.7 MiB | 63.4 ms<br>13.7 MiB | timeout | n/a |

### dense-mul.ZZ-1024b

| n | Poly Mojo | FLINT | SymPy | Singular |
|---|---|---|---|---|
| 100 | < 1 ms<br>1.08 MiB | < 1 ms<br>< 1 MiB | 4.49 ms<br>< 1 MiB | < 1 ms<br>1.73 MiB |
| 1000 | 6.45 ms<br>4.13 MiB | 4.69 ms<br>2.26 MiB | 440 ms<br>1.1 MiB | 5.23 ms<br>3.69 MiB |
| 10000 | 82.7 ms<br>44.4 MiB | 168 ms<br>34.5 MiB | 49.2 s<br>11.3 MiB | 145 ms<br>39.3 MiB |
| 100000 | 1.25 s<br>450 MiB | 1.86 s<br>357 MiB | timeout | 1.83 s<br>396 MiB |

### dense-mul.ZZ-64b

| n | Poly Mojo | FLINT | SymPy | Singular |
|---|---|---|---|---|
| 100 | < 1 ms<br>< 1 MiB | < 1 ms<br>< 1 MiB | 2.11 ms<br>< 1 MiB | < 1 ms<br>1.42 MiB |
| 1000 | < 1 ms<br>< 1 MiB | < 1 ms<br>< 1 MiB | 236 ms<br>< 1 MiB | < 1 ms<br>1.69 MiB |
| 10000 | 5.26 ms<br>4.03 MiB | 6.94 ms<br>1.89 MiB | 25 s<br>5.61 MiB | 10.8 ms<br>5.39 MiB |
| 100000 | 64.9 ms<br>31.2 MiB | 94.9 ms<br>23.9 MiB | timeout | 168 ms<br>42.1 MiB |

### dense-mul.ZZ-8b

| n | Poly Mojo | FLINT | SymPy | Singular |
|---|---|---|---|---|
| 100 | < 1 ms<br>< 1 MiB | < 1 ms<br>< 1 MiB | 2.05 ms<br>< 1 MiB | < 1 ms<br>1.29 MiB |
| 1000 | < 1 ms<br>< 1 MiB | < 1 ms<br>< 1 MiB | 234 ms<br>< 1 MiB | < 1 ms<br>1.5 MiB |
| 10000 | 1.06 ms<br>1.49 MiB | 1.68 ms<br>1.1 MiB | 24.7 s<br>5.26 MiB | 3.33 ms<br>2.67 MiB |
| 100000 | 12 ms<br>12.8 MiB | 21.1 ms<br>6.59 MiB | timeout | 41.6 ms<br>14.7 MiB |

### dense-sqr.GF251

| n | Poly Mojo | FLINT | SymPy | Singular |
|---|---|---|---|---|
| 100 | < 1 ms<br>< 1 MiB | < 1 ms<br>< 1 MiB | 11.2 ms<br>< 1 MiB | < 1 ms<br>1.21 MiB |
| 1000 | < 1 ms<br>< 1 MiB | < 1 ms<br>< 1 MiB | 1.21 s<br>< 1 MiB | < 1 ms<br>1.38 MiB |
| 10000 | < 1 ms<br>< 1 MiB | < 1 ms<br>< 1 MiB | 125 s<br>6.62 MiB | 2.22 ms<br>2.97 MiB |
| 100000 | 4.35 ms<br>6.55 MiB | 13 ms<br>4.97 MiB | timeout | 28.8 ms<br>19.4 MiB |

### dense-sqr.GF4611686018427387847

| n | Poly Mojo | FLINT | SymPy | Singular |
|---|---|---|---|---|
| 100 | < 1 ms<br>< 1 MiB | < 1 ms<br>< 1 MiB | 11.7 ms<br>< 1 MiB | n/a |
| 1000 | < 1 ms<br>< 1 MiB | < 1 ms<br>< 1 MiB | 1.38 s<br>< 1 MiB | n/a |
| 10000 | 1.29 ms<br>1.58 MiB | 3.66 ms<br>1.55 MiB | 131 s<br>6.5 MiB | n/a |
| 100000 | 16 ms<br>4.65 MiB | 49.3 ms<br>11.6 MiB | timeout | n/a |

### dense-sqr.ZZ-1024b

| n | Poly Mojo | FLINT | SymPy | Singular |
|---|---|---|---|---|
| 100 | < 1 ms<br>< 1 MiB | < 1 ms<br>< 1 MiB | 4.27 ms<br>< 1 MiB | < 1 ms<br>1.76 MiB |
| 1000 | 4.92 ms<br>3.57 MiB | 3.28 ms<br>1.59 MiB | 465 ms<br>1.09 MiB | 5.5 ms<br>3.94 MiB |
| 10000 | 69.4 ms<br>39.4 MiB | 101 ms<br>19.9 MiB | 53.6 s<br>11.3 MiB | 162 ms<br>41.4 MiB |
| 100000 | 945 ms<br>411 MiB | 1.22 s<br>203 MiB |  | 2.76 s<br>416 MiB |

