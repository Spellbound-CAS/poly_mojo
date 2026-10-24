# Arithmetic results

Host: Intel(R) Xeon(R) Processor @ 2.10GHz, 4 vCPUs, 16 GiB RAM. Runs are pinned to one CPU, capped at 3G and 600 s.

Versions: Poly Mojo d6549d350302; FLINT python-flint 0.9.0, FLINT 3.6.0; SymPy 1.14.0 (gmpy ground types); Singular 4.4.1.

## Summary

Geometric mean of competitor time ÷ Poly Mojo time over cases both solved with a confirmed answer and where the slower took at least 10 ms (above 1 means Poly Mojo was faster).

| Ladder | FLINT | SymPy | Singular |
|---|---|---|---|
| dense-div.GF251 | 2.09 (*n* = 1) | 4.25e+03 (*n* = 3) | 2.73 (*n* = 1) |
| dense-div.GF4611686018427387847 | 2.59 (*n* = 2) | 1.85e+03 (*n* = 3) | — |
| dense-div.ZZ-1024b | 2.93 (*n* = 3) | 97.9 (*n* = 2) | 44.1 (*n* = 2) |
| dense-div.ZZ-64b | 4.45 (*n* = 2) | 705 (*n* = 2) | 289 (*n* = 3) |
| dense-div.ZZ-8b | 2.24 (*n* = 2) | 2.17e+03 (*n* = 2) | 421 (*n* = 2) |
| dense-mul.GF251 | 2.51 (*n* = 1) | 7.56e+03 (*n* = 3) | 4.66 (*n* = 1) |
| dense-mul.GF4611686018427387847 | 3.08 (*n* = 1) | 4.22e+03 (*n* = 3) | — |

## Answer agreement

Every answer that two or more libraries computed agrees.

## Per-ladder results

Each cell gives the median time per call and the incremental peak memory; `·` marks an instance skipped after a smaller one on the same ladder hit a limit; `n/a` marks an unsupported field.

### dense-div.GF251

| n | Poly Mojo | FLINT | SymPy | Singular |
|---|---|---|---|---|
| 100 | < 1 ms<br>< 1 MiB | < 1 ms<br>< 1 MiB | 13.2 ms<br>< 1 MiB | < 1 ms<br>< 1 MiB |
| 1000 | < 1 ms<br>< 1 MiB | < 1 ms<br>< 1 MiB | 1.42 s<br>< 1 MiB | < 1 ms<br>< 1 MiB |
| 10000 | 2.41 ms<br>1.85 MiB | 3.51 ms<br>< 1 MiB | 184 s<br>5.79 MiB | 5.13 ms<br>1.8 MiB |
| 100000 | 27.7 ms<br>10.8 MiB | 58 ms<br>7.16 MiB | timeout | 75.8 ms<br>12 MiB |

### dense-div.GF4611686018427387847

| n | Poly Mojo | FLINT | SymPy | Singular |
|---|---|---|---|---|
| 100 | < 1 ms<br>< 1 MiB | < 1 ms<br>< 1 MiB | 12.9 ms<br>< 1 MiB | n/a |
| 1000 | < 1 ms<br>< 1 MiB | < 1 ms<br>< 1 MiB | 1.74 s<br>< 1 MiB | n/a |
| 10000 | 8.76 ms<br>3.54 MiB | 20.7 ms<br>1.91 MiB | 149 s<br>3.04 MiB | n/a |
| 100000 | 98.4 ms<br>13.2 MiB | 279 ms<br>25.9 MiB | timeout | n/a |

### dense-div.ZZ-1024b

| n | Poly Mojo | FLINT | SymPy | Singular |
|---|---|---|---|---|
| 100 | 1.84 ms<br>< 1 MiB | 1.31 ms<br>< 1 MiB | 5.18 ms<br>< 1 MiB | 2.6 ms<br>< 1 MiB |
| 1000 | 14.6 ms<br>< 1 MiB | 28.3 ms<br>2.61 MiB | 584 ms<br>< 1 MiB | 261 ms<br>1.12 MiB |
| 10000 | 243 ms<br>2.6 MiB | 686 ms<br>41.5 MiB | 58.2 s<br>7.41 MiB | 26.5 s<br>6.36 MiB |
| 100000 | 2.48 s<br>25.2 MiB | 11.4 s<br>384 MiB | timeout | timeout |

### dense-div.ZZ-64b

| n | Poly Mojo | FLINT | SymPy | Singular |
|---|---|---|---|---|
| 100 | < 1 ms<br>< 1 MiB | < 1 ms<br>< 1 MiB | 2.79 ms<br>< 1 MiB | < 1 ms<br>< 1 MiB |
| 1000 | 1.06 ms<br>< 1 MiB | 3.38 ms<br>< 1 MiB | 299 ms<br>< 1 MiB | 42.5 ms<br>< 1 MiB |
| 10000 | 17.9 ms<br>2.68 MiB | 65.2 ms<br>4.34 MiB | 31.4 s<br>3.01 MiB | 4.59 s<br>2.6 MiB |
| 100000 | 189 ms<br>11.4 MiB | 1.02 s<br>44.5 MiB | timeout | 442 s<br>22.1 MiB |

### dense-div.ZZ-8b

| n | Poly Mojo | FLINT | SymPy | Singular |
|---|---|---|---|---|
| 100 | < 1 ms<br>< 1 MiB | < 1 ms<br>< 1 MiB | 2.65 ms<br>< 1 MiB | < 1 ms<br>< 1 MiB |
| 1000 | < 1 ms<br>< 1 MiB | < 1 ms<br>< 1 MiB | 320 ms<br>< 1 MiB | 7.22 ms<br>< 1 MiB |
| 10000 | 5.88 ms<br>3.07 MiB | 10 ms<br>< 1 MiB | 29.9 s<br>3.27 MiB | 744 ms<br>1.37 MiB |
| 100000 | 60.7 ms<br>7.96 MiB | 179 ms<br>4.2 MiB | timeout | 85 s<br>14.4 MiB |

### dense-mul.GF251

| n | Poly Mojo | FLINT | SymPy | Singular |
|---|---|---|---|---|
| 100 | < 1 ms<br>< 1 MiB | < 1 ms<br>< 1 MiB | 11.8 ms<br>< 1 MiB | < 1 ms<br>1.16 MiB |
| 1000 | < 1 ms<br>< 1 MiB | < 1 ms<br>< 1 MiB | 1.3 s<br>< 1 MiB | < 1 ms<br>1.38 MiB |
| 10000 | < 1 ms<br>1.06 MiB | < 1 ms<br>< 1 MiB | 133 s<br>6.62 MiB | 2.23 ms<br>2.7 MiB |
| 100000 | 6.37 ms<br>10.8 MiB | 16 ms<br>5.14 MiB | timeout | 29.7 ms<br>16.4 MiB |

### dense-mul.GF4611686018427387847

| n | Poly Mojo | FLINT | SymPy | Singular |
|---|---|---|---|---|
| 100 | < 1 ms<br>< 1 MiB | < 1 ms<br>< 1 MiB | 11.4 ms<br>< 1 MiB | n/a |
| 1000 | < 1 ms<br>< 1 MiB | < 1 ms<br>< 1 MiB | 1.21 s<br>< 1 MiB | n/a |
| 10000 | 1.79 ms<br>1.46 MiB | 4.82 ms<br>1.62 MiB | 119 s<br>6.63 MiB | n/a |
| 100000 | 21.9 ms<br>10.7 MiB | 67.3 ms<br>13.8 MiB |  | n/a |

