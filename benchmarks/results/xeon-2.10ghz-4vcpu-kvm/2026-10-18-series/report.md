# Series, evaluation, and interpolation results

Host: Intel(R) Xeon(R) Processor @ 2.10GHz, 4 vCPUs, 16 GiB RAM. Runs are pinned to one CPU, capped at 3G and 600 s.

Versions: Poly Mojo fdba30fb4669; FLINT FLINT 3.3.1 (C, conda-forge libflint); SymPy 1.14.0 (gmpy ground types).

## Summary

Geometric mean of competitor time ÷ Poly Mojo time over cases both solved with a confirmed answer and where the slower took at least 10 ms (above 1 means Poly Mojo was faster).

| Ladder | FLINT | SymPy |
|---|---|---|
| evalq.ZZ | 10.7 (*n* = 2) | 2.7e+03 (*n* = 2) |
| interp.GF4611686018427387847 | 2.09 (*n* = 2) | — |
| interp.QQ | 2.75 (*n* = 1) | 1.58e+04 (*n* = 1) |
| multipoint.GF4611686018427387847 | 2.17 (*n* = 2) | 1.69e+03 (*n* = 2) |
| multipoint.ZZ | 0.573 (*n* = 2) | 59.9 (*n* = 1) |
| series-compose.GF4611686018427387847 | 3.42 (*n* = 3) | 2.75e+04 (*n* = 1) |
| series-compose.QQ | 1.78 (*n* = 2) | 147 (*n* = 2) |
| series-exp.GF4611686018427387847 | 2.06 (*n* = 2) | 7.79e+03 (*n* = 2) |
| series-exp.QQ | 0.298 (*n* = 2) | 80.5 (*n* = 2) |
| series-inv.GF4611686018427387847 | 3 (*n* = 2) | 6.1e+03 (*n* = 2) |
| series-inv.QQ | 1.05 (*n* = 2) | 71 (*n* = 2) |
| series-log.GF4611686018427387847 | 2.5 (*n* = 2) | 8.02e+03 (*n* = 2) |
| series-log.QQ | 1.14 (*n* = 2) | 28.6 (*n* = 1) |
| series-mul.GF4611686018427387847 | 2.17 (*n* = 1) | 7.47e+03 (*n* = 2) |
| series-mul.QQ | — | 514 (*n* = 1) |

## Answer agreement

Every answer that two or more libraries computed agrees.

## Per-ladder results

Each cell gives the median time per call and the incremental peak memory; `·` marks an instance skipped after a smaller one on the same ladder hit a limit; `n/a` marks an unsupported operation.

### evalq.ZZ

| n | Poly Mojo | FLINT | SymPy |
|---|---|---|---|
| 1000 | 1 ms<br>< 1 MiB | 9.41 ms<br>< 1 MiB | 775 ms<br>< 1 MiB |
| 10000 | 23.6 ms<br>2.86 MiB | 274 ms<br>8.49 MiB | 222 s<br>1.34 MiB |
| 100000 | 535 ms<br>26.9 MiB | 5.3 s<br>96.6 MiB | timeout |

### interp.GF4611686018427387847

| n | Poly Mojo | FLINT | SymPy |
|---|---|---|---|
| 1000 | 2.71 ms<br>< 1 MiB | 2.16 ms<br>1.36 MiB | n/a |
| 10000 | 44.8 ms<br>9.47 MiB | 74.5 ms<br>3.15 MiB | n/a |
| 100000 | 660 ms<br>59.7 MiB | 1.74 s<br>31.7 MiB | n/a |

### interp.QQ

| n | Poly Mojo | FLINT | SymPy |
|---|---|---|---|
| 100 | 2.85 ms<br>< 1 MiB | < 1 ms<br>1 MiB | 45.1 s<br>14.9 MiB |
| 1000 | 1.26 s<br>6.75 MiB | 3.45 s<br>6.5 MiB | timeout |
| 10000 | timeout | timeout | · |

### multipoint.GF4611686018427387847

| n | Poly Mojo | FLINT | SymPy |
|---|---|---|---|
| 1000 | 1.85 ms<br>< 1 MiB | 1.12 ms<br>1.16 MiB | 1.09 s<br>< 1 MiB |
| 10000 | 32 ms<br>8.5 MiB | 51.3 ms<br>2.83 MiB | 155 s<br>1.42 MiB |
| 100000 | 430 ms<br>45.7 MiB | 1.27 s<br>29.9 MiB | timeout |

### multipoint.ZZ

| n | Poly Mojo | FLINT | SymPy |
|---|---|---|---|
| 100 | < 1 ms<br>< 1 MiB | < 1 ms<br>< 1 MiB | 6.5 ms<br>< 1 MiB |
| 1000 | 156 ms<br>3.92 MiB | 76.7 ms<br>4.25 MiB | 9.36 s<br>3.95 MiB |
| 10000 | 33.1 s<br>356 MiB | 22.2 s<br>362 MiB | timeout |

### series-compose.GF4611686018427387847

| n | Poly Mojo | FLINT | SymPy |
|---|---|---|---|
| 1000 | 5.71 ms<br>1.6 MiB | 23.1 ms<br>1.77 MiB | 157 s<br>126 MiB |
| 10000 | 334 ms<br>21.7 MiB | 1.33 s<br>24.3 MiB | timeout |
| 100000 | 21 s<br>545 MiB | 52 s<br>681 MiB | · |

### series-compose.QQ

| n | Poly Mojo | FLINT | SymPy |
|---|---|---|---|
| 100 | 1.07 ms<br>1.05 MiB | 1.09 ms<br>< 1 MiB | 51.6 ms<br>< 1 MiB |
| 1000 | 186 ms<br>20.1 MiB | 281 ms<br>9.85 MiB | 82.7 s<br>244 MiB |
| 10000 | 83.5 s<br>1.16 GiB | 174 s<br>894 MiB | timeout |

### series-exp.GF4611686018427387847

| n | Poly Mojo | FLINT | SymPy |
|---|---|---|---|
| 1000 | < 1 ms<br>< 1 MiB | < 1 ms<br>< 1 MiB | 1.78 s<br>< 1 MiB |
| 10000 | 12.4 ms<br>4.37 MiB | 21 ms<br>2.71 MiB | 175 s<br>9.38 MiB |
| 100000 | 121 ms<br>29.7 MiB | 304 ms<br>21.9 MiB | timeout |

### series-exp.QQ

| n | Poly Mojo | FLINT | SymPy |
|---|---|---|---|
| 100 | 1.24 ms<br>< 1 MiB | < 1 ms<br>< 1 MiB | 23.5 ms<br>< 1 MiB |
| 1000 | 187 ms<br>29.7 MiB | 63.2 ms<br>1.46 MiB | 64.1 s<br>3.71 MiB |
| 10000 | 85 s<br>2.32 GiB | 22.4 s<br>1.35 GiB | timeout |

### series-inv.GF4611686018427387847

| n | Poly Mojo | FLINT | SymPy |
|---|---|---|---|
| 1000 | < 1 ms<br>< 1 MiB | < 1 ms<br>< 1 MiB | 808 ms<br>< 1 MiB |
| 10000 | 4.54 ms<br>4.08 MiB | 11.5 ms<br>2.43 MiB | 82.6 s<br>4.86 MiB |
| 100000 | 51.4 ms<br>25.3 MiB | 184 ms<br>20.2 MiB | timeout |

### series-inv.QQ

| n | Poly Mojo | FLINT | SymPy |
|---|---|---|---|
| 100 | < 1 ms<br>< 1 MiB | < 1 ms<br>< 1 MiB | 2.66 ms<br>< 1 MiB |
| 1000 | 16.7 ms<br>8.34 MiB | 17.1 ms<br>3.24 MiB | 496 ms<br>1.47 MiB |
| 10000 | 2.14 s<br>312 MiB | 2.3 s<br>308 MiB | 361 s<br>67.3 MiB |

### series-log.GF4611686018427387847

| n | Poly Mojo | FLINT | SymPy |
|---|---|---|---|
| 1000 | < 1 ms<br>< 1 MiB | < 1 ms<br>< 1 MiB | 1.45 s<br>< 1 MiB |
| 10000 | 6.3 ms<br>5.33 MiB | 14.2 ms<br>2.43 MiB | 142 s<br>8.02 MiB |
| 100000 | 69.9 ms<br>29.1 MiB | 193 ms<br>22 MiB | timeout |

### series-log.QQ

| n | Poly Mojo | FLINT | SymPy |
|---|---|---|---|
| 100 | < 1 ms<br>< 1 MiB | < 1 ms<br>< 1 MiB | 4.47 ms<br>< 1 MiB |
| 1000 | 31.2 ms<br>17.3 MiB | 34.1 ms<br>6.71 MiB | 891 ms<br>1.79 MiB |
| 10000 | 6.5 s<br>1.06 GiB | 7.72 s<br>1.15 GiB | timeout |

### series-mul.GF4611686018427387847

| n | Poly Mojo | FLINT | SymPy |
|---|---|---|---|
| 1000 | < 1 ms<br>< 1 MiB | < 1 ms<br>< 1 MiB | 601 ms<br>< 1 MiB |
| 10000 | 3.01 ms<br>3.45 MiB | 5.8 ms<br>2.04 MiB | 61.1 s<br>2.99 MiB |
| 100000 | 32.6 ms<br>32.5 MiB | 70.7 ms<br>23.6 MiB | timeout |

### series-mul.QQ

| n | Poly Mojo | FLINT | SymPy |
|---|---|---|---|
| 100 | < 1 ms<br>< 1 MiB | < 1 ms<br>< 1 MiB | 1.5 ms<br>< 1 MiB |
| 1000 | < 1 ms<br>< 1 MiB | < 1 ms<br>< 1 MiB | 170 ms<br>< 1 MiB |
| 10000 | 3.46 ms<br>3.18 MiB ? |  |  |

