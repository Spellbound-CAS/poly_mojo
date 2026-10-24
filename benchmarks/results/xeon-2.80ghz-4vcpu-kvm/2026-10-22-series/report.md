# Series, evaluation, and interpolation results

Host: Intel(R) Xeon(R) Processor @ 2.80GHz, 4 vCPUs, 16 GiB RAM. Runs are pinned to one CPU, capped at 3G and 600 s.

Versions: Poly Mojo d3e021c61ab6; FLINT FLINT 3.3.1 (C, conda-forge libflint); SymPy 1.14.0 (gmpy ground types).

## Summary

Geometric mean of competitor time ÷ Poly Mojo time over cases both solved with a confirmed answer and where the slower took at least 10 ms (above 1 means Poly Mojo was faster).

| Ladder | FLINT | SymPy |
|---|---|---|
| evalq.ZZ | 11.5 (*n* = 3) | 2.94e+03 (*n* = 2) |
| interp.GF4611686018427387847 | 2.3 (*n* = 2) | — |
| interp.QQ | 2.35 (*n* = 1) | 1.78e+04 (*n* = 1) |
| multipoint.GF4611686018427387847 | 1.94 (*n* = 2) | 2.39e+03 (*n* = 2) |
| multipoint.ZZ | 0.724 (*n* = 2) | 76.1 (*n* = 1) |
| series-compose.GF4611686018427387847 | 2.96 (*n* = 3) | 2.77e+04 (*n* = 1) |
| series-compose.QQ | 1.86 (*n* = 2) | 201 (*n* = 2) |
| series-exp.GF4611686018427387847 | 1.65 (*n* = 2) | 4.04e+03 (*n* = 1) |
| series-exp.QQ | 0.349 (*n* = 2) | — |
| series-inv.GF4611686018427387847 | 1.96 (*n* = 2) | — |
| series-inv.QQ | 1.06 (*n* = 2) | — |
| series-log.GF4611686018427387847 | 1.99 (*n* = 2) | — |
| series-log.QQ | 1.82 (*n* = 2) | — |
| series-mul.GF4611686018427387847 | 2.94 (*n* = 1) | — |
| series-mul.QQ | — | — |
| series-pow.GF4611686018427387847 | 1.87 (*n* = 2) | — |
| series-pow.QQ | 0.546 (*n* = 1) | — |
| series-revert.GF4611686018427387847 | 1.69 (*n* = 3) | — |
| series-revert.QQ | 1.98 (*n* = 1) | — |
| series-sqrt.GF4611686018427387847 | 1.36 (*n* = 2) | — |
| series-sqrt.QQ | 0.897 (*n* = 2) | — |

## Answer agreement

Every answer that two or more libraries computed agrees.

## Per-ladder results

Each cell gives the median time per call and the incremental peak memory; `·` marks an instance skipped after a smaller one on the same ladder hit a limit; `n/a` marks an unsupported operation.

### evalq.ZZ

| n | Poly Mojo | FLINT | SymPy |
|---|---|---|---|
| 1000 | < 1 ms<br>< 1 MiB | 10.1 ms<br>< 1 MiB | 798 ms<br>< 1 MiB |
| 10000 | 23.5 ms<br>2.64 MiB | 309 ms<br>8.55 MiB | 243 s<br>1.39 MiB |
| 100000 | 541 ms<br>26.8 MiB | 5.85 s<br>96.5 MiB | timeout |

### interp.GF4611686018427387847

| n | Poly Mojo | FLINT | SymPy |
|---|---|---|---|
| 1000 | 3.74 ms<br>1.43 MiB | 2.62 ms<br>1.29 MiB | n/a |
| 10000 | 59.8 ms<br>9.84 MiB | 142 ms<br>3.21 MiB | n/a |
| 100000 | 1.02 s<br>59 MiB | 2.27 s<br>31.8 MiB | n/a |

### interp.QQ

| n | Poly Mojo | FLINT | SymPy |
|---|---|---|---|
| 100 | 3.56 ms<br>< 1 MiB | 1.49 ms<br>1.19 MiB | 63.5 s<br>14.9 MiB |
| 1000 | 1.7 s<br>6.7 MiB | 3.98 s<br>6.71 MiB | timeout |
| 10000 | timeout | timeout | · |

### multipoint.GF4611686018427387847

| n | Poly Mojo | FLINT | SymPy |
|---|---|---|---|
| 1000 | 2.8 ms<br>1.25 MiB | 1.54 ms<br>1.25 MiB | 2.34 s<br>< 1 MiB |
| 10000 | 40.1 ms<br>8.29 MiB | 56.2 ms<br>2.98 MiB | 275 s<br>1.34 MiB |
| 100000 | 655 ms<br>45.8 MiB | 1.76 s<br>30 MiB | timeout |

### multipoint.ZZ

| n | Poly Mojo | FLINT | SymPy |
|---|---|---|---|
| 100 | < 1 ms<br>< 1 MiB | < 1 ms<br>< 1 MiB | 7.58 ms<br>< 1 MiB |
| 1000 | 136 ms<br>4.31 MiB | 101 ms<br>4.44 MiB | 10.3 s<br>3.93 MiB |
| 10000 | 41.8 s<br>363 MiB | 29.3 s<br>362 MiB | timeout |

### series-compose.GF4611686018427387847

| n | Poly Mojo | FLINT | SymPy |
|---|---|---|---|
| 1000 | 7.93 ms<br>1.17 MiB | 21.9 ms<br>1.96 MiB | 219 s<br>126 MiB |
| 10000 | 488 ms<br>20.9 MiB | 1.6 s<br>24.5 MiB | timeout |
| 100000 | 23.1 s<br>519 MiB | 66.3 s<br>681 MiB | · |

### series-compose.QQ

| n | Poly Mojo | FLINT | SymPy |
|---|---|---|---|
| 100 | < 1 ms<br>< 1 MiB | 1.06 ms<br>1.1 MiB | 66.8 ms<br>< 1 MiB |
| 1000 | 172 ms<br>21.9 MiB | 285 ms<br>9.95 MiB | 102 s<br>245 MiB |
| 10000 | 96.1 s<br>1.21 GiB | 200 s<br>894 MiB | timeout |

### series-exp.GF4611686018427387847

| n | Poly Mojo | FLINT | SymPy |
|---|---|---|---|
| 1000 | < 1 ms<br>< 1 MiB | < 1 ms<br>< 1 MiB | 2.76 s<br>< 1 MiB |
| 10000 | 15.6 ms<br>5.23 MiB | 20.8 ms<br>2.79 MiB |  |
| 100000 | 161 ms<br>28.9 MiB | 332 ms<br>22.1 MiB |  |

### series-exp.QQ

| n | Poly Mojo | FLINT | SymPy |
|---|---|---|---|
| 100 | < 1 ms<br>< 1 MiB | < 1 ms<br>< 1 MiB |  |
| 1000 | 155 ms<br>4.59 MiB | 70.6 ms<br>1.63 MiB |  |
| 10000 | 113 s<br>2.29 GiB | 30.2 s<br>1.35 GiB |  |

### series-inv.GF4611686018427387847

| n | Poly Mojo | FLINT | SymPy |
|---|---|---|---|
| 1000 | < 1 ms<br>< 1 MiB | < 1 ms<br>< 1 MiB |  |
| 10000 | 6.38 ms<br>4.53 MiB | 12.3 ms<br>2.56 MiB |  |
| 100000 | 84.1 ms<br>24.9 MiB | 168 ms<br>20.3 MiB |  |

### series-inv.QQ

| n | Poly Mojo | FLINT | SymPy |
|---|---|---|---|
| 100 | < 1 ms<br>< 1 MiB | < 1 ms<br>< 1 MiB |  |
| 1000 | 19.6 ms<br>8.7 MiB | 18.5 ms<br>3.44 MiB |  |
| 10000 | 2.02 s<br>228 MiB | 2.41 s<br>308 MiB |  |

### series-log.GF4611686018427387847

| n | Poly Mojo | FLINT | SymPy |
|---|---|---|---|
| 1000 | < 1 ms<br>< 1 MiB | < 1 ms<br>< 1 MiB |  |
| 10000 | 11.6 ms<br>5.42 MiB | 21.7 ms<br>2.66 MiB |  |
| 100000 | 121 ms<br>27.4 MiB | 258 ms<br>22.1 MiB |  |

### series-log.QQ

| n | Poly Mojo | FLINT | SymPy |
|---|---|---|---|
| 100 | < 1 ms<br>< 1 MiB | < 1 ms<br>< 1 MiB |  |
| 1000 | 25.9 ms<br>10.1 MiB | 34.4 ms<br>6.85 MiB |  |
| 10000 | 5.22 s<br>577 MiB | 13.1 s<br>1.15 GiB |  |

### series-mul.GF4611686018427387847

| n | Poly Mojo | FLINT | SymPy |
|---|---|---|---|
| 1000 | < 1 ms<br>< 1 MiB | < 1 ms<br>< 1 MiB |  |
| 10000 | 2.97 ms<br>1.2 MiB | 6.63 ms<br>2.26 MiB |  |
| 100000 | 29.7 ms<br>10.1 MiB | 87.5 ms<br>23.7 MiB |  |

### series-mul.QQ

| n | Poly Mojo | FLINT | SymPy |
|---|---|---|---|
| 100 | < 1 ms<br>< 1 MiB | < 1 ms<br>< 1 MiB |  |
| 1000 | < 1 ms<br>< 1 MiB | < 1 ms<br>< 1 MiB |  |
| 10000 | 3.93 ms<br>1.86 MiB | 1.46 ms<br>< 1 MiB |  |

### series-pow.GF4611686018427387847

| n | Poly Mojo | FLINT | SymPy |
|---|---|---|---|
| 1000 | < 1 ms<br>< 1 MiB | < 1 ms<br>< 1 MiB |  |
| 10000 | 8.52 ms<br>2.36 MiB | 15.3 ms<br>2.46 MiB |  |
| 100000 | 99.5 ms<br>18.4 MiB | 194 ms<br>25.2 MiB |  |

### series-pow.QQ

| n | Poly Mojo | FLINT | SymPy |
|---|---|---|---|
| 100 | < 1 ms<br>< 1 MiB | < 1 ms<br>< 1 MiB |  |
| 1000 | < 1 ms<br>< 1 MiB | < 1 ms<br>< 1 MiB |  |
| 10000 | 11.3 ms<br>3.25 MiB | 6.16 ms<br>1.62 MiB |  |

### series-revert.GF4611686018427387847

| n | Poly Mojo | FLINT | SymPy |
|---|---|---|---|
| 1000 | 11.8 ms<br>1.85 MiB | 20.6 ms<br>1.37 MiB |  |
| 10000 | 890 ms<br>26.1 MiB | 1.47 s<br>10.6 MiB |  |
| 100000 | 34 s<br>603 MiB | 57.2 s<br>268 MiB |  |

### series-revert.QQ

| n | Poly Mojo | FLINT | SymPy |
|---|---|---|---|
| 100 | 3.92 ms<br>1.28 MiB | 2.55 ms<br>1.45 MiB |  |
| 1000 | 1.99 s<br>81.3 MiB | 3.94 s<br>66.2 MiB |  |
| 10000 | OOM | OOM |  |

### series-sqrt.GF4611686018427387847

| n | Poly Mojo | FLINT | SymPy |
|---|---|---|---|
| 1000 | < 1 ms<br>< 1 MiB | < 1 ms<br>< 1 MiB |  |
| 10000 | 13.6 ms<br>4.57 MiB | 14.1 ms<br>< 1 MiB |  |
| 100000 | 168 ms<br>38.8 MiB | 298 ms<br>17.2 MiB |  |

### series-sqrt.QQ

| n | Poly Mojo | FLINT | SymPy |
|---|---|---|---|
| 100 | 1.04 ms<br>1.19 MiB | < 1 ms<br>1.1 MiB |  |
| 1000 | 120 ms<br>23.7 MiB | 137 ms<br>16.4 MiB |  |
| 10000 | 35.5 s<br>1.81 GiB | 25 s<br>1.31 GiB |  |

