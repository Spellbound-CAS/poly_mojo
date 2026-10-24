# Series, evaluation, and interpolation results

Host: Intel(R) Xeon(R) Processor @ 2.80GHz, 4 vCPUs, 16 GiB RAM. Runs are pinned to one CPU, capped at 3G and 600 s.

Versions: Poly Mojo d1f992ffe139; FLINT FLINT 3.3.1 (C, conda-forge libflint).

## Summary

Geometric mean of competitor time ÷ Poly Mojo time over cases both solved with a confirmed answer and where the slower took at least 10 ms (above 1 means Poly Mojo was faster).

| Ladder | FLINT |
|---|---|
| evalq.ZZ | 11.7 (*n* = 3) |
| interp.GF4611686018427387847 | 1.84 (*n* = 2) |
| interp.QQ | 2.31 (*n* = 1) |
| multipoint.GF4611686018427387847 | 2.02 (*n* = 2) |
| multipoint.ZZ | 0.657 (*n* = 2) |
| series-compose.GF4611686018427387847 | 2.9 (*n* = 3) |
| series-compose.QQ | 1.75 (*n* = 2) |
| series-exp.GF4611686018427387847 | 1.7 (*n* = 2) |
| series-exp.QQ | 0.33 (*n* = 2) |
| series-inv.GF4611686018427387847 | 2.15 (*n* = 2) |
| series-inv.QQ | 1.17 (*n* = 2) |
| series-log.GF4611686018427387847 | 1.92 (*n* = 2) |
| series-log.QQ | 1.56 (*n* = 2) |
| series-mul.GF4611686018427387847 | 2.42 (*n* = 1) |
| series-mul.QQ | — |
| series-pow.GF4611686018427387847 | 1.87 (*n* = 2) |
| series-pow.QQ | 0.509 (*n* = 1) |
| series-revert.GF4611686018427387847 | 1.76 (*n* = 3) |
| series-revert.QQ | 1.6 (*n* = 1) |
| series-sqrt.GF4611686018427387847 | 1.49 (*n* = 2) |
| series-sqrt.QQ | 0.989 (*n* = 2) |

## Answer agreement

Every answer that two or more libraries computed agrees.

## Per-ladder results

Each cell gives the median time per call and the incremental peak memory; `·` marks an instance skipped after a smaller one on the same ladder hit a limit; `n/a` marks an unsupported operation.

### evalq.ZZ

| n | Poly Mojo | FLINT |
|---|---|---|
| 1000 | < 1 ms<br>< 1 MiB | 10.4 ms<br>< 1 MiB |
| 10000 | 25 ms<br>2.85 MiB | 286 ms<br>8.61 MiB |
| 100000 | 504 ms<br>23.9 MiB | 6.45 s<br>96.5 MiB |

### interp.GF4611686018427387847

| n | Poly Mojo | FLINT |
|---|---|---|
| 1000 | 3.85 ms<br>1.29 MiB | 2.46 ms<br>1.27 MiB |
| 10000 | 62.1 ms<br>9.42 MiB | 91 ms<br>3.24 MiB |
| 100000 | 836 ms<br>60.8 MiB | 1.94 s<br>31.8 MiB |

### interp.QQ

| n | Poly Mojo | FLINT |
|---|---|---|
| 100 | 3.02 ms<br>< 1 MiB | 1.18 ms<br>1.23 MiB |
| 1000 | 1.61 s<br>6.31 MiB | 3.72 s<br>6.79 MiB |
| 10000 | timeout | timeout |

### multipoint.GF4611686018427387847

| n | Poly Mojo | FLINT |
|---|---|---|
| 1000 | 2.65 ms<br>< 1 MiB | 1.31 ms<br>1.24 MiB |
| 10000 | 38.6 ms<br>8.27 MiB | 56 ms<br>2.87 MiB |
| 100000 | 568 ms<br>45.2 MiB | 1.6 s<br>30 MiB |

### multipoint.ZZ

| n | Poly Mojo | FLINT |
|---|---|---|
| 100 | < 1 ms<br>< 1 MiB | < 1 ms<br>< 1 MiB |
| 1000 | 125 ms<br>4.29 MiB | 86.3 ms<br>4.35 MiB |
| 10000 | 39 s<br>363 MiB | 24.4 s<br>362 MiB |

### series-compose.GF4611686018427387847

| n | Poly Mojo | FLINT |
|---|---|---|
| 1000 | 7.58 ms<br>1.51 MiB | 21.1 ms<br>1.96 MiB |
| 10000 | 372 ms<br>20.5 MiB | 1.17 s<br>24.6 MiB |
| 100000 | 20.2 s<br>519 MiB | 56.4 s<br>681 MiB |

### series-compose.QQ

| n | Poly Mojo | FLINT |
|---|---|---|
| 100 | < 1 ms<br>< 1 MiB | 1.02 ms<br>1.09 MiB |
| 1000 | 178 ms<br>22.5 MiB | 267 ms<br>9.96 MiB |
| 10000 | 81.8 s<br>1.22 GiB | 167 s<br>894 MiB |

### series-exp.GF4611686018427387847

| n | Poly Mojo | FLINT |
|---|---|---|
| 1000 | < 1 ms<br>< 1 MiB | < 1 ms<br>< 1 MiB |
| 10000 | 14.7 ms<br>4.87 MiB | 19.9 ms<br>2.81 MiB |
| 100000 | 148 ms<br>28.8 MiB | 315 ms<br>22.1 MiB |

### series-exp.QQ

| n | Poly Mojo | FLINT |
|---|---|---|
| 100 | < 1 ms<br>< 1 MiB | < 1 ms<br>< 1 MiB |
| 1000 | 138 ms<br>4.65 MiB | 64.3 ms<br>1.61 MiB |
| 10000 | 91.3 s<br>2.35 GiB | 21.4 s<br>1.35 GiB |

### series-inv.GF4611686018427387847

| n | Poly Mojo | FLINT |
|---|---|---|
| 1000 | < 1 ms<br>< 1 MiB | < 1 ms<br>< 1 MiB |
| 10000 | 6.04 ms<br>4.46 MiB | 11.2 ms<br>2.56 MiB |
| 100000 | 66.8 ms<br>25.6 MiB | 167 ms<br>20.3 MiB |

### series-inv.QQ

| n | Poly Mojo | FLINT |
|---|---|---|
| 100 | < 1 ms<br>< 1 MiB | < 1 ms<br>< 1 MiB |
| 1000 | 14.9 ms<br>7.53 MiB | 17.1 ms<br>3.43 MiB |
| 10000 | 1.89 s<br>228 MiB | 2.24 s<br>308 MiB |

### series-log.GF4611686018427387847

| n | Poly Mojo | FLINT |
|---|---|---|
| 1000 | < 1 ms<br>< 1 MiB | < 1 ms<br>< 1 MiB |
| 10000 | 8.64 ms<br>5.71 MiB | 15 ms<br>2.66 MiB |
| 100000 | 95.7 ms<br>27.4 MiB | 203 ms<br>22.1 MiB |

### series-log.QQ

| n | Poly Mojo | FLINT |
|---|---|---|
| 100 | < 1 ms<br>< 1 MiB | < 1 ms<br>< 1 MiB |
| 1000 | 23.1 ms<br>10 MiB | 33.7 ms<br>6.81 MiB |
| 10000 | 4.7 s<br>594 MiB | 7.86 s<br>1.15 GiB |

### series-mul.GF4611686018427387847

| n | Poly Mojo | FLINT |
|---|---|---|
| 1000 | < 1 ms<br>< 1 MiB | < 1 ms<br>< 1 MiB |
| 10000 | 2.48 ms<br>< 1 MiB | 5.89 ms<br>2.32 MiB |
| 100000 | 30.4 ms<br>10.1 MiB | 73.4 ms<br>23.7 MiB |

### series-mul.QQ

| n | Poly Mojo | FLINT |
|---|---|---|
| 100 | < 1 ms<br>< 1 MiB | < 1 ms<br>< 1 MiB |
| 1000 | < 1 ms<br>< 1 MiB | < 1 ms<br>< 1 MiB |
| 10000 | 3.57 ms<br>1.88 MiB | 1.54 ms<br>< 1 MiB |

### series-pow.GF4611686018427387847

| n | Poly Mojo | FLINT |
|---|---|---|
| 1000 | < 1 ms<br>< 1 MiB | < 1 ms<br>< 1 MiB |
| 10000 | 7.55 ms<br>2.1 MiB | 14 ms<br>2.4 MiB |
| 100000 | 89.4 ms<br>18.4 MiB | 168 ms<br>25.2 MiB |

### series-pow.QQ

| n | Poly Mojo | FLINT |
|---|---|---|
| 100 | < 1 ms<br>< 1 MiB | < 1 ms<br>< 1 MiB |
| 1000 | < 1 ms<br>< 1 MiB | < 1 ms<br>< 1 MiB |
| 10000 | 11.7 ms<br>4.03 MiB | 5.94 ms<br>1.62 MiB |

### series-revert.GF4611686018427387847

| n | Poly Mojo | FLINT |
|---|---|---|
| 1000 | 11.2 ms<br>1.82 MiB | 21.1 ms<br>1.37 MiB |
| 10000 | 784 ms<br>26.3 MiB | 1.3 s<br>10.6 MiB |
| 100000 | 31 s<br>610 MiB | 53.6 s<br>268 MiB |

### series-revert.QQ

| n | Poly Mojo | FLINT |
|---|---|---|
| 100 | 5.44 ms<br>1.27 MiB | 2.56 ms<br>1.46 MiB |
| 1000 | 1.95 s<br>81.3 MiB | 3.13 s<br>66.3 MiB |
| 10000 | OOM | OOM |

### series-sqrt.GF4611686018427387847

| n | Poly Mojo | FLINT |
|---|---|---|
| 1000 | < 1 ms<br>< 1 MiB | < 1 ms<br>< 1 MiB |
| 10000 | 12.4 ms<br>4.66 MiB | 13.9 ms<br>< 1 MiB |
| 100000 | 142 ms<br>39 MiB | 280 ms<br>17.1 MiB |

### series-sqrt.QQ

| n | Poly Mojo | FLINT |
|---|---|---|
| 100 | 1.03 ms<br>< 1 MiB | < 1 ms<br>1.08 MiB |
| 1000 | 111 ms<br>24 MiB | 135 ms<br>16.4 MiB |
| 10000 | 25.1 s<br>1.75 GiB | 20.1 s<br>1.31 GiB |

