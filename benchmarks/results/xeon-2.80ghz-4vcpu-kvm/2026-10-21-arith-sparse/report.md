# Arithmetic results

Host: Intel(R) Xeon(R) Processor @ 2.80GHz, 4 vCPUs, 16 GiB RAM. Runs are pinned to one CPU, capped at 3G and 600 s.

Versions: Poly Mojo d1f992ffe139; FLINT python-flint 0.9.0, FLINT 3.6.0; Singular 4.4.1.

## Summary

Geometric mean of competitor time ÷ Poly Mojo time over cases both solved with a confirmed answer and where the slower took at least 10 ms (above 1 means Poly Mojo was faster).

| Ladder | FLINT | Singular |
|---|---|---|
| fateman-div.ZZ | 4.57 (*n* = 3) | 4.51 (*n* = 4) |
| fateman-mul.ZZ | 1.28 (*n* = 3) | 1.37 (*n* = 3) |
| fateman-pow.ZZ | — | 8.51 (*n* = 1) |
| mp-div.ZZ | 1.47 (*n* = 1) | 2.15 (*n* = 1) |
| mp-mul.ZZ | 1.02 (*n* = 1) | 1.56 (*n* = 1) |
| mp-pow.ZZ | — | — |

## Answer agreement

Every answer that two or more libraries computed agrees.

## Per-ladder results

Each cell gives the median time per call and the incremental peak memory; `·` marks an instance skipped after a smaller one on the same ladder hit a limit; `n/a` marks an unsupported field.

### fateman-div.ZZ

| n | Poly Mojo | FLINT | Singular |
|---|---|---|---|
| 10 | 2.81 ms<br>< 1 MiB | 9.69 ms<br>< 1 MiB | 10.9 ms<br>< 1 MiB |
| 15 | 32.2 ms<br>< 1 MiB | 148 ms<br>< 1 MiB | 150 ms<br>1.36 MiB |
| 20 | 248 ms<br>< 1 MiB | 1.2 s<br>< 1 MiB | 1.33 s<br>7.38 MiB |
| 30 | 5.95 s<br>4.29 MiB | 25.7 s<br>2.59 MiB | 25.3 s<br>27.8 MiB |

### fateman-mul.ZZ

| n | Poly Mojo | FLINT | Singular |
|---|---|---|---|
| 10 | 1.63 ms<br>1.97 MiB | 1.08 ms<br>< 1 MiB | 1.8 ms<br>1.72 MiB |
| 15 | 15.8 ms<br>7.07 MiB | 14.5 ms<br>1.64 MiB | 18.6 ms<br>5.3 MiB |
| 20 | 192 ms<br>25.6 MiB | 130 ms<br>7.59 MiB | 168 ms<br>18.1 MiB |
| 30 | 4.72 s<br>126 MiB | 15.8 s<br>1.46 GiB | 11.9 s<br>1.46 GiB |

### fateman-pow.ZZ

| n | Poly Mojo | FLINT | Singular |
|---|---|---|---|
| 10 | < 1 ms<br>< 1 MiB | < 1 ms<br>< 1 MiB | < 1 ms<br>< 1 MiB |
| 15 | < 1 ms<br>< 1 MiB | < 1 ms<br>< 1 MiB | < 1 ms<br>< 1 MiB |
| 20 | 1.63 ms<br>2.21 MiB | 1.48 ms<br>1.32 MiB | 4.37 ms<br>1.34 MiB |
| 30 | 8.83 ms<br>7.86 MiB | 7.29 ms<br>2.55 MiB | 75.1 ms<br>6.64 MiB |

### mp-div.ZZ

| n | Poly Mojo | FLINT | Singular |
|---|---|---|---|
| 12 | 803 ms<br>< 1 MiB | 1.18 s<br>< 1 MiB | 1.73 s<br>140 MiB |

### mp-mul.ZZ

| n | Poly Mojo | FLINT | Singular |
|---|---|---|---|
| 12 | 1.19 s<br>762 MiB | 1.22 s<br>188 MiB | 1.86 s<br>679 MiB |

### mp-pow.ZZ

| n | Poly Mojo | FLINT | Singular |
|---|---|---|---|
| 12 | < 1 ms<br>1.75 MiB | < 1 ms<br>1.07 MiB | 1.61 ms<br>< 1 MiB |

