# Arithmetic results

Host: Intel(R) Xeon(R) Processor @ 2.10GHz, 4 vCPUs, 16 GiB RAM. Runs are pinned to one CPU, capped at 3G and 600 s.

Versions: Poly Mojo d6549d350302; FLINT python-flint 0.9.0, FLINT 3.6.0; Singular 4.4.1.

## Summary

Geometric mean of competitor time ÷ Poly Mojo time over cases both solved with a confirmed answer and where the slower took at least 10 ms (above 1 means Poly Mojo was faster).

| Ladder | FLINT | Singular |
|---|---|---|
| fateman-div.ZZ | 4.55 (*n* = 3) | 4.56 (*n* = 3) |
| fateman-mul.ZZ | 1.25 (*n* = 3) | 1.32 (*n* = 3) |
| fateman-pow.ZZ | — | 8.46 (*n* = 1) |
| mp-div.ZZ | 1.46 (*n* = 1) | 2.23 (*n* = 1) |
| mp-mul.ZZ | 1.03 (*n* = 1) | 1.58 (*n* = 1) |
| mp-pow.ZZ | — | — |

## Answer agreement

Every answer that two or more libraries computed agrees.

## Per-ladder results

Each cell gives the median time per call and the incremental peak memory; `·` marks an instance skipped after a smaller one on the same ladder hit a limit; `n/a` marks an unsupported field.

### fateman-div.ZZ

| n | Poly Mojo | FLINT | Singular |
|---|---|---|---|
| 10 | 2.49 ms<br>< 1 MiB | 8.54 ms<br>< 1 MiB | 9.28 ms<br>< 1 MiB |
| 15 | 27.9 ms<br>< 1 MiB | 127 ms<br>< 1 MiB | 129 ms<br>1.3 MiB |
| 20 | 211 ms<br>< 1 MiB | 1.11 s<br>< 1 MiB | 1.07 s<br>7.32 MiB |
| 30 | 5.59 s<br>4.3 MiB | 22 s<br>2.7 MiB | 22.6 s<br>27.9 MiB |

### fateman-mul.ZZ

| n | Poly Mojo | FLINT | Singular |
|---|---|---|---|
| 10 | 1.47 ms<br>1.67 MiB | < 1 ms<br>< 1 MiB | 1.73 ms<br>1.72 MiB |
| 15 | 14.5 ms<br>7.09 MiB | 11.8 ms<br>1.64 MiB | 15.7 ms<br>5.36 MiB |
| 20 | 167 ms<br>25.6 MiB | 128 ms<br>7.52 MiB | 136 ms<br>18.1 MiB |
| 30 | 4.28 s<br>126 MiB | 13.2 s<br>1.46 GiB | 11.1 s<br>1.46 GiB |

### fateman-pow.ZZ

| n | Poly Mojo | FLINT | Singular |
|---|---|---|---|
| 10 | < 1 ms<br>< 1 MiB | < 1 ms<br>< 1 MiB | < 1 ms<br>< 1 MiB |
| 15 | < 1 ms<br>< 1 MiB | < 1 ms<br>< 1 MiB | < 1 ms<br>< 1 MiB |
| 20 | 1.68 ms<br>1.76 MiB | 1.43 ms<br>1.22 MiB | 4.35 ms<br>1.35 MiB |
| 30 | 8.74 ms<br>7.87 MiB | 6.94 ms<br>2.5 MiB | 73.9 ms<br>6.64 MiB |

### mp-div.ZZ

| n | Poly Mojo | FLINT | Singular |
|---|---|---|---|
| 12 | 772 ms<br>< 1 MiB | 1.12 s<br>< 1 MiB | 1.72 s<br>140 MiB |

### mp-mul.ZZ

| n | Poly Mojo | FLINT | Singular |
|---|---|---|---|
| 12 | 1.23 s<br>762 MiB | 1.26 s<br>188 MiB | 1.95 s<br>679 MiB |

### mp-pow.ZZ

| n | Poly Mojo | FLINT | Singular |
|---|---|---|---|
| 12 | 1.22 ms<br>1.75 MiB | < 1 ms<br>< 1 MiB | 1.64 ms<br>< 1 MiB |

