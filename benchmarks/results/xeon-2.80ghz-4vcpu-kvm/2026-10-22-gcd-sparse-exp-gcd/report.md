# Gcd, resultant, and discriminant results

Host: Intel(R) Xeon(R) Processor @ 2.80GHz, 4 vCPUs. One pinned CPU per run, 12G and 60 s per run.

Versions: Poly Mojo df9f7ad5f26b; Poly Mojo (modular) df9f7ad5f26b; FLINT python-flint 0.9.0; Singular 4.4.1; SymPy 1.14.0 (gmpy ground types).

Cells give the median time and incremental peak memory; ✗ marks an answer that disagrees with a confirmed one, ? an answer no other library produced. Ratios are competitor time ÷ Poly Mojo time (above 1: Poly Mojo faster) over cases both solved with a confirmed answer and at least 10 ms for the slower tool.

## Summary

| Ladder | Poly Mojo (modular) | FLINT | Singular | SymPy |
|---|---|---|---|---|
| mgcd.dense | 15.5 (5) | 1.05 (1) | 1.24 (1) | 25.3 (5) |
| mgcd.fateman-f1 | 315 (4) |  |  | 155 (2) |
| mgcd.fateman-f2 | 51.6 (5) |  |  | 23 (4) |
| mgcd.fateman-f3 | 101 (1) |  |  | 20 (4) |
| mgcd.sparse | 761 (6) |  |  | 18.6 (5) |
| mres.2v |  |  |  | 81.4 (2) |
| mres.3v |  | 15.6 (2) | 33.1 (3) | 300 (3) |

## mgcd.dense

| Size | Poly Mojo | Poly Mojo (modular) | FLINT | Singular | SymPy |
|---|---|---|---|---|---|
| 5 | < 1 ms<br>< 1 MiB | 7.67 ms<br>< 1 MiB | < 1 ms<br>1.95 MiB | < 1 ms<br>2.72 MiB | 7.02 ms<br>< 1 MiB |
| 7 | 1.56 ms<br>< 1 MiB | 20.8 ms<br>< 1 MiB | 1.07 ms<br>2.04 MiB | 1.28 ms<br>2.8 MiB | 21.1 ms<br>< 1 MiB |
| 9 | 3.13 ms<br>< 1 MiB | 45.6 ms<br>2.06 MiB | 2.04 ms<br>2.02 MiB | 2.64 ms<br>2.75 MiB | 55.9 ms<br>< 1 MiB |
| 11 | 5.48 ms<br>< 1 MiB | 94.5 ms<br>3.17 MiB | 3.59 ms<br>2.05 MiB | 4.5 ms<br>3.21 MiB | 153 ms<br>1.16 MiB |
| 13 | 9.72 ms<br>1.62 MiB | 161 ms<br>4.17 MiB | 5.84 ms<br>2.03 MiB | 8.31 ms<br>3.66 MiB | 354 ms<br>1.78 MiB |
| 15 | 16.8 ms<br>1.61 MiB | 267 ms<br>5.74 MiB | 17.7 ms<br>2.28 MiB | 20.9 ms<br>4.07 MiB | 716 ms<br>2.33 MiB |

## mgcd.fateman-f1

| Size | Poly Mojo | Poly Mojo (modular) | FLINT | Singular | SymPy |
|---|---|---|---|---|---|
| 2 | < 1 ms<br>< 1 MiB | 1.3 ms<br>< 1 MiB | < 1 ms<br>< 1 MiB | < 1 ms<br>1.62 MiB | < 1 ms<br>< 1 MiB |
| 3 | < 1 ms<br>< 1 MiB | 2.39 ms<br>< 1 MiB | < 1 ms<br>< 1 MiB | < 1 ms<br>1.62 MiB | 1.33 ms<br>< 1 MiB |
| 4 | < 1 ms<br>< 1 MiB | 4.91 ms<br>1.77 MiB | < 1 ms<br>< 1 MiB | < 1 ms<br>1.62 MiB | 2.59 ms<br>< 1 MiB |
| 5 | < 1 ms<br>< 1 MiB | 11 ms<br>3.46 MiB | < 1 ms<br>< 1 MiB | < 1 ms<br>1.62 MiB | 4.62 ms<br>< 1 MiB |
| 6 | < 1 ms<br>< 1 MiB | 21.2 ms<br>5.31 MiB | < 1 ms<br>< 1 MiB | < 1 ms<br>1.62 MiB | 7.86 ms<br>< 1 MiB |
| 7 | < 1 ms<br>< 1 MiB | 47.4 ms<br>9.4 MiB | < 1 ms<br>< 1 MiB | < 1 ms<br>1.62 MiB | 12.7 ms<br>< 1 MiB |
| 8 | < 1 ms<br>< 1 MiB | 89.5 ms<br>16.4 MiB | < 1 ms<br>< 1 MiB | < 1 ms<br>1.62 MiB | 18.9 ms<br>< 1 MiB |

## mgcd.fateman-f2

| Size | Poly Mojo | Poly Mojo (modular) | FLINT | Singular | SymPy |
|---|---|---|---|---|---|
| 2 | < 1 ms<br>< 1 MiB | 1.3 ms<br>< 1 MiB | < 1 ms<br>2.17 MiB | < 1 ms<br>3.36 MiB | < 1 ms<br>< 1 MiB |
| 3 | < 1 ms<br>< 1 MiB | 3.6 ms<br>< 1 MiB | < 1 ms<br>2.19 MiB | < 1 ms<br>3.35 MiB | 2.68 ms<br>< 1 MiB |
| 4 | < 1 ms<br>< 1 MiB | 11.4 ms<br>< 1 MiB | 1.01 ms<br>2.38 MiB | 1.11 ms<br>3.36 MiB | 7.67 ms<br>< 1 MiB |
| 5 | 2.18 ms<br>< 1 MiB | 34.4 ms<br>1.21 MiB | 2.39 ms<br>2.46 MiB | 2.07 ms<br>3.38 MiB | 18.3 ms<br>< 1 MiB |
| 6 | 1.41 ms<br>< 1 MiB | 94 ms<br>1.64 MiB | 2.16 ms<br>2.48 MiB | 2.22 ms<br>3.41 MiB | 30.7 ms<br>< 1 MiB |
| 7 | 2.54 ms<br>1.77 MiB | 298 ms<br>3.89 MiB | 4.05 ms<br>2.42 MiB | 3.33 ms<br>3.5 MiB | 89 ms<br>1.27 MiB |
| 8 | 4.43 ms<br>1.91 MiB | 957 ms<br>4.44 MiB | 4.93 ms<br>2.52 MiB | 4.99 ms<br>3.59 MiB | 194 ms<br>3.62 MiB |

## mgcd.fateman-f3

| Size | Poly Mojo | Poly Mojo (modular) | FLINT | Singular | SymPy |
|---|---|---|---|---|---|
| 2 | < 1 ms<br>< 1 MiB | 4.33 ms<br>< 1 MiB | < 1 ms<br>2.23 MiB | < 1 ms<br>3.41 MiB | 1.13 ms<br>< 1 MiB |
| 3 | < 1 ms<br>< 1 MiB | 44.9 ms<br>2.48 MiB | < 1 ms<br>2.19 MiB | < 1 ms<br>3.43 MiB | 2.92 ms<br>< 1 MiB |
| 4 | < 1 ms<br>< 1 MiB | timeout | 1.05 ms<br>2.39 MiB | 1.09 ms<br>3.45 MiB | 7.7 ms<br>< 1 MiB |
| 5 | 2.1 ms<br>< 1 MiB | · | 1.74 ms<br>2.41 MiB | 1.56 ms<br>3.47 MiB | 16.8 ms<br>< 1 MiB |
| 6 | 2.01 ms<br>< 1 MiB | · | 2.39 ms<br>2.48 MiB | 2.45 ms<br>3.52 MiB | 35.6 ms<br>< 1 MiB |
| 7 | 3.16 ms<br>1.54 MiB | · | 3.22 ms<br>2.45 MiB | 3.38 ms<br>3.55 MiB | 79.4 ms<br>1.51 MiB |
| 8 | 4.8 ms<br>1.88 MiB | · | 5.01 ms<br>2.52 MiB | 5.22 ms<br>3.56 MiB | 213 ms<br>3.79 MiB |

## mgcd.sparse

| Size | Poly Mojo | Poly Mojo (modular) | FLINT | Singular | SymPy |
|---|---|---|---|---|---|
| 3 | < 1 ms<br>< 1 MiB | 298 ms<br>3.2 MiB | < 1 ms<br>1.92 MiB | < 1 ms<br>2.64 MiB | 3.42 ms<br>< 1 MiB |
| 4 | 1.17 ms<br>< 1 MiB | 9.7 s<br>21.1 MiB | 1.75 ms<br>2.21 MiB | 1.91 ms<br>3.41 MiB | 10.4 ms<br>< 1 MiB |
| 5 | 1.27 ms<br>< 1 MiB | 19.7 s<br>65.9 MiB | 1.54 ms<br>2.41 MiB | 1.66 ms<br>3.61 MiB | 21.2 ms<br>< 1 MiB |
| 6 | 1.63 ms<br>< 1 MiB | 179 ms<br>27 MiB | 1.4 ms<br>2.2 MiB | 1.52 ms<br>3.43 MiB | 30.8 ms<br>< 1 MiB |
| 7 | 1.81 ms<br>1.14 MiB | 220 ms<br>31.3 MiB | 2.08 ms<br>2.52 MiB | 2.27 ms<br>3.62 MiB | 45.5 ms<br>< 1 MiB |
| 8 | 2.68 ms<br>< 1 MiB | 601 ms<br>3.96 MiB | 2.1 ms<br>2.45 MiB | 2.13 ms<br>3.56 MiB | 84 ms<br>1.02 MiB |

## mres.2v

| Size | Poly Mojo | Poly Mojo (modular) | FLINT | Singular | SymPy |
|---|---|---|---|---|---|
| 2 | < 1 ms<br>< 1 MiB | n/a | < 1 ms<br>< 1 MiB | < 1 ms<br>< 1 MiB | < 1 ms<br>< 1 MiB |
| 4 | < 1 ms<br>< 1 MiB | n/a | < 1 ms<br>< 1 MiB | < 1 ms<br>2.17 MiB | < 1 ms<br>< 1 MiB |
| 6 | < 1 ms<br>< 1 MiB | n/a | < 1 ms<br>< 1 MiB | < 1 ms<br>2.39 MiB | 5.29 ms<br>< 1 MiB |
| 8 | < 1 ms<br>< 1 MiB | n/a | 1.29 ms<br>1.03 MiB | 2.92 ms<br>2.63 MiB | 22.9 ms<br>< 1 MiB |
| 10 | < 1 ms<br>< 1 MiB | n/a | 5.86 ms<br>1.15 MiB | 8.83 ms<br>2.82 MiB | 77 ms<br>< 1 MiB |

## mres.3v

| Size | Poly Mojo | Poly Mojo (modular) | FLINT | Singular | SymPy |
|---|---|---|---|---|---|
| 2 | < 1 ms<br>< 1 MiB | n/a | < 1 ms<br>< 1 MiB | < 1 ms<br>< 1 MiB | < 1 ms<br>< 1 MiB |
| 4 | < 1 ms<br>< 1 MiB | n/a | < 1 ms<br>< 1 MiB | 1.71 ms<br>1.97 MiB | 46.5 ms<br>< 1 MiB |
| 6 | 4.02 ms<br>< 1 MiB | n/a | 8.02 ms<br>1.43 MiB | 88.5 ms<br>3.3 MiB | 1.31 s<br>1.71 MiB |
| 8 | 27.1 ms<br>1.18 MiB | n/a | 361 ms<br>4.09 MiB | 978 ms<br>10.5 MiB | 16.6 s<br>7.65 MiB |
| 10 | 131 ms<br>2.09 MiB | n/a | 2.41 s<br>12.6 MiB | 5.98 s<br>35.6 MiB | timeout |
