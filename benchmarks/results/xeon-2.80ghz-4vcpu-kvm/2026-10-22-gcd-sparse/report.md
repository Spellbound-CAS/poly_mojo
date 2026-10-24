# Gcd, resultant, and discriminant results

Host: Intel(R) Xeon(R) Processor @ 2.80GHz, 4 vCPUs. One pinned CPU per run, 12G and 60 s per run.

Versions: Poly Mojo da10b2c8da4f; Poly Mojo (modular) da10b2c8da4f; FLINT python-flint 0.9.0; Singular 4.4.1; SymPy 1.14.0 (gmpy ground types).

Cells give the median time and incremental peak memory; ✗ marks an answer that disagrees with a confirmed one, ? an answer no other library produced. Ratios are competitor time ÷ Poly Mojo time (above 1: Poly Mojo faster) over cases both solved with a confirmed answer and at least 10 ms for the slower tool.

## Summary

| Ladder | Poly Mojo (modular) | FLINT | Singular | SymPy |
|---|---|---|---|---|
| mgcd.dense | 8.95 (5) | 0.564 (2) | 0.675 (2) | 18.8 (5) |
| mgcd.fateman-f1 | 90.5 (3) |  |  | 36.5 (1) |
| mgcd.fateman-f2 | 51.3 (4) |  |  | 20 (4) |
| mgcd.fateman-f3 | 2.92e+03 (2) |  |  | 20.6 (4) |
| mgcd.sparse | 558 (6) |  |  | 17.3 (4) |
| mres.2v |  |  |  | 21.2 (2) |
| mres.3v |  | 1.27 (3) | 5 (3) | 39.6 (3) |

## mgcd.dense

| Size | Poly Mojo | Poly Mojo (modular) | FLINT | Singular | SymPy |
|---|---|---|---|---|---|
| 5 | < 1 ms<br>< 1 MiB | 4.71 ms<br>< 1 MiB | < 1 ms<br>1.82 MiB | < 1 ms<br>2.44 MiB | 5.18 ms<br>< 1 MiB |
| 7 | 1.23 ms<br>< 1 MiB | 12.6 ms<br>< 1 MiB | < 1 ms<br>1.86 MiB | < 1 ms<br>2.59 MiB | 16.1 ms<br>< 1 MiB |
| 9 | 3.22 ms<br>< 1 MiB | 27.4 ms<br>1.95 MiB | 1.61 ms<br>1.98 MiB | 1.94 ms<br>2.47 MiB | 43.7 ms<br>< 1 MiB |
| 11 | 5.87 ms<br>1.12 MiB | 51.3 ms<br>2.78 MiB | 2.84 ms<br>2.01 MiB | 3.43 ms<br>3.06 MiB | 113 ms<br>1.1 MiB |
| 13 | 10.5 ms<br>2.4 MiB | 90.8 ms<br>4 MiB | 4.55 ms<br>1.91 MiB | 5.58 ms<br>3.38 MiB | 240 ms<br>1.71 MiB |
| 15 | 17.6 ms<br>1.9 MiB | 153 ms<br>5.68 MiB | 13 ms<br>2.23 MiB | 15.2 ms<br>3.91 MiB | 533 ms<br>2.32 MiB |

## mgcd.fateman-f1

| Size | Poly Mojo | Poly Mojo (modular) | FLINT | Singular | SymPy |
|---|---|---|---|---|---|
| 2 | < 1 ms<br>< 1 MiB | < 1 ms<br>< 1 MiB | < 1 ms<br>< 1 MiB | < 1 ms<br>1.55 MiB | < 1 ms<br>< 1 MiB |
| 3 | < 1 ms<br>< 1 MiB | 1.76 ms<br>1.09 MiB | < 1 ms<br>< 1 MiB | < 1 ms<br>1.55 MiB | 1.03 ms<br>< 1 MiB |
| 4 | < 1 ms<br>< 1 MiB | 3.47 ms<br>1.77 MiB | < 1 ms<br>< 1 MiB | < 1 ms<br>1.55 MiB | 1.99 ms<br>< 1 MiB |
| 5 | < 1 ms<br>< 1 MiB | 7.1 ms<br>3.25 MiB | < 1 ms<br>< 1 MiB | < 1 ms<br>1.55 MiB | 3.57 ms<br>< 1 MiB |
| 6 | < 1 ms<br>< 1 MiB | 13.8 ms<br>5.36 MiB | < 1 ms<br>< 1 MiB | < 1 ms<br>1.54 MiB | 5.92 ms<br>< 1 MiB |
| 7 | < 1 ms<br>< 1 MiB | 29 ms<br>9.49 MiB | < 1 ms<br>< 1 MiB | < 1 ms<br>1.56 MiB | 9.55 ms<br>< 1 MiB |
| 8 | < 1 ms<br>< 1 MiB | 58 ms<br>16.3 MiB | < 1 ms<br>< 1 MiB | < 1 ms<br>1.56 MiB | 15 ms<br>< 1 MiB |

## mgcd.fateman-f2

| Size | Poly Mojo | Poly Mojo (modular) | FLINT | Singular | SymPy |
|---|---|---|---|---|---|
| 2 | < 1 ms<br>< 1 MiB | < 1 ms<br>< 1 MiB | < 1 ms<br>2.13 MiB | < 1 ms<br>3.14 MiB | < 1 ms<br>< 1 MiB |
| 3 | < 1 ms<br>< 1 MiB | 2.24 ms<br>< 1 MiB | < 1 ms<br>2.12 MiB | < 1 ms<br>3.07 MiB | 1.98 ms<br>< 1 MiB |
| 4 | < 1 ms<br>< 1 MiB | 6.39 ms<br>< 1 MiB | < 1 ms<br>2.31 MiB | < 1 ms<br>3.18 MiB | 4.97 ms<br>< 1 MiB |
| 5 | 1.48 ms<br>< 1 MiB | 18.6 ms<br>1.02 MiB | 1.11 ms<br>2.42 MiB | 1.19 ms<br>3.17 MiB | 10.6 ms<br>< 1 MiB |
| 6 | 1.17 ms<br>< 1 MiB | 54.7 ms<br>1.68 MiB | 1.71 ms<br>2.44 MiB | 1.76 ms<br>3.3 MiB | 23.3 ms<br>< 1 MiB |
| 7 | 2.07 ms<br>2.14 MiB | 170 ms<br>3.95 MiB | 2.36 ms<br>2.38 MiB | 2.45 ms<br>3.34 MiB | 55.9 ms<br>1.26 MiB |
| 8 | 3.55 ms<br>1.98 MiB | 507 ms<br>4.26 MiB | 3.69 ms<br>2.46 MiB | 3.79 ms<br>3.4 MiB | 147 ms<br>3.52 MiB |

## mgcd.fateman-f3

| Size | Poly Mojo | Poly Mojo (modular) | FLINT | Singular | SymPy |
|---|---|---|---|---|---|
| 2 | < 1 ms<br>< 1 MiB | 2.66 ms<br>< 1 MiB | < 1 ms<br>2.19 MiB | < 1 ms<br>3.27 MiB | < 1 ms<br>< 1 MiB |
| 3 | < 1 ms<br>< 1 MiB | 26.1 ms<br>2.39 MiB | < 1 ms<br>2.23 MiB | < 1 ms<br>3.29 MiB | 2.31 ms<br>< 1 MiB |
| 4 | < 1 ms<br>< 1 MiB | 44 s<br>57 MiB | < 1 ms<br>2.36 MiB | < 1 ms<br>3.29 MiB | 5.89 ms<br>< 1 MiB |
| 5 | 1.6 ms<br>< 1 MiB | timeout | 1.14 ms<br>2.37 MiB | 1.22 ms<br>3.12 MiB | 12.4 ms<br>< 1 MiB |
| 6 | 1.31 ms<br>1.09 MiB | · | 1.77 ms<br>2.38 MiB | 1.81 ms<br>3.37 MiB | 27.5 ms<br>< 1 MiB |
| 7 | 2.31 ms<br>1.61 MiB | · | 2.42 ms<br>2.38 MiB | 2.51 ms<br>3.29 MiB | 59.6 ms<br>1.4 MiB |
| 8 | 3.83 ms<br>1.99 MiB | · | 3.75 ms<br>2.48 MiB | 3.85 ms<br>3.47 MiB | 163 ms<br>3.85 MiB |

## mgcd.sparse

| Size | Poly Mojo | Poly Mojo (modular) | FLINT | Singular | SymPy |
|---|---|---|---|---|---|
| 3 | < 1 ms<br>< 1 MiB | 230 ms<br>3.44 MiB | < 1 ms<br>1.89 MiB | < 1 ms<br>2.58 MiB | 2.61 ms<br>< 1 MiB |
| 4 | < 1 ms<br>< 1 MiB | 6.78 s<br>21.3 MiB | 1.35 ms<br>2.19 MiB | 1.44 ms<br>3.28 MiB | 7.46 ms<br>< 1 MiB |
| 5 | 1.14 ms<br>< 1 MiB | 14 s<br>65.8 MiB | 1.21 ms<br>2.37 MiB | 1.31 ms<br>3.52 MiB | 13.1 ms<br>< 1 MiB |
| 6 | 1.49 ms<br>< 1 MiB | 100 ms<br>27.3 MiB | 1.11 ms<br>2.2 MiB | 1.17 ms<br>3.37 MiB | 25.5 ms<br>< 1 MiB |
| 7 | 1.74 ms<br>< 1 MiB | 131 ms<br>31 MiB | 1.68 ms<br>2.46 MiB | 1.78 ms<br>3.39 MiB | 32.6 ms<br>< 1 MiB |
| 8 | 2.51 ms<br>< 1 MiB | 345 ms<br>4.05 MiB | 1.57 ms<br>2.45 MiB | 1.69 ms<br>3.53 MiB | 61.2 ms<br>1 MiB |

## mres.2v

| Size | Poly Mojo | Poly Mojo (modular) | FLINT | Singular | SymPy |
|---|---|---|---|---|---|
| 2 | < 1 ms<br>< 1 MiB | n/a | < 1 ms<br>< 1 MiB | < 1 ms<br>< 1 MiB | < 1 ms<br>< 1 MiB |
| 4 | < 1 ms<br>< 1 MiB | n/a | < 1 ms<br>< 1 MiB | < 1 ms<br>2.23 MiB | < 1 ms<br>< 1 MiB |
| 6 | < 1 ms<br>< 1 MiB | n/a | < 1 ms<br>< 1 MiB | < 1 ms<br>2.42 MiB | 4.04 ms<br>< 1 MiB |
| 8 | < 1 ms<br>< 1 MiB | n/a | 1 ms<br>1.04 MiB | 2.31 ms<br>2.61 MiB | 17 ms<br>< 1 MiB |
| 10 | 2.4 ms<br>< 1 MiB | n/a | 4.54 ms<br>1.14 MiB | 6.64 ms<br>2.84 MiB | 54.7 ms<br>< 1 MiB |

## mres.3v

| Size | Poly Mojo | Poly Mojo (modular) | FLINT | Singular | SymPy |
|---|---|---|---|---|---|
| 2 | < 1 ms<br>< 1 MiB | n/a | < 1 ms<br>< 1 MiB | < 1 ms<br>< 1 MiB | < 1 ms<br>< 1 MiB |
| 4 | 2.65 ms<br>< 1 MiB | n/a | < 1 ms<br>< 1 MiB | 1.32 ms<br>2.05 MiB | 35.6 ms<br>< 1 MiB |
| 6 | 19.3 ms<br>< 1 MiB | n/a | 6.08 ms<br>1.43 MiB | 61.7 ms<br>3.34 MiB | 908 ms<br>1.76 MiB |
| 8 | 119 ms<br>1.77 MiB | n/a | 260 ms<br>4.14 MiB | 696 ms<br>10.5 MiB | 11.7 s<br>8 MiB |
| 10 | 609 ms<br>2.92 MiB | n/a | 1.81 s<br>12.6 MiB | 4.09 s<br>35.5 MiB | timeout |
