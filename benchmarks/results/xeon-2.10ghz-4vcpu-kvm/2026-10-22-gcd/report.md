# Gcd, resultant, and discriminant results

Host: Intel(R) Xeon(R) Processor @ 2.10GHz, 4 vCPUs. One pinned CPU per run, 12G and 60 s per run.

Versions: Poly Mojo d6549d350302; Poly Mojo (modular) d6549d350302; FLINT python-flint 0.9.0; Singular 4.4.1; SymPy 1.14.0 (gmpy ground types).

Cells give the median time and incremental peak memory; ✗ marks an answer that disagrees with a confirmed one, ? an answer no other library produced. Ratios are competitor time ÷ Poly Mojo time (above 1: Poly Mojo faster) over cases both solved with a confirmed answer and at least 10 ms for the slower tool.

## Summary

| Ladder | Poly Mojo (modular) | FLINT | Singular | SymPy |
|---|---|---|---|---|
| mgcd.dense | 9.11 (5) | 0.722 (1) | 0.831 (1) | 20.4 (5) |
| mgcd.fateman-f1 | 88.7 (3) |  |  | 32.9 (1) |
| mgcd.fateman-f2 | 48.6 (4) |  |  | 26.4 (3) |
| mgcd.fateman-f3 | 2.86e+03 (2) |  |  | 20.3 (4) |
| mgcd.sparse | 528 (6) |  |  | 16.8 (4) |
| mres.2v |  |  |  | 22.3 (2) |
| mres.3v |  | 1.48 (3) | 5.97 (3) | 38.4 (3) |
| udisc.gf |  |  |  | 413 (2) |
| udisc.zz |  | 1.87 (1) | 24 (2) | 23.3 (2) |
| ugcd.gf |  |  |  | 2.27e+03 (5) |
| ugcd.zz.128 |  | 1.78 (1) | 1.81 (1) | 30.8 (3) |
| ugcd.zz.32 |  |  |  | 69.7 (3) |
| ugcd.zz.512 |  | 1.52 (2) | 1.88 (3) | 13.4 (4) |
| ures.gf |  |  |  | 400 (2) |
| ures.zz |  | 1.89 (1) | 24.8 (2) | 24.8 (2) |

## mgcd.dense

| Size | Poly Mojo | Poly Mojo (modular) | FLINT | Singular | SymPy |
|---|---|---|---|---|---|
| 5 | < 1 ms<br>< 1 MiB | 3.98 ms<br>< 1 MiB | < 1 ms<br>1.85 MiB | < 1 ms<br>2.4 MiB | 4.48 ms<br>< 1 MiB |
| 7 | 1.06 ms<br>< 1 MiB | 10.8 ms<br>1.2 MiB | < 1 ms<br>1.95 MiB | < 1 ms<br>2.54 MiB | 13.3 ms<br>< 1 MiB |
| 9 | 2.64 ms<br>< 1 MiB | 21.8 ms<br>1.9 MiB | 1.23 ms<br>1.95 MiB | 1.6 ms<br>2.47 MiB | 39.8 ms<br>< 1 MiB |
| 11 | 4.75 ms<br>1.03 MiB | 42.7 ms<br>2.66 MiB | 2.15 ms<br>2.02 MiB | 2.75 ms<br>3 MiB | 95.4 ms<br>1.16 MiB |
| 13 | 7.86 ms<br>1.95 MiB | 73.8 ms<br>4.59 MiB | 3.57 ms<br>1.96 MiB | 4.67 ms<br>3.38 MiB | 225 ms<br>1.72 MiB |
| 15 | 13.7 ms<br>2.58 MiB | 121 ms<br>5.83 MiB | 9.89 ms<br>2.12 MiB | 11.4 ms<br>3.91 MiB | 441 ms<br>2.39 MiB |

## mgcd.fateman-f1

| Size | Poly Mojo | Poly Mojo (modular) | FLINT | Singular | SymPy |
|---|---|---|---|---|---|
| 2 | < 1 ms<br>< 1 MiB | < 1 ms<br>< 1 MiB | < 1 ms<br>< 1 MiB | < 1 ms<br>1.54 MiB | < 1 ms<br>< 1 MiB |
| 3 | < 1 ms<br>< 1 MiB | 1.51 ms<br>< 1 MiB | < 1 ms<br>< 1 MiB | < 1 ms<br>1.54 MiB | < 1 ms<br>< 1 MiB |
| 4 | < 1 ms<br>< 1 MiB | 3.04 ms<br>2.02 MiB | < 1 ms<br>< 1 MiB | < 1 ms<br>1.54 MiB | 1.68 ms<br>< 1 MiB |
| 5 | < 1 ms<br>< 1 MiB | 6.05 ms<br>3.11 MiB | < 1 ms<br>< 1 MiB | < 1 ms<br>1.52 MiB | 3.39 ms<br>< 1 MiB |
| 6 | < 1 ms<br>< 1 MiB | 12.4 ms<br>5.31 MiB | < 1 ms<br>< 1 MiB | < 1 ms<br>1.54 MiB | 5.12 ms<br>< 1 MiB |
| 7 | < 1 ms<br>< 1 MiB | 23.8 ms<br>9.25 MiB | < 1 ms<br>< 1 MiB | < 1 ms<br>1.53 MiB | 8.43 ms<br>< 1 MiB |
| 8 | < 1 ms<br>< 1 MiB | 51.8 ms<br>16 MiB | < 1 ms<br>< 1 MiB | < 1 ms<br>1.54 MiB | 12.2 ms<br>< 1 MiB |

## mgcd.fateman-f2

| Size | Poly Mojo | Poly Mojo (modular) | FLINT | Singular | SymPy |
|---|---|---|---|---|---|
| 2 | < 1 ms<br>< 1 MiB | < 1 ms<br>< 1 MiB | < 1 ms<br>2.05 MiB | < 1 ms<br>3.05 MiB | < 1 ms<br>< 1 MiB |
| 3 | < 1 ms<br>< 1 MiB | 1.97 ms<br>< 1 MiB | < 1 ms<br>2.09 MiB | < 1 ms<br>3.09 MiB | 1.82 ms<br>< 1 MiB |
| 4 | < 1 ms<br>< 1 MiB | 5.88 ms<br>< 1 MiB | < 1 ms<br>2.26 MiB | < 1 ms<br>3.04 MiB | 4.13 ms<br>< 1 MiB |
| 5 | 1.28 ms<br>< 1 MiB | 15.8 ms<br>< 1 MiB | < 1 ms<br>2.37 MiB | 1 ms<br>3.22 MiB | 8.95 ms<br>< 1 MiB |
| 6 | 1.06 ms<br>< 1 MiB | 49.8 ms<br>1.46 MiB | 1.48 ms<br>2.39 MiB | 1.55 ms<br>3.3 MiB | 18.8 ms<br>< 1 MiB |
| 7 | 1.96 ms<br>2.03 MiB | 147 ms<br>3.98 MiB | 2.14 ms<br>2.43 MiB | 2.22 ms<br>3.32 MiB | 48.8 ms<br>1.21 MiB |
| 8 | 3.19 ms<br>2.02 MiB | 409 ms<br>3.8 MiB | 2.99 ms<br>2.44 MiB | 3.27 ms<br>3.38 MiB | 133 ms<br>3.52 MiB |

## mgcd.fateman-f3

| Size | Poly Mojo | Poly Mojo (modular) | FLINT | Singular | SymPy |
|---|---|---|---|---|---|
| 2 | < 1 ms<br>< 1 MiB | 2.34 ms<br>< 1 MiB | < 1 ms<br>2.07 MiB | < 1 ms<br>3.19 MiB | < 1 ms<br>< 1 MiB |
| 3 | < 1 ms<br>< 1 MiB | 22.8 ms<br>2.48 MiB | < 1 ms<br>2.14 MiB | < 1 ms<br>3.18 MiB | 2.04 ms<br>< 1 MiB |
| 4 | < 1 ms<br>< 1 MiB | 36.8 s<br>56.8 MiB | < 1 ms<br>2.25 MiB | < 1 ms<br>3.23 MiB | 4.95 ms<br>< 1 MiB |
| 5 | 1.36 ms<br>< 1 MiB | timeout | < 1 ms<br>2.36 MiB | 1.07 ms<br>3.19 MiB | 10.3 ms<br>< 1 MiB |
| 6 | 1.08 ms<br>< 1 MiB | · | 1.49 ms<br>2.43 MiB | 1.59 ms<br>3.36 MiB | 21.9 ms<br>< 1 MiB |
| 7 | 1.9 ms<br>1.84 MiB | · | 2.19 ms<br>2.41 MiB | 2.25 ms<br>3.39 MiB | 53.8 ms<br>1.4 MiB |
| 8 | 3.56 ms<br>2.38 MiB | · | 3.06 ms<br>2.48 MiB | 3.36 ms<br>3.43 MiB | 138 ms<br>3.79 MiB |

## mgcd.sparse

| Size | Poly Mojo | Poly Mojo (modular) | FLINT | Singular | SymPy |
|---|---|---|---|---|---|
| 3 | < 1 ms<br>< 1 MiB | 191 ms<br>3.44 MiB | < 1 ms<br>1.92 MiB | < 1 ms<br>2.51 MiB | 2.42 ms<br>< 1 MiB |
| 4 | < 1 ms<br>< 1 MiB | 5.95 s<br>21 MiB | 1.27 ms<br>2.15 MiB | 1.24 ms<br>3.19 MiB | 6.29 ms<br>< 1 MiB |
| 5 | < 1 ms<br>< 1 MiB | 12.3 s<br>65.9 MiB | 1.1 ms<br>2.39 MiB | 1.21 ms<br>3.51 MiB | 12.2 ms<br>< 1 MiB |
| 6 | 1.37 ms<br>< 1 MiB | 80.5 ms<br>27.3 MiB | < 1 ms<br>2.22 MiB | 1.04 ms<br>3.3 MiB | 21.4 ms<br>< 1 MiB |
| 7 | 1.47 ms<br>1.38 MiB | 103 ms<br>31.1 MiB | 1.59 ms<br>2.45 MiB | 1.55 ms<br>3.54 MiB | 27.7 ms<br>< 1 MiB |
| 8 | 2.45 ms<br>< 1 MiB | 295 ms<br>4.12 MiB | 1.44 ms<br>2.41 MiB | 1.55 ms<br>3.52 MiB | 53.3 ms<br>1.08 MiB |

## mres.2v

| Size | Poly Mojo | Poly Mojo (modular) | FLINT | Singular | SymPy |
|---|---|---|---|---|---|
| 2 | < 1 ms<br>< 1 MiB | n/a | < 1 ms<br>< 1 MiB | < 1 ms<br>< 1 MiB | < 1 ms<br>< 1 MiB |
| 4 | < 1 ms<br>< 1 MiB | n/a | < 1 ms<br>< 1 MiB | < 1 ms<br>2.16 MiB | < 1 ms<br>< 1 MiB |
| 6 | < 1 ms<br>< 1 MiB | n/a | < 1 ms<br>< 1 MiB | < 1 ms<br>2.38 MiB | 3.41 ms<br>< 1 MiB |
| 8 | < 1 ms<br>< 1 MiB | n/a | < 1 ms<br>1.03 MiB | 2.3 ms<br>2.59 MiB | 13.4 ms<br>< 1 MiB |
| 10 | 1.93 ms<br>< 1 MiB | n/a | 4.61 ms<br>1.07 MiB | 7.17 ms<br>2.75 MiB | 48.8 ms<br>< 1 MiB |

## mres.3v

| Size | Poly Mojo | Poly Mojo (modular) | FLINT | Singular | SymPy |
|---|---|---|---|---|---|
| 2 | < 1 ms<br>< 1 MiB | n/a | < 1 ms<br>< 1 MiB | < 1 ms<br>< 1 MiB | < 1 ms<br>< 1 MiB |
| 4 | 2.13 ms<br>< 1 MiB | n/a | < 1 ms<br>< 1 MiB | 1.34 ms<br>1.98 MiB | 29.4 ms<br>< 1 MiB |
| 6 | 16 ms<br>< 1 MiB | n/a | 6.64 ms<br>1.37 MiB | 62 ms<br>3.33 MiB | 725 ms<br>1.7 MiB |
| 8 | 100 ms<br>1.73 MiB | n/a | 237 ms<br>4.16 MiB | 668 ms<br>10.4 MiB | 9.03 s<br>7.7 MiB |
| 10 | 491 ms<br>2.81 MiB | n/a | 1.61 s<br>12.6 MiB | 4.03 s<br>35.4 MiB | timeout |

## udisc.gf

| Size | Poly Mojo | Poly Mojo (modular) | FLINT | Singular | SymPy |
|---|---|---|---|---|---|
| 10 | < 1 ms<br>< 1 MiB | n/a | timeout | < 1 ms<br>< 1 MiB | < 1 ms<br>< 1 MiB |
| 20 | < 1 ms<br>< 1 MiB | n/a | · | < 1 ms<br>< 1 MiB | 1 ms<br>< 1 MiB |
| 50 | < 1 ms<br>< 1 MiB | n/a | · | < 1 ms<br>< 1 MiB | 5.25 ms<br>< 1 MiB |
| 100 | < 1 ms<br>< 1 MiB | n/a | · | < 1 ms<br>< 1 MiB | 20.4 ms<br>< 1 MiB |
| 200 | < 1 ms<br>< 1 MiB | n/a | · | 3.15 ms<br>1 MiB | 83.6 ms<br>2.79 MiB |

## udisc.zz

| Size | Poly Mojo | Poly Mojo (modular) | FLINT | Singular | SymPy |
|---|---|---|---|---|---|
| 10 | < 1 ms<br>< 1 MiB | n/a | < 1 ms<br>< 1 MiB | < 1 ms<br>< 1 MiB | < 1 ms<br>< 1 MiB |
| 20 | < 1 ms<br>< 1 MiB | n/a | < 1 ms<br>< 1 MiB | < 1 ms<br>< 1 MiB | < 1 ms<br>< 1 MiB |
| 50 | < 1 ms<br>< 1 MiB | n/a | 1.18 ms<br>1.46 MiB | 3.83 ms<br>< 1 MiB | 4.79 ms<br>< 1 MiB |
| 100 | 2.92 ms<br>< 1 MiB | n/a | 5.63 ms<br>1.78 MiB | 44.2 ms<br>2.79 MiB | 44.4 ms<br>2.21 MiB |
| 200 | 13.8 ms<br>< 1 MiB | n/a | 25.8 ms<br>1.79 MiB | 522 ms<br>14.6 MiB | 493 ms<br>13.8 MiB |

## ugcd.gf

| Size | Poly Mojo | Poly Mojo (modular) | FLINT | Singular | SymPy |
|---|---|---|---|---|---|
| 50 | < 1 ms<br>< 1 MiB | n/a | < 1 ms<br>< 1 MiB | < 1 ms<br>1.21 MiB | 5.49 ms<br>< 1 MiB |
| 100 | < 1 ms<br>< 1 MiB | n/a | < 1 ms<br>< 1 MiB | < 1 ms<br>1.19 MiB | 20.5 ms<br>< 1 MiB |
| 200 | < 1 ms<br>< 1 MiB | n/a | < 1 ms<br>< 1 MiB | < 1 ms<br>1.23 MiB | 83.8 ms<br>2.07 MiB |
| 500 | < 1 ms<br>< 1 MiB | n/a | < 1 ms<br>< 1 MiB | < 1 ms<br>2.08 MiB | 543 ms<br>13.7 MiB |
| 1000 | < 1 ms<br>< 1 MiB | n/a | < 1 ms<br>< 1 MiB | < 1 ms<br>2.19 MiB | 2.31 s<br>55 MiB |
| 2000 | 1.11 ms<br>< 1 MiB | n/a | 1.45 ms<br>1.37 MiB | 1.61 ms<br>2.46 MiB | 9.58 s<br>220 MiB |

## ugcd.zz.128

| Size | Poly Mojo | Poly Mojo (modular) | FLINT | Singular | SymPy |
|---|---|---|---|---|---|
| 50 | < 1 ms<br>< 1 MiB | n/a | < 1 ms<br>< 1 MiB | < 1 ms<br>2.04 MiB | < 1 ms<br>< 1 MiB |
| 100 | < 1 ms<br>< 1 MiB | n/a | < 1 ms<br>< 1 MiB | < 1 ms<br>2.13 MiB | 1.33 ms<br>< 1 MiB |
| 200 | < 1 ms<br>< 1 MiB | n/a | < 1 ms<br>< 1 MiB | < 1 ms<br>2.2 MiB | 5.11 ms<br>< 1 MiB |
| 500 | 1.48 ms<br>< 1 MiB | n/a | 2.3 ms<br>< 1 MiB | 3.55 ms<br>2.98 MiB | 33.8 ms<br>< 1 MiB |
| 1000 | 4.24 ms<br>< 1 MiB | n/a | 6.23 ms<br>1.21 MiB | 9.49 ms<br>3.49 MiB | 130 ms<br>< 1 MiB |
| 2000 | 13.6 ms<br>1.2 MiB | n/a | 24.2 ms<br>2.3 MiB | 24.6 ms<br>4.55 MiB | 571 ms<br>1.13 MiB |

## ugcd.zz.32

| Size | Poly Mojo | Poly Mojo (modular) | FLINT | Singular | SymPy |
|---|---|---|---|---|---|
| 50 | < 1 ms<br>< 1 MiB | n/a | < 1 ms<br>< 1 MiB | < 1 ms<br>1.93 MiB | < 1 ms<br>< 1 MiB |
| 100 | < 1 ms<br>< 1 MiB | n/a | < 1 ms<br>< 1 MiB | < 1 ms<br>1.89 MiB | 1.06 ms<br>< 1 MiB |
| 200 | < 1 ms<br>< 1 MiB | n/a | < 1 ms<br>< 1 MiB | < 1 ms<br>1.92 MiB | 4.16 ms<br>< 1 MiB |
| 500 | < 1 ms<br>< 1 MiB | n/a | 1 ms<br>< 1 MiB | 1.62 ms<br>2.75 MiB | 28.6 ms<br>< 1 MiB |
| 1000 | 1.56 ms<br>< 1 MiB | n/a | 2.61 ms<br>1.04 MiB | 3.94 ms<br>2.95 MiB | 116 ms<br>< 1 MiB |
| 2000 | 4.82 ms<br>< 1 MiB | n/a | 9.07 ms<br>1.78 MiB | 9.79 ms<br>3.71 MiB | 447 ms<br>< 1 MiB |

## ugcd.zz.512

| Size | Poly Mojo | Poly Mojo (modular) | FLINT | Singular | SymPy |
|---|---|---|---|---|---|
| 50 | < 1 ms<br>< 1 MiB | n/a | < 1 ms<br>1.04 MiB | < 1 ms<br>2.24 MiB | 1.05 ms<br>< 1 MiB |
| 100 | < 1 ms<br>< 1 MiB | n/a | < 1 ms<br>1.05 MiB | < 1 ms<br>2.3 MiB | 3.36 ms<br>< 1 MiB |
| 200 | 1.83 ms<br>< 1 MiB | n/a | 2.1 ms<br>1.16 MiB | 2.23 ms<br>2.47 MiB | 11.6 ms<br>< 1 MiB |
| 500 | 5.15 ms<br>< 1 MiB | n/a | 7.07 ms<br>1.4 MiB | 10.4 ms<br>3.73 MiB | 67.7 ms<br>< 1 MiB |
| 1000 | 14.2 ms<br>< 1 MiB | n/a | 19.2 ms<br>2.04 MiB | 27.5 ms<br>4.7 MiB | 256 ms<br>1.22 MiB |
| 2000 | 43.6 ms<br>1.72 MiB | n/a | 74.3 ms<br>3.7 MiB | 73.5 ms<br>7.13 MiB | 943 ms<br>2.13 MiB |

## ures.gf

| Size | Poly Mojo | Poly Mojo (modular) | FLINT | Singular | SymPy |
|---|---|---|---|---|---|
| 10 | < 1 ms<br>< 1 MiB | n/a | timeout | < 1 ms<br>< 1 MiB | < 1 ms<br>< 1 MiB |
| 20 | < 1 ms<br>< 1 MiB | n/a | · | < 1 ms<br>< 1 MiB | < 1 ms<br>< 1 MiB |
| 50 | < 1 ms<br>< 1 MiB | n/a | · | < 1 ms<br>< 1 MiB | 5.12 ms<br>< 1 MiB |
| 100 | < 1 ms<br>< 1 MiB | n/a | · | < 1 ms<br>< 1 MiB | 19.7 ms<br>< 1 MiB |
| 200 | < 1 ms<br>< 1 MiB | n/a | · | 3.1 ms<br>< 1 MiB | 81 ms<br>2.86 MiB |

## ures.zz

| Size | Poly Mojo | Poly Mojo (modular) | FLINT | Singular | SymPy |
|---|---|---|---|---|---|
| 10 | < 1 ms<br>< 1 MiB | n/a | < 1 ms<br>< 1 MiB | < 1 ms<br>< 1 MiB | < 1 ms<br>< 1 MiB |
| 20 | < 1 ms<br>< 1 MiB | n/a | < 1 ms<br>< 1 MiB | < 1 ms<br>< 1 MiB | < 1 ms<br>< 1 MiB |
| 50 | < 1 ms<br>< 1 MiB | n/a | 1.12 ms<br>1.59 MiB | 3.86 ms<br>< 1 MiB | 4.75 ms<br>< 1 MiB |
| 100 | 2.59 ms<br>< 1 MiB | n/a | 5.2 ms<br>1.79 MiB | 41.1 ms<br>2.71 MiB | 43.4 ms<br>2.18 MiB |
| 200 | 12.7 ms<br>< 1 MiB | n/a | 24.1 ms<br>1.82 MiB | 493 ms<br>13.8 MiB | 467 ms<br>13.3 MiB |
