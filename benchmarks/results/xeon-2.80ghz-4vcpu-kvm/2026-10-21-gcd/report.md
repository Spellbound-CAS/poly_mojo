# Gcd, resultant, and discriminant results

Host: Intel(R) Xeon(R) Processor @ 2.80GHz, 4 vCPUs. One pinned CPU per run, 12G and 60 s per run.

Versions: Poly Mojo e88c776db002; Poly Mojo (modular) e88c776db002; FLINT python-flint 0.9.0; Singular 4.4.1; SymPy 1.14.0 (gmpy ground types).

Cells give the median time and incremental peak memory; ✗ marks an answer that disagrees with a confirmed one, ? an answer no other library produced. Ratios are competitor time ÷ Poly Mojo time (above 1: Poly Mojo faster) over cases both solved with a confirmed answer and at least 10 ms for the slower tool.

## Summary

| Ladder | Poly Mojo (modular) | FLINT | Singular | SymPy |
|---|---|---|---|---|
| mgcd.dense | 8.86 (5) | 0.57 (2) | 0.691 (2) | 18.9 (5) |
| mgcd.fateman-f1 | 91.3 (3) |  |  | 36.3 (1) |
| mgcd.fateman-f2 | 49.6 (4) |  |  | 20.1 (4) |
| mgcd.fateman-f3 | 2.93e+03 (2) |  |  | 20.5 (4) |
| mgcd.sparse | 552 (6) |  |  | 17.6 (4) |
| mres.2v |  |  |  | 21.7 (2) |
| mres.3v |  | 1.24 (3) | 4.93 (3) | 39.1 (3) |
| udisc.gf |  |  |  | 610 (2) |
| udisc.zz |  | 1.32 (1) | 15.8 (2) | 16.7 (2) |
| ugcd.gf |  |  |  | 2.87e+03 (5) |
| ugcd.zz.128 |  | 1.71 (1) | 1.89 (2) | 31.4 (3) |
| ugcd.zz.32 |  | 1.74 (1) | 1.91 (1) | 74 (3) |
| ugcd.zz.512 |  | 1.5 (2) | 1.72 (3) | 13 (4) |
| ures.gf |  |  |  | 587 (2) |
| ures.zz |  | 1.39 (1) | 16.7 (2) | 17.1 (2) |

## mgcd.dense

| Size | Poly Mojo | Poly Mojo (modular) | FLINT | Singular | SymPy |
|---|---|---|---|---|---|
| 5 | < 1 ms<br>< 1 MiB | 4.76 ms<br>< 1 MiB | < 1 ms<br>1.87 MiB | < 1 ms<br>2.41 MiB | 5.17 ms<br>< 1 MiB |
| 7 | 1.23 ms<br>< 1 MiB | 12.7 ms<br>1.27 MiB | < 1 ms<br>1.94 MiB | < 1 ms<br>2.58 MiB | 16 ms<br>< 1 MiB |
| 9 | 3.32 ms<br>1.12 MiB | 27.1 ms<br>1.94 MiB | 1.61 ms<br>1.96 MiB | 1.94 ms<br>2.47 MiB | 44.1 ms<br>< 1 MiB |
| 11 | 6.06 ms<br>1.16 MiB | 50 ms<br>3.22 MiB | 2.85 ms<br>1.93 MiB | 3.39 ms<br>3.06 MiB | 113 ms<br>1.18 MiB |
| 13 | 10.3 ms<br>2.14 MiB | 94.7 ms<br>3.81 MiB | 4.57 ms<br>1.9 MiB | 5.72 ms<br>3.31 MiB | 270 ms<br>1.8 MiB |
| 15 | 17.6 ms<br>1.91 MiB | 152 ms<br>5.41 MiB | 12.9 ms<br>2.21 MiB | 15.3 ms<br>3.91 MiB | 513 ms<br>2.5 MiB |

## mgcd.fateman-f1

| Size | Poly Mojo | Poly Mojo (modular) | FLINT | Singular | SymPy |
|---|---|---|---|---|---|
| 2 | < 1 ms<br>< 1 MiB | < 1 ms<br>< 1 MiB | < 1 ms<br>< 1 MiB | < 1 ms<br>1.55 MiB | < 1 ms<br>< 1 MiB |
| 3 | < 1 ms<br>< 1 MiB | 1.8 ms<br>1.16 MiB | < 1 ms<br>< 1 MiB | < 1 ms<br>1.55 MiB | 1.03 ms<br>< 1 MiB |
| 4 | < 1 ms<br>< 1 MiB | 3.56 ms<br>1.82 MiB | < 1 ms<br>< 1 MiB | < 1 ms<br>1.55 MiB | 1.99 ms<br>< 1 MiB |
| 5 | < 1 ms<br>< 1 MiB | 7.34 ms<br>3.46 MiB | < 1 ms<br>< 1 MiB | < 1 ms<br>1.55 MiB | 3.63 ms<br>< 1 MiB |
| 6 | < 1 ms<br>< 1 MiB | 14.7 ms<br>5.38 MiB | < 1 ms<br>< 1 MiB | < 1 ms<br>1.49 MiB | 6.08 ms<br>< 1 MiB |
| 7 | < 1 ms<br>< 1 MiB | 29.8 ms<br>9.51 MiB | < 1 ms<br>< 1 MiB | < 1 ms<br>1.55 MiB | 9.52 ms<br>< 1 MiB |
| 8 | < 1 ms<br>< 1 MiB | 58.8 ms<br>15.8 MiB | < 1 ms<br>< 1 MiB | < 1 ms<br>1.56 MiB | 15.1 ms<br>< 1 MiB |

## mgcd.fateman-f2

| Size | Poly Mojo | Poly Mojo (modular) | FLINT | Singular | SymPy |
|---|---|---|---|---|---|
| 2 | < 1 ms<br>< 1 MiB | < 1 ms<br>< 1 MiB | < 1 ms<br>2.13 MiB | < 1 ms<br>3.13 MiB | < 1 ms<br>< 1 MiB |
| 3 | < 1 ms<br>< 1 MiB | 2.28 ms<br>< 1 MiB | < 1 ms<br>2.18 MiB | < 1 ms<br>3.17 MiB | 2 ms<br>< 1 MiB |
| 4 | < 1 ms<br>< 1 MiB | 6.58 ms<br>1.01 MiB | < 1 ms<br>2.26 MiB | < 1 ms<br>3.19 MiB | 4.95 ms<br>< 1 MiB |
| 5 | 1.52 ms<br>< 1 MiB | 19.2 ms<br>< 1 MiB | 1.13 ms<br>2.41 MiB | 1.19 ms<br>3.27 MiB | 11.1 ms<br>< 1 MiB |
| 6 | 1.17 ms<br>< 1 MiB | 54.6 ms<br>1.6 MiB | 1.73 ms<br>2.39 MiB | 1.78 ms<br>3.29 MiB | 24 ms<br>< 1 MiB |
| 7 | 2.11 ms<br>1.52 MiB | 169 ms<br>3.84 MiB | 2.39 ms<br>2.38 MiB | 2.48 ms<br>3.32 MiB | 55.1 ms<br>1.16 MiB |
| 8 | 3.76 ms<br>2.24 MiB | 489 ms<br>4.25 MiB | 3.65 ms<br>2.48 MiB | 3.83 ms<br>3.43 MiB | 157 ms<br>3.51 MiB |

## mgcd.fateman-f3

| Size | Poly Mojo | Poly Mojo (modular) | FLINT | Singular | SymPy |
|---|---|---|---|---|---|
| 2 | < 1 ms<br>< 1 MiB | 2.82 ms<br>< 1 MiB | < 1 ms<br>2.13 MiB | < 1 ms<br>3.25 MiB | < 1 ms<br>< 1 MiB |
| 3 | < 1 ms<br>< 1 MiB | 27.1 ms<br>2.91 MiB | < 1 ms<br>2.12 MiB | < 1 ms<br>3.28 MiB | 2.32 ms<br>< 1 MiB |
| 4 | < 1 ms<br>< 1 MiB | 44.9 s<br>57 MiB | < 1 ms<br>2.34 MiB | < 1 ms<br>3.3 MiB | 5.89 ms<br>< 1 MiB |
| 5 | 1.61 ms<br>< 1 MiB | timeout | 1.17 ms<br>2.44 MiB | 1.21 ms<br>3.23 MiB | 12.9 ms<br>< 1 MiB |
| 6 | 1.36 ms<br>< 1 MiB | · | 1.77 ms<br>2.43 MiB | 1.82 ms<br>3.36 MiB | 27.7 ms<br>< 1 MiB |
| 7 | 2.34 ms<br>1.64 MiB | · | 2.42 ms<br>2.36 MiB | 2.6 ms<br>3.4 MiB | 62.1 ms<br>1.41 MiB |
| 8 | 3.93 ms<br>1.82 MiB | · | 3.85 ms<br>2.49 MiB | 4.1 ms<br>3.39 MiB | 161 ms<br>3.75 MiB |

## mgcd.sparse

| Size | Poly Mojo | Poly Mojo (modular) | FLINT | Singular | SymPy |
|---|---|---|---|---|---|
| 3 | < 1 ms<br>< 1 MiB | 223 ms<br>3.18 MiB | < 1 ms<br>1.89 MiB | < 1 ms<br>2.58 MiB | 2.63 ms<br>< 1 MiB |
| 4 | < 1 ms<br>< 1 MiB | 6.88 s<br>21 MiB | 1.37 ms<br>2.26 MiB | 1.45 ms<br>3.18 MiB | 8.14 ms<br>< 1 MiB |
| 5 | 1.16 ms<br>< 1 MiB | 14 s<br>66.1 MiB | 1.27 ms<br>2.39 MiB | 1.34 ms<br>3.52 MiB | 13.7 ms<br>< 1 MiB |
| 6 | 1.48 ms<br>< 1 MiB | 102 ms<br>27.5 MiB | 1.09 ms<br>2.2 MiB | 1.2 ms<br>3.22 MiB | 25.6 ms<br>< 1 MiB |
| 7 | 1.78 ms<br>< 1 MiB | 128 ms<br>31 MiB | 1.65 ms<br>2.46 MiB | 1.79 ms<br>3.51 MiB | 33.8 ms<br>< 1 MiB |
| 8 | 2.6 ms<br>< 1 MiB | 340 ms<br>3.72 MiB | 1.6 ms<br>2.35 MiB | 1.73 ms<br>3.53 MiB | 64.9 ms<br>1.01 MiB |

## mres.2v

| Size | Poly Mojo | Poly Mojo (modular) | FLINT | Singular | SymPy |
|---|---|---|---|---|---|
| 2 | < 1 ms<br>< 1 MiB | n/a | < 1 ms<br>< 1 MiB | < 1 ms<br>< 1 MiB | < 1 ms<br>< 1 MiB |
| 4 | < 1 ms<br>< 1 MiB | n/a | < 1 ms<br>< 1 MiB | < 1 ms<br>2.17 MiB | < 1 ms<br>< 1 MiB |
| 6 | < 1 ms<br>< 1 MiB | n/a | < 1 ms<br>< 1 MiB | < 1 ms<br>2.41 MiB | 4.14 ms<br>< 1 MiB |
| 8 | < 1 ms<br>< 1 MiB | n/a | 1.02 ms<br>1.04 MiB | 2.38 ms<br>2.62 MiB | 17.4 ms<br>< 1 MiB |
| 10 | 2.52 ms<br>< 1 MiB | n/a | 4.7 ms<br>1.14 MiB | 7 ms<br>2.83 MiB | 60.6 ms<br>< 1 MiB |

## mres.3v

| Size | Poly Mojo | Poly Mojo (modular) | FLINT | Singular | SymPy |
|---|---|---|---|---|---|
| 2 | < 1 ms<br>< 1 MiB | n/a | < 1 ms<br>< 1 MiB | < 1 ms<br>< 1 MiB | < 1 ms<br>< 1 MiB |
| 4 | 2.86 ms<br>< 1 MiB | n/a | < 1 ms<br>< 1 MiB | 1.36 ms<br>1.96 MiB | 39.2 ms<br>< 1 MiB |
| 6 | 20.6 ms<br>< 1 MiB | n/a | 6.22 ms<br>1.43 MiB | 62.3 ms<br>3.36 MiB | 927 ms<br>1.68 MiB |
| 8 | 125 ms<br>1.91 MiB | n/a | 261 ms<br>4.03 MiB | 708 ms<br>10.4 MiB | 12.1 s<br>7.57 MiB |
| 10 | 606 ms<br>2.93 MiB | n/a | 1.83 s<br>12.6 MiB | 4.24 s<br>35.8 MiB | timeout |

## udisc.gf

| Size | Poly Mojo | Poly Mojo (modular) | FLINT | Singular | SymPy |
|---|---|---|---|---|---|
| 10 | < 1 ms<br>< 1 MiB | n/a | timeout | < 1 ms<br>< 1 MiB | < 1 ms<br>< 1 MiB |
| 20 | < 1 ms<br>< 1 MiB | n/a | · | < 1 ms<br>< 1 MiB | 1.32 ms<br>< 1 MiB |
| 50 | < 1 ms<br>< 1 MiB | n/a | · | < 1 ms<br>< 1 MiB | 8.18 ms<br>< 1 MiB |
| 100 | < 1 ms<br>< 1 MiB | n/a | · | < 1 ms<br>< 1 MiB | 28.8 ms<br>< 1 MiB |
| 200 | < 1 ms<br>< 1 MiB | n/a | · | 3.53 ms<br>1 MiB | 129 ms<br>2.79 MiB |

## udisc.zz

| Size | Poly Mojo | Poly Mojo (modular) | FLINT | Singular | SymPy |
|---|---|---|---|---|---|
| 10 | < 1 ms<br>< 1 MiB | n/a | < 1 ms<br>< 1 MiB | < 1 ms<br>< 1 MiB | < 1 ms<br>< 1 MiB |
| 20 | < 1 ms<br>< 1 MiB | n/a | < 1 ms<br>< 1 MiB | < 1 ms<br>< 1 MiB | < 1 ms<br>< 1 MiB |
| 50 | 1.23 ms<br>< 1 MiB | n/a | 1.31 ms<br>1.52 MiB | 4.29 ms<br>< 1 MiB | 5.7 ms<br>< 1 MiB |
| 100 | 4.98 ms<br>< 1 MiB | n/a | 6.57 ms<br>1.85 MiB | 49.3 ms<br>2.8 MiB | 55.4 ms<br>2.27 MiB |
| 200 | 24.4 ms<br>< 1 MiB | n/a | 32.1 ms<br>1.86 MiB | 615 ms<br>14.6 MiB | 609 ms<br>13.8 MiB |

## ugcd.gf

| Size | Poly Mojo | Poly Mojo (modular) | FLINT | Singular | SymPy |
|---|---|---|---|---|---|
| 50 | < 1 ms<br>< 1 MiB | n/a | timeout | < 1 ms<br>1.22 MiB | 7.69 ms<br>< 1 MiB |
| 100 | < 1 ms<br>< 1 MiB | n/a | · | < 1 ms<br>1.23 MiB | 33.3 ms<br>< 1 MiB |
| 200 | < 1 ms<br>< 1 MiB | n/a | · | < 1 ms<br>1.19 MiB | 118 ms<br>2.16 MiB |
| 500 | < 1 ms<br>< 1 MiB | n/a | · | < 1 ms<br>2.09 MiB | 844 ms<br>13.6 MiB |
| 1000 | < 1 ms<br>< 1 MiB | n/a | · | < 1 ms<br>2.2 MiB | 3.41 s<br>55 MiB |
| 2000 | 1.51 ms<br>< 1 MiB | n/a | · | 2.05 ms<br>2.51 MiB | 13.4 s<br>220 MiB |

## ugcd.zz.128

| Size | Poly Mojo | Poly Mojo (modular) | FLINT | Singular | SymPy |
|---|---|---|---|---|---|
| 50 | < 1 ms<br>< 1 MiB | n/a | < 1 ms<br>< 1 MiB | < 1 ms<br>2.04 MiB | < 1 ms<br>< 1 MiB |
| 100 | < 1 ms<br>< 1 MiB | n/a | < 1 ms<br>< 1 MiB | < 1 ms<br>2.14 MiB | 1.63 ms<br>< 1 MiB |
| 200 | < 1 ms<br>< 1 MiB | n/a | < 1 ms<br>< 1 MiB | < 1 ms<br>2.11 MiB | 6.6 ms<br>< 1 MiB |
| 500 | 1.91 ms<br>< 1 MiB | n/a | 2.79 ms<br>< 1 MiB | 4.08 ms<br>2.98 MiB | 42.7 ms<br>< 1 MiB |
| 1000 | 5.27 ms<br>< 1 MiB | n/a | 7.77 ms<br>1.21 MiB | 10.9 ms<br>3.48 MiB | 171 ms<br>< 1 MiB |
| 2000 | 16.6 ms<br>< 1 MiB | n/a | 28.5 ms<br>2.3 MiB | 28.9 ms<br>4.54 MiB | 711 ms<br>1.06 MiB |

## ugcd.zz.32

| Size | Poly Mojo | Poly Mojo (modular) | FLINT | Singular | SymPy |
|---|---|---|---|---|---|
| 50 | < 1 ms<br>< 1 MiB | n/a | < 1 ms<br>< 1 MiB | < 1 ms<br>1.94 MiB | < 1 ms<br>< 1 MiB |
| 100 | < 1 ms<br>< 1 MiB | n/a | < 1 ms<br>< 1 MiB | < 1 ms<br>1.95 MiB | 1.33 ms<br>< 1 MiB |
| 200 | < 1 ms<br>< 1 MiB | n/a | < 1 ms<br>< 1 MiB | < 1 ms<br>1.93 MiB | 5.17 ms<br>< 1 MiB |
| 500 | < 1 ms<br>< 1 MiB | n/a | 1.14 ms<br>< 1 MiB | 1.85 ms<br>2.75 MiB | 37 ms<br>< 1 MiB |
| 1000 | 1.97 ms<br>< 1 MiB | n/a | 3.14 ms<br>1.04 MiB | 4.48 ms<br>2.98 MiB | 144 ms<br>< 1 MiB |
| 2000 | 5.86 ms<br>< 1 MiB | n/a | 10.2 ms<br>1.72 MiB | 11.2 ms<br>3.72 MiB | 621 ms<br>< 1 MiB |

## ugcd.zz.512

| Size | Poly Mojo | Poly Mojo (modular) | FLINT | Singular | SymPy |
|---|---|---|---|---|---|
| 50 | < 1 ms<br>< 1 MiB | n/a | < 1 ms<br>1.02 MiB | < 1 ms<br>2.23 MiB | 1.19 ms<br>< 1 MiB |
| 100 | 1.12 ms<br>< 1 MiB | n/a | < 1 ms<br>1.08 MiB | 1.05 ms<br>2.32 MiB | 3.89 ms<br>< 1 MiB |
| 200 | 2.4 ms<br>< 1 MiB | n/a | 2.4 ms<br>1.14 MiB | 2.61 ms<br>2.56 MiB | 14.3 ms<br>< 1 MiB |
| 500 | 6.63 ms<br>< 1 MiB | n/a | 8.67 ms<br>1.39 MiB | 12 ms<br>3.79 MiB | 79.8 ms<br>< 1 MiB |
| 1000 | 18 ms<br>< 1 MiB | n/a | 24.3 ms<br>2.05 MiB | 31.5 ms<br>4.79 MiB | 302 ms<br>1.15 MiB |
| 2000 | 53.8 ms<br>1.63 MiB | n/a | 89.9 ms<br>3.7 MiB | 87.3 ms<br>7.13 MiB | 1.28 s<br>2.14 MiB |

## ures.gf

| Size | Poly Mojo | Poly Mojo (modular) | FLINT | Singular | SymPy |
|---|---|---|---|---|---|
| 10 | < 1 ms<br>< 1 MiB | n/a | timeout | < 1 ms<br>< 1 MiB | < 1 ms<br>< 1 MiB |
| 20 | < 1 ms<br>< 1 MiB | n/a | · | < 1 ms<br>< 1 MiB | 1.39 ms<br>< 1 MiB |
| 50 | < 1 ms<br>< 1 MiB | n/a | · | < 1 ms<br>< 1 MiB | 7.35 ms<br>< 1 MiB |
| 100 | < 1 ms<br>< 1 MiB | n/a | · | < 1 ms<br>< 1 MiB | 29.5 ms<br>< 1 MiB |
| 200 | < 1 ms<br>< 1 MiB | n/a | · | 3.57 ms<br>< 1 MiB | 117 ms<br>2.8 MiB |

## ures.zz

| Size | Poly Mojo | Poly Mojo (modular) | FLINT | Singular | SymPy |
|---|---|---|---|---|---|
| 10 | < 1 ms<br>< 1 MiB | n/a | < 1 ms<br>< 1 MiB | < 1 ms<br>< 1 MiB | < 1 ms<br>< 1 MiB |
| 20 | < 1 ms<br>< 1 MiB | n/a | < 1 ms<br>< 1 MiB | < 1 ms<br>< 1 MiB | < 1 ms<br>< 1 MiB |
| 50 | 1.22 ms<br>< 1 MiB | n/a | 1.25 ms<br>1.65 MiB | 4.28 ms<br>< 1 MiB | 5.67 ms<br>< 1 MiB |
| 100 | 4.58 ms<br>< 1 MiB | n/a | 5.97 ms<br>1.85 MiB | 47.2 ms<br>2.71 MiB | 51.8 ms<br>2.22 MiB |
| 200 | 21.4 ms<br>< 1 MiB | n/a | 29.7 ms<br>1.91 MiB | 581 ms<br>13.9 MiB | 555 ms<br>13.3 MiB |
