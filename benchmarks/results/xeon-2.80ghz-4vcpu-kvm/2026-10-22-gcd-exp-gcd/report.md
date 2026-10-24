# Gcd, resultant, and discriminant results

Host: Intel(R) Xeon(R) Processor @ 2.80GHz, 4 vCPUs. One pinned CPU per run, 12G and 60 s per run.

Versions: Poly Mojo df9f7ad5f26b; Poly Mojo (modular) df9f7ad5f26b; FLINT python-flint 0.9.0; Singular 4.4.1; SymPy 1.14.0 (gmpy ground types).

Cells give the median time and incremental peak memory; ✗ marks an answer that disagrees with a confirmed one, ? an answer no other library produced. Ratios are competitor time ÷ Poly Mojo time (above 1: Poly Mojo faster) over cases both solved with a confirmed answer and at least 10 ms for the slower tool.

## Summary

| Ladder | Poly Mojo (modular) | FLINT | Singular | SymPy |
|---|---|---|---|---|
| mgcd.dense | 16.4 (6) | 1.04 (1) | 1.24 (1) | 25.7 (6) |
| mgcd.fateman-f1 | 335 (4) |  |  | 140 (3) |
| mgcd.fateman-f2 | 50.4 (5) |  |  | 22.3 (4) |
| mgcd.fateman-f3 | 111 (1) |  |  | 21.1 (4) |
| mgcd.sparse | 752 (6) |  |  | 21.1 (4) |
| mres.2v |  |  |  | 91.9 (2) |
| mres.3v |  | 15 (2) | 29.5 (3) | 271 (3) |
| udisc.gf |  |  |  | 645 (2) |
| udisc.zz |  | 1.24 (1) | 14.1 (2) | 14.1 (2) |
| ugcd.gf |  |  |  | 1.71e+03 (6) |
| ugcd.zz.128 |  | 1.63 (1) | 1.71 (2) | 31.3 (3) |
| ugcd.zz.32 |  | 1.62 (1) | 1.8 (1) | 77.8 (3) |
| ugcd.zz.512 |  | 1.34 (3) | 1.56 (3) | 10.2 (4) |
| ures.gf |  |  |  | 778 (2) |
| ures.zz |  | 1.12 (1) | 15 (2) | 15.4 (2) |

## mgcd.dense

| Size | Poly Mojo | Poly Mojo (modular) | FLINT | Singular | SymPy |
|---|---|---|---|---|---|
| 5 | < 1 ms<br>< 1 MiB | 12 ms<br>< 1 MiB | < 1 ms<br>1.95 MiB | < 1 ms<br>2.71 MiB | 12.1 ms<br>< 1 MiB |
| 7 | 1.54 ms<br>< 1 MiB | 22.9 ms<br>1.14 MiB | 1.05 ms<br>1.98 MiB | 1.34 ms<br>2.8 MiB | 22.1 ms<br>< 1 MiB |
| 9 | 2.82 ms<br>1.02 MiB | 43.3 ms<br>1.84 MiB | 2.01 ms<br>2.09 MiB | 2.58 ms<br>2.75 MiB | 61.9 ms<br>< 1 MiB |
| 11 | 5.27 ms<br>< 1 MiB | 87.6 ms<br>2.84 MiB | 3.78 ms<br>2.11 MiB | 4.73 ms<br>3.22 MiB | 149 ms<br>1.17 MiB |
| 13 | 9.45 ms<br>1.47 MiB | 155 ms<br>4.19 MiB | 5.84 ms<br>2.02 MiB | 7.85 ms<br>3.66 MiB | 370 ms<br>1.82 MiB |
| 15 | 16.3 ms<br>1.23 MiB | 276 ms<br>5.48 MiB | 16.9 ms<br>2.32 MiB | 20.2 ms<br>4.07 MiB | 729 ms<br>2.39 MiB |

## mgcd.fateman-f1

| Size | Poly Mojo | Poly Mojo (modular) | FLINT | Singular | SymPy |
|---|---|---|---|---|---|
| 2 | < 1 ms<br>< 1 MiB | 1.36 ms<br>< 1 MiB | < 1 ms<br>< 1 MiB | < 1 ms<br>1.62 MiB | < 1 ms<br>< 1 MiB |
| 3 | < 1 ms<br>< 1 MiB | 2.38 ms<br>< 1 MiB | < 1 ms<br>< 1 MiB | < 1 ms<br>1.62 MiB | 1.34 ms<br>< 1 MiB |
| 4 | < 1 ms<br>< 1 MiB | 4.89 ms<br>1.75 MiB | < 1 ms<br>< 1 MiB | < 1 ms<br>1.62 MiB | 2.76 ms<br>< 1 MiB |
| 5 | < 1 ms<br>< 1 MiB | 11.2 ms<br>3.32 MiB | < 1 ms<br>< 1 MiB | < 1 ms<br>1.62 MiB | 4.58 ms<br>< 1 MiB |
| 6 | < 1 ms<br>< 1 MiB | 21.5 ms<br>5.31 MiB | < 1 ms<br>< 1 MiB | < 1 ms<br>1.62 MiB | 10.3 ms<br>< 1 MiB |
| 7 | < 1 ms<br>< 1 MiB | 53.5 ms<br>9.44 MiB | < 1 ms<br>< 1 MiB | < 1 ms<br>1.62 MiB | 13.4 ms<br>< 1 MiB |
| 8 | < 1 ms<br>< 1 MiB | 97.4 ms<br>16.1 MiB | < 1 ms<br>< 1 MiB | < 1 ms<br>1.62 MiB | 19.7 ms<br>< 1 MiB |

## mgcd.fateman-f2

| Size | Poly Mojo | Poly Mojo (modular) | FLINT | Singular | SymPy |
|---|---|---|---|---|---|
| 2 | < 1 ms<br>< 1 MiB | 1.41 ms<br>< 1 MiB | < 1 ms<br>2.17 MiB | < 1 ms<br>3.35 MiB | 1.01 ms<br>< 1 MiB |
| 3 | < 1 ms<br>< 1 MiB | 3.57 ms<br>< 1 MiB | < 1 ms<br>2.19 MiB | < 1 ms<br>3.27 MiB | 2.72 ms<br>< 1 MiB |
| 4 | < 1 ms<br>< 1 MiB | 10.6 ms<br>< 1 MiB | 1 ms<br>2.32 MiB | 1.04 ms<br>3.4 MiB | 6.39 ms<br>< 1 MiB |
| 5 | 1.94 ms<br>< 1 MiB | 30.2 ms<br>< 1 MiB | 1.55 ms<br>2.46 MiB | 1.54 ms<br>3.43 MiB | 15.1 ms<br>< 1 MiB |
| 6 | 1.4 ms<br>1.18 MiB | 98.9 ms<br>1.31 MiB | 2.22 ms<br>2.48 MiB | 2.27 ms<br>3.36 MiB | 30.4 ms<br>< 1 MiB |
| 7 | 2.55 ms<br>2.14 MiB | 280 ms<br>3.53 MiB | 3.09 ms<br>2.42 MiB | 3.18 ms<br>3.39 MiB | 77.7 ms<br>1.23 MiB |
| 8 | 4.48 ms<br>2.11 MiB | 930 ms<br>4.32 MiB | 4.78 ms<br>2.52 MiB | 5.09 ms<br>3.58 MiB | 216 ms<br>3.52 MiB |

## mgcd.fateman-f3

| Size | Poly Mojo | Poly Mojo (modular) | FLINT | Singular | SymPy |
|---|---|---|---|---|---|
| 2 | < 1 ms<br>< 1 MiB | 4.38 ms<br>< 1 MiB | < 1 ms<br>2.16 MiB | < 1 ms<br>3.43 MiB | 1.13 ms<br>< 1 MiB |
| 3 | < 1 ms<br>< 1 MiB | 47.1 ms<br>2.68 MiB | < 1 ms<br>2.21 MiB | < 1 ms<br>3.33 MiB | 3.16 ms<br>< 1 MiB |
| 4 | < 1 ms<br>< 1 MiB | timeout | 1.04 ms<br>2.32 MiB | 1.2 ms<br>3.41 MiB | 8.44 ms<br>< 1 MiB |
| 5 | 2.11 ms<br>< 1 MiB | · | 1.58 ms<br>2.46 MiB | 1.57 ms<br>3.27 MiB | 17 ms<br>< 1 MiB |
| 6 | 1.71 ms<br>< 1 MiB | · | 2.3 ms<br>2.42 MiB | 2.29 ms<br>3.44 MiB | 37 ms<br>< 1 MiB |
| 7 | 3.06 ms<br>1.77 MiB | · | 3.25 ms<br>2.42 MiB | 3.29 ms<br>3.54 MiB | 80.2 ms<br>1.48 MiB |
| 8 | 4.89 ms<br>1.64 MiB | · | 4.88 ms<br>2.52 MiB | 5.23 ms<br>3.55 MiB | 211 ms<br>3.79 MiB |

## mgcd.sparse

| Size | Poly Mojo | Poly Mojo (modular) | FLINT | Singular | SymPy |
|---|---|---|---|---|---|
| 3 | < 1 ms<br>< 1 MiB | 283 ms<br>3.11 MiB | < 1 ms<br>1.93 MiB | < 1 ms<br>2.69 MiB | 3.34 ms<br>< 1 MiB |
| 4 | 1.12 ms<br>< 1 MiB | 9.32 s<br>21 MiB | 1.76 ms<br>2.2 MiB | 1.88 ms<br>3.33 MiB | 9.87 ms<br>< 1 MiB |
| 5 | 1.23 ms<br>< 1 MiB | 19.7 s<br>63.8 MiB | 1.58 ms<br>2.38 MiB | 1.66 ms<br>3.6 MiB | 17 ms<br>< 1 MiB |
| 6 | 1.49 ms<br>< 1 MiB | 168 ms<br>27.2 MiB | 1.45 ms<br>2.28 MiB | 1.58 ms<br>3.34 MiB | 32 ms<br>< 1 MiB |
| 7 | 1.83 ms<br>1 MiB | 183 ms<br>31.4 MiB | 2.11 ms<br>2.47 MiB | 2.17 ms<br>3.61 MiB | 40.1 ms<br>< 1 MiB |
| 8 | 2.71 ms<br>< 1 MiB | 655 ms<br>3.8 MiB | 2 ms<br>2.39 MiB | 2.18 ms<br>3.61 MiB | 83.5 ms<br>1.07 MiB |

## mres.2v

| Size | Poly Mojo | Poly Mojo (modular) | FLINT | Singular | SymPy |
|---|---|---|---|---|---|
| 2 | < 1 ms<br>< 1 MiB | n/a | < 1 ms<br>< 1 MiB | < 1 ms<br>< 1 MiB | < 1 ms<br>< 1 MiB |
| 4 | < 1 ms<br>< 1 MiB | n/a | < 1 ms<br>< 1 MiB | < 1 ms<br>2.21 MiB | 1.02 ms<br>< 1 MiB |
| 6 | < 1 ms<br>< 1 MiB | n/a | < 1 ms<br>< 1 MiB | < 1 ms<br>2.41 MiB | 5.42 ms<br>< 1 MiB |
| 8 | < 1 ms<br>< 1 MiB | n/a | 1.23 ms<br>1.04 MiB | 3 ms<br>2.64 MiB | 30.5 ms<br>< 1 MiB |
| 10 | < 1 ms<br>< 1 MiB | n/a | 5.95 ms<br>1.15 MiB | 8.49 ms<br>2.83 MiB | 75.1 ms<br>< 1 MiB |

## mres.3v

| Size | Poly Mojo | Poly Mojo (modular) | FLINT | Singular | SymPy |
|---|---|---|---|---|---|
| 2 | < 1 ms<br>< 1 MiB | n/a | < 1 ms<br>< 1 MiB | < 1 ms<br>< 1 MiB | < 1 ms<br>< 1 MiB |
| 4 | < 1 ms<br>< 1 MiB | n/a | < 1 ms<br>< 1 MiB | 1.82 ms<br>2.01 MiB | 45.5 ms<br>< 1 MiB |
| 6 | 3.99 ms<br>< 1 MiB | n/a | 7.79 ms<br>1.44 MiB | 78.8 ms<br>3.34 MiB | 1.31 s<br>1.7 MiB |
| 8 | 33.3 ms<br>1.17 MiB | n/a | 375 ms<br>4.18 MiB | 939 ms<br>10.4 MiB | 15.5 s<br>7.7 MiB |
| 10 | 133 ms<br>1.73 MiB | n/a | 2.66 s<br>12.6 MiB | 6.14 s<br>35.6 MiB | timeout |

## udisc.gf

| Size | Poly Mojo | Poly Mojo (modular) | FLINT | Singular | SymPy |
|---|---|---|---|---|---|
| 10 | < 1 ms<br>< 1 MiB | n/a | timeout | < 1 ms<br>< 1 MiB | < 1 ms<br>< 1 MiB |
| 20 | < 1 ms<br>< 1 MiB | n/a | · | < 1 ms<br>< 1 MiB | 1.59 ms<br>< 1 MiB |
| 50 | < 1 ms<br>< 1 MiB | n/a | · | < 1 ms<br>< 1 MiB | 8.82 ms<br>< 1 MiB |
| 100 | < 1 ms<br>< 1 MiB | n/a | · | 1.67 ms<br>< 1 MiB | 35.5 ms<br>< 1 MiB |
| 200 | < 1 ms<br>< 1 MiB | n/a | · | 6.01 ms<br>1 MiB | 150 ms<br>2.78 MiB |

## udisc.zz

| Size | Poly Mojo | Poly Mojo (modular) | FLINT | Singular | SymPy |
|---|---|---|---|---|---|
| 10 | < 1 ms<br>< 1 MiB | n/a | < 1 ms<br>< 1 MiB | < 1 ms<br>< 1 MiB | < 1 ms<br>< 1 MiB |
| 20 | < 1 ms<br>< 1 MiB | n/a | < 1 ms<br>< 1 MiB | < 1 ms<br>< 1 MiB | < 1 ms<br>< 1 MiB |
| 50 | 2.61 ms<br>< 1 MiB | n/a | 2.18 ms<br>1.52 MiB | 4.8 ms<br>< 1 MiB | 7.46 ms<br>< 1 MiB |
| 100 | 6.97 ms<br>< 1 MiB | n/a | 8.87 ms<br>1.84 MiB | 60 ms<br>2.8 MiB | 62.5 ms<br>2.27 MiB |
| 200 | 32.6 ms<br>< 1 MiB | n/a | 40.3 ms<br>1.88 MiB | 753 ms<br>14.6 MiB | 720 ms<br>13.9 MiB |

## ugcd.gf

| Size | Poly Mojo | Poly Mojo (modular) | FLINT | Singular | SymPy |
|---|---|---|---|---|---|
| 50 | < 1 ms<br>< 1 MiB | n/a | < 1 ms<br>< 1 MiB | < 1 ms<br>1.16 MiB | 10.8 ms<br>< 1 MiB |
| 100 | < 1 ms<br>< 1 MiB | n/a | < 1 ms<br>< 1 MiB | < 1 ms<br>1.23 MiB | 57.7 ms<br>< 1 MiB |
| 200 | < 1 ms<br>< 1 MiB | n/a | < 1 ms<br>< 1 MiB | < 1 ms<br>1.16 MiB | 148 ms<br>2.06 MiB |
| 500 | < 1 ms<br>< 1 MiB | n/a | < 1 ms<br>< 1 MiB | < 1 ms<br>2.09 MiB | 1.05 s<br>13.8 MiB |
| 1000 | < 1 ms<br>< 1 MiB | n/a | < 1 ms<br>< 1 MiB | 1.07 ms<br>2.18 MiB | 4.37 s<br>55.1 MiB |
| 2000 | 2.6 ms<br>< 1 MiB | n/a | 2.39 ms<br>1.37 MiB | 2.57 ms<br>2.58 MiB | 18.2 s<br>220 MiB |

## ugcd.zz.128

| Size | Poly Mojo | Poly Mojo (modular) | FLINT | Singular | SymPy |
|---|---|---|---|---|---|
| 50 | < 1 ms<br>< 1 MiB | n/a | < 1 ms<br>< 1 MiB | < 1 ms<br>2.04 MiB | < 1 ms<br>< 1 MiB |
| 100 | < 1 ms<br>< 1 MiB | n/a | < 1 ms<br>< 1 MiB | < 1 ms<br>2.14 MiB | 1.99 ms<br>< 1 MiB |
| 200 | 1.19 ms<br>< 1 MiB | n/a | < 1 ms<br>< 1 MiB | 1.07 ms<br>2.14 MiB | 7.87 ms<br>< 1 MiB |
| 500 | 2.81 ms<br>< 1 MiB | n/a | 3.15 ms<br>< 1 MiB | 5.11 ms<br>2.98 MiB | 62.6 ms<br>< 1 MiB |
| 1000 | 7.31 ms<br>< 1 MiB | n/a | 9.07 ms<br>1.22 MiB | 13.2 ms<br>3.46 MiB | 226 ms<br>< 1 MiB |
| 2000 | 21.4 ms<br>< 1 MiB | n/a | 34.9 ms<br>2.32 MiB | 34.6 ms<br>4.49 MiB | 950 ms<br>1.12 MiB |

## ugcd.zz.32

| Size | Poly Mojo | Poly Mojo (modular) | FLINT | Singular | SymPy |
|---|---|---|---|---|---|
| 50 | < 1 ms<br>< 1 MiB | n/a | < 1 ms<br>< 1 MiB | < 1 ms<br>2 MiB | < 1 ms<br>< 1 MiB |
| 100 | < 1 ms<br>< 1 MiB | n/a | < 1 ms<br>< 1 MiB | < 1 ms<br>1.96 MiB | 2.24 ms<br>< 1 MiB |
| 200 | < 1 ms<br>< 1 MiB | n/a | < 1 ms<br>< 1 MiB | < 1 ms<br>2 MiB | 6.55 ms<br>< 1 MiB |
| 500 | 1.13 ms<br>< 1 MiB | n/a | 1.35 ms<br>< 1 MiB | 2.19 ms<br>2.67 MiB | 50.8 ms<br>< 1 MiB |
| 1000 | 2.61 ms<br>< 1 MiB | n/a | 3.62 ms<br>1.04 MiB | 5.55 ms<br>3.02 MiB | 208 ms<br>< 1 MiB |
| 2000 | 7.69 ms<br>< 1 MiB | n/a | 12.5 ms<br>1.82 MiB | 13.9 ms<br>3.72 MiB | 1.01 s<br>< 1 MiB |

## ugcd.zz.512

| Size | Poly Mojo | Poly Mojo (modular) | FLINT | Singular | SymPy |
|---|---|---|---|---|---|
| 50 | < 1 ms<br>< 1 MiB | n/a | < 1 ms<br>1.06 MiB | < 1 ms<br>2.29 MiB | 1.41 ms<br>< 1 MiB |
| 100 | 2.26 ms<br>< 1 MiB | n/a | 1.26 ms<br>1.05 MiB | 1.35 ms<br>2.32 MiB | 4.73 ms<br>< 1 MiB |
| 200 | 5.47 ms<br>< 1 MiB | n/a | 2.88 ms<br>1.14 MiB | 3.49 ms<br>2.5 MiB | 16.3 ms<br>< 1 MiB |
| 500 | 10.3 ms<br>< 1 MiB | n/a | 11.3 ms<br>1.49 MiB | 16.1 ms<br>3.67 MiB | 103 ms<br>< 1 MiB |
| 1000 | 24.9 ms<br>< 1 MiB | n/a | 34.4 ms<br>2.06 MiB | 40.4 ms<br>4.73 MiB | 415 ms<br>1.16 MiB |
| 2000 | 75.4 ms<br>1.83 MiB | n/a | 118 ms<br>3.7 MiB | 113 ms<br>7.12 MiB | 1.63 s<br>2.14 MiB |

## ures.gf

| Size | Poly Mojo | Poly Mojo (modular) | FLINT | Singular | SymPy |
|---|---|---|---|---|---|
| 10 | < 1 ms<br>< 1 MiB | n/a | timeout | < 1 ms<br>< 1 MiB | < 1 ms<br>< 1 MiB |
| 20 | < 1 ms<br>< 1 MiB | n/a | · | < 1 ms<br>< 1 MiB | 1.64 ms<br>< 1 MiB |
| 50 | < 1 ms<br>< 1 MiB | n/a | · | < 1 ms<br>< 1 MiB | 9.32 ms<br>< 1 MiB |
| 100 | < 1 ms<br>< 1 MiB | n/a | · | 1.58 ms<br>< 1 MiB | 35.3 ms<br>< 1 MiB |
| 200 | < 1 ms<br>< 1 MiB | n/a | · | 6.12 ms<br>< 1 MiB | 172 ms<br>2.79 MiB |

## ures.zz

| Size | Poly Mojo | Poly Mojo (modular) | FLINT | Singular | SymPy |
|---|---|---|---|---|---|
| 10 | < 1 ms<br>< 1 MiB | n/a | < 1 ms<br>< 1 MiB | < 1 ms<br>< 1 MiB | < 1 ms<br>< 1 MiB |
| 20 | < 1 ms<br>< 1 MiB | n/a | < 1 ms<br>< 1 MiB | < 1 ms<br>< 1 MiB | < 1 ms<br>< 1 MiB |
| 50 | 2.04 ms<br>< 1 MiB | n/a | 1.63 ms<br>1.63 MiB | 4.87 ms<br>< 1 MiB | 6.78 ms<br>< 1 MiB |
| 100 | 6.55 ms<br>< 1 MiB | n/a | 9.42 ms<br>1.87 MiB | 61.7 ms<br>2.77 MiB | 64.2 ms<br>2.18 MiB |
| 200 | 31.8 ms<br>< 1 MiB | n/a | 35.8 ms<br>1.91 MiB | 764 ms<br>13.9 MiB | 771 ms<br>13.3 MiB |
