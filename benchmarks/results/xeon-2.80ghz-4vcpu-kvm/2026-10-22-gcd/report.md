# Gcd, resultant, and discriminant results

Host: Intel(R) Xeon(R) Processor @ 2.80GHz, 4 vCPUs. One pinned CPU per run, 12G and 60 s per run.

Versions: Poly Mojo 319a0ee35646; Poly Mojo (modular) 319a0ee35646; FLINT python-flint 0.9.0; Singular 4.4.1; SymPy 1.14.0 (gmpy ground types).

Cells give the median time and incremental peak memory; ✗ marks an answer that disagrees with a confirmed one, ? an answer no other library produced. Ratios are competitor time ÷ Poly Mojo time (above 1: Poly Mojo faster) over cases both solved with a confirmed answer and at least 10 ms for the slower tool.

## Summary

| Ladder | Poly Mojo (modular) | FLINT | Singular | SymPy |
|---|---|---|---|---|
| mgcd.dense | 8.98 (5) | 0.573 (2) | 0.689 (2) | 18.9 (5) |
| mgcd.fateman-f1 | 89.7 (3) |  |  | 35.7 (1) |
| mgcd.fateman-f2 | 51.2 (4) |  |  | 20.3 (4) |
| mgcd.fateman-f3 | 3e+03 (2) |  |  | 19.9 (4) |
| mgcd.sparse | 571 (6) |  |  | 17.5 (4) |
| mres.2v |  |  |  | 22.1 (2) |
| mres.3v |  | 1.28 (3) | 4.8 (3) | 38.9 (3) |
| udisc.gf |  |  |  | 612 (2) |
| udisc.zz |  | 1.33 (1) | 16.3 (2) | 16.2 (2) |
| ugcd.gf |  |  |  | 2.83e+03 (5) |
| ugcd.zz.128 |  | 1.76 (1) | 1.92 (2) | 32.2 (3) |
| ugcd.zz.32 |  | 1.8 (1) | 1.94 (1) | 77.1 (3) |
| ugcd.zz.512 |  | 1.5 (2) | 1.89 (3) | 13.8 (4) |
| ures.gf |  |  |  | 613 (2) |
| ures.zz |  | 1.34 (1) | 16.4 (2) | 17.8 (2) |

## mgcd.dense

| Size | Poly Mojo | Poly Mojo (modular) | FLINT | Singular | SymPy |
|---|---|---|---|---|---|
| 5 | < 1 ms<br>< 1 MiB | 4.79 ms<br>< 1 MiB | < 1 ms<br>2.02 MiB | < 1 ms<br>3.02 MiB | 5.15 ms<br>< 1 MiB |
| 7 | 1.22 ms<br>< 1 MiB | 12.1 ms<br>1.15 MiB | < 1 ms<br>2.12 MiB | < 1 ms<br>3.13 MiB | 16.1 ms<br>< 1 MiB |
| 9 | 3.16 ms<br>1.26 MiB | 27.7 ms<br>1.87 MiB | 1.61 ms<br>2.08 MiB | 1.94 ms<br>3.02 MiB | 43.9 ms<br>< 1 MiB |
| 11 | 5.9 ms<br>1.15 MiB | 51.8 ms<br>2.91 MiB | 2.81 ms<br>2.11 MiB | 3.43 ms<br>3.61 MiB | 113 ms<br>1.17 MiB |
| 13 | 10.4 ms<br>1.89 MiB | 89.2 ms<br>4.27 MiB | 4.56 ms<br>2.09 MiB | 5.57 ms<br>3.99 MiB | 243 ms<br>1.79 MiB |
| 15 | 17.2 ms<br>1.95 MiB | 153 ms<br>6.08 MiB | 12.9 ms<br>2.38 MiB | 15.2 ms<br>4.4 MiB | 506 ms<br>2.39 MiB |

## mgcd.fateman-f1

| Size | Poly Mojo | Poly Mojo (modular) | FLINT | Singular | SymPy |
|---|---|---|---|---|---|
| 2 | < 1 ms<br>< 1 MiB | < 1 ms<br>< 1 MiB | < 1 ms<br>< 1 MiB | < 1 ms<br>1.65 MiB | < 1 ms<br>< 1 MiB |
| 3 | < 1 ms<br>< 1 MiB | 1.77 ms<br>< 1 MiB | < 1 ms<br>< 1 MiB | < 1 ms<br>1.73 MiB | 1.03 ms<br>< 1 MiB |
| 4 | < 1 ms<br>< 1 MiB | 3.48 ms<br>1.8 MiB | < 1 ms<br>< 1 MiB | < 1 ms<br>1.73 MiB | 2 ms<br>< 1 MiB |
| 5 | < 1 ms<br>< 1 MiB | 7.31 ms<br>3.14 MiB | < 1 ms<br>< 1 MiB | < 1 ms<br>1.73 MiB | 3.5 ms<br>< 1 MiB |
| 6 | < 1 ms<br>< 1 MiB | 14.3 ms<br>5.82 MiB | < 1 ms<br>< 1 MiB | < 1 ms<br>1.73 MiB | 5.99 ms<br>< 1 MiB |
| 7 | < 1 ms<br>< 1 MiB | 28.3 ms<br>9.6 MiB | < 1 ms<br>< 1 MiB | < 1 ms<br>1.74 MiB | 9.81 ms<br>< 1 MiB |
| 8 | < 1 ms<br>< 1 MiB | 56.7 ms<br>16.2 MiB | < 1 ms<br>< 1 MiB | < 1 ms<br>1.68 MiB | 14.9 ms<br>< 1 MiB |

## mgcd.fateman-f2

| Size | Poly Mojo | Poly Mojo (modular) | FLINT | Singular | SymPy |
|---|---|---|---|---|---|
| 2 | < 1 ms<br>< 1 MiB | < 1 ms<br>< 1 MiB | < 1 ms<br>2.37 MiB | < 1 ms<br>3.68 MiB | < 1 ms<br>< 1 MiB |
| 3 | < 1 ms<br>< 1 MiB | 2.28 ms<br>< 1 MiB | < 1 ms<br>2.4 MiB | < 1 ms<br>3.71 MiB | 2.02 ms<br>< 1 MiB |
| 4 | < 1 ms<br>< 1 MiB | 6.64 ms<br>< 1 MiB | < 1 ms<br>2.59 MiB | < 1 ms<br>3.77 MiB | 5 ms<br>< 1 MiB |
| 5 | 1.51 ms<br>< 1 MiB | 19.4 ms<br>< 1 MiB | 1.14 ms<br>2.68 MiB | 1.2 ms<br>3.76 MiB | 11 ms<br>< 1 MiB |
| 6 | 1.17 ms<br>< 1 MiB | 54.8 ms<br>1.79 MiB | 1.75 ms<br>2.63 MiB | 1.77 ms<br>3.76 MiB | 24.4 ms<br>< 1 MiB |
| 7 | 2.1 ms<br>2.15 MiB | 167 ms<br>3.25 MiB | 2.4 ms<br>2.63 MiB | 2.46 ms<br>3.78 MiB | 55.5 ms<br>1.23 MiB |
| 8 | 3.58 ms<br>1.69 MiB | 510 ms<br>4.3 MiB | 3.7 ms<br>2.73 MiB | 3.82 ms<br>3.92 MiB | 151 ms<br>3.51 MiB |

## mgcd.fateman-f3

| Size | Poly Mojo | Poly Mojo (modular) | FLINT | Singular | SymPy |
|---|---|---|---|---|---|
| 2 | < 1 ms<br>< 1 MiB | 2.7 ms<br>< 1 MiB | < 1 ms<br>2.38 MiB | < 1 ms<br>3.76 MiB | < 1 ms<br>< 1 MiB |
| 3 | < 1 ms<br>< 1 MiB | 26.5 ms<br>2.45 MiB | < 1 ms<br>2.46 MiB | < 1 ms<br>3.76 MiB | 2.3 ms<br>< 1 MiB |
| 4 | < 1 ms<br>< 1 MiB | 45.2 s<br>57.1 MiB | < 1 ms<br>2.59 MiB | < 1 ms<br>3.83 MiB | 6 ms<br>< 1 MiB |
| 5 | 1.58 ms<br>< 1 MiB | timeout | 1.15 ms<br>2.62 MiB | 1.23 ms<br>3.89 MiB | 12.8 ms<br>< 1 MiB |
| 6 | 1.36 ms<br>< 1 MiB | · | 1.76 ms<br>2.65 MiB | 1.83 ms<br>3.9 MiB | 27.1 ms<br>< 1 MiB |
| 7 | 2.44 ms<br>2.03 MiB | · | 2.45 ms<br>2.6 MiB | 2.55 ms<br>3.92 MiB | 62.1 ms<br>1.39 MiB |
| 8 | 3.99 ms<br>1.8 MiB | · | 3.77 ms<br>2.73 MiB | 4.01 ms<br>3.93 MiB | 151 ms<br>3.78 MiB |

## mgcd.sparse

| Size | Poly Mojo | Poly Mojo (modular) | FLINT | Singular | SymPy |
|---|---|---|---|---|---|
| 3 | < 1 ms<br>< 1 MiB | 239 ms<br>3.45 MiB | < 1 ms<br>2.07 MiB | < 1 ms<br>2.98 MiB | 2.61 ms<br>< 1 MiB |
| 4 | 1.01 ms<br>< 1 MiB | 6.96 s<br>21.2 MiB | 1.36 ms<br>2.5 MiB | 1.46 ms<br>3.74 MiB | 7.58 ms<br>< 1 MiB |
| 5 | 1.15 ms<br>< 1 MiB | 14.5 s<br>65.8 MiB | 1.22 ms<br>2.68 MiB | 1.35 ms<br>3.95 MiB | 13.6 ms<br>< 1 MiB |
| 6 | 1.5 ms<br>< 1 MiB | 110 ms<br>27 MiB | 1.11 ms<br>2.45 MiB | 1.19 ms<br>3.67 MiB | 24.9 ms<br>< 1 MiB |
| 7 | 1.69 ms<br>1.04 MiB | 127 ms<br>31.2 MiB | 1.64 ms<br>2.74 MiB | 1.74 ms<br>3.96 MiB | 33.3 ms<br>< 1 MiB |
| 8 | 2.56 ms<br>< 1 MiB | 341 ms<br>3.71 MiB | 1.6 ms<br>2.62 MiB | 1.67 ms<br>3.92 MiB | 61.4 ms<br>1.02 MiB |

## mres.2v

| Size | Poly Mojo | Poly Mojo (modular) | FLINT | Singular | SymPy |
|---|---|---|---|---|---|
| 2 | < 1 ms<br>< 1 MiB | n/a | < 1 ms<br>< 1 MiB | < 1 ms<br>< 1 MiB | < 1 ms<br>< 1 MiB |
| 4 | < 1 ms<br>< 1 MiB | n/a | < 1 ms<br>< 1 MiB | < 1 ms<br>2.39 MiB | < 1 ms<br>< 1 MiB |
| 6 | < 1 ms<br>< 1 MiB | n/a | < 1 ms<br>< 1 MiB | < 1 ms<br>2.56 MiB | 4.1 ms<br>< 1 MiB |
| 8 | < 1 ms<br>< 1 MiB | n/a | 1.01 ms<br>1.08 MiB | 2.32 ms<br>2.82 MiB | 17.2 ms<br>< 1 MiB |
| 10 | 2.36 ms<br>< 1 MiB | n/a | 4.61 ms<br>1.18 MiB | 6.8 ms<br>2.98 MiB | 58.1 ms<br>< 1 MiB |

## mres.3v

| Size | Poly Mojo | Poly Mojo (modular) | FLINT | Singular | SymPy |
|---|---|---|---|---|---|
| 2 | < 1 ms<br>< 1 MiB | n/a | < 1 ms<br>< 1 MiB | < 1 ms<br>< 1 MiB | < 1 ms<br>< 1 MiB |
| 4 | 2.67 ms<br>< 1 MiB | n/a | < 1 ms<br>< 1 MiB | 1.3 ms<br>2.16 MiB | 35.1 ms<br>< 1 MiB |
| 6 | 19.9 ms<br>< 1 MiB | n/a | 6.2 ms<br>1.5 MiB | 63.7 ms<br>3.47 MiB | 910 ms<br>1.73 MiB |
| 8 | 123 ms<br>1.86 MiB | n/a | 271 ms<br>4.24 MiB | 678 ms<br>10.7 MiB | 12.1 s<br>7.87 MiB |
| 10 | 620 ms<br>2.48 MiB | n/a | 1.89 s<br>12.8 MiB | 3.9 s<br>35.7 MiB | timeout |

## udisc.gf

| Size | Poly Mojo | Poly Mojo (modular) | FLINT | Singular | SymPy |
|---|---|---|---|---|---|
| 10 | < 1 ms<br>< 1 MiB | n/a | timeout | < 1 ms<br>< 1 MiB | < 1 ms<br>< 1 MiB |
| 20 | < 1 ms<br>< 1 MiB | n/a | · | < 1 ms<br>< 1 MiB | 1.32 ms<br>< 1 MiB |
| 50 | < 1 ms<br>< 1 MiB | n/a | · | < 1 ms<br>< 1 MiB | 7.5 ms<br>< 1 MiB |
| 100 | < 1 ms<br>< 1 MiB | n/a | · | < 1 ms<br>< 1 MiB | 28 ms<br>< 1 MiB |
| 200 | < 1 ms<br>< 1 MiB | n/a | · | 3.77 ms<br>1.07 MiB | 133 ms<br>2.85 MiB |

## udisc.zz

| Size | Poly Mojo | Poly Mojo (modular) | FLINT | Singular | SymPy |
|---|---|---|---|---|---|
| 10 | < 1 ms<br>< 1 MiB | n/a | < 1 ms<br>< 1 MiB | < 1 ms<br>< 1 MiB | < 1 ms<br>< 1 MiB |
| 20 | < 1 ms<br>< 1 MiB | n/a | < 1 ms<br>1.07 MiB | < 1 ms<br>< 1 MiB | < 1 ms<br>< 1 MiB |
| 50 | 1.17 ms<br>< 1 MiB | n/a | 1.35 ms<br>1.61 MiB | 4.29 ms<br>1.04 MiB | 5.75 ms<br>< 1 MiB |
| 100 | 5.13 ms<br>< 1 MiB | n/a | 6.53 ms<br>1.94 MiB | 50.5 ms<br>2.84 MiB | 54.2 ms<br>2.26 MiB |
| 200 | 25.2 ms<br>< 1 MiB | n/a | 33.5 ms<br>1.94 MiB | 678 ms<br>14.7 MiB | 626 ms<br>14 MiB |

## ugcd.gf

| Size | Poly Mojo | Poly Mojo (modular) | FLINT | Singular | SymPy |
|---|---|---|---|---|---|
| 50 | < 1 ms<br>< 1 MiB | n/a | timeout | < 1 ms<br>1.36 MiB | 7.73 ms<br>< 1 MiB |
| 100 | < 1 ms<br>< 1 MiB | n/a | · | < 1 ms<br>1.31 MiB | 29.9 ms<br>< 1 MiB |
| 200 | < 1 ms<br>< 1 MiB | n/a | · | < 1 ms<br>1.39 MiB | 128 ms<br>2.06 MiB |
| 500 | < 1 ms<br>< 1 MiB | n/a | · | < 1 ms<br>2.32 MiB | 803 ms<br>13.8 MiB |
| 1000 | < 1 ms<br>< 1 MiB | n/a | · | < 1 ms<br>2.4 MiB | 3.37 s<br>55.1 MiB |
| 2000 | 1.5 ms<br>< 1 MiB | n/a | · | 1.96 ms<br>2.71 MiB | 14.1 s<br>220 MiB |

## ugcd.zz.128

| Size | Poly Mojo | Poly Mojo (modular) | FLINT | Singular | SymPy |
|---|---|---|---|---|---|
| 50 | < 1 ms<br>< 1 MiB | n/a | < 1 ms<br>< 1 MiB | < 1 ms<br>2.24 MiB | < 1 ms<br>< 1 MiB |
| 100 | < 1 ms<br>< 1 MiB | n/a | < 1 ms<br>< 1 MiB | < 1 ms<br>2.33 MiB | 1.62 ms<br>< 1 MiB |
| 200 | < 1 ms<br>< 1 MiB | n/a | < 1 ms<br>< 1 MiB | < 1 ms<br>2.33 MiB | 6.49 ms<br>< 1 MiB |
| 500 | 1.88 ms<br>< 1 MiB | n/a | 2.72 ms<br>< 1 MiB | 4.06 ms<br>3.3 MiB | 42.2 ms<br>< 1 MiB |
| 1000 | 5.25 ms<br>< 1 MiB | n/a | 7.82 ms<br>1.22 MiB | 10.8 ms<br>3.73 MiB | 179 ms<br>< 1 MiB |
| 2000 | 16.6 ms<br>< 1 MiB | n/a | 29.2 ms<br>2.39 MiB | 29.5 ms<br>4.74 MiB | 721 ms<br>1.05 MiB |

## ugcd.zz.32

| Size | Poly Mojo | Poly Mojo (modular) | FLINT | Singular | SymPy |
|---|---|---|---|---|---|
| 50 | < 1 ms<br>< 1 MiB | n/a | < 1 ms<br>< 1 MiB | < 1 ms<br>2.05 MiB | < 1 ms<br>< 1 MiB |
| 100 | < 1 ms<br>< 1 MiB | n/a | < 1 ms<br>< 1 MiB | < 1 ms<br>2.03 MiB | 1.35 ms<br>< 1 MiB |
| 200 | < 1 ms<br>< 1 MiB | n/a | < 1 ms<br>< 1 MiB | < 1 ms<br>2.08 MiB | 5.29 ms<br>< 1 MiB |
| 500 | < 1 ms<br>< 1 MiB | n/a | 1.15 ms<br>< 1 MiB | 1.83 ms<br>2.96 MiB | 37.7 ms<br>< 1 MiB |
| 1000 | 1.92 ms<br>< 1 MiB | n/a | 3.12 ms<br>1.05 MiB | 4.48 ms<br>3.22 MiB | 144 ms<br>< 1 MiB |
| 2000 | 5.77 ms<br>< 1 MiB | n/a | 10.4 ms<br>1.83 MiB | 11.2 ms<br>3.93 MiB | 664 ms<br>< 1 MiB |

## ugcd.zz.512

| Size | Poly Mojo | Poly Mojo (modular) | FLINT | Singular | SymPy |
|---|---|---|---|---|---|
| 50 | < 1 ms<br>< 1 MiB | n/a | < 1 ms<br>1.03 MiB | < 1 ms<br>2.38 MiB | 1.24 ms<br>< 1 MiB |
| 100 | 1.11 ms<br>< 1 MiB | n/a | < 1 ms<br>1.09 MiB | 1.05 ms<br>2.45 MiB | 4.03 ms<br>< 1 MiB |
| 200 | 2.36 ms<br>< 1 MiB | n/a | 2.43 ms<br>1.15 MiB | 2.61 ms<br>2.65 MiB | 14.2 ms<br>< 1 MiB |
| 500 | 6.46 ms<br>< 1 MiB | n/a | 8.72 ms<br>1.43 MiB | 12.2 ms<br>3.93 MiB | 82.6 ms<br>< 1 MiB |
| 1000 | 17.5 ms<br>< 1 MiB | n/a | 24.1 ms<br>2.05 MiB | 32.5 ms<br>5 MiB | 297 ms<br>1.25 MiB |
| 2000 | 52.8 ms<br>1.82 MiB | n/a | 86.9 ms<br>3.82 MiB | 102 ms<br>7.34 MiB | 1.47 s<br>2.12 MiB |

## ures.gf

| Size | Poly Mojo | Poly Mojo (modular) | FLINT | Singular | SymPy |
|---|---|---|---|---|---|
| 10 | < 1 ms<br>< 1 MiB | n/a | timeout | < 1 ms<br>< 1 MiB | < 1 ms<br>< 1 MiB |
| 20 | < 1 ms<br>< 1 MiB | n/a | · | < 1 ms<br>< 1 MiB | 1.36 ms<br>< 1 MiB |
| 50 | < 1 ms<br>< 1 MiB | n/a | · | < 1 ms<br>< 1 MiB | 7.66 ms<br>< 1 MiB |
| 100 | < 1 ms<br>< 1 MiB | n/a | · | < 1 ms<br>< 1 MiB | 32.2 ms<br>< 1 MiB |
| 200 | < 1 ms<br>< 1 MiB | n/a | · | 3.56 ms<br>1.03 MiB | 117 ms<br>2.79 MiB |

## ures.zz

| Size | Poly Mojo | Poly Mojo (modular) | FLINT | Singular | SymPy |
|---|---|---|---|---|---|
| 10 | < 1 ms<br>< 1 MiB | n/a | < 1 ms<br>< 1 MiB | < 1 ms<br>< 1 MiB | < 1 ms<br>< 1 MiB |
| 20 | < 1 ms<br>< 1 MiB | n/a | < 1 ms<br>< 1 MiB | < 1 ms<br>< 1 MiB | < 1 ms<br>< 1 MiB |
| 50 | 1.18 ms<br>< 1 MiB | n/a | 1.25 ms<br>1.73 MiB | 4.28 ms<br>1.07 MiB | 5.56 ms<br>< 1 MiB |
| 100 | 4.56 ms<br>< 1 MiB | n/a | 6.28 ms<br>1.94 MiB | 47.5 ms<br>2.83 MiB | 52 ms<br>2.16 MiB |
| 200 | 22 ms<br>< 1 MiB | n/a | 29.3 ms<br>2.01 MiB | 569 ms<br>13.9 MiB | 610 ms<br>13.3 MiB |
