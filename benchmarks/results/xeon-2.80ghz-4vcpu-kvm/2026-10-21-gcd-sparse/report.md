# Gcd, resultant, and discriminant results

Host: Intel(R) Xeon(R) Processor @ 2.80GHz, 4 vCPUs. One pinned CPU per run, 12G and 60 s per run.

Versions: Poly Mojo d1f992ffe139; Poly Mojo (modular) d1f992ffe139; FLINT python-flint 0.9.0; Singular 4.4.1; SymPy 1.14.0 (gmpy ground types).

Cells give the median time and incremental peak memory; ✗ marks an answer that disagrees with a confirmed one, ? an answer no other library produced. Ratios are competitor time ÷ Poly Mojo time (above 1: Poly Mojo faster) over cases both solved with a confirmed answer and at least 10 ms for the slower tool.

## Summary

| Ladder | Poly Mojo (modular) | FLINT | Singular | SymPy |
|---|---|---|---|---|
| mgcd.dense | 8.93 (5) | 0.566 (2) | 0.672 (2) | 18.9 (5) |
| mgcd.fateman-f1 | 89 (3) |  |  | 36.4 (1) |
| mgcd.fateman-f2 | 53 (4) |  |  | 20.2 (4) |
| mgcd.fateman-f3 | 2.98e+03 (2) |  |  | 20.3 (4) |
| mgcd.sparse | 557 (6) |  |  | 17.2 (4) |
| mres.2v |  |  |  | 21.6 (2) |
| mres.3v |  | 1.31 (3) | 5.03 (3) | 39.4 (3) |

## mgcd.dense

| Size | Poly Mojo | Poly Mojo (modular) | FLINT | Singular | SymPy |
|---|---|---|---|---|---|
| 5 | < 1 ms<br>< 1 MiB | 4.65 ms<br>1.12 MiB | < 1 ms<br>1.82 MiB | < 1 ms<br>2.38 MiB | 5.24 ms<br>< 1 MiB |
| 7 | 1.25 ms<br>< 1 MiB | 12.4 ms<br>1.16 MiB | < 1 ms<br>1.86 MiB | < 1 ms<br>2.58 MiB | 16.4 ms<br>< 1 MiB |
| 9 | 3.23 ms<br>< 1 MiB | 27.2 ms<br>1.54 MiB | 1.63 ms<br>1.96 MiB | 1.99 ms<br>2.53 MiB | 45.9 ms<br>< 1 MiB |
| 11 | 5.98 ms<br>< 1 MiB | 51.5 ms<br>2.44 MiB | 2.81 ms<br>1.93 MiB | 3.46 ms<br>3 MiB | 112 ms<br>1.18 MiB |
| 13 | 10.3 ms<br>1.92 MiB | 94.5 ms<br>3.91 MiB | 4.64 ms<br>1.9 MiB | 5.63 ms<br>3.44 MiB | 252 ms<br>1.82 MiB |
| 15 | 18.4 ms<br>2.32 MiB | 158 ms<br>6.07 MiB | 13.1 ms<br>2.14 MiB | 15.2 ms<br>3.85 MiB | 527 ms<br>2.46 MiB |

## mgcd.fateman-f1

| Size | Poly Mojo | Poly Mojo (modular) | FLINT | Singular | SymPy |
|---|---|---|---|---|---|
| 2 | < 1 ms<br>< 1 MiB | < 1 ms<br>< 1 MiB | < 1 ms<br>< 1 MiB | < 1 ms<br>1.48 MiB | < 1 ms<br>< 1 MiB |
| 3 | < 1 ms<br>< 1 MiB | 1.8 ms<br>< 1 MiB | < 1 ms<br>< 1 MiB | < 1 ms<br>1.55 MiB | 1.07 ms<br>< 1 MiB |
| 4 | < 1 ms<br>< 1 MiB | 3.5 ms<br>1.75 MiB | < 1 ms<br>< 1 MiB | < 1 ms<br>1.49 MiB | 2 ms<br>< 1 MiB |
| 5 | < 1 ms<br>< 1 MiB | 7.33 ms<br>3.52 MiB | < 1 ms<br>< 1 MiB | < 1 ms<br>1.55 MiB | 3.64 ms<br>< 1 MiB |
| 6 | < 1 ms<br>< 1 MiB | 14.3 ms<br>5.63 MiB | < 1 ms<br>< 1 MiB | < 1 ms<br>1.55 MiB | 6.04 ms<br>< 1 MiB |
| 7 | < 1 ms<br>< 1 MiB | 29 ms<br>9.28 MiB | < 1 ms<br>< 1 MiB | < 1 ms<br>1.55 MiB | 9.55 ms<br>< 1 MiB |
| 8 | < 1 ms<br>< 1 MiB | 57.5 ms<br>16.1 MiB | < 1 ms<br>< 1 MiB | < 1 ms<br>1.56 MiB | 15.2 ms<br>< 1 MiB |

## mgcd.fateman-f2

| Size | Poly Mojo | Poly Mojo (modular) | FLINT | Singular | SymPy |
|---|---|---|---|---|---|
| 2 | < 1 ms<br>< 1 MiB | < 1 ms<br>< 1 MiB | < 1 ms<br>2.14 MiB | < 1 ms<br>3.1 MiB | < 1 ms<br>< 1 MiB |
| 3 | < 1 ms<br>< 1 MiB | 2.37 ms<br>< 1 MiB | < 1 ms<br>2.12 MiB | < 1 ms<br>3.13 MiB | 1.99 ms<br>< 1 MiB |
| 4 | < 1 ms<br>< 1 MiB | 6.6 ms<br>1.07 MiB | < 1 ms<br>2.25 MiB | < 1 ms<br>3.1 MiB | 4.85 ms<br>< 1 MiB |
| 5 | 1.51 ms<br>< 1 MiB | 20 ms<br>< 1 MiB | 1.14 ms<br>2.36 MiB | 1.2 ms<br>3.23 MiB | 11.5 ms<br>< 1 MiB |
| 6 | 1.17 ms<br>1.18 MiB | 55.4 ms<br>1.56 MiB | 1.71 ms<br>2.38 MiB | 1.78 ms<br>3.25 MiB | 24.2 ms<br>< 1 MiB |
| 7 | 2.17 ms<br>1.5 MiB | 175 ms<br>3.26 MiB | 2.37 ms<br>2.33 MiB | 2.48 ms<br>3.2 MiB | 57.5 ms<br>1.22 MiB |
| 8 | 3.66 ms<br>1.98 MiB | 564 ms<br>4.39 MiB | 3.75 ms<br>2.48 MiB | 3.89 ms<br>3.39 MiB | 144 ms<br>3.52 MiB |

## mgcd.fateman-f3

| Size | Poly Mojo | Poly Mojo (modular) | FLINT | Singular | SymPy |
|---|---|---|---|---|---|
| 2 | < 1 ms<br>< 1 MiB | 2.74 ms<br>< 1 MiB | < 1 ms<br>2.21 MiB | < 1 ms<br>3.22 MiB | < 1 ms<br>< 1 MiB |
| 3 | < 1 ms<br>< 1 MiB | 28 ms<br>2.46 MiB | < 1 ms<br>2.17 MiB | < 1 ms<br>3.2 MiB | 2.3 ms<br>< 1 MiB |
| 4 | < 1 ms<br>< 1 MiB | 45.1 s<br>56.9 MiB | < 1 ms<br>2.24 MiB | < 1 ms<br>3.18 MiB | 6 ms<br>< 1 MiB |
| 5 | 1.57 ms<br>< 1 MiB | timeout | 1.16 ms<br>2.37 MiB | 1.23 ms<br>3.2 MiB | 12.5 ms<br>< 1 MiB |
| 6 | 1.34 ms<br>< 1 MiB | · | 1.81 ms<br>2.44 MiB | 1.8 ms<br>3.33 MiB | 27.2 ms<br>< 1 MiB |
| 7 | 2.36 ms<br>1.9 MiB | · | 2.43 ms<br>2.44 MiB | 2.53 ms<br>3.27 MiB | 60.1 ms<br>1.4 MiB |
| 8 | 3.89 ms<br>2.18 MiB | · | 3.76 ms<br>2.46 MiB | 3.96 ms<br>3.4 MiB | 159 ms<br>3.73 MiB |

## mgcd.sparse

| Size | Poly Mojo | Poly Mojo (modular) | FLINT | Singular | SymPy |
|---|---|---|---|---|---|
| 3 | < 1 ms<br>< 1 MiB | 228 ms<br>3.12 MiB | < 1 ms<br>1.89 MiB | < 1 ms<br>2.56 MiB | 2.63 ms<br>< 1 MiB |
| 4 | 1.03 ms<br>< 1 MiB | 7.38 s<br>20.8 MiB | 1.37 ms<br>2.25 MiB | 1.47 ms<br>3.18 MiB | 7.63 ms<br>< 1 MiB |
| 5 | 1.15 ms<br>< 1 MiB | 14.6 s<br>63.9 MiB | 1.25 ms<br>2.39 MiB | 1.33 ms<br>3.48 MiB | 13.2 ms<br>< 1 MiB |
| 6 | 1.51 ms<br>1.06 MiB | 99.7 ms<br>27.4 MiB | 1.12 ms<br>2.26 MiB | 1.2 ms<br>3.34 MiB | 25.6 ms<br>< 1 MiB |
| 7 | 1.75 ms<br>1.02 MiB | 126 ms<br>31.4 MiB | 1.65 ms<br>2.54 MiB | 1.74 ms<br>3.52 MiB | 31.2 ms<br>< 1 MiB |
| 8 | 2.52 ms<br>< 1 MiB | 338 ms<br>3.32 MiB | 1.6 ms<br>2.44 MiB | 1.69 ms<br>3.5 MiB | 63.2 ms<br>1 MiB |

## mres.2v

| Size | Poly Mojo | Poly Mojo (modular) | FLINT | Singular | SymPy |
|---|---|---|---|---|---|
| 2 | < 1 ms<br>< 1 MiB | n/a | < 1 ms<br>< 1 MiB | < 1 ms<br>< 1 MiB | < 1 ms<br>< 1 MiB |
| 4 | < 1 ms<br>< 1 MiB | n/a | < 1 ms<br>< 1 MiB | < 1 ms<br>2.18 MiB | < 1 ms<br>< 1 MiB |
| 6 | < 1 ms<br>< 1 MiB | n/a | < 1 ms<br>< 1 MiB | < 1 ms<br>2.35 MiB | 4.09 ms<br>< 1 MiB |
| 8 | < 1 ms<br>< 1 MiB | n/a | 1.03 ms<br>1.04 MiB | 2.32 ms<br>2.63 MiB | 17.6 ms<br>< 1 MiB |
| 10 | 2.41 ms<br>< 1 MiB | n/a | 4.58 ms<br>1.14 MiB | 6.78 ms<br>2.84 MiB | 56.8 ms<br>< 1 MiB |

## mres.3v

| Size | Poly Mojo | Poly Mojo (modular) | FLINT | Singular | SymPy |
|---|---|---|---|---|---|
| 2 | < 1 ms<br>< 1 MiB | n/a | < 1 ms<br>< 1 MiB | < 1 ms<br>< 1 MiB | < 1 ms<br>< 1 MiB |
| 4 | 2.68 ms<br>< 1 MiB | n/a | < 1 ms<br>< 1 MiB | 1.32 ms<br>2 MiB | 35.4 ms<br>< 1 MiB |
| 6 | 19.4 ms<br>< 1 MiB | n/a | 6.2 ms<br>1.48 MiB | 60.5 ms<br>3.36 MiB | 937 ms<br>1.71 MiB |
| 8 | 123 ms<br>1.84 MiB | n/a | 280 ms<br>4.11 MiB | 714 ms<br>10.4 MiB | 11.8 s<br>7.55 MiB |
| 10 | 592 ms<br>2.8 MiB | n/a | 1.84 s<br>12.6 MiB | 4.15 s<br>35.9 MiB | timeout |
