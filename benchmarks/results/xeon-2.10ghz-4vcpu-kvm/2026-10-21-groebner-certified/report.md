# Gröbner basis pilot results

> **Certified-mode run on claude/exp-qq-trace c5bbd8d, 2026-10-21** (wave 6 task 5). Scope: grevlex ladders over QQ and GF(65521) only; tools Poly Mojo (comparable), Poly Mojo (certified) and msolve; 60 s limit and 9G cap for every case; one process per case, at least 3 repetitions and 0.5 s. Compare with 2026-10-21-groebner (nightly d1f992f), where certified mode solved 25/44 QQ and 34/44 GF grevlex cases.
> - Certified mode now builds the QQ basis from its trace (exact replay, then inputs, S-pairs and reducedness checked), and `check_groebner` replays the trace and calls `check_groebner_basis`.
> - msolve is timed as a whole process less a baseline, so its sub-20 ms times are noisy; with both tools at 10 ms or more, msolve/certified over QQ is 0.506 (n = 20), and 0.265 (n = 9) with both at 0.1 s or more.

Host: Intel(R) Xeon(R) Processor @ 2.10GHz, 4 vCPUs (KVM), 16 GiB RAM, Ubuntu 24.04.5 LTS, Linux 6.18.44-fc-v140. Every run is pinned to one CPU and single-threaded, capped at 9G of memory and 60 s.

Versions: Poly Mojo c5bbd8d66f1d; Poly Mojo (certified) c5bbd8d66f1d; msolve 0.11.0.

Records: 264 (229 ok, 15 skipped, 20 timeout).

## Summary

Ratios cover cases that both tools solved with a confirmed answer and on which the slower of the two took at least 10 ms (memory: at least 1 MiB). Poly Mojo runs in comparable mode (no witnesses).

**Time** (geometric mean of competitor ÷ Poly Mojo; above 1 means Poly Mojo used less; *n* = cases compared)

| Field, order | msolve |
|---|---|
| GF(65521), grevlex | 7.2 (*n* = 31) |
| QQ, grevlex | 1.47 (*n* = 30) |

**Incremental peak memory** (geometric mean of competitor ÷ Poly Mojo; above 1 means Poly Mojo used less; *n* = cases compared)

| Field, order | msolve |
|---|---|
| GF(65521), grevlex | 0.353 (*n* = 26) |
| QQ, grevlex | 0.351 (*n* = 29) |

**Solved cases** (finished within the limits with an answer no other library contradicts):

| Field, order | Poly Mojo | Poly Mojo (certified) | msolve |
|---|---|---|---|
| GF(65521), grevlex | 42/44 | 40/44 | 42/44 |
| QQ, grevlex | 35/44 | 32/44 | 38/44 |

## Answer agreement

Every answer that two or more libraries computed agrees.

Answers computed by only one library (marked ? in the tables): cyclic-8.QQ.grevlex, katsura-10.QQ.grevlex, randquad-8.QQ.grevlex.

## Certification overhead

Comparable mode computes without witnesses. Certified mode retains cofactor rows and runs `sparse.check_groebner`, which verifies both ideal containments and the Buchberger criterion.

| Case | Comparable | Certified compute | Check | Certified peak ÷ comparable peak |
|---|---|---|---|---|
| cyclic-7.GF65521.grevlex | 43.4 ms | 122 ms | 127 ms | 2.29 |
| cyclic-8.GF65521.grevlex | 762 ms | 2.54 s | 3.01 s | 3.98 |
| cyclic-6.QQ.grevlex | 15.2 ms | 50.8 ms | 18.8 ms | 1.35 |
| cyclic-7.QQ.grevlex | 1.26 s | 15.7 s | 10.4 s | 2.04 |
| eco-9.GF65521.grevlex | 9.73 ms | 18.6 ms | 25.5 ms | 2.32 |
| eco-10.GF65521.grevlex | 38.6 ms | 113 ms | 175 ms | 2.56 |
| eco-11.GF65521.grevlex | 212 ms | 704 ms | 922 ms | 3.25 |
| eco-12.GF65521.grevlex | 1.46 s | 5.38 s | 7.11 s | 6.7 |
| eco-7.QQ.grevlex | 6.32 ms | 13.1 ms | 5.39 ms | 1.43 |
| eco-8.QQ.grevlex | 28.8 ms | 46.6 ms | 23.3 ms | 1.09 |
| eco-9.QQ.grevlex | 116 ms | 253 ms | 119 ms | 1.36 |
| eco-10.QQ.grevlex | 802 ms | 1.5 s | 931 ms | 1.56 |
| eco-11.QQ.grevlex | 5.41 s | 11.1 s | 8.16 s | 1.62 |
| katsura-7.GF65521.grevlex | 4.74 ms | 10.6 ms | 14.4 ms | 1.72 |
| katsura-8.GF65521.grevlex | 22.9 ms | 73.9 ms | 85.8 ms | 2.13 |
| katsura-9.GF65521.grevlex | 139 ms | 526 ms | 534 ms | 2.47 |
| katsura-10.GF65521.grevlex | 1.1 s | 3.7 s | 4.25 s | 3.2 |
| katsura-6.QQ.grevlex | 13.1 ms | 34.1 ms | 18.6 ms | 1.61 |
| katsura-7.QQ.grevlex | 80.3 ms | 197 ms | 115 ms | 1.37 |
| katsura-8.QQ.grevlex | 591 ms | 1.26 s | 793 ms | 1.77 |
| katsura-9.QQ.grevlex | 4.78 s | 10.1 s | 7.2 s | 2 |
| noon-6.GF65521.grevlex | 14.9 ms | 17.9 ms | 20.5 ms | 1.37 |
| noon-7.GF65521.grevlex | 94.6 ms | 126 ms | 147 ms | 1.3 |
| noon-8.GF65521.grevlex | 750 ms | 991 ms | 1.23 s | 2.05 |
| noon-9.GF65521.grevlex | 7.04 s | 7.79 s | 8.89 s | 2.59 |
| noon-5.QQ.grevlex | 16.4 ms | 20.3 ms | 8.44 ms | 1.18 |
| noon-6.QQ.grevlex | 55.2 ms | 135 ms | 47.5 ms | 1.52 |
| noon-7.QQ.grevlex | 345 ms | 820 ms | 301 ms | 1.71 |
| noon-8.QQ.grevlex | 2.89 s | 6.81 s | 2.11 s | 1.46 |
| randquad-8.GF65521.grevlex | 16.4 ms | 59.7 ms | 105 ms | 1.71 |
| randquad-9.GF65521.grevlex | 95.8 ms | 316 ms | 454 ms | 2.01 |
| randquad-10.GF65521.grevlex | 746 ms | 2.47 s | 3.28 s | 2.34 |
| randquad-4.QQ.grevlex | 8 ms | 12.8 ms | 6.35 ms | 1.21 |
| randquad-5.QQ.grevlex | 55 ms | 97.9 ms | 48.6 ms | 1.08 |
| randquad-6.QQ.grevlex | 527 ms | 1.84 s | 754 ms | 1.25 |
| reimer-6.GF65521.grevlex | 25 ms | 39.2 ms | 54 ms | 1.33 |
| reimer-5.QQ.grevlex | 70.2 ms | 36.3 ms | 21.6 ms | 0.48 |
| reimer-6.QQ.grevlex | 6.53 s | 378 ms | 247 ms | 0.173 |

## Per-ladder results

Each cell gives the median time per call and the incremental peak memory. `·` marks an instance skipped after a smaller one on the same ladder hit a limit; `n/a` marks an unsupported order or field.

#### cyclic-n over GF(65521), grevlex

| n | vars | Poly Mojo | Poly Mojo (certified) | msolve |
|---|---|---|---|---|
| 4 | 4 | < 1 ms<br>< 1 MiB | < 1 ms<br>< 1 MiB | < 1 ms<br>< 1 MiB |
| 5 | 5 | < 1 ms<br>< 1 MiB | < 1 ms<br>< 1 MiB | < 1 ms<br>< 1 MiB |
| 6 | 6 | 1.73 ms<br>1.14 MiB | 2.42 ms<br>2.11 MiB | 3.27 ms<br>< 1 MiB |
| 7 | 7 | 43.4 ms<br>11.4 MiB | 122 ms<br>26.1 MiB | 127 ms<br>3.36 MiB |
| 8 | 8 | 762 ms<br>54.1 MiB | 2.54 s<br>215 MiB | 2.1 s<br>25.4 MiB |
| 9 | 9 | timeout | timeout | timeout |

#### eco-n over GF(65521), grevlex

| n | vars | Poly Mojo | Poly Mojo (certified) | msolve |
|---|---|---|---|---|
| 5 | 5 | < 1 ms<br>< 1 MiB | < 1 ms<br>< 1 MiB | 5.08 ms<br>< 1 MiB |
| 6 | 6 | < 1 ms<br>< 1 MiB | < 1 ms<br>< 1 MiB | 18.8 ms<br>< 1 MiB |
| 7 | 7 | < 1 ms<br>< 1 MiB | 1.02 ms<br>< 1 MiB | < 1 ms<br>< 1 MiB |
| 8 | 8 | 2.93 ms<br>1.61 MiB | 3.96 ms<br>2.63 MiB | 22 ms<br>< 1 MiB |
| 9 | 9 | 9.73 ms<br>4.53 MiB | 18.6 ms<br>10.5 MiB | 39.7 ms<br>1.25 MiB |
| 10 | 10 | 38.6 ms<br>11.1 MiB | 113 ms<br>28.4 MiB | 121 ms<br>2.76 MiB |
| 11 | 11 | 212 ms<br>30.4 MiB | 704 ms<br>99 MiB | 737 ms<br>10.4 MiB |
| 12 | 12 | 1.46 s<br>86.6 MiB | 5.38 s<br>581 MiB | 5.12 s<br>36.5 MiB |
| 13 | 13 | 9.58 s<br>283 MiB | timeout | 36 s<br>162 MiB |

#### katsura-n over GF(65521), grevlex

| n | vars | Poly Mojo | Poly Mojo (certified) | msolve |
|---|---|---|---|---|
| 3 | 4 | < 1 ms<br>< 1 MiB | < 1 ms<br>< 1 MiB | 9.93 ms<br>< 1 MiB |
| 4 | 5 | < 1 ms<br>< 1 MiB | < 1 ms<br>< 1 MiB | 21.5 ms<br>< 1 MiB |
| 5 | 6 | < 1 ms<br>< 1 MiB | < 1 ms<br>< 1 MiB | 9.96 ms<br>< 1 MiB |
| 6 | 7 | 1.3 ms<br>1.11 MiB | 2.09 ms<br>1.92 MiB | < 1 ms<br>< 1 MiB |
| 7 | 8 | 4.74 ms<br>2.62 MiB | 10.6 ms<br>4.5 MiB | 36.3 ms<br>1 MiB |
| 8 | 9 | 22.9 ms<br>8.62 MiB | 73.9 ms<br>18.3 MiB | 78.2 ms<br>2.25 MiB |
| 9 | 10 | 139 ms<br>25.2 MiB | 526 ms<br>62.3 MiB | 383 ms<br>5.88 MiB |
| 10 | 11 | 1.1 s<br>81.2 MiB | 3.7 s<br>260 MiB | 2.65 s<br>22.4 MiB |
| 11 | 12 | 7.4 s<br>259 MiB | timeout | 20.6 s<br>91.3 MiB |
| 12 | 13 | timeout | · | timeout |

#### noon-n over GF(65521), grevlex

| n | vars | Poly Mojo | Poly Mojo (certified) | msolve |
|---|---|---|---|---|
| 3 | 3 | < 1 ms<br>< 1 MiB | < 1 ms<br>< 1 MiB | 2.22 ms<br>< 1 MiB |
| 4 | 4 | < 1 ms<br>< 1 MiB | < 1 ms<br>< 1 MiB | 19.3 ms<br>< 1 MiB |
| 5 | 5 | 2.46 ms<br>1.68 MiB | 2.95 ms<br>1.74 MiB | 25.4 ms<br>< 1 MiB |
| 6 | 6 | 14.9 ms<br>4.84 MiB | 17.9 ms<br>6.65 MiB | 55.3 ms<br>1.5 MiB |
| 7 | 7 | 94.6 ms<br>20.8 MiB | 126 ms<br>26.9 MiB | 205 ms<br>6.04 MiB |
| 8 | 8 | 750 ms<br>66.9 MiB | 991 ms<br>137 MiB | 1.45 s<br>35 MiB |
| 9 | 9 | 7.04 s<br>426 MiB | 7.79 s<br>1.08 GiB | 12 s<br>188 MiB |

#### randquad-n over GF(65521), grevlex

| n | vars | Poly Mojo | Poly Mojo (certified) | msolve |
|---|---|---|---|---|
| 3 | 3 | < 1 ms<br>< 1 MiB | < 1 ms<br>< 1 MiB | < 1 ms<br>< 1 MiB |
| 4 | 4 | < 1 ms<br>< 1 MiB | < 1 ms<br>< 1 MiB | 25.4 ms<br>< 1 MiB |
| 5 | 5 | < 1 ms<br>< 1 MiB | < 1 ms<br>< 1 MiB | 20.9 ms<br>< 1 MiB |
| 6 | 6 | 1.12 ms<br>< 1 MiB | 1.72 ms<br>1.57 MiB | 24.2 ms<br>< 1 MiB |
| 7 | 7 | 3.82 ms<br>2.39 MiB | 7.83 ms<br>3.98 MiB | 41.2 ms<br>< 1 MiB |
| 8 | 8 | 16.4 ms<br>6.52 MiB | 59.7 ms<br>11.1 MiB | 77 ms<br>1.75 MiB |
| 9 | 9 | 95.8 ms<br>22 MiB | 316 ms<br>44.1 MiB | 313 ms<br>5.37 MiB |
| 10 | 10 | 746 ms<br>68.2 MiB | 2.47 s<br>159 MiB | 1.93 s<br>18.1 MiB |

#### reimer-n over GF(65521), grevlex

| n | vars | Poly Mojo | Poly Mojo (certified) | msolve |
|---|---|---|---|---|
| 3 | 3 | < 1 ms<br>< 1 MiB | < 1 ms<br>< 1 MiB | 8.16 ms<br>< 1 MiB |
| 4 | 4 | < 1 ms<br>< 1 MiB | < 1 ms<br>1.02 MiB | 24.1 ms<br>< 1 MiB |
| 5 | 5 | 3.19 ms<br>1.51 MiB | 3.54 ms<br>2.18 MiB | 25.1 ms<br>< 1 MiB |
| 6 | 6 | 25 ms<br>6.65 MiB | 39.2 ms<br>8.83 MiB | 78.5 ms<br>2.5 MiB |

#### cyclic-n over QQ, grevlex

| n | vars | Poly Mojo | Poly Mojo (certified) | msolve |
|---|---|---|---|---|
| 4 | 4 | 1.1 ms<br>< 1 MiB | 1.63 ms<br>< 1 MiB | < 1 ms<br>< 1 MiB |
| 5 | 5 | 3.27 ms<br>1.46 MiB | 7.13 ms<br>1.65 MiB | < 1 ms<br>< 1 MiB |
| 6 | 6 | 15.2 ms<br>3.14 MiB | 50.8 ms<br>4.23 MiB | 20.6 ms<br>1.5 MiB |
| 7 | 7 | 1.26 s<br>53.3 MiB | 15.7 s<br>109 MiB | 840 ms<br>10.2 MiB |
| 8 | 8 | timeout | timeout | 20.1 s<br>63.8 MiB ? |
| 9 | 9 | · | · | timeout |

#### eco-n over QQ, grevlex

| n | vars | Poly Mojo | Poly Mojo (certified) | msolve |
|---|---|---|---|---|
| 5 | 5 | 1.69 ms<br>< 1 MiB | 2.59 ms<br>1.12 MiB | 16.2 ms<br>< 1 MiB |
| 6 | 6 | 2.84 ms<br>1.07 MiB | 4.9 ms<br>1.48 MiB | 21.1 ms<br>< 1 MiB |
| 7 | 7 | 6.32 ms<br>1.64 MiB | 13.1 ms<br>2.34 MiB | 18.4 ms<br>< 1 MiB |
| 8 | 8 | 28.8 ms<br>4.63 MiB | 46.6 ms<br>5.07 MiB | 32.9 ms<br>1.94 MiB |
| 9 | 9 | 116 ms<br>12.9 MiB | 253 ms<br>17.5 MiB | 79.1 ms<br>3.41 MiB |
| 10 | 10 | 802 ms<br>36.4 MiB | 1.5 s<br>56.8 MiB | 390 ms<br>9.1 MiB |
| 11 | 11 | 5.41 s<br>125 MiB | 11.1 s<br>202 MiB | 2.69 s<br>37.9 MiB |
| 12 | 12 | 37.9 s<br>450 MiB | timeout | 18.5 s<br>162 MiB |
| 13 | 13 | timeout | · | timeout |

#### katsura-n over QQ, grevlex

| n | vars | Poly Mojo | Poly Mojo (certified) | msolve |
|---|---|---|---|---|
| 3 | 4 | 1.23 ms<br>< 1 MiB | 1.98 ms<br>< 1 MiB | 15.7 ms<br>< 1 MiB |
| 4 | 5 | 2.03 ms<br>< 1 MiB | 3.85 ms<br>1.58 MiB | 26.4 ms<br>< 1 MiB |
| 5 | 6 | 4.21 ms<br>1.37 MiB | 8.73 ms<br>2.1 MiB | 19.6 ms<br>< 1 MiB |
| 6 | 7 | 13.1 ms<br>2.42 MiB | 34.1 ms<br>3.89 MiB | 38.9 ms<br>< 1 MiB |
| 7 | 8 | 80.3 ms<br>7.55 MiB | 197 ms<br>10.4 MiB | 61 ms<br>2.69 MiB |
| 8 | 9 | 591 ms<br>25.5 MiB | 1.26 s<br>45.1 MiB | 325 ms<br>8.41 MiB |
| 9 | 10 | 4.78 s<br>88.6 MiB | 10.1 s<br>177 MiB | 1.86 s<br>37.4 MiB |
| 10 | 11 | timeout | timeout | 13.9 s<br>145 MiB ? |
| 11 | 12 | · | · | timeout |
| 12 | 13 | · | · | · |

#### noon-n over QQ, grevlex

| n | vars | Poly Mojo | Poly Mojo (certified) | msolve |
|---|---|---|---|---|
| 3 | 3 | 1.39 ms<br>1 MiB | 1.92 ms<br>1.24 MiB | < 1 ms<br>< 1 MiB |
| 4 | 4 | 3.14 ms<br>1.38 MiB | 4.77 ms<br>1.32 MiB | 41 ms<br>< 1 MiB |
| 5 | 5 | 16.4 ms<br>2.85 MiB | 20.3 ms<br>3.37 MiB | 34.3 ms<br>< 1 MiB |
| 6 | 6 | 55.2 ms<br>8.75 MiB | 135 ms<br>13.3 MiB | 88.8 ms<br>4.16 MiB |
| 7 | 7 | 345 ms<br>33.2 MiB | 820 ms<br>56.9 MiB | 474 ms<br>21.2 MiB |
| 8 | 8 | 2.89 s<br>195 MiB | 6.81 s<br>285 MiB | 3.13 s<br>130 MiB |
| 9 | 9 | 22.1 s<br>1 GiB | timeout | 22.1 s<br>642 MiB |

#### randquad-n over QQ, grevlex

| n | vars | Poly Mojo | Poly Mojo (certified) | msolve |
|---|---|---|---|---|
| 3 | 3 | 2.84 ms<br>< 1 MiB | 2.78 ms<br>< 1 MiB | 5.71 ms<br>< 1 MiB |
| 4 | 4 | 8 ms<br>1.38 MiB | 12.8 ms<br>1.66 MiB | 35.1 ms<br>< 1 MiB |
| 5 | 5 | 55 ms<br>3.78 MiB | 97.9 ms<br>4.07 MiB | 108 ms<br>1.43 MiB |
| 6 | 6 | 527 ms<br>12.5 MiB | 1.84 s<br>15.7 MiB | 462 ms<br>3.32 MiB |
| 7 | 7 | 9.12 s<br>47.1 MiB | timeout | 3.22 s<br>19.2 MiB |
| 8 | 8 | timeout | · | 35.7 s<br>165 MiB ? |
| 9 | 9 | · | · | timeout |
| 10 | 10 | · | · | · |

#### reimer-n over QQ, grevlex

| n | vars | Poly Mojo | Poly Mojo (certified) | msolve |
|---|---|---|---|---|
| 3 | 3 | 1.59 ms<br>< 1 MiB | 2.14 ms<br>1.19 MiB | < 1 ms<br>< 1 MiB |
| 4 | 4 | 5.7 ms<br>1.57 MiB | 6.89 ms<br>1.77 MiB | 25.5 ms<br>< 1 MiB |
| 5 | 5 | 70.2 ms<br>9.07 MiB | 36.3 ms<br>4.35 MiB | 45.5 ms<br>1.41 MiB |
| 6 | 6 | 6.53 s<br>139 MiB | 378 ms<br>24 MiB | 227 ms<br>5.79 MiB |

## Plots

![plot-cyclic.GF65521.grevlex.svg](plot-cyclic.GF65521.grevlex.svg)
![plot-cyclic.QQ.grevlex.svg](plot-cyclic.QQ.grevlex.svg)
![plot-eco.GF65521.grevlex.svg](plot-eco.GF65521.grevlex.svg)
![plot-eco.QQ.grevlex.svg](plot-eco.QQ.grevlex.svg)
![plot-katsura.GF65521.grevlex.svg](plot-katsura.GF65521.grevlex.svg)
![plot-katsura.QQ.grevlex.svg](plot-katsura.QQ.grevlex.svg)
![plot-noon.GF65521.grevlex.svg](plot-noon.GF65521.grevlex.svg)
![plot-noon.QQ.grevlex.svg](plot-noon.QQ.grevlex.svg)
![plot-randquad.GF65521.grevlex.svg](plot-randquad.GF65521.grevlex.svg)
![plot-randquad.QQ.grevlex.svg](plot-randquad.QQ.grevlex.svg)
![plot-reimer.GF65521.grevlex.svg](plot-reimer.GF65521.grevlex.svg)
![plot-reimer.QQ.grevlex.svg](plot-reimer.QQ.grevlex.svg)
