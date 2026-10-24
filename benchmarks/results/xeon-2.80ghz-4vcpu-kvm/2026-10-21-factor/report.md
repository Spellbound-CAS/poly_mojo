# Factorization

Records: `2026-10-21-factor`. Times are medians of in-process repetitions on one pinned core. Ratio is the other tool's time over Poly Mojo's (above 1: Poly Mojo is faster). "agree" means every tool that finished returned the same normalized factor multiset and its product equals the input.

| case | polymojo | flint | singular | sympy | flint ratio | singular ratio | sympy ratio | agree |
|---|---|---|---|---|---|---|---|---|
| zz-cyclotomic-360 | 0.07 ms | 81.94 ms | 4.13 s | 140.16 ms | 1.23e+03 | 6.21e+04 | 2.11e+03 | yes |
| zz-cyclotomic-2520 | 0.43 ms | 52.30 s | timeout | 12.14 s | 1.21e+05 |  | 2.82e+04 | yes |
| zz-cyclotomic-5040 | 0.87 ms | timeout | skipped | 52.16 s |  |  | 6e+04 | yes |
| zz-product-k2-d10-b16 | 0.20 ms | 0.35 ms | 0.44 ms | 7.24 ms | 1.78 | 2.25 | 36.6 | yes |
| zz-product-k2-d25-b16 | 2.01 ms | 0.90 ms | 1.08 ms | 51.88 ms | 0.449 | 0.536 | 25.8 | yes |
| zz-product-k2-d50-b16 | 7.60 ms | 2.50 ms | 3.11 ms | 277.26 ms | 0.329 | 0.41 | 36.5 | yes |
| zz-product-k2-d100-b16 | 41.32 ms | 13.70 ms | 15.92 ms | 2.04 s | 0.331 | 0.385 | 49.4 | yes |
| zz-product-k2-d10-b256 | 2.28 ms | 1.33 ms | 1.32 ms | 27.18 ms | 0.583 | 0.579 | 11.9 | yes |
| zz-product-k2-d25-b256 | 5.91 ms | 3.67 ms | 3.69 ms | 222.34 ms | 0.621 | 0.624 | 37.6 | yes |
| zz-product-k2-d50-b256 | 24.16 ms | 8.17 ms | 8.55 ms | 777.83 ms | 0.338 | 0.354 | 32.2 | yes |
| zz-product-k2-d100-b256 | 84.39 ms | 37.25 ms | 36.88 ms | 7.62 s | 0.441 | 0.437 | 90.4 | yes |
| zz-product-k8-d10-b16 | 18.46 ms | 6.36 ms | 6.76 ms | 1.11 s | 0.345 | 0.366 | 60.2 | yes |
| zz-product-k8-d25-b16 | 96.59 ms | 17.63 ms | 20.97 ms | 11.52 s | 0.183 | 0.217 | 119 | yes |
| zz-product-k8-d50-b16 | 244.00 ms | 56.64 ms | 70.21 ms | 65.83 s | 0.232 | 0.288 | 270 | yes |
| zz-product-k8-d100-b16 | 739.09 ms | 137.33 ms | 167.19 ms | timeout | 0.186 | 0.226 |  | yes |
| zz-product-k8-d10-b256 | 54.34 ms | 21.67 ms | 22.44 ms | 3.89 s | 0.399 | 0.413 | 71.7 | yes |
| zz-product-k8-d25-b256 | 225.57 ms | 84.86 ms | 86.83 ms | 60.58 s | 0.376 | 0.385 | 269 | yes |
| zz-product-k8-d50-b256 | 601.36 ms | 264.49 ms | 262.13 ms | timeout | 0.44 | 0.436 |  | yes |
| zz-product-k8-d100-b256 | 1.17 s | 615.08 ms | 555.80 ms | skipped | 0.526 | 0.476 |  | yes |
| zz-repeated-d10 | 0.25 ms | 0.16 ms | 0.21 ms | 11.81 ms | 0.667 | 0.851 | 47.8 | yes |
| zz-repeated-d20 | 0.75 ms | 0.41 ms | 0.54 ms | 39.27 ms | 0.549 | 0.724 | 52.2 | yes |
| zz-repeated-d40 | 2.43 ms | 1.62 ms | 1.89 ms | 199.30 ms | 0.666 | 0.778 | 82.1 | yes |
| zz-swinnerton-dyer-4 | 0.30 ms | 0.42 ms | 0.47 ms | 7.20 ms | 1.41 | 1.58 | 24.3 | yes |
| zz-swinnerton-dyer-5 | 5.97 ms | 2.10 ms | 2.63 ms | 194.31 ms | 0.352 | 0.44 | 32.5 | yes |
| zz-swinnerton-dyer-6 | 13.20 ms | 9.24 ms | 12.14 ms | timeout | 0.7 | 0.92 |  | yes |
| zz-swinnerton-dyer-7 | 92.33 ms | 69.53 ms | 198.97 ms | skipped | 0.753 | 2.16 |  | yes |
| gf-p3-d100 | 0.25 ms | 0.35 ms | 0.53 ms | 102.80 ms | 1.39 | 2.11 | 409 | yes |
| gf-p3-d250 | 2.36 ms | 1.83 ms | 2.66 ms | 1.07 s | 0.776 | 1.13 | 453 | yes |
| gf-p3-d500 | 5.25 ms | 16.26 ms | 26.60 ms | 24.98 s | 3.1 | 5.07 | 4.76e+03 | yes |
| gf-p3-d1000 | 40.11 ms | 38.25 ms | 75.14 ms | 70.67 s | 0.954 | 1.87 | 1.76e+03 | yes |
| gf-p3-d2000 | 181.04 ms | 264.42 ms | 296.90 ms | timeout | 1.46 | 1.64 |  | yes |
| gf-p61-d100 | 6.36 ms | 4.58 ms | unsupported | 1.53 s | 0.72 |  | 241 | yes |
| gf-p61-d250 | 45.02 ms | 69.46 ms | unsupported | 24.00 s | 1.54 |  | 533 | yes |
| gf-p61-d500 | 195.96 ms | 233.79 ms | unsupported | 199.44 s | 1.19 |  | 1.02e+03 | yes |
| gf-p61-d1000 | 844.34 ms | 905.07 ms | unsupported | timeout | 1.07 |  |  | yes |
| gf-p61-d2000 | 2.72 s | 5.83 s | unsupported | skipped | 2.14 |  |  | yes |
| gf-p65521-d100 | 1.20 ms | 1.45 ms | 1.58 ms | 1.19 s | 1.2 | 1.31 | 992 | yes |
| gf-p65521-d250 | 12.73 ms | 11.60 ms | 12.96 ms | 24.65 s | 0.911 | 1.02 | 1.94e+03 | yes |
| gf-p65521-d500 | 45.57 ms | 59.56 ms | 56.62 ms | 190.15 s | 1.31 | 1.24 | 4.17e+03 | yes |
| gf-p65521-d1000 | 216.10 ms | 235.45 ms | 203.56 ms | timeout | 1.09 | 0.942 |  | yes |
| gf-p65521-d2000 | 901.76 ms | 1.36 s | 1.06 s | skipped | 1.5 | 1.17 |  | yes |
| mzz-dense-k5 | 0.23 ms | 0.50 ms | 701.62 ms | 16.01 ms | 2.17 | 3.07e+03 | 70.1 | yes |
| mzz-dense-k7 | 0.38 ms | 0.86 ms | 1.42 ms | 38.09 ms | 2.27 | 3.77 | 101 | yes |
| mzz-dense-k10 | 0.93 ms | 1.27 ms | 2.69 ms | 135.74 ms | 1.37 | 2.91 | 146 | yes |
| mzz-dense-k12 | 1.43 ms | 1.90 ms | 4.42 ms | 296.86 ms | 1.32 | 3.08 | 207 | yes |
| mzz-dense-k15 | 2.49 ms | 3.20 ms | 8.41 ms | 946.61 ms | 1.28 | 3.38 | 380 | yes |
| mzz-sparse-n3-t10 | 2.47 ms | 1.88 ms | 1.80 ms | 33.04 ms | 0.761 | 0.728 | 13.4 | yes |
| mzz-sparse-n3-t30 | 4.81 ms | 2.51 ms | 3.03 ms | 57.24 ms | 0.523 | 0.631 | 11.9 | yes |
| mzz-sparse-n4-t20 | 10.13 ms | 4.85 ms | 5.81 ms | 227.75 ms | 0.479 | 0.573 | 22.5 | yes |
| mzz-sparse-n5-t30 | 31.37 ms | 34.41 ms | 15.19 ms | 1.59 s | 1.1 | 0.484 | 50.6 | yes |
| mzz-sparse-n6-t30 | 51.94 ms | 128.38 ms | 65.04 ms | 5.86 s | 2.47 | 1.25 | 113 | yes |
| mzz-sparse-n8-t30 | 93.61 ms | 781.18 ms | 456.41 ms | 87.93 s | 8.34 | 4.88 | 939 | yes |
| mzz-sparse-n8-t100 | 834.72 ms | 42.56 s | 7.68 s | 64.40 s | 51 | 9.2 | 77.2 | yes |
| mzz-wang-f6 | 0.04 ms | 0.02 ms | 0.03 ms | 1.28 ms | 0.529 | 0.772 | 31.8 | yes |
| mzz-wang-f1 | 0.12 ms | 0.03 ms | 0.06 ms | 1.70 ms | 0.256 | 0.483 | 13.9 | yes |
| mzz-wang-f7 | 0.80 ms | 0.44 ms | 0.68 ms | 10.49 ms | 0.55 | 0.854 | 13.1 | yes |
| mzz-wang-f3 | 0.55 ms | 0.44 ms | 0.53 ms | 7.32 ms | 0.805 | 0.968 | 13.4 | yes |
| mzz-wang-f4 | 0.57 ms | 1.27 ms | 1.29 ms | 7.94 ms | 2.22 | 2.26 | 13.9 | yes |
| mzz-wang-w2 | 0.93 ms | 0.47 ms | 0.49 ms | 11.48 ms | 0.511 | 0.529 | 12.4 | yes |
| mzz-wang-f2 | 0.48 ms | 0.42 ms | 0.45 ms | 5.60 ms | 0.883 | 0.953 | 11.8 | yes |
| mzz-wang-f5 | 1.47 ms | 1.28 ms | 1.51 ms | 88.99 ms | 0.873 | 1.03 | 60.7 | yes |
| mzz-wang-w1 | 1.20 ms | 1.38 ms | 1.45 ms | 22.94 ms | 1.15 | 1.21 | 19.2 | yes |
