# Factorization

Records: `2026-10-22-factor`. Times are medians of in-process repetitions on one pinned core. Ratio is the other tool's time over Poly Mojo's (above 1: Poly Mojo is faster). "agree" means every tool that finished returned the same normalized factor multiset and its product equals the input.

| case | polymojo | flint | singular | sympy | flint ratio | singular ratio | sympy ratio | agree |
|---|---|---|---|---|---|---|---|---|
| zz-cyclotomic-360 | 0.06 ms | 88.78 ms | 3.97 s | 140.09 ms | 1.37e+03 | 6.11e+04 | 2.16e+03 | yes |
| zz-cyclotomic-2520 | 0.42 ms | 48.45 s | timeout | 11.97 s | 1.15e+05 |  | 2.85e+04 | yes |
| zz-cyclotomic-5040 | 0.80 ms | timeout | skipped | 49.01 s |  |  | 6.09e+04 | yes |
| zz-product-k2-d10-b16 | 0.18 ms | 0.60 ms | 0.82 ms | 7.16 ms | 3.26 | 4.45 | 39 | yes |
| zz-product-k2-d25-b16 | 1.32 ms | 0.99 ms | 1.11 ms | 45.30 ms | 0.75 | 0.843 | 34.4 | yes |
| zz-product-k2-d50-b16 | 5.50 ms | 2.71 ms | 3.22 ms | 361.91 ms | 0.493 | 0.586 | 65.9 | yes |
| zz-product-k2-d100-b16 | 22.87 ms | 13.96 ms | 23.35 ms | 2.15 s | 0.61 | 1.02 | 94.2 | yes |
| zz-product-k2-d10-b256 | 2.54 ms | 1.36 ms | 1.45 ms | 30.46 ms | 0.538 | 0.57 | 12 | yes |
| zz-product-k2-d25-b256 | 5.16 ms | 3.77 ms | 3.96 ms | 223.59 ms | 0.73 | 0.768 | 43.3 | yes |
| zz-product-k2-d50-b256 | 19.13 ms | 9.12 ms | 10.29 ms | 938.23 ms | 0.477 | 0.538 | 49 | yes |
| zz-product-k2-d100-b256 | 61.71 ms | 42.65 ms | 44.80 ms | 8.38 s | 0.691 | 0.726 | 136 | yes |
| zz-product-k8-d10-b16 | 8.17 ms | 6.36 ms | 7.17 ms | 1.20 s | 0.779 | 0.877 | 147 | yes |
| zz-product-k8-d25-b16 | 48.65 ms | 17.69 ms | 21.15 ms | 9.77 s | 0.364 | 0.435 | 201 | yes |
| zz-product-k8-d50-b16 | 105.43 ms | 57.08 ms | 73.26 ms | 66.22 s | 0.541 | 0.695 | 628 | yes |
| zz-product-k8-d100-b16 | 261.06 ms | 141.85 ms | 200.17 ms | timeout | 0.543 | 0.767 |  | yes |
| zz-product-k8-d10-b256 | 48.84 ms | 31.14 ms | 25.30 ms | 4.44 s | 0.638 | 0.518 | 91 | yes |
| zz-product-k8-d25-b256 | 218.62 ms | 100.28 ms | 119.27 ms | 63.60 s | 0.459 | 0.546 | 291 | yes |
| zz-product-k8-d50-b256 | 363.96 ms | 222.30 ms | 297.45 ms | timeout | 0.611 | 0.817 |  | yes |
| zz-product-k8-d100-b256 | 888.27 ms | 706.34 ms | 635.06 ms | skipped | 0.795 | 0.715 |  | yes |
| zz-repeated-d10 | 0.22 ms | 0.18 ms | 0.22 ms | 12.04 ms | 0.796 | 0.984 | 53.6 | yes |
| zz-repeated-d20 | 0.61 ms | 0.49 ms | 0.57 ms | 40.67 ms | 0.802 | 0.925 | 66.2 | yes |
| zz-repeated-d40 | 1.81 ms | 1.76 ms | 2.12 ms | 198.19 ms | 0.973 | 1.17 | 110 | yes |
| zz-swinnerton-dyer-4 | 0.30 ms | 0.54 ms | 0.52 ms | 7.68 ms | 1.81 | 1.74 | 25.5 | yes |
| zz-swinnerton-dyer-5 | 3.01 ms | 3.57 ms | 2.81 ms | 186.89 ms | 1.19 | 0.932 | 62.1 | yes |
| zz-swinnerton-dyer-6 | 10.24 ms | 8.96 ms | 12.52 ms | timeout | 0.875 | 1.22 |  | yes |
| zz-swinnerton-dyer-7 | 83.47 ms | 68.43 ms | 215.68 ms | skipped | 0.82 | 2.58 |  | yes |
| gf-p3-d100 | 0.21 ms | 0.33 ms | 0.51 ms | 101.39 ms | 1.54 | 2.37 | 472 | yes |
| gf-p3-d250 | 1.73 ms | 1.82 ms | 2.64 ms | 1.05 s | 1.05 | 1.52 | 607 | yes |
| gf-p3-d500 | 5.48 ms | 17.35 ms | 26.18 ms | 24.83 s | 3.16 | 4.78 | 4.53e+03 | yes |
| gf-p3-d1000 | 46.33 ms | 37.82 ms | 65.97 ms | 67.79 s | 0.816 | 1.42 | 1.46e+03 | yes |
| gf-p3-d2000 | 233.21 ms | 240.85 ms | 276.72 ms | timeout | 1.03 | 1.19 |  | yes |
| gf-p61-d100 | 5.38 ms | 6.76 ms | unsupported | 1.29 s | 1.26 |  | 240 | yes |
| gf-p61-d250 | 37.15 ms | 69.54 ms | unsupported | 21.61 s | 1.87 |  | 582 | yes |
| gf-p61-d500 | 190.95 ms | 292.73 ms | unsupported | 182.67 s | 1.53 |  | 957 | yes |
| gf-p61-d1000 | 730.77 ms | 1.22 s | unsupported | timeout | 1.66 |  |  | yes |
| gf-p61-d2000 | 2.39 s | 6.66 s | unsupported | skipped | 2.79 |  |  | yes |
| gf-p65521-d100 | 0.98 ms | 1.40 ms | 1.56 ms | 1.14 s | 1.43 | 1.6 | 1.16e+03 | yes |
| gf-p65521-d250 | 10.05 ms | 12.62 ms | 12.62 ms | 21.27 s | 1.26 | 1.26 | 2.12e+03 | yes |
| gf-p65521-d500 | 48.37 ms | 65.13 ms | 57.73 ms | 181.30 s | 1.35 | 1.19 | 3.75e+03 | yes |
| gf-p65521-d1000 | 265.59 ms | 427.64 ms | 331.98 ms | timeout | 1.61 | 1.25 |  | yes |
| gf-p65521-d2000 | 941.87 ms | 1.68 s | 1.05 s | skipped | 1.78 | 1.11 |  | yes |
| mzz-dense-k5 | 0.21 ms | 0.51 ms | 609.29 ms | 14.83 ms | 2.4 | 2.87e+03 | 69.8 | yes |
| mzz-dense-k7 | 0.34 ms | 1.03 ms | 1.53 ms | 36.14 ms | 3.01 | 4.45 | 105 | yes |
| mzz-dense-k10 | 1.27 ms | 1.22 ms | 2.91 ms | 133.92 ms | 0.965 | 2.3 | 106 | yes |
| mzz-dense-k12 | 1.40 ms | 1.79 ms | 5.09 ms | 352.59 ms | 1.28 | 3.64 | 252 | yes |
| mzz-dense-k15 | 2.86 ms | 3.11 ms | 9.22 ms | 1.18 s | 1.09 | 3.22 | 412 | yes |
| mzz-sparse-n3-t10 | 2.35 ms | 2.91 ms | 2.29 ms | 37.39 ms | 1.24 | 0.972 | 15.9 | yes |
| mzz-sparse-n3-t30 | 4.60 ms | 2.83 ms | 3.28 ms | 71.73 ms | 0.616 | 0.713 | 15.6 | yes |
| mzz-sparse-n4-t20 | 9.83 ms | 7.55 ms | 6.83 ms | 290.56 ms | 0.768 | 0.695 | 29.6 | yes |
| mzz-sparse-n5-t30 | 29.91 ms | 52.42 ms | 16.82 ms | 1.61 s | 1.75 | 0.562 | 53.8 | yes |
| mzz-sparse-n6-t30 | 65.23 ms | 170.76 ms | 107.27 ms | 5.62 s | 2.62 | 1.64 | 86.2 | yes |
| mzz-sparse-n8-t30 | 74.06 ms | 733.64 ms | 488.78 ms | 89.36 s | 9.91 | 6.6 | 1.21e+03 | yes |
| mzz-sparse-n8-t100 | 781.73 ms | 42.05 s | 8.74 s | 215.93 s | 53.8 | 11.2 | 276 | yes |
| mzz-wang-f6 | 0.04 ms | 0.02 ms | 0.03 ms | 1.25 ms | 0.576 | 0.802 | 33.3 | yes |
| mzz-wang-f1 | 0.11 ms | 0.03 ms | 0.06 ms | 2.81 ms | 0.291 | 0.536 | 26.5 | yes |
| mzz-wang-f7 | 0.68 ms | 0.70 ms | 1.07 ms | 10.93 ms | 1.04 | 1.58 | 16.2 | yes |
| mzz-wang-f4 | 0.60 ms | 1.95 ms | 1.51 ms | 8.74 ms | 3.28 | 2.53 | 14.7 | yes |
| mzz-wang-f3 | 0.51 ms | 0.77 ms | 0.69 ms | 9.16 ms | 1.52 | 1.36 | 18 | yes |
| mzz-wang-w2 | 0.78 ms | 0.78 ms | 0.86 ms | 13.31 ms | 1 | 1.1 | 17.1 | yes |
| mzz-wang-f2 | 0.51 ms | 0.66 ms | 0.85 ms | 7.69 ms | 1.29 | 1.67 | 15.1 | yes |
| mzz-wang-f5 | 1.74 ms | 1.43 ms | 2.45 ms | 77.40 ms | 0.821 | 1.41 | 44.5 | yes |
| mzz-wang-w1 | 1.14 ms | 1.53 ms | 1.68 ms | 21.74 ms | 1.35 | 1.48 | 19.1 | yes |
