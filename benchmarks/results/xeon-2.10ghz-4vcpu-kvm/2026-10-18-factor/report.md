# Factorization

Records: `2026-10-18-factor`. Times are medians of in-process repetitions on one pinned core. Ratio is the other tool's time over Poly Mojo's (above 1: Poly Mojo is faster). "agree" means every tool that finished returned the same normalized factor multiset and its product equals the input.

| case | polymojo | flint | singular | sympy | flint ratio | singular ratio | sympy ratio | agree |
|---|---|---|---|---|---|---|---|---|
| zz-cyclotomic-360 | 0.06 ms | 73.06 ms | 3.12 s | 115.52 ms | 1.32e+03 | 5.65e+04 | 2.09e+03 | yes |
| zz-cyclotomic-2520 | 0.35 ms | 43.30 s | timeout | 9.75 s | 1.25e+05 |  | 2.81e+04 | yes |
| zz-cyclotomic-5040 | 0.65 ms | timeout | skipped | 43.56 s |  |  | 6.66e+04 | yes |
| zz-product-k2-d10-b16 | 0.29 ms | 0.30 ms | 0.37 ms | 6.27 ms | 1.05 | 1.3 | 21.8 | yes |
| zz-product-k2-d25-b16 | 1.93 ms | 0.76 ms | 0.93 ms | 39.54 ms | 0.395 | 0.482 | 20.5 | yes |
| zz-product-k2-d50-b16 | 7.60 ms | 2.09 ms | 2.73 ms | 231.49 ms | 0.274 | 0.359 | 30.4 | yes |
| zz-product-k2-d100-b16 | 40.66 ms | 12.47 ms | 13.91 ms | 1.75 s | 0.307 | 0.342 | 43.1 | yes |
| zz-product-k2-d10-b256 | 3.16 ms | 1.19 ms | 1.24 ms | 23.34 ms | 0.377 | 0.394 | 7.4 | yes |
| zz-product-k2-d25-b256 | 11.49 ms | 3.17 ms | 3.28 ms | 154.65 ms | 0.276 | 0.285 | 13.5 | yes |
| zz-product-k2-d50-b256 | 46.82 ms | 7.36 ms | 7.68 ms | 755.01 ms | 0.157 | 0.164 | 16.1 | yes |
| zz-product-k2-d100-b256 | 107.20 ms | 31.35 ms | 33.37 ms | 7.13 s | 0.292 | 0.311 | 66.5 | yes |
| zz-product-k8-d10-b16 | 17.03 ms | 6.96 ms | 5.65 ms | 752.55 ms | 0.408 | 0.332 | 44.2 | yes |
| zz-product-k8-d25-b16 | 140.61 ms | 15.39 ms | 17.35 ms | 8.68 s | 0.109 | 0.123 | 61.7 | yes |
| zz-product-k8-d50-b16 | 386.00 ms | 49.05 ms | 58.39 ms | 54.28 s | 0.127 | 0.151 | 141 | yes |
| zz-product-k8-d100-b16 | 1.53 s | 117.99 ms | 146.76 ms | timeout | 0.0772 | 0.096 |  | yes |
| zz-product-k8-d10-b256 | 131.99 ms | 15.68 ms | 19.76 ms | 3.95 s | 0.119 | 0.15 | 30 | yes |
| zz-product-k8-d25-b256 | 386.22 ms | 80.32 ms | 107.12 ms | 54.86 s | 0.208 | 0.277 | 142 | yes |
| zz-product-k8-d50-b256 | 1.52 s | 204.52 ms | 268.07 ms | timeout | 0.135 | 0.177 |  | yes |
| zz-product-k8-d100-b256 | 4.63 s | 533.18 ms | 516.93 ms | skipped | 0.115 | 0.112 |  | yes |
| zz-repeated-d10 | 0.28 ms | 0.15 ms | 0.17 ms | 7.66 ms | 0.517 | 0.597 | 27.2 | yes |
| zz-repeated-d20 | 0.64 ms | 0.36 ms | 0.47 ms | 30.19 ms | 0.564 | 0.743 | 47.3 | yes |
| zz-repeated-d40 | 2.09 ms | 1.41 ms | 1.77 ms | 189.07 ms | 0.676 | 0.847 | 90.5 | yes |
| zz-swinnerton-dyer-4 | 0.25 ms | 0.35 ms | 0.39 ms | 5.72 ms | 1.39 | 1.54 | 22.5 | yes |
| zz-swinnerton-dyer-5 | 4.40 ms | 1.78 ms | 2.08 ms | 148.35 ms | 0.404 | 0.474 | 33.7 | yes |
| zz-swinnerton-dyer-6 | 21.19 ms | 8.34 ms | 10.13 ms | timeout | 0.394 | 0.478 |  | yes |
| zz-swinnerton-dyer-7 | 447.59 ms | 58.30 ms | 197.40 ms | skipped | 0.13 | 0.441 |  | yes |
| gf-p3-d100 | 0.22 ms | 0.31 ms | 0.46 ms | 83.56 ms | 1.42 | 2.12 | 388 | yes |
| gf-p3-d250 | 1.70 ms | 1.41 ms | 2.32 ms | 939.64 ms | 0.83 | 1.37 | 554 | yes |
| gf-p3-d500 | 6.83 ms | 14.40 ms | 27.92 ms | 19.79 s | 2.11 | 4.09 | 2.9e+03 | yes |
| gf-p3-d1000 | 34.60 ms | 33.87 ms | 61.17 ms | 57.27 s | 0.979 | 1.77 | 1.66e+03 | yes |
| gf-p3-d2000 | 213.75 ms | 185.92 ms | 248.24 ms | timeout | 0.87 | 1.16 |  | yes |
| gf-p61-d100 | 4.74 ms | 3.84 ms | unsupported | 1.11 s | 0.811 |  | 233 | yes |
| gf-p61-d250 | 34.93 ms | 58.97 ms | unsupported | 18.91 s | 1.69 |  | 541 | yes |
| gf-p61-d500 | 170.34 ms | 264.04 ms | unsupported | 156.33 s | 1.55 |  | 918 | yes |
| gf-p61-d1000 | 666.75 ms | 960.50 ms | unsupported | timeout | 1.44 |  |  | yes |
| gf-p61-d2000 | 2.24 s | 5.55 s | unsupported | skipped | 2.48 |  |  | yes |
| gf-p65521-d100 | 0.98 ms | 1.22 ms | 1.36 ms | 959.54 ms | 1.24 | 1.39 | 976 | yes |
| gf-p65521-d250 | 8.57 ms | 10.81 ms | 11.67 ms | 17.45 s | 1.26 | 1.36 | 2.04e+03 | yes |
| gf-p65521-d500 | 33.26 ms | 44.64 ms | 48.86 ms | 145.19 s | 1.34 | 1.47 | 4.37e+03 | yes |
| gf-p65521-d1000 | 173.14 ms | 239.42 ms | 185.42 ms | timeout | 1.38 | 1.07 |  | yes |
| gf-p65521-d2000 | 761.92 ms | 1.61 s | 892.35 ms | skipped | 2.11 | 1.17 |  | yes |
| mzz-dense-k5 | 0.33 ms | 0.48 ms | 0.69 ms | 12.50 ms | 1.45 | 2.11 | 38.1 | yes |
| mzz-dense-k7 | 0.39 ms | 0.75 ms | 1.35 ms | 32.80 ms | 1.92 | 3.46 | 84 | yes |
| mzz-dense-k10 | 0.82 ms | 1.06 ms | 2.36 ms | 122.21 ms | 1.29 | 2.87 | 149 | yes |
| mzz-dense-k12 | 1.38 ms | 1.57 ms | 3.91 ms | 250.79 ms | 1.14 | 2.84 | 182 | yes |
| mzz-dense-k15 | 2.63 ms | 2.78 ms | 7.40 ms | 855.17 ms | 1.06 | 2.81 | 325 | yes |
| mzz-sparse-n3-t10 | 2.20 ms | 1.62 ms | 1.69 ms | 30.08 ms | 0.735 | 0.769 | 13.7 | yes |
| mzz-sparse-n3-t30 | 5.17 ms | 2.22 ms | 2.79 ms | 48.12 ms | 0.43 | 0.54 | 9.31 | yes |
| mzz-sparse-n4-t20 | 9.50 ms | 5.77 ms | 5.15 ms | 188.61 ms | 0.607 | 0.542 | 19.8 | yes |
| mzz-sparse-n5-t30 | 30.59 ms | 44.03 ms | 20.91 ms | 1.29 s | 1.44 | 0.683 | 42 | yes |
| mzz-sparse-n6-t30 | 43.26 ms | 104.46 ms | 63.97 ms | 5.18 s | 2.41 | 1.48 | 120 | yes |
| mzz-sparse-n8-t30 | 88.03 ms | 667.28 ms | 455.38 ms | 76.87 s | 7.58 | 5.17 | 873 | yes |
| mzz-sparse-n8-t100 | 849.12 ms | 42.22 s | 7.01 s | 80.89 s | 49.7 | 8.26 | 95.3 | yes |
| mzz-wang-f6 | 0.06 ms | 0.03 ms | 0.03 ms | 1.06 ms | 0.461 | 0.396 | 16.7 | yes |
| mzz-wang-f1 | 0.11 ms | 0.03 ms | 0.08 ms | 1.42 ms | 0.263 | 0.712 | 12.8 | yes |
| mzz-wang-f7 | 0.77 ms | 0.38 ms | 0.73 ms | 9.42 ms | 0.49 | 0.953 | 12.3 | yes |
| mzz-wang-f4 | 0.53 ms | 1.10 ms | 1.20 ms | 6.26 ms | 2.09 | 2.29 | 11.9 | yes |
| mzz-wang-f3 | 0.64 ms | 0.44 ms | 0.48 ms | 5.79 ms | 0.689 | 0.755 | 9.12 | yes |
| mzz-wang-w2 | 0.81 ms | 0.51 ms | 0.47 ms | 8.73 ms | 0.636 | 0.582 | 10.8 | yes |
| mzz-wang-f2 | 0.47 ms | 0.35 ms | 0.39 ms | 4.61 ms | 0.743 | 0.824 | 9.83 | yes |
| mzz-wang-f5 | 1.82 ms | 1.19 ms | 1.09 ms | 67.25 ms | 0.651 | 0.599 | 36.9 | yes |
| mzz-wang-w1 | 1.21 ms | 1.29 ms | 1.42 ms | 21.03 ms | 1.07 | 1.17 | 17.4 | yes |
