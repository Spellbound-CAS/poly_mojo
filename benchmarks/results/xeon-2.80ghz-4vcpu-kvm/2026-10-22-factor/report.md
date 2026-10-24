# Factorization

Records: `2026-10-22-factor`. Times are medians of in-process repetitions on one pinned core. Ratio is the other tool's time over Poly Mojo's (above 1: Poly Mojo is faster). "agree" means every tool that finished returned the same normalized factor multiset and its product equals the input.

| case | polymojo | flint | singular | sympy | flint ratio | singular ratio | sympy ratio | agree |
|---|---|---|---|---|---|---|---|---|
| zz-cyclotomic-360 | 0.07 ms | 82.15 ms | 3.52 s | 167.26 ms | 1.23e+03 | 5.28e+04 | 2.51e+03 | yes |
| zz-cyclotomic-2520 | 0.43 ms | 52.42 s | timeout | 12.32 s | 1.23e+05 |  | 2.88e+04 | yes |
| zz-cyclotomic-5040 | 0.83 ms | timeout | skipped | 53.96 s |  |  | 6.51e+04 | yes |
| zz-product-k2-d10-b16 | 0.19 ms | 0.35 ms | 0.44 ms | 7.39 ms | 1.81 | 2.29 | 38.3 | yes |
| zz-product-k2-d25-b16 | 2.07 ms | 0.88 ms | 1.05 ms | 49.48 ms | 0.427 | 0.506 | 23.9 | yes |
| zz-product-k2-d50-b16 | 8.80 ms | 2.53 ms | 3.12 ms | 279.17 ms | 0.287 | 0.354 | 31.7 | yes |
| zz-product-k2-d100-b16 | 41.41 ms | 13.84 ms | 15.87 ms | 2.20 s | 0.334 | 0.383 | 53.1 | yes |
| zz-product-k2-d10-b256 | 2.27 ms | 1.30 ms | 1.35 ms | 27.22 ms | 0.575 | 0.595 | 12 | yes |
| zz-product-k2-d25-b256 | 5.72 ms | 3.72 ms | 3.71 ms | 216.97 ms | 0.651 | 0.648 | 37.9 | yes |
| zz-product-k2-d50-b256 | 25.55 ms | 8.57 ms | 8.51 ms | 808.76 ms | 0.335 | 0.333 | 31.6 | yes |
| zz-product-k2-d100-b256 | 73.79 ms | 37.51 ms | 37.20 ms | 8.08 s | 0.508 | 0.504 | 109 | yes |
| zz-product-k8-d10-b16 | 17.74 ms | 6.05 ms | 6.69 ms | 1.06 s | 0.341 | 0.377 | 59.9 | yes |
| zz-product-k8-d25-b16 | 92.73 ms | 17.14 ms | 19.74 ms | 10.79 s | 0.185 | 0.213 | 116 | yes |
| zz-product-k8-d50-b16 | 234.50 ms | 55.88 ms | 69.42 ms | 68.08 s | 0.238 | 0.296 | 290 | yes |
| zz-product-k8-d100-b16 | 771.25 ms | 126.95 ms | 179.00 ms | timeout | 0.165 | 0.232 |  | yes |
| zz-product-k8-d10-b256 | 54.74 ms | 21.83 ms | 22.10 ms | 3.97 s | 0.399 | 0.404 | 72.5 | yes |
| zz-product-k8-d25-b256 | 206.99 ms | 81.73 ms | 83.21 ms | 57.84 s | 0.395 | 0.402 | 279 | yes |
| zz-product-k8-d50-b256 | 488.79 ms | 217.52 ms | 279.41 ms | timeout | 0.445 | 0.572 |  | yes |
| zz-product-k8-d100-b256 | 1.16 s | 566.13 ms | 568.60 ms | skipped | 0.488 | 0.49 |  | yes |
| zz-repeated-d10 | 0.27 ms | 0.17 ms | 0.21 ms | 11.26 ms | 0.606 | 0.76 | 41.3 | yes |
| zz-repeated-d20 | 0.69 ms | 0.40 ms | 0.53 ms | 37.41 ms | 0.584 | 0.764 | 53.9 | yes |
| zz-repeated-d40 | 2.11 ms | 1.56 ms | 1.86 ms | 210.88 ms | 0.74 | 0.879 | 99.7 | yes |
| zz-swinnerton-dyer-4 | 0.30 ms | 0.42 ms | 0.47 ms | 6.84 ms | 1.38 | 1.55 | 22.5 | yes |
| zz-swinnerton-dyer-5 | 5.70 ms | 2.10 ms | 2.57 ms | 159.03 ms | 0.368 | 0.45 | 27.9 | yes |
| zz-swinnerton-dyer-6 | 12.73 ms | 9.38 ms | 12.17 ms | timeout | 0.737 | 0.956 |  | yes |
| zz-swinnerton-dyer-7 | 93.78 ms | 69.25 ms | 207.84 ms | skipped | 0.738 | 2.22 |  | yes |
| gf-p3-d100 | 0.38 ms | 0.36 ms | 0.52 ms | 103.56 ms | 0.925 | 1.35 | 269 | yes |
| gf-p3-d250 | 2.46 ms | 1.78 ms | 2.65 ms | 991.14 ms | 0.723 | 1.08 | 403 | yes |
| gf-p3-d500 | 5.15 ms | 16.64 ms | 27.74 ms | 22.83 s | 3.23 | 5.39 | 4.43e+03 | yes |
| gf-p3-d1000 | 41.85 ms | 38.53 ms | 74.57 ms | 77.05 s | 0.921 | 1.78 | 1.84e+03 | yes |
| gf-p3-d2000 | 176.12 ms | 239.19 ms | 309.40 ms | timeout | 1.36 | 1.76 |  | yes |
| gf-p61-d100 | 6.05 ms | 4.44 ms | unsupported | 1.48 s | 0.735 |  | 245 | yes |
| gf-p61-d250 | 45.61 ms | 76.32 ms | unsupported | 23.94 s | 1.67 |  | 525 | yes |
| gf-p61-d500 | 217.41 ms | 234.41 ms | unsupported | 201.39 s | 1.08 |  | 926 | yes |
| gf-p61-d1000 | 791.01 ms | 1.12 s | unsupported | timeout | 1.42 |  |  | yes |
| gf-p61-d2000 | 3.00 s | 6.01 s | unsupported | skipped | 2 |  |  | yes |
| gf-p65521-d100 | 1.22 ms | 1.40 ms | 1.57 ms | 1.22 s | 1.15 | 1.29 | 999 | yes |
| gf-p65521-d250 | 12.76 ms | 11.93 ms | 12.87 ms | 25.12 s | 0.935 | 1.01 | 1.97e+03 | yes |
| gf-p65521-d500 | 47.04 ms | 59.59 ms | 54.93 ms | 198.62 s | 1.27 | 1.17 | 4.22e+03 | yes |
| gf-p65521-d1000 | 232.44 ms | 231.63 ms | 203.25 ms | timeout | 0.997 | 0.874 |  | yes |
| gf-p65521-d2000 | 941.16 ms | 1.41 s | 1.12 s | skipped | 1.49 | 1.19 |  | yes |
| mzz-dense-k5 | 0.23 ms | 0.50 ms | 0.75 ms | 15.19 ms | 2.23 | 3.33 | 67.1 | yes |
| mzz-dense-k7 | 0.40 ms | 0.88 ms | 1.47 ms | 38.01 ms | 2.19 | 3.68 | 95.1 | yes |
| mzz-dense-k10 | 0.92 ms | 1.25 ms | 2.72 ms | 140.17 ms | 1.36 | 2.96 | 153 | yes |
| mzz-dense-k12 | 1.53 ms | 1.91 ms | 4.33 ms | 320.16 ms | 1.24 | 2.82 | 209 | yes |
| mzz-dense-k15 | 2.52 ms | 3.28 ms | 8.54 ms | 819.30 ms | 1.3 | 3.39 | 325 | yes |
| mzz-sparse-n3-t10 | 2.52 ms | 1.86 ms | 1.83 ms | 33.06 ms | 0.737 | 0.724 | 13.1 | yes |
| mzz-sparse-n3-t30 | 4.72 ms | 2.52 ms | 3.05 ms | 55.63 ms | 0.534 | 0.646 | 11.8 | yes |
| mzz-sparse-n4-t20 | 10.23 ms | 4.83 ms | 5.92 ms | 204.11 ms | 0.472 | 0.578 | 19.9 | yes |
| mzz-sparse-n5-t30 | 28.79 ms | 32.92 ms | 15.27 ms | 1.32 s | 1.14 | 0.531 | 45.7 | yes |
| mzz-sparse-n6-t30 | 43.96 ms | 118.27 ms | 66.07 ms | 5.76 s | 2.69 | 1.5 | 131 | yes |
| mzz-sparse-n8-t30 | 82.37 ms | 772.20 ms | 445.89 ms | 85.92 s | 9.37 | 5.41 | 1.04e+03 | yes |
| mzz-sparse-n8-t100 | 802.69 ms | 42.30 s | 8.07 s | 144.50 s | 52.7 | 10.1 | 180 | yes |
| mzz-wang-f6 | 0.04 ms | 0.02 ms | 0.03 ms | 1.27 ms | 0.521 | 0.788 | 31.3 | yes |
| mzz-wang-f1 | 0.12 ms | 0.03 ms | 0.06 ms | 1.77 ms | 0.25 | 0.475 | 14.2 | yes |
| mzz-wang-f7 | 0.79 ms | 0.43 ms | 0.73 ms | 10.59 ms | 0.54 | 0.934 | 13.5 | yes |
| mzz-wang-f4 | 0.57 ms | 1.21 ms | 1.27 ms | 7.89 ms | 2.1 | 2.22 | 13.7 | yes |
| mzz-wang-f3 | 0.54 ms | 0.50 ms | 0.52 ms | 7.24 ms | 0.915 | 0.961 | 13.4 | yes |
| mzz-wang-w2 | 0.84 ms | 0.51 ms | 0.51 ms | 10.92 ms | 0.612 | 0.602 | 13 | yes |
| mzz-wang-f2 | 0.48 ms | 0.40 ms | 0.45 ms | 5.52 ms | 0.844 | 0.942 | 11.6 | yes |
| mzz-wang-f5 | 1.62 ms | 1.29 ms | 1.51 ms | 85.00 ms | 0.798 | 0.929 | 52.5 | yes |
| mzz-wang-w1 | 1.18 ms | 1.37 ms | 1.44 ms | 22.15 ms | 1.16 | 1.21 | 18.7 | yes |
