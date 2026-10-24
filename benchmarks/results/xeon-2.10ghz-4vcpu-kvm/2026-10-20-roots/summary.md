# Real roots, 2026-10-20 (nightly d6f6b2a; Poly Mojo rerun rows at 55bf6bb)

Pinned single core, 120 s case timeout; ladders stop after a timeout (skipped). records.jsonl is the main run;
polymojo-rerun.jsonl redoes Poly Mojo rows the main run lost to a mid-run revision change (random1 refine53, random64,
random1024, repeated); sympy.jsonl redoes SymPy with CONDA_PREFIX set. Times in ms; ratio = ours / fastest peer
(msolve's whole-process 0.00 rows excluded).

```
case (ms)                       ours     flint    msolve     sympy  ours/best
mignotte-10.count               0.03      0.00unsupporte      0.22     22.11
mignotte-20.count               0.06      0.00unsupporte      0.34     20.66
mignotte-50.count               0.23      0.01unsupporte      0.88     24.76
mignotte-100.count              0.81      0.03unsupporte      2.09     28.72
mignotte-200.count              2.59      0.09unsupporte      6.82     27.84
mignotte-400.count              9.76      0.35unsupporte     25.10     28.03
mignotte-10.exact            timeoutunsupporteunsupporte     error         -
mignotte-20.exact            skippedunsupporteunsupporte   timeout         -
mignotte-50.exact            skippedunsupporteunsupporte   skipped         -
mignotte-100.exact           skippedunsupporteunsupporte   skipped         -
mignotte-200.exact           skippedunsupporteunsupporte   skipped         -
mignotte-400.exact           skippedunsupporteunsupporte   skipped         -
mignotte-10.isolate             0.35      0.80      0.00      0.53      0.66
mignotte-20.isolate             0.66      5.37      0.00      0.89      0.74
mignotte-50.isolate             2.29    139.55     12.81      3.49      0.65
mignotte-100.isolate            5.48   2412.50    468.69     11.12      0.49
mignotte-200.isolate           19.45  51271.25  14497.07     42.08      0.46
mignotte-400.isolate           62.57   timeout   timeout    174.34      0.36
mignotte-10.refine1000          9.56      3.24      0.00   1231.55      2.95
mignotte-20.refine1000         10.63     11.89      0.00   timeout      0.89
mignotte-50.refine1000        193.53    160.20     32.57   skipped      5.94
mignotte-100.refine1000       215.45   2285.92    505.30   skipped      0.43
mignotte-200.refine1000       354.63  50913.97  14986.42   skipped      0.02
mignotte-400.refine1000        72.72   timeout   timeout   skipped         -
mignotte-10.refine53            0.76      0.81      0.00   1058.83      0.94
mignotte-20.refine53            1.00      5.31      0.00      9.32      0.19
mignotte-50.refine53            2.73    144.73     10.13     54.05      0.27
mignotte-100.refine53           6.50   2359.43    487.13    174.97      0.04
mignotte-200.refine53          20.90  51179.74  14884.80   1188.38      0.02
mignotte-400.refine53          65.91   timeout   timeout  14065.83      0.00
wilkinson-20.count              0.21      0.04unsupporte      0.88      5.62
wilkinson-40.count              1.87      1.85unsupporte      2.74      1.01
wilkinson-80.count             12.80     90.92unsupporte     11.80      1.08
wilkinson-160.count            90.16   4104.85unsupporte     75.07      1.20
wilkinson-320.count           688.73 115734.47unsupporte    645.23      1.07
wilkinson-640.count          4955.74   timeoutunsupporte   6581.84      0.75
wilkinson-20.exact              3.65unsupporteunsupporte     27.45      0.13
wilkinson-40.exact             21.61unsupporteunsupporte     81.91      0.26
wilkinson-80.exact            118.08unsupporteunsupporte    309.75      0.38
wilkinson-160.exact          1450.25unsupporteunsupporte   1561.85      0.93
wilkinson-320.exact         48898.66unsupporteunsupporte  10315.44      4.74
wilkinson-640.exact          timeoutunsupporteunsupporte  89048.17         -
wilkinson-20.isolate            0.91      2.61      0.00      0.94      0.96
wilkinson-40.isolate            4.56     23.65      0.00      4.44      1.03
wilkinson-80.isolate           21.52    534.50     17.97     28.02      1.20
wilkinson-160.isolate         155.73   7448.94     90.17    223.37      1.73
wilkinson-320.isolate        1488.72   timeout    414.69   2194.03      3.59
wilkinson-640.isolate       20389.04   skipped   2947.51  23497.67      6.92
wilkinson-20.refine1000        51.96     11.95      0.00      0.91     57.38
wilkinson-40.refine1000       134.24     55.90      0.00      4.62     29.03
wilkinson-80.refine1000       384.81    678.80     12.99     27.67     29.62
wilkinson-160.refine1000     1908.59   8059.76     87.03    216.79     21.93
wilkinson-320.refine1000    18061.15   timeout    439.64   2255.35     41.08
wilkinson-640.refine1000     timeout   skipped   3144.00  22872.41         -
wilkinson-20.refine53           6.26      3.10      0.00      0.98      6.39
wilkinson-40.refine53          22.78     26.96      0.00      4.66      4.89
wilkinson-80.refine53         109.59    553.55     10.36     29.19     10.58
wilkinson-160.refine53        774.23   7746.71     81.65    225.53      9.48
wilkinson-320.refine53       9640.19   timeout    410.50   2274.63     23.48
wilkinson-640.refine53       timeout   skipped   2855.08  22726.56         -
chebyshev-100.count             2.27     33.31unsupporte      5.35      0.42
chebyshev-200.count            10.20   1067.37unsupporte     22.29      0.46
chebyshev-400.count            46.56  30912.21unsupporte     96.12      0.48
chebyshev-700.count           163.29   timeoutunsupporte    359.54      0.45
chebyshev-1000.count          384.93   skippedunsupporte    849.73      0.45
chebyshev-100.exact          timeoutunsupporteunsupporte     error         -
chebyshev-200.exact          skippedunsupporteunsupporte     error         -
chebyshev-400.exact          skippedunsupporteunsupporte     error         -
chebyshev-700.exact          skippedunsupporteunsupporte   timeout         -
chebyshev-1000.exact         skippedunsupporteunsupporte   skipped         -
chebyshev-100.isolate          10.74     79.33     36.01    361.01      0.30
chebyshev-200.isolate          49.71    879.72    211.70   3367.91      0.23
chebyshev-400.isolate         268.18  11916.93   2798.90  36887.28      0.10
chebyshev-700.isolate        1053.09   timeout  31400.53   timeout      0.03
chebyshev-1000.isolate       2638.03   skipped   timeout   skipped         -
chebyshev-100.refine1000      362.01    128.91   1194.05   timeout      2.81
chebyshev-200.refine1000     1279.24   1029.42   9013.36   skipped      1.24
chebyshev-400.refine1000     3494.05  12458.02  72909.42   skipped      0.28
chebyshev-700.refine1000    13759.66   timeout   timeout   skipped         -
chebyshev-1000.refine1000   30714.60   skipped   skipped   skipped         -
chebyshev-100.refine53         59.23     80.59     41.19   7989.64      1.44
chebyshev-200.refine53        249.84    874.00    235.00  75138.00      1.06
chebyshev-400.refine53       1150.87  12187.13   2777.77   timeout      0.41
chebyshev-700.refine53       4794.80   timeout  31406.04   skipped      0.15
chebyshev-1000.refine53     12960.45   skipped   timeout   skipped         -
random1-100.count             168.18      1.65unsupporte   1506.21    101.69
random1-200.count            3242.55     23.53unsupporte  65682.43    137.80
random1-500.count            timeout    858.06unsupporte   timeout         -
random1-1000.count           skipped  12760.18unsupporte   skipped         -
random1-2000.count           skipped   timeoutunsupporte   skipped         -
random1-100.exact            timeoutunsupporteunsupporte     error         -
random1-200.exact            skippedunsupporteunsupporte     error         -
random1-500.exact            skippedunsupporteunsupporte     error         -
random1-1000.exact           skippedunsupporteunsupporte   timeout         -
random1-2000.exact           skippedunsupporteunsupporte   skipped         -
random1-100.isolate           109.34     28.07      3.24      8.10     33.72
random1-200.isolate          2798.83    185.14     21.22     31.72    131.89
random1-500.isolate          timeout   1812.27    144.17    243.10         -
random1-1000.isolate         skipped   8516.41    711.94    592.32         -
random1-2000.isolate         skipped  64090.11  12167.68   6086.33         -
random1-100.refine1000        189.44    207.42     86.99  36771.71      2.18
random1-200.refine1000       5836.81    868.59    106.58   timeout     54.77
random1-500.refine1000       timeout   5062.17   1161.75   skipped         -
random1-1000.refine1000      skipped  20892.36   3889.95   skipped         -
random1-2000.refine1000      skipped  98964.79  25907.32   skipped         -
random1-100.refine53          180.74     36.36      1.88    291.88     96.00
random1-200.refine53         5690.53    220.75     24.12    784.86    235.97
random1-500.refine53         timeout   1849.34    132.03  41397.31         -
random1-1000.refine53        skipped  11150.00    653.84  51100.94         -
random1-2000.refine53        skipped  62219.78  11927.99   timeout         -
random1024-100.count         timeout   5837.59unsupporte   timeout         -
random1024-200.count         skipped  53816.70unsupporte   skipped         -
random1024-500.count         skipped   timeoutunsupporte   skipped         -
random1024-1000.count        skipped   skippedunsupporte   skipped         -
random1024-2000.count        skipped   skippedunsupporte   skipped         -
random1024-100.exact         timeoutunsupporteunsupporte     error         -
random1024-200.exact         skippedunsupporteunsupporte     error         -
random1024-500.exact         skippedunsupporteunsupporte   timeout         -
random1024-1000.exact        skippedunsupporteunsupporte   skipped         -
random1024-2000.exact        skippedunsupporteunsupporte   skipped         -
random1024-100.isolate       timeout     32.96    136.16     20.24         -
random1024-200.isolate       skipped    174.64    193.34     68.30         -
random1024-500.isolate       skipped   1799.43    555.21    360.45         -
random1024-1000.isolate      skipped   6959.10   2430.53   2444.44         -
random1024-2000.isolate      skipped  59500.03  40092.58  36800.56         -
random1024-100.refine1000    timeout    270.98    203.65  65998.40         -
random1024-200.refine1000    skipped   1061.15    388.14   timeout         -
random1024-500.refine1000    skipped   6764.25   1675.60   skipped         -
random1024-1000.refine1000   skipped  24571.10   7223.60   skipped         -
random1024-2000.refine1000   skipped 108798.73  99017.63   skipped         -
random1024-100.refine53      timeout     39.42    132.35    570.36         -
random1024-200.refine53      skipped    203.53    184.62   1781.88         -
random1024-500.refine53      skipped   1813.74    556.51  16404.52         -
random1024-1000.refine53     skipped  11924.29   2343.91  89498.07         -
random1024-2000.refine53     skipped 109290.25  38692.86   timeout         -
random64-100.count          12077.98    104.50unsupporte 105209.87    115.58
random64-200.count           timeout   1220.22unsupporte   timeout         -
random64-500.count           skipped  31446.60unsupporte   skipped         -
random64-1000.count          skipped   timeoutunsupporte   skipped         -
random64-2000.count          skipped   skippedunsupporte   skipped         -
random64-100.exact          23331.45unsupporteunsupporte    336.79     69.28
random64-200.exact           timeoutunsupporteunsupporte     error         -
random64-500.exact           skippedunsupporteunsupporte     error         -
random64-1000.exact          skippedunsupporteunsupporte   timeout         -
random64-2000.exact          skippedunsupporteunsupporte   skipped         -
random64-100.isolate        11684.48     41.49      9.13      6.76   1727.22
random64-200.isolate         timeout    200.22     33.64     38.67         -
random64-500.isolate         skipped   1985.59    187.38    413.00         -
random64-1000.isolate        skipped   8929.23   1639.57   2210.10         -
random64-2000.isolate        skipped  48323.27  13330.58   4644.28         -
random64-100.refine1000     26411.36    218.69     11.67      6.84   3860.67
random64-200.refine1000      timeout    857.44    116.97   timeout         -
random64-500.refine1000      skipped   5126.28   1519.66   skipped         -
random64-1000.refine1000     skipped  20726.95  10463.05   skipped         -
random64-2000.refine1000     skipped  84254.27  39726.72   skipped         -
random64-100.refine53       26238.04     50.07      8.01      6.99   3756.10
random64-200.refine53        timeout    245.47     30.80    779.14         -
random64-500.refine53        skipped   1942.00    176.49  17003.54         -
random64-1000.refine53       skipped  12026.79   1654.77   timeout         -
random64-2000.refine53       skipped  48458.10  13393.11   skipped         -
repeated-4.count                1.09      0.03unsupporte      0.56     43.16
repeated-8.count               11.47      0.39unsupporte      3.62     29.31
repeated-16.count             104.74     14.88unsupporte    140.15      7.04
repeated-32.count            5711.87    581.23unsupporte  12399.33      9.83
repeated-64.count            timeout  20231.19unsupporte   timeout         -
repeated-4.exact              391.11unsupporteunsupporte     error         -
repeated-8.exact             timeoutunsupporteunsupporte     error         -
repeated-16.exact            skippedunsupporteunsupporte     error         -
repeated-32.exact            skippedunsupporteunsupporte     error         -
repeated-64.exact            skippedunsupporteunsupporte     error         -
repeated-4.isolate              0.87     error      1.35     wrong      0.65
repeated-8.isolate              5.94     error      3.42     wrong      1.74
repeated-16.isolate           132.38     error     19.50     wrong      6.79
repeated-32.isolate          6939.82     error     34.82     wrong    199.32
repeated-64.isolate          timeout     error    123.69     wrong         -
repeated-4.refine1000          17.67     error     16.76     91.53      1.05
repeated-8.refine1000          40.32     error     14.69    435.11      2.75
repeated-16.refine1000        352.80     error     40.65   2082.48      8.68
repeated-32.refine1000      21543.89     error    119.40  19029.21    180.43
repeated-64.refine1000       timeout     error    797.81   timeout         -
repeated-4.refine53             2.60     error      4.33      3.53      0.74
repeated-8.refine53            11.12     error     11.04     13.86      1.01
repeated-16.refine53          263.87     error     21.71     47.38     12.16
repeated-32.refine53        21449.05     error     31.33    270.76    684.51
repeated-64.refine53         timeout     error    109.18   2054.62         -
```
