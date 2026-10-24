# Real roots, 2026-10-21 (nightly 841d966: Descartes real roots by default)

Clean worktree at 841d966, same harness settings as 2026-10-20-roots: pinned CPU 1, 120 s case timeout, 12G cap,
ladders stop after a timeout (skipped), CONDA_PREFIX set for every tool.

- records.jsonl: groups 1-2 (mignotte, wilkinson, chebyshev, random1) with all five tools; group 3 (random64,
  random1024, repeated) without SymPy, at the coordinator's request.
- rerun-mignotte-wilkinson.jsonl: group 1 again without SymPy, because a documentation build overlapped its first pass.
  It overrides records.jsonl for those rows. Group-1 SymPy ratios therefore divide first-pass SymPy times by rerun times.

Ratios are peer time / Poly Mojo time inside one run, so host speed cancels (>1: Poly Mojo faster).
Old = 2026-10-20-roots (records + polymojo-rerun + sympy files). msolve whole-process zero rows are excluded.
No Poly Mojo answer is wrong; outcomes: 142 ok, 11 timeout, 32 skipped.

## Geometric-mean competitor ratios, old -> new

| mode | FLINT/ours | msolve/ours | SymPy/ours |
|---|---|---|---|
| isolate | 7.26 -> 114 (15) | 0.994 -> 17 (14) | 1.84 -> 24.4 (17) |
| refine53 | 2.83 -> 74.7 (15) | 0.279 -> 15.8 (14) | 5.98 -> 94.8 (15) |
| refine1000 | 0.951 -> 33.6 (15) | 0.359 -> 7.68 (13) | 0.55 -> 52.8 (7) |
| count (Sturm) | 0.35 -> 0.375 (21) | unsupported | 2.74 -> 4.02 (19) |
| exact | unsupported | unsupported | 1.76 -> 0.21 (5) |

Rows where a peer still leads (isolate/refine): mignotte-50.refine1000 (msolve 0.25, FLINT 0.77),
mignotte-20.refine1000 (msolve 0.21), wilkinson-20.refine1000 (msolve 0.40), repeated-16.refine1000 (msolve 0.87),
random64-500.isolate (msolve 0.94), repeated-32.isolate (msolve 0.99).

The exact-mode regression (wilkinson-20.exact 3.6 -> 103 ms) is in algebraic_roots, not isolation: on Wilkinson 20,
algebraic_roots takes 3.3 ms at d6f6b2a, 70 ms at 11df858 (CAD revision 9, before Descartes) and 76 ms at 841d966,
while isolate_real_roots takes 0.58, 0.54 and 0.21 ms.

## Per case

```
case                         ours old ms ours new ms FLINT/ours old FLINT/ours new msolve/ours old msolve/ours new SymPy/ours old SymPy/ours new
chebyshev-100.count                2.271       1.814           14.7             16    unsupported    unsupported           2.35           2.79
chebyshev-200.count               10.201       8.402            105            110    unsupported    unsupported           2.18           2.39
chebyshev-400.count               46.561      38.958            664            661    unsupported    unsupported           2.06           2.16
chebyshev-700.count              163.291     143.459        timeout        timeout    unsupported    unsupported            2.2           2.19
chebyshev-1000.count             384.935     335.606        skipped        skipped    unsupported    unsupported           2.21           2.33
mignotte-10.count                  0.025       0.022         0.0452         0.0724    unsupported    unsupported           8.84           8.48
mignotte-20.count                  0.055       0.040         0.0484         0.0528    unsupported    unsupported           6.21           7.57
mignotte-50.count                  0.228       0.142         0.0404         0.0503    unsupported    unsupported           3.86           5.79
mignotte-100.count                 0.807       0.526         0.0348         0.0406    unsupported    unsupported           2.59           3.76
mignotte-200.count                 2.590       1.705         0.0359         0.0418    unsupported    unsupported           2.63           8.34
mignotte-400.count                 9.760       6.848         0.0357         0.0386    unsupported    unsupported           2.57           6.85
random1-100.count                168.175     148.291        0.00983        0.00957    unsupported    unsupported           8.96            8.7
random1-200.count               3242.548    2735.861        0.00726        0.00742    unsupported    unsupported           20.3           20.9
random1-500.count                timeout     timeout             ok             ok    unsupported    unsupported        timeout        timeout
random1-1000.count               skipped     skipped             ok             ok    unsupported    unsupported        skipped        skipped
random1-2000.count               skipped     skipped        timeout        timeout    unsupported    unsupported        skipped        skipped
random1024-100.count             timeout     timeout             ok             ok    unsupported    unsupported        timeout        timeout
random1024-200.count             skipped     skipped             ok             ok    unsupported    unsupported        skipped        skipped
random1024-500.count             skipped     skipped        timeout        timeout    unsupported    unsupported        skipped              -
random1024-1000.count            skipped     skipped        skipped        skipped    unsupported    unsupported        skipped              -
random1024-2000.count            skipped     skipped        skipped        skipped    unsupported    unsupported        skipped              -
random64-100.count             12077.978   23172.123        0.00865        0.00862    unsupported    unsupported           8.71              -
random64-200.count               timeout     timeout             ok             ok    unsupported    unsupported        timeout              -
random64-500.count               skipped     skipped             ok             ok    unsupported    unsupported        skipped              -
random64-1000.count              skipped     skipped        timeout        timeout    unsupported    unsupported        skipped              -
random64-2000.count              skipped     skipped        skipped        skipped    unsupported    unsupported        skipped              -
repeated-4.count                   1.086       0.975         0.0232         0.0228    unsupported    unsupported          0.513              -
repeated-8.count                  11.466       9.268         0.0341         0.0373    unsupported    unsupported          0.316              -
repeated-16.count                104.736      93.161          0.142          0.146    unsupported    unsupported           1.34              -
repeated-32.count               5711.871    4995.648          0.102          0.105    unsupported    unsupported           2.17              -
repeated-64.count                timeout     timeout             ok             ok    unsupported    unsupported        timeout              -
wilkinson-20.count                 0.207       0.170          0.178          0.184    unsupported    unsupported           4.28           4.46
wilkinson-40.count                 1.874       1.594          0.986           1.03    unsupported    unsupported           1.46           4.03
wilkinson-80.count                12.804      10.610            7.1           7.39    unsupported    unsupported          0.922            2.1
wilkinson-160.count               90.157      79.160           45.5           41.2    unsupported    unsupported          0.833            1.8
wilkinson-320.count              688.728     621.166            168            170    unsupported    unsupported          0.937            1.8
wilkinson-640.count             4955.745    4635.357        timeout        timeout    unsupported    unsupported           1.33           2.38
chebyshev-100.exact              timeout     timeout    unsupported    unsupported    unsupported    unsupported          error          error
chebyshev-200.exact              skipped     skipped    unsupported    unsupported    unsupported    unsupported          error          error
chebyshev-400.exact              skipped     skipped    unsupported    unsupported    unsupported    unsupported          error          error
chebyshev-700.exact              skipped     skipped    unsupported    unsupported    unsupported    unsupported        timeout        timeout
chebyshev-1000.exact             skipped     skipped    unsupported    unsupported    unsupported    unsupported        skipped        skipped
mignotte-10.exact                timeout     timeout    unsupported    unsupported    unsupported    unsupported          error          error
mignotte-20.exact                skipped     skipped    unsupported    unsupported    unsupported    unsupported        timeout        timeout
mignotte-50.exact                skipped     skipped    unsupported    unsupported    unsupported    unsupported        skipped        skipped
mignotte-100.exact               skipped     skipped    unsupported    unsupported    unsupported    unsupported        skipped        skipped
mignotte-200.exact               skipped     skipped    unsupported    unsupported    unsupported    unsupported        skipped        skipped
mignotte-400.exact               skipped     skipped    unsupported    unsupported    unsupported    unsupported        skipped        skipped
random1-100.exact                timeout     timeout    unsupported    unsupported    unsupported    unsupported          error          error
random1-200.exact                skipped     skipped    unsupported    unsupported    unsupported    unsupported          error          error
random1-500.exact                skipped     skipped    unsupported    unsupported    unsupported    unsupported          error          error
random1-1000.exact               skipped     skipped    unsupported    unsupported    unsupported    unsupported        timeout        timeout
random1-2000.exact               skipped     skipped    unsupported    unsupported    unsupported    unsupported        skipped        skipped
random1024-100.exact             timeout     timeout    unsupported    unsupported    unsupported    unsupported          error              -
random1024-200.exact             skipped     skipped    unsupported    unsupported    unsupported    unsupported          error              -
random1024-500.exact             skipped     skipped    unsupported    unsupported    unsupported    unsupported        timeout              -
random1024-1000.exact            skipped     skipped    unsupported    unsupported    unsupported    unsupported        skipped              -
random1024-2000.exact            skipped     skipped    unsupported    unsupported    unsupported    unsupported        skipped              -
random64-100.exact             23331.449   11288.361    unsupported    unsupported    unsupported    unsupported         0.0144              -
random64-200.exact               timeout     timeout    unsupported    unsupported    unsupported    unsupported          error              -
random64-500.exact               skipped     skipped    unsupported    unsupported    unsupported    unsupported          error              -
random64-1000.exact              skipped     skipped    unsupported    unsupported    unsupported    unsupported        timeout              -
random64-2000.exact              skipped     skipped    unsupported    unsupported    unsupported    unsupported        skipped              -
repeated-4.exact                 391.106    2268.326    unsupported    unsupported    unsupported    unsupported          error              -
repeated-8.exact                 timeout     timeout    unsupported    unsupported    unsupported    unsupported          error              -
repeated-16.exact                skipped     skipped    unsupported    unsupported    unsupported    unsupported          error              -
repeated-32.exact                skipped     skipped    unsupported    unsupported    unsupported    unsupported          error              -
repeated-64.exact                skipped     skipped    unsupported    unsupported    unsupported    unsupported          error              -
wilkinson-20.exact                 3.645     102.590    unsupported    unsupported    unsupported    unsupported           7.53          0.529
wilkinson-40.exact                21.613     640.151    unsupported    unsupported    unsupported    unsupported           3.79          0.176
wilkinson-80.exact               118.080    3265.331    unsupported    unsupported    unsupported    unsupported           2.62          0.173
wilkinson-160.exact             1450.250   16130.232    unsupported    unsupported    unsupported    unsupported           1.08          0.156
wilkinson-320.exact            48898.659  101168.972    unsupported    unsupported    unsupported    unsupported          0.211          0.162
wilkinson-640.exact              timeout     timeout    unsupported    unsupported    unsupported    unsupported             ok        timeout
chebyshev-100.isolate             10.739       2.460           7.39             28           3.35           17.7           33.6            139
chebyshev-200.isolate             49.707      19.554           17.7           39.6           4.26           8.26           67.8            144
chebyshev-400.isolate            268.179     158.667           44.4           62.5           10.4           15.3            138            192
chebyshev-700.isolate           1053.085    1729.892        timeout        timeout           29.8             16        timeout        timeout
chebyshev-1000.isolate          2638.032    4925.075        skipped        skipped        timeout        timeout        skipped        skipped
mignotte-10.isolate                0.352       0.066           2.28           10.8             ok           23.3           1.52           7.24
mignotte-20.isolate                0.658       0.057           8.16             86             ok             ok           1.35           14.2
mignotte-50.isolate                2.286       0.193           61.1            683           5.61            186           1.53           36.3
mignotte-100.isolate               5.477       0.659            441           3257           85.6            848           2.03           33.1
mignotte-200.isolate              19.452       3.335           2636          14022            745           4765           2.16           23.4
mignotte-400.isolate              62.568      17.536        timeout        timeout        timeout        timeout           2.79           19.1
random1-100.isolate              109.338       0.544          0.257           45.9         0.0297             ok         0.0741           13.9
random1-200.isolate             2798.832       2.550         0.0662           64.3        0.00758           7.84         0.0113             11
random1-500.isolate              timeout      23.269             ok           66.8             ok           5.18             ok           9.46
random1-1000.isolate             skipped     198.717             ok           36.8             ok           2.97             ok           2.69
random1-2000.isolate             skipped    1673.968             ok           32.8             ok           6.23             ok           3.38
random1024-100.isolate           timeout       2.681             ok           23.2             ok            123             ok              -
random1024-200.isolate           skipped      15.400             ok           21.4             ok           26.6             ok              -
random1024-500.isolate           skipped     106.180             ok           31.4             ok           10.6             ok              -
random1024-1000.isolate          skipped     758.167             ok           16.2             ok              6             ok              -
random1024-2000.isolate          skipped   11365.150             ok           9.79             ok           6.52             ok              -
random64-100.isolate           11684.484       1.569        0.00355           25.1       0.000782           3.05       0.000579              -
random64-200.isolate             timeout       6.551             ok           27.8             ok           1.38             ok              -
random64-500.isolate             skipped     173.003             ok             11             ok          0.943             ok              -
random64-1000.isolate            skipped     296.531             ok           28.6             ok           5.27             ok              -
random64-2000.isolate            skipped    1907.003             ok           24.9             ok            6.8             ok              -
repeated-4.isolate                 0.873       0.090          error          error           1.54             ok          wrong              -
repeated-8.isolate                 5.937       0.248          error          error          0.576             ok          wrong              -
repeated-16.isolate              132.381       1.332          error          error          0.147             ok          wrong              -
repeated-32.isolate             6939.818       6.122          error          error        0.00502          0.991          wrong              -
repeated-64.isolate              timeout      35.661          error          error             ok            2.4          wrong              -
wilkinson-20.isolate               0.908       0.120           2.88           19.2             ok           33.9           1.04           6.59
wilkinson-40.isolate               4.559       0.490           5.19             43             ok           40.4          0.975             16
wilkinson-80.isolate              21.525       2.585           24.8            180          0.835           13.3            1.3           18.8
wilkinson-160.isolate            155.733      19.731           47.8            339          0.579           5.36           1.43           18.6
wilkinson-320.isolate           1488.723     189.110        timeout            599          0.279           2.22           1.47           19.1
wilkinson-640.isolate          20389.041    1842.745        skipped        timeout          0.145           1.64           1.15             17
chebyshev-100.refine1000         362.010      58.929          0.356           1.88            3.3           16.9        timeout        timeout
chebyshev-200.refine1000        1279.242     234.930          0.805            3.9           7.05           32.2        skipped        skipped
chebyshev-400.refine1000        3494.050     937.408           3.57           10.8           20.9           65.7        skipped        skipped
chebyshev-700.refine1000       13759.663    3685.705        timeout        timeout        timeout        timeout        skipped        skipped
chebyshev-1000.refine1000      30714.596    8950.497        skipped        skipped        skipped        skipped        skipped        skipped
mignotte-10.refine1000             9.561       1.316          0.338           2.18             ok           4.24            129           1691
mignotte-20.refine1000            10.627       2.199           1.12            4.8             ok          0.207        timeout        timeout
mignotte-50.refine1000           193.535     183.618          0.828           0.77          0.168           0.25        skipped        skipped
mignotte-100.refine1000          215.452     323.312           10.6           6.79           2.35           1.53        skipped        skipped
mignotte-200.refine1000          354.633     449.802            144            101           42.3           34.8        skipped        skipped
mignotte-400.refine1000           72.723     116.011        timeout        timeout        timeout        timeout        skipped        skipped
random1-100.refine1000           189.444       2.681           1.09           71.2          0.459           24.7            194          12018
random1-200.refine1000          5836.812       2.673          0.149            305         0.0183             32        timeout          40114
random1-500.refine1000           timeout      11.369             ok            385             ok           95.1        skipped        timeout
random1-1000.refine1000          skipped      20.533             ok            946             ok            164        skipped        skipped
random1-2000.refine1000          skipped      97.041             ok            913             ok            246        skipped        skipped
random1024-100.refine1000        timeout      16.683             ok           29.7             ok           22.2             ok              -
random1024-200.refine1000        skipped      22.671             ok           84.1             ok           33.2        timeout              -
random1024-500.refine1000        skipped      64.187             ok            178             ok           43.6        skipped              -
random1024-1000.refine1000       skipped     128.243             ok            354             ok            103        skipped              -
random1024-2000.refine1000       skipped     756.251             ok        timeout             ok        timeout        skipped              -
random64-100.refine1000        26411.356       0.031        0.00828           9084       0.000442             ok       0.000259              -
random64-200.refine1000          timeout       3.408             ok            253             ok           26.3        timeout              -
random64-500.refine1000          skipped      20.765             ok            260             ok           69.9        skipped              -
random64-1000.refine1000         skipped      50.300             ok            412             ok            206        skipped              -
random64-2000.refine1000         skipped     105.525             ok            787             ok            363        skipped              -
repeated-4.refine1000             17.666       1.362          error          error          0.949             ok           5.18              -
repeated-8.refine1000             40.323       2.955          error          error          0.364             ok           10.8              -
repeated-16.refine1000           352.801       6.356          error          error          0.115          0.873            5.9              -
repeated-32.refine1000         21543.890      18.434          error          error        0.00554           4.89          0.883              -
repeated-64.refine1000           timeout      69.690          error          error             ok             10        timeout              -
wilkinson-20.refine1000           51.959       0.109           0.23           99.6             ok          0.397         0.0174           7.15
wilkinson-40.refine1000          134.243       0.483          0.416            105             ok           40.2         0.0345           8.09
wilkinson-80.refine1000          384.810       2.582           1.76            241         0.0338           15.2         0.0719            9.4
wilkinson-160.refine1000        1908.593      19.904           4.22            361         0.0456           5.38          0.114           9.41
wilkinson-320.refine1000       18061.155     162.520        timeout            722         0.0243           2.46          0.125             11
wilkinson-640.refine1000         timeout    1249.666        skipped        timeout             ok            2.4             ok           14.8
chebyshev-100.refine53            59.231      12.797           1.36           5.46          0.695           2.95            135            548
chebyshev-200.refine53           249.838      51.411            3.5           14.7          0.941            3.3            301           1260
chebyshev-400.refine53          1150.867     274.671           10.6           37.5           2.41           9.04        timeout        timeout
chebyshev-700.refine53          4794.802    1239.411        timeout        timeout           6.55           22.2        skipped        skipped
chebyshev-1000.refine53        12960.451    2969.372        skipped        skipped        timeout        timeout        skipped        skipped
mignotte-10.refine53               0.758       0.460           1.06           1.47             ok           3.58           1396           4017
mignotte-20.refine53               1.004       0.685           5.29            6.8             ok           1.54           9.28           23.2
mignotte-50.refine53               2.728       1.635             53           76.8           3.71           20.9           19.8             59
mignotte-100.refine53              6.496       3.338            363            650             75            142           26.9           92.8
mignotte-200.refine53             20.905      12.102           2448           3689            712           1252           56.8            171
mignotte-400.refine53             65.914      61.505        timeout        timeout        timeout        timeout            213            419
random1-100.refine53             180.738       0.659          0.201           48.7         0.0104           6.77           1.61            372
random1-200.refine53            5690.532       0.740         0.0388            286        0.00424           25.7          0.138            977
random1-500.refine53             timeout       3.254             ok            488             ok           37.4             ok          11089
random1-1000.refine53            skipped       7.014             ok           1437             ok             89             ok           6941
random1-2000.refine53            skipped      42.245             ok           1345             ok            262        timeout        timeout
random1024-100.refine53          timeout       7.914             ok           10.1             ok           31.5             ok              -
random1024-200.refine53          skipped       9.228             ok           43.8             ok           42.5             ok              -
random1024-500.refine53          skipped      31.701             ok            107             ok           35.3             ok              -
random1024-1000.refine53         skipped      74.195             ok            293             ok           63.6             ok              -
random1024-2000.refine53         skipped     353.982             ok        timeout             ok            214        timeout              -
random64-100.refine53          26238.036       0.023        0.00191           2102       0.000305            170       0.000266              -
random64-200.refine53            timeout       0.972             ok            273             ok           16.3             ok              -
random64-500.refine53            skipped       9.304             ok            221             ok           15.6             ok              -
random64-1000.refine53           skipped      22.271             ok            513             ok             76        timeout              -
random64-2000.refine53           skipped      51.540             ok            884             ok            248        skipped              -
repeated-4.refine53                2.605       0.193          error          error           1.66             ok           1.36              -
repeated-8.refine53               11.116       0.630          error          error          0.993             ok           1.25              -
repeated-16.refine53             263.868       1.301          error          error         0.0823             ok           0.18              -
repeated-32.refine53           21449.047       4.884          error          error        0.00146           1.32         0.0126              -
repeated-64.refine53             timeout      24.377          error          error             ok           3.62             ok              -
wilkinson-20.refine53              6.264       0.109          0.495           25.4             ok             ok          0.156           7.13
wilkinson-40.refine53             22.779       0.467           1.18           50.7             ok           44.5          0.204            8.2
wilkinson-80.refine53            109.586       2.924           5.05            174         0.0945           13.6          0.266           7.73
wilkinson-160.refine53           774.225      19.819             10            345          0.105           5.27          0.291           9.23
wilkinson-320.refine53          9640.193     148.532        timeout            775         0.0426           2.85          0.236           11.8
wilkinson-640.refine53           timeout    1403.987        skipped        timeout             ok           2.07             ok           12.7

```
