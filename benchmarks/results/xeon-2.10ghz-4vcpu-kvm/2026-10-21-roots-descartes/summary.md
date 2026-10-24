# Real roots, 2026-10-21, step 3 (claude/exp-roots 86297d6: exact real algebraic values on Descartes)

Clean worktree at 86297d6 (code identical to 000b7b7), same harness settings as 2026-10-21-roots: pinned CPU 1,
120 s case timeout, 12G cap, ladders stop after a timeout, CONDA_PREFIX set, every tool in one pass.

Ratios are peer time / Poly Mojo time inside one run, so host speed cancels (>1: Poly Mojo faster).
Old = 2026-10-21-roots (nightly 841d966, step 2). No Poly Mojo answer is wrong; outcomes: 170 ok, 5 timeout, 10 skipped.

## Geometric-mean competitor ratios, old -> new

| mode | FLINT/ours | msolve/ours | SymPy/ours |
|---|---|---|---|
| exact | unsupported | unsupported | 0.21 -> 15.2 (5) |
| isolate | 57.5 -> 58.7 (27) | 10.3 -> 10.8 (27) | 18.9 -> 14.4 (20) |
| refine53 | 126 -> 123 (26) | 24.8 -> 23.9 (28) | 140 -> 114 (18) |
| refine1000 | 83.3 -> 81.7 (26) | 20.6 -> 20.3 (26) | 45 -> 42.8 (8) |
| count (Sturm, unchanged) | 0.276 -> 0.285 (20) | unsupported | 4.02 -> 2.98 (19) |

Exact mode now solves every row the old run timed out on except chebyshev-700/1000, e.g. mignotte-10 timeout -> 0.26 ms,
random1-100 timeout -> 8.3 ms, chebyshev-100 timeout -> 91 ms, repeated-4 2268 -> 0.76 ms, wilkinson-640 timeout -> 4.1 s.
Isolate and refine are unchanged (step 3 does not touch them); their SymPy ratios move with SymPy's own run-to-run time.

## Per case

```
case                         ours old ms ours new ms FLINT/ours old FLINT/ours new msolve/ours old msolve/ours new SymPy/ours old SymPy/ours new
chebyshev-100.count                1.814       2.253             16             15    unsupported    unsupported           2.79           2.63
chebyshev-200.count                8.402      10.487            110            108    unsupported    unsupported           2.39           2.21
chebyshev-400.count               38.958      47.534            661            689    unsupported    unsupported           2.16           2.24
chebyshev-700.count              143.459     175.782        timeout        timeout    unsupported    unsupported           2.19           2.17
chebyshev-1000.count             335.606     399.322        skipped        skipped    unsupported    unsupported           2.33           2.53
mignotte-10.count                  0.022       0.026         0.0724         0.0749    unsupported    unsupported           8.48            9.2
mignotte-20.count                  0.040       0.047         0.0528         0.0568    unsupported    unsupported           7.57           8.28
mignotte-50.count                  0.142       0.171         0.0503         0.0537    unsupported    unsupported           5.79           5.35
mignotte-100.count                 0.526       0.633         0.0406         0.0446    unsupported    unsupported           3.76           3.46
mignotte-200.count                 1.705       1.975         0.0418         0.0486    unsupported    unsupported           8.34           3.58
mignotte-400.count                 6.848       7.992         0.0386         0.0446    unsupported    unsupported           6.85           3.33
random1-100.count                148.291     196.842        0.00957        0.00903    unsupported    unsupported            8.7           8.43
random1-200.count               2735.861    3424.006        0.00742        0.00758    unsupported    unsupported           20.9           20.6
random1-500.count                timeout     timeout             ok             ok    unsupported    unsupported        timeout        timeout
random1-1000.count               skipped     skipped             ok             ok    unsupported    unsupported        skipped        skipped
random1-2000.count               skipped     skipped        timeout        timeout    unsupported    unsupported        skipped        skipped
random1024-100.count             timeout     timeout             ok             ok    unsupported    unsupported        timeout        timeout
random1024-200.count             skipped     skipped             ok             ok    unsupported    unsupported        skipped        skipped
random1024-500.count             skipped     skipped        timeout        timeout    unsupported    unsupported              -        skipped
random1024-1000.count            skipped     skipped        skipped        skipped    unsupported    unsupported              -        skipped
random1024-2000.count            skipped     skipped        skipped        skipped    unsupported    unsupported              -        skipped
random64-100.count             23172.123   13387.311        0.00862        0.00817    unsupported    unsupported              -           8.34
random64-200.count               timeout     timeout             ok             ok    unsupported    unsupported              -        timeout
random64-500.count               skipped     skipped             ok             ok    unsupported    unsupported              -        skipped
random64-1000.count              skipped     skipped        timeout        timeout    unsupported    unsupported              -        skipped
random64-2000.count              skipped     skipped        skipped        skipped    unsupported    unsupported              -        skipped
repeated-4.count                   0.975       1.138         0.0228         0.0226    unsupported    unsupported              -          0.525
repeated-8.count                   9.268      11.771         0.0373         0.0349    unsupported    unsupported              -          0.315
repeated-16.count                 93.161     105.748          0.146          0.146    unsupported    unsupported              -           1.22
repeated-32.count               4995.648    5706.988          0.105          0.107    unsupported    unsupported              -           2.23
repeated-64.count                timeout     timeout             ok             ok    unsupported    unsupported              -        timeout
wilkinson-20.count                 0.170       0.216          0.184          0.171    unsupported    unsupported           4.46           4.27
wilkinson-40.count                 1.594       1.931           1.03           1.48    unsupported    unsupported           4.03           1.42
wilkinson-80.count                10.610      13.306           7.39           6.79    unsupported    unsupported            2.1           0.88
wilkinson-160.count               79.160      93.118           41.2           44.1    unsupported    unsupported            1.8          0.856
wilkinson-320.count              621.166     739.415            170        timeout    unsupported    unsupported            1.8          0.891
wilkinson-640.count             4635.357    5308.474        timeout        skipped    unsupported    unsupported           2.38           1.24
chebyshev-100.exact              timeout      91.412    unsupported    unsupported    unsupported    unsupported          error          error
chebyshev-200.exact              skipped     790.435    unsupported    unsupported    unsupported    unsupported          error          error
chebyshev-400.exact              skipped    7666.519    unsupported    unsupported    unsupported    unsupported          error        timeout
chebyshev-700.exact              skipped     timeout    unsupported    unsupported    unsupported    unsupported        timeout        skipped
chebyshev-1000.exact             skipped     skipped    unsupported    unsupported    unsupported    unsupported        skipped        skipped
mignotte-10.exact                timeout       0.258    unsupported    unsupported    unsupported    unsupported          error          error
mignotte-20.exact                skipped       0.420    unsupported    unsupported    unsupported    unsupported        timeout        timeout
mignotte-50.exact                skipped       2.966    unsupported    unsupported    unsupported    unsupported        skipped        skipped
mignotte-100.exact               skipped      36.775    unsupported    unsupported    unsupported    unsupported        skipped        skipped
mignotte-200.exact               skipped     852.386    unsupported    unsupported    unsupported    unsupported        skipped        skipped
mignotte-400.exact               skipped   21617.900    unsupported    unsupported    unsupported    unsupported        skipped        skipped
random1-100.exact                timeout       8.340    unsupported    unsupported    unsupported    unsupported          error          error
random1-200.exact                skipped      32.181    unsupported    unsupported    unsupported    unsupported          error          error
random1-500.exact                skipped     274.131    unsupported    unsupported    unsupported    unsupported          error          error
random1-1000.exact               skipped     942.773    unsupported    unsupported    unsupported    unsupported        timeout        timeout
random1-2000.exact               skipped   12106.705    unsupported    unsupported    unsupported    unsupported        skipped        skipped
random1024-100.exact             timeout      15.652    unsupported    unsupported    unsupported    unsupported              -          error
random1024-200.exact             skipped      97.115    unsupported    unsupported    unsupported    unsupported              -          error
random1024-500.exact             skipped     462.894    unsupported    unsupported    unsupported    unsupported              -        timeout
random1024-1000.exact            skipped    6095.067    unsupported    unsupported    unsupported    unsupported              -        skipped
random1024-2000.exact            skipped   30282.520    unsupported    unsupported    unsupported    unsupported              -        skipped
random64-100.exact             11288.361      12.614    unsupported    unsupported    unsupported    unsupported              -           31.8
random64-200.exact               timeout      52.661    unsupported    unsupported    unsupported    unsupported              -          error
random64-500.exact               skipped     608.812    unsupported    unsupported    unsupported    unsupported              -          error
random64-1000.exact              skipped    5706.838    unsupported    unsupported    unsupported    unsupported              -        timeout
random64-2000.exact              skipped   10257.304    unsupported    unsupported    unsupported    unsupported              -        skipped
repeated-4.exact                2268.326       0.763    unsupported    unsupported    unsupported    unsupported              -          error
repeated-8.exact                 timeout       1.984    unsupported    unsupported    unsupported    unsupported              -          error
repeated-16.exact                skipped      13.466    unsupported    unsupported    unsupported    unsupported              -          error
repeated-32.exact                skipped      75.509    unsupported    unsupported    unsupported    unsupported              -          error
repeated-64.exact                skipped     596.927    unsupported    unsupported    unsupported    unsupported              -          error
wilkinson-20.exact               102.590       1.401    unsupported    unsupported    unsupported    unsupported          0.529           20.4
wilkinson-40.exact               640.151       5.848    unsupported    unsupported    unsupported    unsupported          0.176           14.6
wilkinson-80.exact              3265.331      24.647    unsupported    unsupported    unsupported    unsupported          0.173           12.2
wilkinson-160.exact            16130.232     118.455    unsupported    unsupported    unsupported    unsupported          0.156           13.3
wilkinson-320.exact           101168.972     630.310    unsupported    unsupported    unsupported    unsupported          0.162           16.7
wilkinson-640.exact              timeout    4121.322    unsupported    unsupported    unsupported    unsupported        timeout           21.9
chebyshev-100.isolate              2.460       2.890             28             29           17.7           14.3            139            144
chebyshev-200.isolate             19.554      24.883           39.6           39.7           8.26           9.43            144            140
chebyshev-400.isolate            158.667     199.531           62.5           64.7           15.3           14.5            192            192
chebyshev-700.isolate           1729.892    1962.362        timeout        timeout             16           16.9        timeout        timeout
chebyshev-1000.isolate          4925.075    5627.789        skipped        skipped        timeout        timeout        skipped        skipped
mignotte-10.isolate                0.066       0.063           10.8             13           23.3             ok           7.24           8.89
mignotte-20.isolate                0.057       0.064             86           86.1             ok             ok           14.2           14.7
mignotte-50.isolate                0.193       0.236            683            622            186           57.8           36.3           15.4
mignotte-100.isolate               0.659       0.756           3257           3543            848            694           33.1           16.2
mignotte-200.isolate               3.335       3.347          14022          16324           4765           4602           23.4           13.4
mignotte-400.isolate              17.536      19.087        timeout        timeout        timeout        timeout           19.1           9.73
random1-100.isolate                0.544       0.645           45.9           44.5             ok           0.97           13.9           13.3
random1-200.isolate                2.550       2.987           64.3           64.7           7.84            6.8             11           12.2
random1-500.isolate               23.269      28.094           66.8           67.3           5.18           5.17           9.46           9.48
random1-1000.isolate             198.717     217.957           36.8           39.8           2.97           3.25           2.69           2.89
random1-2000.isolate            1673.968    2012.424           32.8           33.9           6.23            6.3           3.38           3.16
random1024-100.isolate             2.681       2.230           23.2             15            123           66.1              -             11
random1024-200.isolate            15.400       9.246           21.4           21.6           26.6           22.2              -           7.85
random1024-500.isolate           106.180      61.737           31.4           31.8           10.6           9.81              -           6.07
random1024-1000.isolate          758.167     446.078           16.2           16.9              6           5.96              -            6.2
random1024-2000.isolate        11365.150    5893.650           9.79           11.4           6.52           7.41              -           6.88
random64-100.isolate               1.569       1.787           25.1             26           3.05           13.5              -           3.87
random64-200.isolate               6.551       7.617           27.8           29.5           1.38           4.64              -           5.48
random64-500.isolate             173.003     185.940             11           11.4          0.943           1.03              -           2.41
random64-1000.isolate            296.531     332.186           28.6           29.3           5.27           5.16              -           7.17
random64-2000.isolate           1907.003    1982.916           24.9           26.3            6.8           7.38              -           2.51
repeated-4.isolate                 0.090       0.084          error          error             ok             ok              -          wrong
repeated-8.isolate                 0.248       0.292          error          error             ok           31.2              -          wrong
repeated-16.isolate                1.332       1.474          error          error             ok           14.8              -          wrong
repeated-32.isolate                6.122       6.983          error          error          0.991           4.43              -          wrong
repeated-64.isolate               35.661      47.052          error          error            2.4           2.28              -          wrong
wilkinson-20.isolate               0.120       0.138           19.2           20.2           33.9             ok           6.59           6.84
wilkinson-40.isolate               0.490       0.578             43           44.2           40.4             ok             16           8.29
wilkinson-80.isolate               2.585       3.177            180            178           13.3           6.63           18.8           8.91
wilkinson-160.isolate             19.731      23.236            339            339           5.36           4.31           18.6           9.82
wilkinson-320.isolate            189.110     206.212            599        timeout           2.22           2.46           19.1           11.6
wilkinson-640.isolate           1842.745    2049.151        timeout        skipped           1.64           1.52             17           11.5
chebyshev-100.refine1000          58.929      79.006           1.88           1.65           16.9           15.5        timeout        timeout
chebyshev-200.refine1000         234.930     273.985            3.9           4.35           32.2           35.4        skipped        skipped
chebyshev-400.refine1000         937.408    1259.849           10.8           10.9           65.7             62        skipped        skipped
chebyshev-700.refine1000        3685.705    4871.330        timeout        timeout        timeout        timeout        skipped        skipped
chebyshev-1000.refine1000       8950.497   10371.296        skipped        skipped        skipped        skipped        skipped        skipped
mignotte-10.refine1000             1.316       1.639           2.18              2           4.24             ok           1691            839
mignotte-20.refine1000             2.199       2.758            4.8           4.57          0.207             ok        timeout        timeout
mignotte-50.refine1000           183.618     223.792           0.77          0.741           0.25          0.082        skipped        skipped
mignotte-100.refine1000          323.312     381.590           6.79           6.51           1.53           1.46        skipped        skipped
mignotte-200.refine1000          449.802     489.407            101            107           34.8           30.8        skipped        skipped
mignotte-400.refine1000          116.011     133.539        timeout        timeout        timeout        timeout        skipped        skipped
random1-100.refine1000             2.681       3.184           71.2           73.1           24.7             20          12018          12354
random1-200.refine1000             2.673       3.255            305            303             32           34.9          40114        timeout
random1-500.refine1000            11.369      13.733            385            386           95.1           90.6        timeout        skipped
random1-1000.refine1000           20.533      24.752            946            887            164            161        skipped        skipped
random1-2000.refine1000           97.041     114.098            913            927            246            249        skipped        skipped
random1024-100.refine1000         16.683      10.426           29.7             27           22.2           21.7              -           6534
random1024-200.refine1000         22.671      11.820           84.1           89.3           33.2             34              -        timeout
random1024-500.refine1000         64.187      33.687            178            203           43.6             48              -        skipped
random1024-1000.refine1000       128.243      67.355            354            401            103            119              -        skipped
random1024-2000.refine1000       756.251     448.097        timeout        timeout        timeout            247              -        skipped
random64-100.refine1000            0.031       0.030           9084           8334             ok            661              -            231
random64-200.refine1000            3.408       3.640            253            266           26.3           34.7              -        timeout
random64-500.refine1000           20.765      25.943            260            220           69.9           71.8              -        skipped
random64-1000.refine1000          50.300      56.365            412            399            206            215              -        skipped
random64-2000.refine1000         105.525     112.615            787            822            363            384              -        skipped
repeated-4.refine1000              1.362       1.391          error          error             ok           10.8              -           66.2
repeated-8.refine1000              2.955       3.493          error          error             ok           5.12              -            113
repeated-16.refine1000             6.356       7.257          error          error          0.873           4.42              -            287
repeated-32.refine1000            18.434      21.514          error          error           4.89           7.97              -            961
repeated-64.refine1000            69.690     106.259          error             10             10           7.72              -        timeout
wilkinson-20.refine1000            0.109       0.135           99.6             89          0.397             ok           7.15           7.47
wilkinson-40.refine1000            0.483       0.589            105            101           40.2             ok           8.09           8.03
wilkinson-80.refine1000            2.582       3.533            241            198           15.2           6.59            9.4            8.2
wilkinson-160.refine1000          19.904      23.491            361            358           5.38           3.51           9.41             10
wilkinson-320.refine1000         162.520     175.271            722        timeout           2.46           2.44             11           12.9
wilkinson-640.refine1000        1249.666    1375.112        timeout        skipped            2.4           2.32           14.8             17
chebyshev-100.refine53            12.797      15.820           5.46           5.07           2.95           3.04            548            570
chebyshev-200.refine53            51.411      69.094           14.7           13.6            3.3           3.24           1260           1177
chebyshev-400.refine53           274.671     324.998           37.5           38.6           9.04           8.87        timeout        timeout
chebyshev-700.refine53          1239.411    1514.498        timeout        timeout           22.2           21.6        skipped        skipped
chebyshev-1000.refine53         2969.372    3636.189        skipped        skipped        timeout        timeout        skipped        skipped
mignotte-10.refine53               0.460       0.692           1.47           1.31           3.58             ok           4017           1809
mignotte-20.refine53               0.685       0.837            6.8           6.68           1.54             ok           23.2           12.2
mignotte-50.refine53               1.635       1.990           76.8           70.6           20.9           5.37             59           31.7
mignotte-100.refine53              3.338       3.956            650            671            142            117           92.8           47.7
mignotte-200.refine53             12.102      13.401           3689           3998           1252           1121            171           92.4
mignotte-400.refine53             61.505      69.719        timeout        timeout        timeout        timeout            419            216
random1-100.refine53               0.659       0.842           48.7           44.7           6.77           1.82            372            377
random1-200.refine53               0.740       0.940            286            258           25.7           31.7            977            953
random1-500.refine53               3.254       4.172            488            460           37.4           36.4          11089          10787
random1-1000.refine53              7.014       8.570           1437           1435             89           82.1           6941           6566
random1-2000.refine53             42.245      46.547           1345           1493            262            280        timeout        timeout
random1024-100.refine53            7.914       4.500           10.1           9.75           31.5           32.8              -            144
random1024-200.refine53            9.228       5.927           43.8           41.4           42.5           34.6              -            365
random1024-500.refine53           31.701      16.703            107            124           35.3           37.7              -           1074
random1024-1000.refine53          74.195      36.086            293            354           63.6           69.3              -           2691
random1024-2000.refine53         353.982     204.674        timeout        timeout            214            210              -        timeout
random64-100.refine53              0.023       0.030           2102           1862            170            727              -            230
random64-200.refine53              0.972       1.077            273            269           16.3           30.4              -            777
random64-500.refine53              9.304      15.183            221            153           15.6           12.1              -           1204
random64-1000.refine53            22.271      25.637            513            515             76             66              -        timeout
random64-2000.refine53            51.540      51.503            884            997            248            289              -        skipped
repeated-4.refine53                0.193       0.248          error          error             ok             ok              -           13.6
repeated-8.refine53                0.630       0.760          error          error             ok           19.3              -           18.8
repeated-16.refine53               1.301       1.631          error          error             ok           12.8              -           30.9
repeated-32.refine53               4.884       5.692          error          error           1.32           6.59              -           46.1
repeated-64.refine53              24.377      29.019          error          error           3.62           3.69              -           74.8
wilkinson-20.refine53              0.109       0.131           25.4           24.3             ok             ok           7.13           7.43
wilkinson-40.refine53              0.467       0.581           50.7           46.3           44.5             ok            8.2           7.91
wilkinson-80.refine53              2.924       3.685            174            162           13.6           3.12           7.73            7.7
wilkinson-160.refine53            19.819      23.173            345            344           5.27           4.63           9.23           10.4
wilkinson-320.refine53           148.532     169.542            775        timeout           2.85           2.58           11.8           12.8
wilkinson-640.refine53          1403.987    1376.173        timeout        skipped           2.07           2.18           12.7           16.6

```
