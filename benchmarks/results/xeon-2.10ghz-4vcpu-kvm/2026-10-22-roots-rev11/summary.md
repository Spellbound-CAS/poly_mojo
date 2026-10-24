# Real roots, revision 11 (f8b820d, library build 85bd9c9f)

Full roots family against FLINT 3.6, msolve 0.11 and SymPy, every mode, requested and matched widths, pinned
to one CPU, 120 s timeout, with build/ made by scripts/build.py from f8b820d (recorded in environment.json).
Record `version` labels vary because results-only commits moved HEAD during the run; the sources are
identical to f8b820d throughout.

Revision 11 counts roots by Descartes' rule: the count mode times `dense.count_roots`, the call a user makes,
which now returns Descartes evidence. `polymojo-roots-sturm` times the explicit Sturm route
(`sturm_sequence` + `roots.count_roots`), the previous count path.

Outcomes: Poly Mojo 369 ok, 5 timeout, 9 skipped, 0 wrong across its tools; every count case is solved.

## Counting

FLINT time / Poly Mojo time, geometric mean: 0.373 (21 cases solved on 2026-10-22-roots) -> 2.42 on the same
21 cases and 8.46 on all 28 cases both solve now. SymPy: 2.49 -> 10.4 (24 cases).

| Case | Before (Sturm) | Descartes count | FLINT |
|---|---|---|---|
| random1-500 | timeout | 26 ms | 954 ms |
| random64-200 | timeout | 7.4 ms | 1.28 s |
| random1024-100 | timeout | 2.0 ms | 6.14 s |
| repeated-64 | timeout | 39 ms | 19.9 s |
| random64-100 | 10.4 s | 1.7 ms | 92 ms |
| random1-200 | 3.02 s | 2.6 ms | 30 ms |
| repeated-32 | 4.83 s | 6.9 ms | 688 ms |
| wilkinson-320 | 0.66 s | 196 ms | 118 s |

Still below parity: mignotte (FLINT 2-330 us, 0.02-0.04x), wilkinson-20 (0.22x), repeated-4 (0.31x).
Slower than the Sturm route itself: high-degree Chebyshev counts (chebyshev-400 175 vs 45 ms, -700 1.9 vs
0.17 s, -1000 5.6 vs 0.37 s), where every root is real and Descartes isolates them all; FLINT times out
from chebyshev-700 and takes 30.5 s on chebyshev-400.

## All modes against 2026-10-22-roots

```
Requested width: geometric mean of peer time / Poly Mojo time over cases with both (n), 10-21 -> 10-22
(the 10-22 column uses the same case set as 10-21 where both runs have both times; "all" = every 10-22 case)
  isolate    FLINT: 57.5 -> 55.9 (n=27); all 55.9 (n=27); msolve: 7.95 -> 8.89 (n=26); all 8.9 (n=30); SymPy: 9.97 -> 9.83 (n=30); all 9.83 (n=30)
  refine53   FLINT: 127 -> 136 (n=27); all 136 (n=27); msolve: 21.9 -> 23.1 (n=27); all 21 (n=29); SymPy: 134 -> 137 (n=30); all 137 (n=30)
  refine1000 FLINT: 84.2 -> 87.5 (n=27); all 87.5 (n=27); msolve: 18.8 -> 19.4 (n=28); all 19.9 (n=30); SymPy: 106 -> 111 (n=14); all 111 (n=14)
  exact      SymPy: 16.7 -> 17.4 (n=7); all 17.4 (n=7)
  count      FLINT: 0.373 -> 2.42 (n=21); all 8.46 (n=28); SymPy: 2.49 -> 10.4 (n=24); all 10.4 (n=24)

Matched width (10-22): geometric mean of peer time / Poly Mojo time at that peer's achieved width, next to the requested-width ratio on the same cases
  isolate    FLINT   matched 11.6  requested 57.8  (n=26); matched range 0.48-180
  isolate    msolve  matched 5.13  requested 10.8  (n=25); matched range 0.652-559
  refine53   FLINT   matched 38.3  requested 124  (n=26); matched range 0.819-574
  refine53   msolve  matched 14.7  requested 26.3  (n=24); matched range 1.55-649
  refine1000 FLINT   matched 42.1  requested 73.7  (n=26); matched range 0.192-673
  refine1000 msolve  matched 12.2  requested 24.6  (n=25); matched range 0.0313-170

Achieved widths (10-22): log2 of width over the requested 2^-bits (isolate: log2 width); widest over cases / median of case medians; longest endpoint bits
  isolate
    polymojo-roots                 cases  36  widest     -1.0  median     -2.7  narrowest-widest    -11.0  endpoint bits max 4523
    flint-roots                    cases  26  widest     -8.9  median   -471.7  narrowest-widest  -3429.7  endpoint bits max 8236
    msolve-roots                   cases  28  widest     -8.0  median   -135.0  narrowest-widest  -1542.0  endpoint bits max 3350
    sympy-roots                    cases  23  widest      0.0  median      0.0  narrowest-widest     -6.1  endpoint bits max 28
    polymojo-roots@flint-roots     cases  22  widest    -56.4  median   -499.8  narrowest-widest  -3430.1  endpoint bits max 18017
    polymojo-roots@msolve-roots    cases  28  widest     -9.0  median   -134.2  narrowest-widest  -1543.2  endpoint bits max 3347
  refine53
    polymojo-roots                 cases  30  widest     -0.1  median     -1.5  narrowest-widest     -1.9  endpoint bits max 4661
    flint-roots                    cases  26  widest    -20.0  median   -419.1  narrowest-widest  -3376.7  endpoint bits max 8236
    msolve-roots                   cases  28  widest    -53.0  median    -83.0  narrowest-widest  -1489.0  endpoint bits max 3350
    sympy-roots                    cases  23  widest     -0.0  median     -1.3  narrowest-widest     -3.1  endpoint bits max 2650
    polymojo-roots@flint-roots     cases  22  widest    -67.4  median   -448.7  narrowest-widest  -3377.1  endpoint bits max 18017
    polymojo-roots@msolve-roots    cases  28  widest    -54.0  median    -82.8  narrowest-widest  -1490.4  endpoint bits max 4125
  refine1000
    polymojo-roots                 cases  30  widest     -0.4  median     -1.4  narrowest-widest     -1.9  endpoint bits max 6849
    flint-roots                    cases  26  widest   -542.4  median  -1032.3  narrowest-widest  -2429.7  endpoint bits max 8236
    msolve-roots                   cases  27  widest  -1047.0  median  -1058.0  narrowest-widest  -2078.0  endpoint bits max 6422
    sympy-roots                    cases   7  widest     -0.1  median     -0.8  narrowest-widest     -0.6  endpoint bits max 1010
    polymojo-roots@flint-roots     cases  22  widest   -544.0  median  -1033.7  narrowest-widest  -2430.1  endpoint bits max 18017
    polymojo-roots@msolve-roots    cases  27  widest  -1048.0  median  -1055.0  narrowest-widest  -2079.4  endpoint bits max 11971

Poly Mojo outcomes: {('polymojo-roots', 'ok'): 184, ('polymojo-roots', 'timeout'): 1, ('polymojo-roots-sturm', 'ok'): 24, ('polymojo-roots-sturm', 'skipped'): 9, ('polymojo-roots-sturm', 'timeout'): 4, ('polymojo-roots-sturm', 'unsupported'): 148, ('polymojo-roots@flint-roots', 'ok'): 78, ('polymojo-roots@flint-roots', 'unsupported'): 107, ('polymojo-roots@msolve-roots', 'ok'): 83, ('polymojo-roots@msolve-roots', 'unsupported'): 102}
```
