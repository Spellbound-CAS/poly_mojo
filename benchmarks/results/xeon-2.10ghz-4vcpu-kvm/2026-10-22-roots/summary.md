# Real roots, 2026-10-22 (nightly d6549d3)

Full roots family, every mode, against FLINT 3.6, msolve 0.11 and SymPy, pinned to one CPU, 120 s timeout.
New this run: each interval answer records its achieved widths and endpoint sizes, and Poly Mojo is also
timed at FLINT's and msolve's achieved widths (`polymojo-roots@flint-roots`, `polymojo-roots@msolve-roots`;
method in `report.md` and `benchmarks/harness/real_roots.py`).

Library: POLYMOJO_LIB=/tmp/claude-0/libW, precompiled from the d6549d3 sources (src/ is unchanged since
f231e73). Its three packages are byte-identical (sha256) to a `scripts/build.py` build of the same sources,
build id 0ac23506e95f49e6bcf04dc28e3ea030. This run predates the harness library guard (cb817f4), so
environment.json has no `poly_mojo_libraries` entry. Driver binary pm_real_roots_5883490feb7ad92b, built
after the library.

Outcomes: Poly Mojo 170 ok, 5 timeout, 10 skipped, 0 wrong (as on 10-21). The five `wrong` records are
SymPy's isolation of the repeated-root family, also wrong on 10-21.

## Ratios against 10-21 (2026-10-21-roots-descartes, 86297d6) and at matched width

Geometric means of peer time / Poly Mojo time (above 1 = Poly Mojo faster).

```
Requested width: geometric mean of peer time / Poly Mojo time over cases with both (n), 10-21 -> 10-22
(the 10-22 column uses the same case set as 10-21 where both runs have both times; "all" = every 10-22 case)
  isolate    FLINT: 58.7 -> 57.5 (n=27); all 57.5 (n=27); msolve: 10.7 -> 7.95 (n=26); all 7.74 (n=27); SymPy: 10.4 -> 9.97 (n=30); all 9.97 (n=30)
  refine53   FLINT: 123 -> 121 (n=26); all 127 (n=27); msolve: 26.3 -> 21.9 (n=27); all 20.8 (n=28); SymPy: 136 -> 134 (n=30); all 134 (n=30)
  refine1000 FLINT: 81.7 -> 80.4 (n=26); all 84.2 (n=27); msolve: 21.1 -> 18.8 (n=28); all 18.8 (n=28); SymPy: 109 -> 106 (n=14); all 219 (n=16)
  exact      SymPy: 17.8 -> 16.7 (n=7); all 16.7 (n=7)
  count      FLINT: 0.285 -> 0.274 (n=20); all 0.373 (n=21); SymPy: 2.51 -> 2.49 (n=24); all 2.49 (n=24)

Matched width (10-22): geometric mean of peer time / Poly Mojo time at that peer's achieved width, next to the requested-width ratio on the same cases
  isolate    FLINT   matched 11.7  requested 59.3  (n=26); matched range 0.766-165
  isolate    msolve  matched 4.34  requested 9.36  (n=22); matched range 0.069-547
  refine53   FLINT   matched 35.5  requested 115  (n=26); matched range 0.848-556
  refine53   msolve  matched 17.2  requested 30.8  (n=22); matched range 2.3-611
  refine1000 FLINT   matched 41.1  requested 70.4  (n=26); matched range 0.182-632
  refine1000 msolve  matched 12.5  requested 24.8  (n=24); matched range 0.0305-147

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
    sympy-roots                    cases   9  widest     -0.1  median     -0.7  narrowest-widest     -0.7  endpoint bits max 1010
    polymojo-roots@flint-roots     cases  22  widest   -544.0  median  -1033.7  narrowest-widest  -2430.1  endpoint bits max 18017
    polymojo-roots@msolve-roots    cases  27  widest  -1048.0  median  -1055.0  narrowest-widest  -2079.4  endpoint bits max 11971

Poly Mojo outcomes: {('polymojo-roots', 'ok'): 170, ('polymojo-roots', 'skipped'): 10, ('polymojo-roots', 'timeout'): 5, ('polymojo-roots@flint-roots', 'ok'): 78, ('polymojo-roots@flint-roots', 'unsupported'): 107, ('polymojo-roots@msolve-roots', 'ok'): 83, ('polymojo-roots@msolve-roots', 'unsupported'): 102}
```

Every tool ran 11-14% faster than on 10-21 (Poly Mojo 0.89x, FLINT 0.86x, SymPy 0.85x of its 10-21 times
on the same cases), so the host was faster; msolve, timed as whole processes, sped up more (0.66x isolate,
0.73x refine53, 0.78x refine1000), which accounts for most of the drop in the msolve ratios. The library
code for these modes is unchanged since 10-21 apart from step 3, which does not touch them.

## Widths

At the requested width Poly Mojo and SymPy return intervals just under 2^-bits. FLINT and msolve return
much narrower ones: at refine53 FLINT's median is about 2^-472 and msolve's about 2^-136; at refine1000
about 2^-2032 and 2^-2058; FLINT's isolating intervals have a median width of about 2^-472. At their widths
Poly Mojo is still ahead on the geometric mean (isolate 11.7x vs FLINT, 4.3x vs msolve; refine53 35.5x and
17.2x; refine1000 41.1x and 12.5x), against 59.3x, 9.4x, 115x, 30.8x, 70.4x and 24.8x at the requested width.

Cases where the peer is faster at its own width (peer time / Poly Mojo time):

| Case | Peer | Matched width | Ratio |
|---|---|---|---|
| mignotte-50.refine1000 | msolve | 2^-2084 | 0.031 |
| random1-200.isolate | msolve | 2^-136 | 0.069 |
| mignotte-50.refine1000 | FLINT | 2^-1882 | 0.18 |
| mignotte-100.refine1000 | msolve | 2^-2054 | 0.27 |
| repeated-8.refine1000 | msolve | 2^-2047 | 0.56 |
| mignotte-10.isolate | FLINT | 2^-224 | 0.77 |
| random64-500.isolate | msolve | 2^-388 | 0.80 |
| repeated-32.isolate | msolve | 2^-41 | 0.81 |
| random64-200.isolate | msolve | 2^-135 | 0.82 |
| mignotte-50.refine53 | FLINT | 2^-858 | 0.85 |
| mignotte-10.refine53 | FLINT | 2^-224 | 0.87 |
| mignotte-50.isolate | FLINT | 2^-858 | 0.93 |

High-precision refinement of Mignotte's close pair is the clear weak spot (msolve 33x faster on
mignotte-50 at 2^-2084).
