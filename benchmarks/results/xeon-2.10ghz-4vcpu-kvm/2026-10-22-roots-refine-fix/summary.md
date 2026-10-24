# Refine rows on the close-pair QIR fix (1ba89dc, library build 9879f9e5)

refine53 and refine1000 for every tool, at requested and matched widths, pinned to one CPU, 120 s timeout,
build/ made by scripts/build.py from 1ba89dc. Records before cf2c3c9 carry 1ba89dc; later ones carry results-only
commits with identical sources. Poly Mojo: 181 ok, 0 wrong.

Peer time / Poly Mojo time, geometric mean, against 2026-10-22-roots-rev11 (same cases):

| Mode | Poly Mojo vs | requested width | matched width |
|---|---|---|---|
| refine53 | FLINT | 136 -> 133 | 38.3 -> 54.8 |
| refine53 | msolve | 21.0 -> 30.6 | 14.7 -> 19.4 |
| refine53 | SymPy | 137 -> 139 | |
| refine1000 | FLINT | 87.5 -> 106 | 42.1 -> 67.0 |
| refine1000 | msolve | 19.9 -> 31.4 | 12.2 -> 23.3 |
| refine1000 | SymPy | 111 -> 109 | |

Close pairs below their separation: mignotte-50 refine1000 207 -> 8.8 ms (at msolve's width 1201 -> 10.7 ms;
msolve 49 ms), mignotte-100 361 -> 34 ms (at msolve's width 2007 -> 38 ms; msolve 475 ms), mignotte-200
438 -> 227 ms (at msolve's width 3878 -> 243 ms; msolve 13.7 s).

Regression: mignotte-400 refine1000 126 ms -> 1.79 s. Its pair is 2^-1335 apart, finer than the requested
2^-1000, so the whole refinement runs while the enclosure is wider than the pair's gap. Newton is then linear,
QIR succeeds about half the time on small grids, and each step now pays for the high-precision evaluation the
fix adds; plain bisection, which the old noisy guess fell back to, is cheaper there. The same effect, smaller,
shows when refining mignotte-100/200 to widths above their gaps (for example 2^-300). Peers faster at their own
width: msolve on mignotte-10 refine1000 (0.16) and repeated-4 refine53 (0.35), FLINT on mignotte-10 refine53
(0.88), all sub-millisecond cases.
