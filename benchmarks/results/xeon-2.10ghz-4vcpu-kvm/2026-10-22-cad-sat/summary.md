# CAD conjunctions on nightly f07fb00, against 2026-10-21-cad-sat

Host `xeon-2.10ghz-4vcpu-kvm`, the class of 2026-10-21-cad-sat. Clean worktree at f07fb00, `pixi run --locked build`.
`python -m harness.cad run --tools polymojo-solve,z3,qepcad-sat --cases '^sat-'`, timeout 300 s, min 0.5 s / 3 reps,
12 GiB cap, each tool pinned to CPU 1 one at a time. Z3 5.3.1 (`z3-solver` wheel), QEPCAD B 1.80.

## `cad.solve` / Z3 per ladder (geomean over cases both solve; below 1 means Poly Mojo is faster)

| ladder | 2026-10-21-cad-sat | 2026-10-22 | solve solved 10-21 -> 10-22 |
|---|---|---|---|
| bank | 0.62x (n=18) | 0.44x (n=18) | 18 -> 18 |
| hong | 0.84x (n=6) | 0.55x (n=6) | 6 -> 6 |
| motzkin | 0.29x (n=2) | 0.25x (n=2) | 2 -> 2 |
| sos | 0.06x (n=3) | 0.04x (n=3) | 3 -> 6 |
| sphere-plane | 0.50x (n=8) | 0.39x (n=8) | 8 -> 8 |
| all | 0.49x (n=37, faster on 24) | 0.35x (n=37, faster on 30) | 37 -> 40 |

`cad.solve` now finishes all 40 conjunctions; the sos ladder (sos-d6/d8) is where the others fail.
Against QEPCAD B where it needs at least 0.05 s (n=6): solve/QEPCAD geomean 0.05x.

Verdicts agree across all three tools on every case. Not ok: Z3 times out on sos-d6-unsat and skips both sos-d8
cases; QEPCAD errors on sos-d6-sat, times out on sos-d8-sat, skips sos-d8-unsat.
