# CAD full decomposition on claude/exp-cad 65172dc (wave 7), against 2026-10-22-cad (f07fb00)

Host `xeon-2.10ghz-4vcpu-kvm`, the class of 2026-10-22-cad. Clean worktree at 65172dc, library from
`pixi run --locked build`, drivers built by the harness. `python -m harness.cad run`, timeout 300 s,
min 0.5 s / 3 reps, 12 GiB cap, each Poly Mojo tool pinned to CPU 1, one at a time.
The QEPCAD B records are reused unchanged from 2026-10-22-cad (same host class, same QEPCAD build).

mccallum-trivariate-random-polynomial ran on its own with `--cell-limit 4`: the harness's independent
sympy check of sample points (exact zeros through `minimal_polynomial`) did not finish within two hours
on its high-degree samples. `cad.check` still verified the whole certificate (check_ok true).

## Poly Mojo, f07fb00 against 65172dc, the 26 bank cases

| mode | solved | geomean speedup on cases both solve |
|---|---|---|
| lazard | 21 -> 22 | 1.05x (n=21) |
| mccallum | 15 -> 15 | 1.11x (n=15) |
| stream | 21 -> 22 | 1.03x (n=21) |

Newly solved: mccallum-trivariate-random-polynomial in lazard (97.6 s, check 29 s) and stream (126 s);
QEPCAD takes 7.1 s. McCallum mode still times out on it.

Lifting cases (lazard, f07fb00 -> 65172dc, QEPCAD): off-center-ellipse 7.20 -> 5.20 s (0.61);
collins-and-johnson 1.34 -> 1.24 s (0.19); simplified-edges-square-product 11.10 -> 11.18 s (3.73);
x-axis-ellipse-problem 8.90 -> 8.86 s (0.87); simplified-putnum-example 55.7 -> 50.4 s (0.14).
McCallum mode: collins 1.55 -> 1.02 s, off-center 8.80 -> 6.03 s, putnum 75.8 -> 61.8 s.

Every record has check_ok true; cells and projection factor counts are unchanged.

## Still failing

Timeouts in every mode: cyclic-4, solotareff-3, piano-mover-s-problem-davenport.
edges-square-product: out of memory (lazard), timeout (stream). McCallum refuses the cases it refused before.
