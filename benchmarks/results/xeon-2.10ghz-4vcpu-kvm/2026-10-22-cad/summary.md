# CAD full decomposition on nightly f07fb00, against 2026-10-21-cad

Host `xeon-2.10ghz-4vcpu-kvm` (Intel Xeon @ 2.10GHz, 4 vCPU), the class of 2026-10-21-cad.
Clean worktree at f07fb00, library from `pixi run --locked build`, drivers built by the harness.
`python -m harness.cad run`, timeout 300 s, min 0.5 s / 3 reps, 12 GiB cap. Every tool ran pinned to CPU 1,
one at a time, nothing else running; piano-mover-s-problem-davenport ran last, on its own.
QEPCAD B 1.80 (Tarski source build, only `qepcad`). `report.md` is the harness report, `records.jsonl` every record.

## Poly Mojo, 2026-10-21 (841d966) against 2026-10-22 (f07fb00), the 26 bank cases

| mode | solved 10-21 -> 10-22 | geomean speedup on cases both solve |
|---|---|---|
| lazard | 19 -> 21 | 3.74x (1.69-17.4x, n=19) |
| mccallum | 13 -> 15 | 4.27x (1.98-18.9x, n=13) |
| stream | 19 -> 21 | 3.12x (1.79-20.8x, n=19) |

Newly solved: simplified-edges-square-product and simplified-putnum-example (all three modes where applicable).

## Against QEPCAD B (ratios within each run; QEPCAD time is wall minus a trivial-input start-up run)

On cases where QEPCAD needs at least 0.05 s: lazard/QEPCAD geomean 53.6x (n=4) on 10-21, 10.7x (n=6) now;
stream/QEPCAD 81.1x -> 17.4x. Per case now (lazard/QEPCAD): simplified-yangxia 1.5, simplified-edges-square-product 3.0,
collins-and-johnson 7.2, x-axis-ellipse-problem 10.3, off-center-ellipse 11.8, simplified-putnum-example 392.

Solved by QEPCAD only: cyclic-4, mccallum-trivariate-random-polynomial, solotareff-3. Solved by Poly Mojo only: none.

## Still failing

- Timeouts (300 s) in every mode: cyclic-4 (McCallum refuses it), mccallum-trivariate-random-polynomial, solotareff-3,
  piano-mover-s-problem-davenport (QEPCAD errors on it too).
- edges-square-product: Lazard out of memory at the 12 GiB cap, stream timeout, McCallum refuses.
- McCallum refuses 8 cases with "McCallum nullification requires Lazard mode at this cell" (documented, not a failure).

No record is `wrong`; every Poly Mojo sample point passed the independent sign check and `cad.check`.
