# CAD on nightly 841d966 (final), against 2026-10-18-cad

Host `xeon-2.10ghz-4vcpu-kvm` (Intel Xeon @ 2.10GHz, 4 vCPU), the same class as 2026-10-18-cad.
Clean worktree at 841d966, library from `pixi run build`, drivers built by the harness.
Settings as 2026-10-18-cad: `python -m harness.cad run`, timeout 300 s, min 0.5 s / 3 reps, 12 GiB cap.
Poly Mojo full CAD ran pinned to CPU 1. QEPCAD B, Z3 and `polymojo-solve` ran pinned to CPU 3, partly at the same time.
QEPCAD B 1.80 came from a Tarski source build; only `qepcad` is used. Z3 is 5.3.1, from the `z3-solver` wheel.
`report.md` is the harness report. `records.jsonl` holds every record, including errors and timeouts.

## Poly Mojo, old (2026-10-18, 9a538f5/55c07db) against new (841d966), the 21 cases of 2026-10-18

| mode | solved old → new | geomean speedup on cases both solve | total on those cases |
|---|---|---|---|
| lazard | 9 → 16 | 2.96× (1.66–11.4×, n=9) | 11.08 s → 3.28 s |
| mccallum | 7 → 13 | 2.87× (1.65–10.1×, n=7) | 9.90 s → 3.07 s |
| stream | 9 → 16 | 3.48× (1.57–20.6×, n=9) | 18.53 s → 3.66 s |

These cases were timeouts in the old run and finish now:
- arnon-84, arnon-84-2, collins-and-johnson, off-center-ellipse, positivity-of-quartic;
- simplified-yangxia (lazard and mccallum);
- sphere-and-catastrophe;
- simplified-edges-square-product (stream).

## Competitor ratios (new run only; the old run has no competitor records)

- **Full CAD against QEPCAD B.** QEPCAD's time is its wall time minus a trivial-input start-up run, which is 0 s on most small cases.
  - On the 4 cases where QEPCAD needs at least 0.05 s, Poly Mojo is slower: lazard/QEPCAD has a geomean of 53.6× and stream/QEPCAD 81.1×.
  - Solved by QEPCAD only: cyclic-4, mccallum-trivariate-random-polynomial, simplified-putnum-example and solotareff-3; also simplified-edges-square-product (lazard) and simplified-yangxia (stream).
  - Solved by Poly Mojo only: none. On piano-mover-s-problem-davenport both fail: QEPCAD errors and Poly Mojo times out.
- **Conjunctions, `cad.solve` against Z3 (nlsat).**
  - Over 36 of 40 cases: geomean 16.7×, median 10.1×; Poly Mojo is faster on 2.
  - Neither tool solves sos-d8: Z3 is skipped on both sos-d8 cases, and Poly Mojo times out on sos-d6-sat and is skipped on the rest of the sos ladder.
  - Z3 times out on sos-d6-unsat.
- **Conjunctions against QEPCAD B.** On the 5 cases where QEPCAD needs at least 0.05 s, the geomean is 11.5×.
- **Correctness.** Every satisfiable/unsatisfiable verdict agrees across all three tools. No record is `wrong`, and every Poly Mojo sample point passed the independent sign check.

## Errors

- McCallum mode refuses 8 cases with "McCallum nullification requires Lazard mode at this cell". This is a documented refusal, not a failure.
- Lazard mode runs out of memory (12 GiB cap) on edges-square-product, which is new in the case bank.
