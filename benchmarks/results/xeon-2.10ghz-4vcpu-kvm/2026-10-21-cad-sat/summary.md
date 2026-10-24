# Conjunctions: cad.solve by coverings search (claude/exp-cad-sat 8c7a9a7), against Z3

Host `xeon-2.10ghz-4vcpu-kvm` (Intel Xeon @ 2.10GHz, 4 vCPU), the class of 2026-10-21-cad.
Library from `pixi run --locked build` at 8c7a9a7 (nightly 42a55a7 merged), drivers built by the harness.
`python -m harness.cad run --tools polymojo-solve,z3 --cases '^sat-'`, timeout 300 s, min 0.5 s / 3 reps,
12 GiB cap, both tools pinned to CPU 1 one at a time. Z3 5.3.1 from the `z3-solver` wheel. QEPCAD B was not
run here; its verdicts are in 2026-10-21-cad. `report.md` is the harness report, `records.jsonl` every record.

`cad.solve` now runs the conflict-driven search over cylindrical algebraic coverings (Contract P revision 10)
and returns a witness or a checkable covering tree; no global projection is built.

## Against Z3 (Poly Mojo time / Z3 time, cases both solve)

- 37 cases: median 0.38x, geomean 0.49x; Poly Mojo is faster on 24 (2026-10-21-cad: median 10.1x,
  geomean 16.7x over 36, faster on 2).
- Slowest relative to Z3: solotareff-3 23.7x (0.037 s), edges-square-product 6.0x, sphere-plane-n5-unsat 3.6x,
  simplified-putnum-example 2.5x, sphere-plane-n4-unsat 2.0x.
- sos-d6-sat, a timeout on 2026-10-21, takes 0.4 ms; sos-d4-sat 9.5 s -> 0.3 ms; solotareff-3 68 s -> 0.037 s.
- Neither tool solves sos-d6-unsat within 300 s (Poly Mojo's search takes about 2 minutes and the
  timed region also includes `cad.check`); the sos ladder therefore stops before sos-d8 for both.

## Correctness

- Every satisfiable/unsatisfiable verdict agrees with Z3 here and with Z3 and QEPCAD B in 2026-10-21-cad.
- All 24 Poly Mojo sample points passed the harness's independent exact sign check.
- `cad.check` accepted every Poly Mojo result (witnesses and covering trees).
