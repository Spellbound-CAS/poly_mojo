# CAD conjunctions on claude/exp-cad 65172dc (wave 7), against 2026-10-22-cad-sat (f07fb00)

Host `xeon-2.10ghz-4vcpu-kvm`. Clean worktree at 65172dc, `pixi run --locked build`.
`python -m harness.cad run --tools polymojo-solve --cases '^sat-'`, timeout 300 s, min 0.5 s / 3 reps,
CPU 1. The Z3 and QEPCAD B records are reused unchanged from 2026-10-22-cad-sat.

`cad.solve` finishes all 40 conjunctions with the same verdicts. Against f07fb00 it is 0.96x
(geomean, n=40, 0.89x-1.01x): wave 7 changed lifting arithmetic, not the solve path, so the
difference is within this host's run-to-run spread. solve/Z3 geomean 0.36x (n=37), 0.35x before.
