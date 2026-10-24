"""Native/wide interval producers are accepted only by the exact final checker."""
from apn_mojo import Rational
from polynomials import Budget, QQ, NoContext, DensePoly, MonomialOrder, real_root, PolyErrorKind
from polynomials._algebraic import _Frame
from polynomials._cad_sample_field import _sample_field, _specialize_field
from polynomials._cad_field_stack import _construct_field_stack, _check_field_stack, _ksectors
from polynomials._cad_interval_roots import _interval_field_roots, _IntervalRootStats
from support.corpus import parse_polynomial
from support.check import Checker


def main() raises:
    var work = Budget()
    var c = Checker("interval field root producer")
    var frame = _Frame()
    var native_nodes = 0
    var wide_nodes = 0
    for algebraic in range(2):
        if algebraic:
            var m = DensePoly[QQ].from_coefficients(0, NoContext(), [QQ(Rational(-2)), QQ.zero(NoContext()), QQ.one(NoContext())], work)
            frame.generator = real_root(m, QQ(Rational(1)), QQ(Rational(2)), work)
            frame.variables.append(0)
            frame.representations.append(DensePoly[QQ].from_coefficients(0, NoContext(), [QQ.zero(NoContext()), QQ.one(NoContext())], work))
        var field = _sample_field(frame, work)
        for expr in ["3*x1 - 1", "x1", "-2*x1 + 3", "7", "x1^3 - x1", "3*x1^2 - x1", "x1^2 + 1", "x1^4 - 5*x1^2 + 4", "x1^2 - 100000000000000000000", "x1^2 - 1/100000000000000000000"]:
            var f = parse_polynomial(expr, MonomialOrder(), work)
            var inputs: List[List[DensePoly[QQ]]] = [_specialize_field(f, 1, field, work)]
            var exact = _construct_field_stack(inputs, field, work)
            var stats = _IntervalRootStats()
            var roots = _interval_field_roots(exact.basis, field, stats, work)
            c.check(Bool(roots), "bounded interval attempt completes simple real/nonreal roots")
            if roots:
                var candidate = exact.copy()
                candidate.roots = roots.take()
                candidate.sectors = _ksectors(candidate.roots, candidate.basis, field, work)
                c.check(_check_field_stack(inputs, candidate, field, work), "exact verifier accepts interval final roots")
                c.check(candidate.sectors == exact.sectors and len(candidate.roots) == len(exact.roots), "same root count and canonical sectors as exact producer")
            native_nodes += stats.native_nodes
            wide_nodes += stats.wide_nodes
        if algebraic:
            for expr in ["x1 - x0", "2*x1 + x0", "x1 - 2*x0 + 1", "x1^2 - x0", "x1^2 + x0", "x1^2 - 2*x0*x1 - 1/100000000*x1 + 2 + 1/100000000*x0"]:
                var inputs: List[List[DensePoly[QQ]]] = [_specialize_field(parse_polynomial(expr, MonomialOrder(), work), 1, field, work)]
                var exact = _construct_field_stack(inputs, field, work)
                var stats = _IntervalRootStats()
                var roots = _interval_field_roots(exact.basis, field, stats, work)
                if roots:
                    exact.roots = roots.take()
                    exact.sectors = _ksectors(exact.roots, exact.basis, field, work)
                c.check(_check_field_stack(inputs, exact, field, work), "close algebraic roots use intervals or exact fallback safely")
                native_nodes += stats.native_nodes
                wide_nodes += stats.wide_nodes
    c.check(native_nodes > 0, "native bounded transforms exercised")
    c.check(wide_nodes > 0, "wide transforms handle native overflow and precision escalation")
    var field = _sample_field(_Frame(), work)
    var inputs: List[List[DensePoly[QQ]]] = [_specialize_field(parse_polynomial("x1^2 - 2", MonomialOrder(), work), 1, field, work)]
    var stack = _construct_field_stack(inputs, field, work)
    var stats = _IntervalRootStats()
    var limited = Budget(10000)
    c.check(not _interval_field_roots(stack.basis, field, stats, limited, cap=10), "local cap requests exact fallback")
    c.check(limited.remaining() < 10000, "failed interval attempt charges caller")
    limited = Budget(0)
    var exhausted = False
    try:
        _ = _interval_field_roots(stack.basis, field, stats, limited)
    except e:
        exhausted = e.kind == PolyErrorKind.LIMIT
    c.check(exhausted, "global limit is not converted into local fallback or capacity")
    var fallback = _construct_field_stack(inputs, field, work, interval=True, interval_cap=0)
    c.check(fallback.isolation.exact_fallbacks == 1 and _check_field_stack(inputs, fallback, field, work), "capped interval producer falls back to a complete exact stack")
    c.check(fallback.sectors == stack.sectors, "exact fallback retains canonical samples")
    c.finish()
