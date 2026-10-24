"""Differential root counts for division-free Sturm chains in ordered fields."""
from apn_mojo import Rational, Integer
from polynomials import Budget, QQ, NoContext, DensePoly, MonomialOrder, real_root, PolyErrorKind
from polynomials._dyadic import _Dyadic, _add, _mul
from polynomials._algebraic import _Frame
from polynomials._algebraic_filter import _constant, _gcd_over_field
from polynomials._cad_sample_field import _SampleField, _sample_field, _specialize_field, _kcopy, _kderivative, _kdivmod, _kscale, _ksign, _kevaluate_sign, _kevaluate
from polynomials._cad_field_stack import _ksequence, _kvariation, _ktotal, _kgcd
from support.check import Checker
from support.corpus import parse_polynomial


def _euclidean(f: List[DensePoly[QQ]], field: _SampleField, mut work: Budget) raises -> List[List[DensePoly[QQ]]]:
    var chain: List[List[DensePoly[QQ]]] = [_kcopy(f, work), _kderivative(f, field, work)]
    while len(chain[len(chain) - 1]):
        var n = len(chain)
        var division = _kdivmod(chain[n - 2], chain[n - 1], field, work)
        if not len(division[1]):
            break
        chain.append(_kscale(division[1], _constant(QQ(Rational(-1)), work), field, work))
    return chain^


def main() raises:
    var work = Budget()
    var c = Checker("primitive field Sturm")
    var order = MonomialOrder()
    var large = Integer(1) << 128
    var box = _Dyadic(large.copy(), large.copy(), 0)
    for multiply in range(2):
        # Pay for integer arithmetic but exhaust the budget in _rounded.
        # LIMIT must reach the interval producer's exact-fallback handler.
        var limited = Budget(UInt64(17 if multiply else 5))
        var exhausted = False
        try:
            if multiply:
                _ = _mul(box, box, 64, limited)
            else:
                _ = _add(box, box, 64, limited)
        except e:
            exhausted = e.kind == PolyErrorKind.LIMIT
        c.check(exhausted, "dyadic rounding preserves work-limit failures")
    for degree in [Int(1), Int(2), Int(4)]:
        var frame = _Frame()
        if degree > 1:
            var coefficients = List[QQ](length=degree+1, fill=QQ.zero(NoContext()))
            coefficients[0] = QQ(Rational(-2))
            coefficients[degree] = QQ.one(NoContext())
            var m = DensePoly[QQ].from_coefficients(0, NoContext(), coefficients^, work)
            frame.generator = real_root(m, QQ(Rational(1)), QQ(Rational(2)), work)
        frame.variables.append(0)
        frame.representations.append(DensePoly[QQ].from_coefficients(0, NoContext(), [QQ.zero(NoContext()), QQ.one(NoContext())], work) if degree > 1 else _constant(QQ(Rational(3, 2)), work))
        var field = _sample_field(frame, work)
        var product = _specialize_field(parse_polynomial("x1^2 + x1 - x0*x1 - x0", order, work), 1, field, work)
        var linear = _specialize_field(parse_polynomial("x1 - x0", order, work), 1, field, work)
        var other = _specialize_field(parse_polynomial("x1 - x0 + 1", order, work), 1, field, work)
        var unit: List[DensePoly[QQ]] = [_constant(QQ(Rational(3)), work)]
        c.check(_kgcd(product, linear, field, work) == linear, "exact evaluation detects a shared linear factor")
        for prime_limit in [0, 1, 4096]:
            c.check(_gcd_over_field(product, linear, field.modulus, work, prime_limit=prime_limit) == linear, "modular GCD allowance and exact fallback preserve the common factor")
        c.check(len(_kgcd(product, other, field, work)) == 1, "linear field factor is coprime when its root is not shared")
        c.check(len(_kgcd(product, unit, field, work)) == 1, "nonzero constants are units over every sample field")
        var divided = _kdivmod(product, unit, field, work)
        c.check(len(divided[1]) == 0 and _kscale(divided[0], unit[0], field, work) == product, "constant division preserves the exact polynomial")
        if degree == 2:
            # Pell approximants cause severe cancellation. The exact norm
            # p^2-2*q^2 = (-1)^n fixes the sign without an interval oracle.
            var p = Integer(1)
            var q = Integer(1)
            for n in range(1, 49):
                var value = DensePoly[QQ].from_coefficients(0, NoContext(), [QQ(Rational(p)), QQ(Rational(-q))], work)
                c.check(_ksign(value, field, work) == (-1 if n % 2 else 1), "outward-rounded signs resolve Pell cancellation")
                var next_p = p + Integer(2)*q
                q = p + q
                p = next_p^
        for expression in [
            "x1^2 - x0", "-x1^2 + x0", "x1^4 - x0",
            "-3*x1^4 + 2*x0*x1^2 + 1", "x1^5 - 5*x1^3 + 4*x1",
            "x1^4 - 2*x0*x1^2 + x0^2", "x0*x1^3 + x1^2 - 2*x1 + 1",
            "-x0*x1^5 + x1^3 - 3*x1 + 2", "x1^3 - 3*x1 + 2",
            "x0*x1 + 1", "-x0*x1 - 1"
        ]:
            var f = _specialize_field(parse_polynomial(expression, order, work), 1, field, work)
            var exact = _euclidean(f, field, work)
            var primitive = _ksequence(f, field, work)
            c.check(_ktotal(primitive, field, work) == _ktotal(exact, field, work), "total roots agree with ordinary field Euclid")
            for numerator in range(-6, 7):
                var point = QQ(Rational(numerator, 2))
                c.check(_kevaluate_sign(f, point, field, work) == _ksign(_kevaluate(f, point, field, work), field, work), "endpoint filter agrees with exact field evaluation, including zero")
                c.check(_kvariation(primitive, point, field, work) == _kvariation(exact, point, field, work), "endpoint variations agree, including zeros and repeated roots")
    c.finish()
