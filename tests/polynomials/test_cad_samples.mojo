"""Canonical CAD sector samples: the simplest rational of each gap between roots."""
from apn_mojo import Integer, Rational
from polynomials import Budget, QQ, SparsePoly, RealAlgebraic, MonomialOrder, real_rational, cad, check_cad, algebraic_compare
from support.check import Checker
from support.corpus import parse_polynomial


def between(value: QQ, left: Optional[RealAlgebraic], right: Optional[RealAlgebraic], mut work: Budget) raises -> Bool:
    var point = real_rational(value, work)
    if left and algebraic_compare(left.value(), point, work).relation >= 0:
        return False
    return not right or algebraic_compare(point, right.value(), work).relation < 0


def simplest_in_gap(sample: QQ, left: Optional[RealAlgebraic], right: Optional[RealAlgebraic], mut work: Budget) raises -> Bool:
    """No rational of smaller denominator, nor integer of smaller magnitude, lies strictly in the gap."""
    if not between(sample, left, right, work):
        return False
    # These gaps are wide: their simplest rationals have small denominators.
    if sample.value.denominator().magnitude_bit_length() > 20 or sample.value.numerator().magnitude_bit_length() > 40:
        return False
    var p = Int(sample.value.numerator())
    var q = Int(sample.value.denominator())
    if q == 1:
        for n in range(-abs(p) + 1, abs(p)):
            if between(QQ(Rational(n)), left, right, work):
                return False
        return True
    # A bounded gap: candidates p'/d near the gap, from the roots' outer bounds.
    var low = left.value().low.value
    var high = right.value().high.value
    for d in range(1, q):
        var start = Int((low * Rational(d)).floor())
        var stop = Int((high * Rational(d)).ceil())
        for n in range(start, stop + 1):
            if between(QQ(Rational(Integer(n), Integer(d))), left, right, work):
                return False
    return True


def main() raises:
    var work = Budget()
    var c = Checker("CAD samples")
    var order = MonomialOrder()
    var systems: List[List[String]] = [["x0^2 + x1^2 - 1"], ["x1^2 - x0^3 - 2", "x0^2 + 3*x1 - 5"], ["3*x0^2 - 7", "x0*x1^2 - 11"]]
    for system in systems:
        var inputs = List[SparsePoly[QQ]]()
        for text in system:
            inputs.append(parse_polynomial(text, order, work))
        for mode in range(2):
            var decomposition = cad(inputs, [Int64(0), Int64(1)], work, lazard=Bool(mode))
            c.check(check_cad(inputs, decomposition, work), "certified decomposition")
            var canonical = True
            for stack in decomposition.stacks:
                for i in range(len(stack.samples)):
                    ref sample = stack.samples[i]
                    if i % 2 == 1:
                        continue
                    var gap = i // 2
                    var left: Optional[RealAlgebraic] = None
                    var right: Optional[RealAlgebraic] = None
                    if gap > 0:
                        left = stack.roots[gap - 1].copy()
                    if gap < len(stack.roots):
                        right = stack.roots[gap].copy()
                    if not sample.rational or not simplest_in_gap(sample.low, left, right, work):
                        canonical = False
            c.check(canonical, "sector samples are the simplest rationals of their gaps")
    c.finish()
