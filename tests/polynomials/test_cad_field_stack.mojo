"""Shared-field stacks: embedding, completeness, nullification and corruptions."""
from apn_mojo import Rational, Integer
from polynomials import Budget, QQ, NoContext, DensePoly, RealAlgebraic, MonomialOrder, real_root, real_rational, algebraic_roots, algebraic_compare, PolyErrorKind
from polynomials._algebraic import _Frame
from polynomials._algebraic_filter import _constant, _child_frame
from polynomials._cad_sample_field import _sample_field, _check_sample_field, _ksign, _specialize_field, _lazard_field, _kevaluate
from polynomials._cad_field_stack import _construct_field_stack, _check_field_stack, _krefine
from polynomials._cad_lifting import _sector
from support.check import Checker
from support.corpus import parse_polynomial


def main() raises:
    var work = Budget()
    var c = Checker("field-relative final stack")
    var order = MonomialOrder()
    var zero = QQ.zero(NoContext())
    var one = QQ.one(NoContext())
    var m = DensePoly[QQ].from_coefficients(0, NoContext(), [QQ(Rational(-2)), zero.copy(), one.copy()], work)
    var gamma = real_root(m, one, QQ(Rational(2)), work)
    var frame = _Frame()
    frame.generator = gamma.copy()
    frame.variables.append(0)
    var z = DensePoly[QQ].from_coefficients(0, NoContext(), [zero.copy(), one.copy()], work)
    frame.representations.append(z.copy())
    var field = _sample_field(frame, work)
    c.check(_check_sample_field(field, work), "one irreducible field and real embedding")
    c.check(_ksign(z, field, work) == 1, "positive embedding sign")
    var other = field.copy()
    other.embedding = real_root(m, QQ(Rational(-2)), QQ(Rational(-1)), work)
    c.check(_check_sample_field(other, work) and _ksign(z, other, work) == -1, "same modulus with different real embedding")
    c.check(_ksign(m, field, work) == 0, "zero is exact reduction modulo the shared modulus")

    # Cancelling coefficients can occur before and between surviving ones.
    # Specialization must preserve every surviving degree and borrowed input.
    var cancellation = parse_polynomial("x0^2*x1^4 - 2*x1^4 + x0*x1^3 + x0^2*x1^2 - 2*x1^2 + 3*x1 + x0^2 - 2", order, work)
    var original_cancellation = cancellation.copy()
    var cancellation_rows = _specialize_field(cancellation, 1, field, work)
    c.check(len(cancellation_rows) == 4 and cancellation_rows[0].is_zero() and cancellation_rows[2].is_zero(), "cancelled specialization degrees remain zero")
    c.check(cancellation_rows[1] == _constant(QQ(Rational(3)), work) and cancellation_rows[3] == z, "surviving specialization coefficients keep their degrees")
    c.check(cancellation == original_cancellation and field.coordinates[0] == z, "specialization preserves borrowed polynomial and sample coordinates")

    # y^2 - sqrt(2): two relative fourth roots, no absolute root construction.
    var f = parse_polynomial("x1^2 - x0", order, work)
    var rows = _specialize_field(f, 1, field, work)
    var inputs: List[List[DensePoly[QQ]]] = [rows.copy()]
    var stack = _construct_field_stack(inputs, field, work)
    c.check(len(stack.roots) == 2 and _check_field_stack(inputs, stack, field, work), "complete field-relative fourth-root stack")
    var minus = real_root(m, QQ(Rational(-2)), QQ(Rational(-1)), work)
    var absolute = algebraic_roots([minus.copy(), real_rational(zero, work), real_rational(one, work)], work)
    c.check(len(absolute.roots) == len(stack.roots), "root count agrees with old algebraic path")
    for i in range(len(stack.roots)):
        var lo = real_rational(stack.roots[i].low, work)
        var hi = real_rational(stack.roots[i].high, work)
        c.check(algebraic_compare(lo, absolute.roots[i], work).relation < 0 and algebraic_compare(absolute.roots[i], hi, work).relation < 0, "relative bracket contains the corresponding old root")
    for i in range(len(stack.sectors)):
        var expected = _sector(absolute.roots, i, work)
        c.check(stack.sectors[i] == expected.low, "canonical simplest rational agrees")

    for damage in range(11):
        var bad = stack.copy()
        if damage == 0:
            _ = bad.roots.pop()
            _ = bad.sectors.pop()
            _ = bad.signs.pop()
            _ = bad.signs.pop()
        elif damage == 1:
            bad.roots[0].basis = 12
        elif damage == 2:
            bad.roots[1] = bad.roots[0].copy()
        elif damage == 3:
            bad.multiplicities[0][0] = 2
        elif damage == 4:
            bad.leading[0] = -1
        elif damage == 5:
            bad.signs[0][0] = -bad.signs[0][0]
        elif damage == 6:
            bad.sectors[1] = bad.roots[1].high.copy()
        elif damage == 7:
            bad.basis[0][0] = _constant(QQ(Rational(-2)), work)
        elif damage == 8:
            bad.roots[0].low = QQ(Rational(-100))
            bad.roots[0].high = QQ(Rational(-99))
            bad.sectors[0] = QQ(Rational(-101))
        elif damage == 9:
            bad.roots[0].high = bad.roots[0].low.copy()
        else:
            bad.signs[0][0] = 2
        c.check(not _check_field_stack(inputs, bad, field, work), "forged final field stack rejected")
    c.check(not _check_field_stack(inputs, stack, other, work), "conjugate embedding cannot reuse positive-field roots")
    var forged = field.copy()
    var extra = DensePoly[QQ].from_coefficients(0, NoContext(), [QQ(Rational(-3)), zero.copy(), one.copy()], work)
    forged.modulus = m.mul(extra, work)
    forged.embedding.value().polynomial = forged.modulus.copy()
    forged.embedding.value().low = QQ(Rational(7, 5))
    forged.embedding.value().high = QQ(Rational(3, 2))
    c.check(not _check_sample_field(forged, work), "isolating an irrational root does not make a reducible modulus a field")
    forged = field.copy()
    forged.coordinates[0] = m.copy()
    c.check(not _check_sample_field(forged, work), "unreduced coordinate map rejected")

    # Shared factors with different multiplicities and exact rational roots.
    inputs = List[List[DensePoly[QQ]]]()
    for expr in ["x1^3 - x0*x1^2", "x1^2 - 2*x0*x1 + 2", "0", "x0"]:
        inputs.append(_specialize_field(parse_polynomial(expr, order, work), 1, field, work))
    stack = _construct_field_stack(inputs, field, work)
    c.check(len(stack.roots) == 2 and _check_field_stack(inputs, stack, field, work), "shared and repeated field roots, zero and constant inputs")
    c.check(stack.signs[0][0] == -1 and stack.signs[2][0] == -1 and stack.signs[4][0] == 1, "even and odd multiplicities propagate signs")
    c.check(stack.roots[0].low.is_zero() and stack.roots[0].high.is_zero(), "rational section retained as an exact root")

    var rational = _sample_field(_Frame(), work)
    inputs = List[List[DensePoly[QQ]]]()
    inputs.append(_specialize_field(parse_polynomial("x0^2 - 2", order, work), 0, rational, work))
    stack = _construct_field_stack(inputs, rational, work)
    c.check(_check_field_stack(inputs, stack, rational, work), "rational base uses the same final-stack obligations")
    inputs = [_specialize_field(parse_polynomial("3*x0 - 1", order, work), 0, rational, work)]
    stack = _construct_field_stack(inputs, rational, work)
    c.check(_check_field_stack(inputs, stack, rational, work), "nondyadic rational root has terminating canonical sampling")
    c.check(stack.sectors == [QQ(Rational(0)), QQ(Rational(1))], "nondyadic boundary preserves simplest samples")
    inputs = [_specialize_field(parse_polynomial("3*x0^2 - x0", order, work), 0, rational, work)]
    stack = _construct_field_stack(inputs, rational, work)
    c.check(_check_field_stack(inputs, stack, rational, work) and stack.sectors[1] == QQ(Rational(1, 4)), "nondyadic rational boundary between two sections terminates exactly")

    var child = _child_frame(frame, 1, parse_polynomial("x1^2 - 3", order, work), real_root(extra, one, QQ(Rational(2)), work), work)
    c.check(child.valid and Bool(child.generator), "a further lift combines two independent irrational coordinates")
    var combined = _sample_field(child, work)
    c.check(_check_sample_field(combined, work) and combined.modulus.degree() == 4, "child sample uses one degree-four field")
    var regular = parse_polynomial("x0*x1*x2^2 - x0 + x1", order, work)
    var regular_lift = _lazard_field(regular, 2, combined, work)
    c.check(regular_lift[0] == [Int(0), Int(0)] and regular_lift[1] == _specialize_field(regular, 2, combined, work), "nonzero full specialization gives zero prefix orders")
    inputs = List[List[DensePoly[QQ]]]()
    for expr in ["x2 - 1", "2*x2 - 2", "-3*x2 + 3", "x2 + 1", "0"]:
        inputs.append(_specialize_field(parse_polynomial(expr, order, work), 2, combined, work))
    stack = _construct_field_stack(inputs, combined, work, interval=True)
    c.check(len(stack.basis) == 2 and len(stack.roots) == 2 and _check_field_stack(inputs, stack, combined, work), "duplicate linear factors share exact rational sections in an algebraic field")
    var nullified = parse_polynomial("x0^2*x2^2 - 2*x2^2 + x1^2*x2 - 3*x2", order, work)
    var lazard = _lazard_field(nullified, 2, combined, work)
    c.check(lazard[0] == [Int(0), Int(1)], "prefix Lazard order distinguishes later nullification")
    c.check(len(_specialize_field(nullified, 2, combined, work)) == 0, "original polynomial is zero despite nonzero Lazard lift")
    c.check(len(lazard[1]) == 2 and _ksign(lazard[1][1], combined, work) == 1, "Lazard derivative has positive algebraic leading coefficient")
    inputs = List[List[DensePoly[QQ]]]()
    inputs.append(_specialize_field(parse_polynomial("x2^2 - x0 - x1", order, work), 2, combined, work))
    stack = _construct_field_stack(inputs, combined, work)
    c.check(len(stack.roots) == 2 and _check_field_stack(inputs, stack, combined, work), "another relative lift over the combined sample field")
    var limited = Budget(0)
    var exhausted = False
    try:
        _ = _check_field_stack(inputs, stack, combined, limited)
    except e:
        exhausted = e.kind == PolyErrorKind.LIMIT
    c.check(exhausted, "final verification propagates global exhaustion")
    c.finish()
