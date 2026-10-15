"""Finite-field fractions retain moduli in all exact workflows."""
from poly import rational, modint, io
from std.testing import assert_true, assert_equal
from polynomials import PolyError, PolyErrorKind

def main() raises:
    var a = modint(1, 5)
    var zero = modint(0, 5)
    var f = rational.asarray([a.copy()], [a.copy(), zero.copy(), a.copy()])
    assert_true(rational.check(f))
    assert_equal(rational.evaluate(f, a).value, UInt64(3))
    assert_equal(rational.evaluate(rational.derivative(f), a).value, UInt64(2))
    var parts = rational.partial_fractions(f)
    assert_true(rational.check(parts))
    var expansion = rational.laurent(f, zero, terms=4)
    assert_true(rational.check(expansion))
    var saved = io.dumps(f)
    var copy = io.loads[type_of(f)](saved)
    var wrong_header = saved.replace(String("RATIONAL(1,GF)\n5\n"), String("RATIONAL(1,GF)\n7\n"))
    assert_true(wrong_header != saved)
    var header_rejected = False
    try:
        _ = io.loads[type_of(f)](wrong_header)
    except e:
        header_rejected = e.kind == PolyErrorKind.ENCODING
    assert_true(header_rejected)
    assert_true(rational.check(copy))
    assert_equal(copy.modulus(), UInt64(5))
    var other = rational.asarray([modint(1, 7)], [modint(1, 7)])
    var failed = False
    try:
        _ = rational.add(f, other)
    except e:
        failed = e.kind == PolyErrorKind.DOMAIN
    assert_true(failed)
    var inseparable = rational.asarray([modint(1, 2)], [modint(1, 2), modint(0, 2), modint(1, 2)])
    assert_true(rational.check(rational.partial_fractions(inseparable)))
    assert_equal(rational.evaluate(rational.derivative(inseparable), modint(0, 2)).value, UInt64(0))
    assert_true(rational.check(rational.laurent(inseparable, modint(0, 2), terms=4)))
    print("completion finite rational checks passed")
