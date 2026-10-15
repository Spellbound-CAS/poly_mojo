"""Exchange representations with core algorithms using explicit variable IDs."""
from poly import series, basis, rational
from apn_mojo import ArithmeticContext, FloatFormat, Rational
from std.testing import assert_equal


def main() raises:
    var context = ArithmeticContext(format=FloatFormat(precision=80))
    var s = series.asarray([Rational(1, 3), Rational(2)], order=5, context=context)
    var core_s = series.to_core(s, variable=-19)
    var restored_s = series.from_core(core_s, context=context)
    assert_equal(restored_s[0].to_json(), s[0].to_json())
    print("Restored series order:", restored_s.order())
    var b = basis.asarray([1, 2, 0], basis="bernstein", context=context)
    var core_b = basis.to_core(b, variable=42)
    assert_equal(len(basis.from_core(core_b, context=context)), 3)
    var f = rational.asarray([1], [1, 1], context=context)
    var ids: List[Int64] = [9]
    var core_f = rational.to_core(f, variables=ids)
    var restored_f = rational.from_core[1](core_f, variables=ids)
    assert_equal(rational.evaluate(restored_f, 1, context=context).to_rational_exact(), Rational(1, 2))
    print("Representation copies preserve stored coefficients")
