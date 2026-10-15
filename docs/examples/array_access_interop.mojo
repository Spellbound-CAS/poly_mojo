"""Keep numerical coefficients unchanged when crossing into core storage."""
from apn_mojo import ArithmeticContext, FloatFormat, Rational
from poly import sparse
from std.testing import assert_equal


def main() raises:
    var context = ArithmeticContext(format=FloatFormat(precision=80))
    var p = sparse.asarray([((2, 0), Rational(1, 3)), ((0, 1), Rational(2))], context=context)
    var variables: List[Int64] = [19, -4]
    var core = sparse.to_core(p, variables=variables)
    var restored = sparse.from_core[ndim=2](core, variables=variables)
    assert_equal(restored.to_list(), p.to_list())
    assert_equal(restored[1][1].to_json(), p[1][1].to_json())
    print("Numerical coefficients preserved across variable maps")
