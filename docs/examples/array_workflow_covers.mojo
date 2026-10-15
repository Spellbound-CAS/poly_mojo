"""Cover an algebraic polynomial directly and verify its streamed evidence."""
from poly import dense, algebraic, roots, coefficients
from polynomials import Budget, NumericOptions
from std.testing import assert_equal, assert_true

def main() raises:
    var a = dense.root_between([-2, 0, 1], 1, 2)
    var p = dense.asarray(coefficients((algebraic.negative(a), 1)))
    var options = NumericOptions(initial_precision=24, max_precision=32,
                                 absolute_bits=3, max_subdivisions=128)
    var cover = dense.real_root_cover(p, options=options, complete_only=True)
    assert_true(roots.check(cover))
    assert_equal(roots.intervals(cover)[0][2], 1)
    var work = Budget()
    var checker = dense.real_root_checker(p, work, options=options)
    var receipt = dense.real_root_cover_to(p, checker, work, options=options,
                                         complete_only=True)
    assert_true(checker.finish(receipt, work))
    print("One root of x - sqrt(2), with conversion and cover evidence checked")
