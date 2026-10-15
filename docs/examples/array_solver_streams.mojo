"""Check root and CAD streams without collecting every witness."""
from poly import dense, cad
from polynomials import Budget, NumericOptions
from std.testing import assert_true


def main() raises:
    var work = Budget()
    var options = NumericOptions(absolute_bits=4)
    var checker = dense.real_root_checker([-2, 0, 1], work, options=options)
    var receipt = dense.real_root_cover_to([-2, 0, 1], checker, work, options=options, complete_only=True)
    assert_true(checker.finish(receipt, work))
    print("Real-root stream checked")

    var inputs = [[(0, -1), (2, 1)]]
    var geometry = cad.checker(inputs, work)
    var traversal = cad.decompose_to(inputs, geometry, work)
    assert_true(geometry.finish(traversal, work))
    print("CAD stream checked")
