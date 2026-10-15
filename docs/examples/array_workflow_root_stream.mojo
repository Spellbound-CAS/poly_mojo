"""Check exact complex roots as they arrive, including multiplicity."""
from poly import dense
from polynomials import Budget
from std.testing import assert_equal, assert_true

def main() raises:
    var p = dense.asarray([1, 0, 2, 0, 1])  # (x*x + 1)**2
    var work = Budget()
    var checker = dense.complex_roots_checker(p, work)
    var receipt = dense.complex_roots_to(p, checker, work)
    assert_true(checker.finish(receipt, work))
    assert_equal(receipt.items, 2)
    assert_equal(receipt.degree, 4)
    print("Distinct complex roots:", receipt.items)
    print("Total multiplicity:", receipt.degree)
