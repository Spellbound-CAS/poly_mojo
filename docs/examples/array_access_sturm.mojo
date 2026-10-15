"""Reuse one checked Sturm chain for several interval queries."""
from poly import dense, roots, io
from std.testing import assert_equal, assert_true


def main() raises:
    var plan = dense.sturm_sequence([-2, 0, 1])
    assert_true(roots.check(plan))
    var negative = roots.count_roots(plan, -2, 0)
    var positive = roots.count_roots(plan, 0, 2)
    assert_true(roots.check(negative))
    assert_true(roots.check(positive))
    assert_equal(roots.sign_changes(plan, 0), 1)
    var chain = roots.polynomials(plan)
    assert_equal(dense.degree(chain[1]), 1)
    var restored = io.loads[roots.SturmPlan](io.dumps(plan))
    assert_equal(roots.count_roots(restored, -2, 2).count(), 2)
    print("Negative roots:", negative.count())
    print("Positive roots:", positive.count())
