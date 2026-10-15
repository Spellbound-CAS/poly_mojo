"""AA specialization certificates bind source coefficients, assignments, and roots."""
from poly import dense, sparse, algebraic, roots, io
from polynomials import Budget, PolyError, PolyErrorKind
from std.testing import assert_equal, assert_true, assert_false

def main() raises:
    var a = dense.root_between([-2, 0, 1], 1, 2)
    var source = sparse.asarray([((1, 0), a.copy()), ((0, 1), a.copy())])
    var assignments = [(1, a.copy())]
    var specialized = sparse.specialize(source, assignments)
    assert_equal(specialized.domain(), "AA")
    assert_equal(specialized.ndim(), 2)
    assert_true(algebraic.check(specialized))
    var coefficients = algebraic.coefficients(specialized, axis=0)
    assert_equal(algebraic.compare(coefficients[0], 2), 0)
    assert_equal(algebraic.compare(coefficients[1], a), 0)
    var result = sparse.roots_at(source, assignments, axis=0)
    assert_true(roots.check(result))
    assert_true(algebraic.check(result.specialization()))
    assert_equal(len(result), 1)
    assert_equal(algebraic.compare(result[0], algebraic.negative(a)), 0)
    assert_true(roots.check(sparse.real_roots(specialized, axis=0)))
    var saved = io.loads[type_of(result)](io.dumps(result))
    assert_true(roots.check(saved))
    var saved_specialization = io.loads[type_of(specialized)](io.dumps(specialized))
    assert_true(algebraic.check(saved_specialization))
    saved._specialization._source = sparse.asarray([((1, 0), algebraic.asreal(1)), ((0, 1), a.copy())])
    assert_false(roots.check(saved))
    var assignment_forgery = io.loads[type_of(result)](io.dumps(result))
    var checking = Budget()
    assignment_forgery._specialization._point[][0][1] = algebraic._copy_witness(algebraic.asreal(0).witness(), checking)
    assert_false(roots.check(assignment_forgery))
    var roots_forgery = io.loads[type_of(result)](io.dumps(result))
    roots_forgery._roots = dense.real_roots([-1, 1])
    assert_false(roots.check(roots_forgery))
    var zero = sparse.specialize([((1, 0), a.copy()), ((0, 1), algebraic.negative(a))], [(0, 1), (1, 1)])
    assert_true(algebraic.check(zero))
    assert_equal(len(sparse.asarray(zero)), 0)
    var unchanged = sparse.specialize(source, List[Tuple[Int, Int]]())
    assert_true(algebraic.check(unchanged))
    for mode in range(3):
        var rejected = False
        try:
            if mode == 0:
                _ = sparse.specialize(source, [(1, 1), (1, 2)])
            elif mode == 1:
                _ = sparse.roots_at(source, [(0, 1)], axis=0)
            else:
                _ = sparse.roots_at(source, List[Tuple[Int, Int]](), axis=0)
        except error:
            rejected = error.kind == PolyErrorKind.DOMAIN
        assert_true(rejected)
    var work = Budget()
    _ = sparse.specialize(source, assignments, work)
    var spent = UInt64.MAX - work.remaining()
    var exact = Budget(spent)
    _ = sparse.specialize(source, assignments, exact)
    assert_equal(exact.remaining(), UInt64(0))
    var limited = False
    try:
        _ = sparse.specialize(source, assignments, max_steps=spent - 1)
    except error:
        limited = error.kind == PolyErrorKind.LIMIT
    assert_true(limited)
    assert_true(algebraic.check(specialized))
    assert_true(roots.check(result))
    print("workflow specialization checks passed")
