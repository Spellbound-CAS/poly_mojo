"""Shared adapters preserve variable maps, factor ownership, and work limits."""
from poly import dense, sparse
from apn_mojo import Integer, Rational, ArithmeticContext, FloatFormat
from polynomials import Budget, PolyErrorKind
from std.testing import assert_equal, assert_true


def variable_maps() raises:
    # 2*x^2*y + 3*x*z^2 - 7 = 167 at (2, 3, 5); the fourth axis is unused.
    var terms = [((2, 1, 0, 0), 2), ((1, 0, 2, 0), 3), ((0, 0, 0, 0), -7)]
    var exact = sparse.asarray(terms)
    var context = ArithmeticContext(format=FloatFormat(precision=53))
    var numeric = sparse.asarray(terms, context=context)
    var ids: List[Int64] = [-9, 4, 17, 30]
    for a in range(3):
        for b in range(3):
            if a == b:
                continue
            var c = 3 - a - b
            var mapping: List[Int64] = [ids[a], ids[b], ids[c], ids[3]]
            var core = sparse.to_core(exact, variables=mapping)
            var numeric_core = sparse.to_core(numeric, variables=mapping)
            assert_equal(core.order.variables, mapping)
            assert_equal(numeric_core.order.variables, mapping)
            var restored = sparse.from_core[ndim=4](core, variables=mapping)
            var numeric_restored = sparse.from_core[ndim=4](numeric_core, variables=mapping)
            assert_equal(restored.to_list(), exact.to_list())
            assert_equal(numeric_restored.to_list(), numeric.to_list())
            assert_equal(sparse.evaluate(restored, (2, 3, 5, 0)), Integer(167))
            assert_equal(sparse.evaluate(numeric_restored, (2, 3, 5, 0), context=context).to_rational_exact(), Rational(167))
            # Copying either result back into mutable core storage cannot alter it.
            core.terms.clear()
            numeric_core.terms.clear()
            assert_equal(len(restored), 3)
            assert_equal(len(numeric_restored), 3)
    var core = sparse.to_core(exact, variables=ids)
    var numeric_core = sparse.to_core(numeric, variables=ids)
    var work = Budget()
    _ = sparse.from_core[ndim=4](core, work, variables=ids)
    var used = UInt64.MAX - work.remaining()
    var full = Budget(used)
    _ = sparse.from_core[ndim=4](core, full, variables=ids)
    assert_equal(full.remaining(), UInt64(0))
    for allowance in [UInt64(0), used // 2, used - 1]:
        var failed = False
        try:
            var limited = Budget(allowance)
            _ = sparse.from_core[ndim=4](core, limited, variables=ids)
        except e:
            failed = e.kind == PolyErrorKind.LIMIT
        assert_true(failed)
        assert_equal(sparse.from_core[ndim=4](core, variables=ids).to_list(), exact.to_list())
    work = Budget()
    _ = sparse.from_core[ndim=4](numeric_core, work, variables=ids)
    used = UInt64.MAX - work.remaining()
    full = Budget(used)
    _ = sparse.from_core[ndim=4](numeric_core, full, variables=ids)
    assert_equal(full.remaining(), UInt64(0))
    for allowance in [UInt64(0), used // 2, used - 1]:
        var failed = False
        try:
            var limited = Budget(allowance)
            _ = sparse.from_core[ndim=4](numeric_core, limited, variables=ids)
        except e:
            failed = e.kind == PolyErrorKind.LIMIT
        assert_true(failed)
        assert_equal(sparse.from_core[ndim=4](numeric_core, variables=ids).to_list(), numeric.to_list())


def factor_outputs() raises:
    # 6*(x - 1)^2*(x + 2); both operations retain the unit and multiplicities.
    var input = dense.asarray([Rational(12), Rational(-18), Rational(0), Rational(6)])
    var terms = sparse.asarray([(0, Rational(12)), (1, Rational(-18)), (3, Rational(6))])
    for squarefree in [False, True]:
        var work = Budget()
        var result = dense.squarefree(input, work) if squarefree else dense.factor(input, work)
        assert_equal(result[0], Rational(6))
        assert_equal(len(result[1]), 2)
        var product = dense.asarray([result[0].copy()])
        for entry in result[1]:
            var multiplicity = 2 if entry[0][0] == Rational(-1) else 1
            assert_equal(entry[0][0], Rational(-1) if multiplicity == 2 else Rational(2))
            assert_equal(entry[0][1], Rational(1))
            assert_equal(entry[1], multiplicity)
            product = dense.multiply(product, dense.power(entry[0], entry[1]))
            var independent = entry[0].to_list()
            independent.clear()
            assert_equal(len(entry[0]), 2)
        assert_equal(product.to_list(), input.to_list())
        var used = UInt64.MAX - work.remaining()
        var full = Budget(used)
        _ = dense.squarefree(input, full) if squarefree else dense.factor(input, full)
        assert_equal(full.remaining(), UInt64(0))
        var failed = False
        try:
            var limited = Budget(used - 1)
            _ = dense.squarefree(input, limited) if squarefree else dense.factor(input, limited)
        except e:
            failed = e.kind == PolyErrorKind.LIMIT
        assert_true(failed)
        assert_equal(input.to_list(), [Rational(12), Rational(-18), Rational(0), Rational(6)])
        work = Budget()
        var sparse_result = sparse.squarefree(terms, work) if squarefree else sparse.factor(terms, work)
        assert_equal(sparse_result[0], Rational(6))
        assert_equal(len(sparse_result[1]), 2)
        var sparse_product = sparse.asarray([(0, sparse_result[0].copy())])
        for entry in sparse_result[1]:
            sparse_product = sparse.multiply(sparse_product, sparse.power(entry[0], entry[1]))
        assert_equal(sparse_product.to_list(), terms.to_list())
        used = UInt64.MAX - work.remaining()
        full = Budget(used)
        _ = sparse.squarefree(terms, full) if squarefree else sparse.factor(terms, full)
        assert_equal(full.remaining(), UInt64(0))
        failed = False
        try:
            var limited = Budget(used - 1)
            _ = sparse.squarefree(terms, limited) if squarefree else sparse.factor(terms, limited)
        except e:
            failed = e.kind == PolyErrorKind.LIMIT
        assert_true(failed)
        assert_equal(terms.to_list(), [(0, Rational(12)), (1, Rational(-18)), (3, Rational(6))])


def main() raises:
    variable_maps()
    factor_outputs()
    print("cleanup adapter checks passed")
