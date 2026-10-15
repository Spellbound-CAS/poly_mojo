"""Domain, budget, nonmutation, and corrupted-evidence solver boundaries."""
from apn_mojo import Rational, ExactComplex
from polynomials import Budget, PolyError, PolyErrorKind, NumericOptions, check_root_count, check_algebraic_roots, check_real_roots_exact, check_complex_roots_exact, check_membership, check_cad, check_cylindrical_covering
from poly import dense, sparse, algebraic, roots, cad, coefficients
from std.memory import ArcPointer
from std.testing import assert_equal, assert_true


def workflow(mut budget: Budget) raises PolyError:
    var p = dense.asarray([-1, 1], budget)
    var result = dense.real_roots(p, budget)
    _ = roots.check(result, budget)
    var basis = sparse.groebner([[(0, -1), (1, 1)]], budget)
    _ = sparse.contains([(0, -1), (1, 1)], basis, budget)
    var cover = dense.real_root_cover(p, budget, options=NumericOptions(absolute_bits=2), complete_only=True)
    _ = roots.check(cover, budget)
    var constraint = cad.zero([(0, -1), (1, 1)], budget)
    var solution = cad.solve([constraint^], budget)
    _ = cad.check(solution, budget)


def main() raises:
    var meter = Budget()
    workflow(meter)
    var cost = UInt64.MAX - meter.remaining()
    var exact = Budget(cost)
    workflow(exact)
    assert_equal(exact.remaining(), UInt64(0))
    var short = Budget(cost - 1)
    var failed = False
    try:
        workflow(short)
    except e:
        failed = e.kind == PolyErrorKind.LIMIT
    assert_true(failed)
    var p = dense.asarray([-2, 0, 1])
    for scenario in range(13):
        failed = False
        try:
            if scenario == 0:
                _ = dense.real_roots(List[Int]())
            elif scenario == 1:
                _ = dense.count_roots([-1, 1], 1, 2)
            elif scenario == 2:
                _ = dense.root_between(p, -2, 2)
            elif scenario == 3:
                _ = algebraic.reciprocal(0)
            elif scenario == 4:
                _ = algebraic.refine(algebraic.asreal(1), 0)
            elif scenario == 5:
                _ = sparse.real_roots([((0, 0), -2), ((2, 0), 1)])
            elif scenario == 6:
                _ = sparse.groebner([[(1, 1)]], order="kernel")
            elif scenario == 7:
                _ = sparse.groebner([[((1, 0), 1)]], axes=[0, 0])
            elif scenario == 8:
                _ = cad.decompose([[((1, 0), 1)]], axes=[1])
            elif scenario == 9:
                _ = dense.real_root_cover(p, options=NumericOptions(absolute_bits=-1))
            elif scenario == 10:
                _ = dense.real_root_cover(p, max_steps=0)
            elif scenario == 11:
                _ = sparse.real_roots([(1000000000, 1)], max_steps=10000)
            else:
                _ = dense.real_roots(p, max_steps=0)
        except e:
            failed = e.kind == (PolyErrorKind.LIMIT if scenario >= 10 else PolyErrorKind.DOMAIN)
        assert_true(failed)
    assert_equal(p.to_list(), dense.asarray([-2, 0, 1]).to_list())
    # Heterogeneous generators promote once across the whole tuple.
    var integer = sparse.asarray([((1, 0), 1)])
    var rational = sparse.asarray([((0, 1), Rational(1, 3))])
    var mixed = sparse.groebner((integer.copy(), rational.copy()))
    assert_equal(mixed[0].domain(), "QQ")
    assert_true(sparse.check_groebner(mixed))
    var imaginary = sparse.asarray([((0, 1), ExactComplex(0, 1))])
    var gaussian = sparse.groebner((integer.copy(), imaginary.copy()))
    assert_equal(gaussian[0].domain(), "QQi")
    assert_true(sparse.check_groebner(gaussian))
    var modular = sparse.asarray([((0, 1), 1)], modulus=5)
    var finite = sparse.groebner((integer.copy(), modular.copy()))
    assert_equal(finite[0].domain(), "GF(5)")
    assert_true(sparse.check_groebner(finite))
    # Corruption is applied only to independent copies of borrowed evidence.
    var counted = dense.count_roots(p, -2, 2)
    var count_evidence = counted.evidence().copy()
    count_evidence.count += 1
    var q = dense.to_core(coefficients((Rational(-2), 0, 1)))
    var work = Budget()
    assert_true(not check_root_count(q, count_evidence, work))
    var exact_roots = dense.real_roots(p)
    var evidence = exact_roots.evidence().copy()
    evidence.accepted.clear()
    assert_true(not check_algebraic_roots(exact_roots._data[].coefficients, evidence, work))
    assert_true(roots.check(exact_roots))
    var wrong_order = exact_roots.copy()
    var swapped = exact_roots._values._values[].copy()
    swapped.swap_elements(0, 1)
    wrong_order._values._values = ArcPointer(swapped^)
    assert_true(not roots.check(wrong_order))
    var options = NumericOptions(absolute_bits=2)
    var cover = dense.real_root_cover(p, options=options, complete_only=True)
    var bad_cover = cover.evidence().copy()
    bad_cover.cells.clear()
    assert_true(not check_real_roots_exact(q, bad_cover, options, work))
    assert_true(roots.check(cover))
    var basis = sparse.groebner([[(1, 1)]])
    var member = sparse.membership([(1, 1)], basis)
    var bad_member = member.value().evidence().copy()
    bad_member.cofactors.clear()
    assert_true(not check_membership(member.value()._data[].polynomial, member.value()._data[].basis[].generators, bad_member, work))
    var inputs = [[(0, -1), (2, 1)]]
    var decomposition = cad.decompose(inputs)
    var bad_cad = decomposition.evidence().copy()
    bad_cad.cells.clear()
    assert_true(not check_cad(decomposition._data[].inputs, bad_cad, work))
    var solution = cad.solve([cad.negative(inputs[0])])
    var bad_solution = solution.evidence().copy()
    bad_solution.witness_cell = 99999
    assert_true(not check_cylindrical_covering(solution._data[].inputs, solution._data[].masks, bad_solution, work))
    print("array solver boundary checks passed")
