"""Algebraic root covers retain conversion evidence across refinement and streams."""
from poly import dense, sparse, algebraic, roots, io, coefficients
from polynomials import Budget, NumericOptions, PolyError, PolyErrorKind, QQ
from apn_mojo import Rational
from std.testing import assert_true, assert_false, assert_equal

def main() raises:
    var a = dense.root_between([-2, 0, 1], 1, 2)
    var p = dense.asarray(coefficients((algebraic.negative(a), 1)))
    var options = NumericOptions(initial_precision=24, max_precision=32, absolute_bits=2, max_subdivisions=128)
    var result = dense.real_root_cover(p, options=options, complete_only=True)
    assert_true(roots.check(result))
    assert_true(result.complete())
    assert_true(result.accurate())
    assert_true(result.isolated())
    var regions = roots.intervals(result)
    assert_equal(len(regions), 1)
    assert_equal(regions[0][2], 1)
    assert_true(regions[0][0] < Rational(3, 2))
    assert_true(regions[0][1] > Rational(7, 5))
    var saved = io.loads[type_of(result)](io.dumps(result))
    assert_true(roots.check(saved))
    var fine = NumericOptions(initial_precision=32, max_precision=48, absolute_bits=4, max_subdivisions=256)
    assert_true(roots.check(roots.refine(result, fine, complete_only=True)))
    var work = Budget()
    var checker = dense.real_root_checker(p, work, options=options)
    var receipt = dense.real_root_cover_to(p, checker, work, options=options, complete_only=True)
    assert_true(checker.finish(receipt, work))
    receipt.items += 1
    assert_false(checker.finish(receipt, work))
    var sp = [(0, algebraic.negative(a)), (1, algebraic.asreal(1))]
    var sc = sparse.real_root_checker(sp, options=options)
    var sr = sparse.real_root_cover_to(sp, sc, options=options, complete_only=True)
    assert_true(sc.finish(sr, work))
    assert_true(roots.check(sparse.real_root_cover(sp, options=options)))
    var zp = [algebraic.negative(algebraic.ascomplex(a, 1)), algebraic.ascomplex(1)]
    var zr = dense.complex_root_cover(zp, options=options, complete_only=True)
    assert_true(roots.check(zr))
    assert_equal(roots.rectangles(zr)[0][4], 1)
    assert_true(roots.check(io.loads[type_of(zr)](io.dumps(zr))))
    var zc = dense.complex_root_checker(zp, options=options)
    var zreceipt = dense.complex_root_cover_to(zp, zc, options=options, complete_only=True)
    assert_true(zc.finish(zreceipt, work))
    var zs = [(0, zp[0].copy()), (1, zp[1].copy())]
    var zsc = sparse.complex_root_checker(zs, options=options)
    var zsr = sparse.complex_root_cover_to(zs, zsc, options=options, complete_only=True)
    assert_true(zsc.finish(zsr, work))
    assert_true(roots.check(sparse.complex_root_cover(zs, options=options)))
    # A forged source or conversion must not be accepted with an old cover.
    saved._source = dense.asarray(coefficients((algebraic.negative(a), 2)))
    assert_false(roots.check(saved))
    var damaged = io.loads[type_of(result)](io.dumps(result))
    damaged._converted[].witnesses[0].real.low = QQ(Rational(100))
    assert_false(roots.check(damaged))
    var limited = False
    try:
        _ = dense.real_root_cover(p, options=options, max_steps=0)
    except error:
        limited = error.kind == PolyErrorKind.LIMIT
    assert_true(limited)
    assert_true(roots.check(result))
    print("workflow algebraic cover checks passed")
