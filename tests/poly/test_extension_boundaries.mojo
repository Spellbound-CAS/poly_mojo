"""Invalid inputs, shared allowances, one-pass samples, and corrupted proofs."""
from std.testing import assert_equal, assert_true, assert_false
from std.memory import ArcPointer
from std.iter import Iterable, Iterator, StopIteration
from apn_mojo import ArithmeticContext, FloatFormat, BallContext, Ball, Float, Rational, Batch, Integer
from polynomials import Budget, PolyError, PolyErrorKind, NumericOptions, NumericOutcome, QQ
from poly import dense, sparse, series, basis, roots, algebraic, evaluation


@fieldwise_init
struct Cursor(Iterator):
    comptime Element = Int
    var state: ArcPointer[Int]
    var calls: ArcPointer[Int]

    def __next__(mut self) raises StopIteration -> Int:
        self.calls[] += 1
        if self.state[] == 3:
            raise StopIteration()
        var value = self.state[]
        self.state[] += 1
        return value


@fieldwise_init
struct Once(Iterable):
    var state: ArcPointer[Int]
    var calls: ArcPointer[Int]
    comptime IteratorType[iterable_mut: Bool, //, iterable_origin: Origin[mut=iterable_mut]]: Iterator = Cursor

    def __iter__(self) -> Cursor:
        return Cursor(self.state.copy(), self.calls.copy())


def workflow(context: ArithmeticContext, bc: BallContext, mut budget: Budget) raises PolyError:
    var p = dense.interpolate([0, 1, 2], [1, 4, 9], budget, context=context)
    _ = dense.evaluate_many(p, [0, 1, 2], budget, context=context).to_list(budget)
    var q = dense.asarray([1, 1], budget, context=context)
    var d = dense.divide_checked(p, q, budget, context=context)
    _ = dense.check_division(d, budget)
    _ = d.quotient(budget)
    var s = series.asarray([1, -1], budget, context=context, order=4)
    _ = series.coefficients(series.reciprocal(s, budget, context=context), budget)
    var b = basis.from_polynomial(p, budget, basis="bernstein", context=context)
    _ = basis.coefficients(basis.subdivide(b, 1, budget, context=context)[0], budget)
    _ = evaluation.check(dense.enclose(p, 2, budget, context=bc), budget)
    _ = evaluation.check(dense.adaptive_evaluate([1, 1], 2, budget, options=NumericOptions(initial_precision=8, max_precision=32, absolute_bits=8)), budget)
    # Keep exact roots in the numerical workflow to catch exponent-limit type regressions.
    var exact = dense.real_roots([algebraic.asreal(-1, budget), algebraic.asreal(1, budget)], budget)
    _ = roots.check(exact, budget)
    var specialized = sparse.roots_at([((1, 0), 1), ((0, 1), -1)], [(1, 2)], budget, axis=0)
    _ = roots.check(specialized, budget)


def main() raises:
    var ctx = ArithmeticContext(format=FloatFormat(precision=80))
    var bc = BallContext(precision=80)
    var meter = Budget()
    workflow(ctx, bc, meter)
    var cost = UInt64.MAX - meter.remaining()
    var exact = Budget(cost)
    workflow(ctx, bc, exact)
    assert_equal(exact.remaining(), UInt64(0))
    var short = Budget(cost - 1)
    var failed = False
    try:
        workflow(ctx, bc, short)
    except e:
        failed = e.kind == PolyErrorKind.LIMIT
    assert_true(failed)
    var source = Once(ArcPointer(Int(0)), ArcPointer(Int(0)))
    var p = dense.asarray([1, 2], context=ctx)
    failed = False
    try:
        _ = dense.evaluate_many(p, source, context=ctx, max_steps=0)
    except e:
        failed = e.kind == PolyErrorKind.LIMIT
    assert_true(failed)
    assert_equal(source.calls[], 0)
    var sampled = dense.evaluate_many(p, source, context=ctx)
    assert_equal(source.calls[], 4)
    assert_equal(sampled[2].to_rational_exact(), Rational(5))
    var nodes = Once(ArcPointer(Int(0)), ArcPointer(Int(0)))
    var samples = Once(ArcPointer(Int(0)), ArcPointer(Int(0)))
    assert_equal(dense.interpolate(nodes, samples, context=ctx)[1].to_rational_exact(), Rational(1))
    assert_equal(nodes.calls[], 4)
    assert_equal(samples.calls[], 4)
    var s = series.asarray([1, -1], order=4, context=ctx)
    var uncertain = dense.asarray([Ball(1, precision=80), Ball(Rational(0), Rational(1), precision=80)], context=bc)
    var sp = sparse.asarray([((1, 0), 1)], context=ctx)
    var unspecified = BallContext()
    var matrix = Batch[Integer]([1, 2, 3, 4]).reshape([2, 2])
    for scenario in range(24):
        failed = False
        try:
            if scenario == 0:
                _ = dense.interpolate([0, 0], [1, 2], context=ctx)
            elif scenario == 1:
                _ = dense.interpolate([0, 1], [1], context=ctx)
            elif scenario == 2:
                _ = dense.divmod(dense.asarray([1], context=bc), uncertain, context=bc)
            elif scenario == 3:
                _ = series.asarray([1], context=ctx)
            elif scenario == 4:
                _ = series.asarray(s, order=5, context=ctx)
            elif scenario == 5:
                _ = series.add(s, series.truncate(s, 2, context=ctx), context=ctx)
            elif scenario == 6:
                _ = series.log(series.asarray([2], order=4, context=ctx), context=ctx)
            elif scenario == 7:
                _ = basis.asarray([1, 2], basis="newton", nodes=[0, 1], context=ctx)
            elif scenario == 8:
                _ = basis.from_polynomial([0, 0, 1], basis="bernstein", degree=1, context=ctx)
            elif scenario == 9:
                _ = basis.asarray(List[Int](), basis="bernstein", context=ctx)
            elif scenario == 10:
                _ = basis.from_polynomial([1], basis="chebyshev", nodes=[0], context=ctx)
            elif scenario == 11:
                _ = sparse.substitute(sp, [(0, sp.copy()), (0, sp.copy())], context=ctx)
            elif scenario == 12:
                _ = sparse.affine_transform(sp, 1, 2, axis=2, context=ctx)
            elif scenario == 13:
                _ = dense.enclose([1], 1, context=unspecified)
            elif scenario == 14:
                _ = dense.adaptive_evaluate([1], 1, options=NumericOptions(initial_precision=0))
            elif scenario == 15:
                _ = sparse.roots_at([((1, 0), 1)], [(0, 1)], axis=0)
            elif scenario == 16:
                _ = sparse.roots_at([((1, 0), 1), ((0, 1), 1)], List[Tuple[Int, Int]](), axis=0)
            elif scenario == 17:
                _ = dense.real_roots([algebraic.asreal(0)])
            elif scenario == 18:
                _ = dense.evaluate_many(p, matrix, context=ctx)
            elif scenario == 19:
                var core = basis.to_core(basis.asarray([1], basis="chebyshev", context=ctx))
                core.kind = 9
                _ = basis.from_core(core, context=ctx)
            elif scenario == 20:
                var core = series.to_core(s)
                core.order = -1
                _ = series.from_core(core, context=ctx)
            elif scenario == 21:
                _ = dense.adaptive_evaluate(p, 1, options=NumericOptions(absolute_bits=0, relative_bits=0))
            elif scenario == 22:
                _ = series.reciprocal(series.asarray([0, 1], order=4, context=ctx), context=ctx)
            else:
                _ = sparse.evaluate_many(sparse.asarray([(1, 1)], context=ctx), matrix, context=ctx)
        except e:
            failed = e.kind == PolyErrorKind.DOMAIN
        assert_true(failed, "scenario " + String(scenario))
    assert_equal(p[0].to_rational_exact(), Rational(1))
    assert_equal(s[0].to_rational_exact(), Rational(1))
    var divided = dense.divide_checked(p, dense.asarray([1], context=ctx), context=ctx)
    var corrupted = divided.evidence().copy()
    corrupted.quotient.coefficients[0] = Float(100, context=ctx)
    var wrong = divided.copy()
    wrong._proof = ArcPointer(corrupted^)
    assert_false(dense.check_division(wrong))
    assert_true(dense.check_division(divided))
    var root_values = dense.real_roots([algebraic.asreal(-1), algebraic.asreal(1)])
    root_values._data[].proof.accepted.clear()
    assert_false(roots.check(root_values))
    var specialized = sparse.roots_at([((1, 0), 1), ((0, 1), -1)], [(1, 2)], axis=0)
    specialized._roots = dense.real_roots([-3, 1])
    assert_false(roots.check(specialized))
    var bad_axis = sparse.roots_at([((1, 0), 1), ((0, 1), -1)], [(1, 2)], axis=0)
    bad_axis._axis = 1
    assert_false(roots.check(bad_axis))
    print("array extension boundaries checks passed")
