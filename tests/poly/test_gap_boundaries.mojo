"""Work limits, invalid shapes, canonical decoding, and corrupted evidence."""
from std.testing import assert_equal, assert_true
from std.memory import ArcPointer
from std.iter import Iterable, Iterator, StopIteration
from apn_mojo import Float, FloatFormat, ArithmeticContext, BallContext, Rational, Integer, Batch
from polynomials import Budget, PolyError, PolyErrorKind, QQ, NumericOptions, PolynomialTextSource, PolynomialTextSink
from poly import dense, sparse, roots, algebraic, batch


@fieldwise_init
struct Source(PolynomialTextSource):
    var text: String
    var position: Int

    def read(mut self, mut budget: Budget) raises PolyError -> Optional[String]:
        budget.charge(1)
        if self.position == self.text.byte_length():
            return None
        var chunk = String(self.text[byte=self.position:self.position + 1])
        self.position += 1
        return chunk^


@fieldwise_init
struct Sink(PolynomialTextSink):
    var text: String

    def write(mut self, var chunk: String, mut budget: Budget) raises PolyError:
        budget.charge(UInt64(chunk.byte_length()) + 1)
        self.text.write(chunk)


@fieldwise_init
struct _Counter(Copyable, Movable):
    var calls: Int
    var position: Int


@fieldwise_init
struct Once(Iterable):
    var state: ArcPointer[_Counter]
    comptime IteratorType[iterable_mut: Bool, //, iterable_origin: Origin[mut=iterable_mut]]: Iterator = OnceIterator

    def __iter__(self) -> OnceIterator:
        return OnceIterator(self.state.copy())


@fieldwise_init
struct OnceIterator(Iterator):
    comptime Element = Int
    var state: ArcPointer[_Counter]

    def __next__(mut self) raises StopIteration -> Int:
        self.state[].calls += 1
        if self.state[].position == 3:
            raise StopIteration()
        self.state[].position += 1
        return self.state[].position


def plus(a: Int, b: Int, mut budget: Budget) raises PolyError -> Int:
    budget.charge(1)
    return a + b


def workflow(context: ArithmeticContext, mut budget: Budget) raises PolyError:
    var p = dense.asarray([1, 2, 3], budget, context=context)
    var d = dense.derivative(p, budget, context=context)
    _ = dense.evaluate(d, 2, budget, context=context)
    var text = dense.dumps(p, budget)
    _ = dense.loads_numeric[Float](text, budget)
    var a = batch.asarray([1, 2], budget)
    var b = batch.scalar(3, budget)
    _ = batch.vmap(plus, a, b, budget)
    _ = batch.accumulate(plus, a, 0, budget)
    var exact = dense.evaluate([1, 2], algebraic.asreal(3, budget), budget)
    _ = algebraic.check(exact, budget)
    _ = roots.check(roots.diagnostics([1, 2], 3, budget), budget)


def main() raises:
    var ctx = ArithmeticContext(format=FloatFormat(precision=80))
    var source = Once(ArcPointer(_Counter(0, 0)))
    var blocked = False
    try:
        _ = dense.asarray(source, context=ctx, max_steps=0)
    except e:
        blocked = e.kind == PolyErrorKind.LIMIT
    assert_true(blocked)
    assert_equal(source.state[].calls, 0)
    var once = dense.asarray(source, context=ctx)
    assert_equal(len(once), 3)
    assert_equal(source.state[].calls, 4)
    var meter = Budget()
    workflow(ctx, meter)
    var cost = UInt64.MAX - meter.remaining()
    var exact = Budget(cost)
    workflow(ctx, exact)
    assert_equal(exact.remaining(), UInt64(0))
    var short = Budget(cost - 1)
    var failed = False
    try:
        workflow(ctx, short)
    except e:
        failed = e.kind == PolyErrorKind.LIMIT
    assert_true(failed)
    var nonfinite = Float.infinity()
    var implicit_ball_context = BallContext()
    for scenario in range(14):
        failed = False
        try:
            if scenario == 0:
                _ = dense.asarray([nonfinite], context=ctx)
            elif scenario == 1:
                _ = dense.asarray([1], context=implicit_ball_context)
            elif scenario == 2:
                _ = batch.asarray([1, 2], shape=List[Int]([3]))
            elif scenario == 3:
                _ = batch.asarray([1, 2], shape=List[Int]([-2]))
            elif scenario == 4:
                _ = batch.broadcast_to(batch.asarray([1, 2]), [3])
            elif scenario == 5:
                _ = batch.vmap(plus, batch.asarray([1, 2]), batch.asarray([3, 4, 5]))
            elif scenario == 6:
                _ = sparse.specialize([((1, 0), 1)], [(2, algebraic.asreal(1))])
            elif scenario == 7:
                _ = sparse.specialize([((1, 0), 1)], [(0, algebraic.asreal(1)), (0, algebraic.asreal(2))])
            elif scenario == 8:
                _ = dense.loads[Integer]("invalid")
            elif scenario == 9:
                _ = dense.loads[Rational](dense.dumps([1, 2]))
            elif scenario == 10:
                _ = dense.asarray([1, 2], context=ctx, max_steps=0)
            elif scenario == 11:
                _ = batch.asarray([1, 2], max_steps=0)
            elif scenario == 12:
                _ = dense.loads[Integer](dense.dumps([1, 2]), max_steps=0)
            else:
                _ = roots.diagnostics([1, 2], 0, max_steps=0)
        except e:
            var expected = PolyErrorKind.LIMIT if scenario >= 10 else PolyErrorKind.ENCODING if scenario >= 8 else PolyErrorKind.DOMAIN
            failed = e.kind == expected
        assert_true(failed)
    var multidimensional = Batch[Integer]([1, 2, 3, 4]).reshape([2, 2])
    for scenario in range(2):
        failed = False
        try:
            if scenario == 0:
                _ = dense.asarray(multidimensional, context=ctx)
            else:
                _ = dense.evaluate(multidimensional, algebraic.asreal(1))
        except e:
            failed = e.kind == PolyErrorKind.DOMAIN
        assert_true(failed)
    var sink = Sink(String())
    dense.dump([1, -2, 3], sink)
    assert_equal(sink.text, dense.dumps([1, -2, 3]))
    assert_equal(dense.loads[Integer](sink.text).to_list(), dense.load[Integer](Source(sink.text.copy(), 0)).to_list())
    var sp = sparse.asarray([((2, 0), 1), ((0, 1), 2)])
    sink = Sink(String())
    sparse.dump(sp, sink)
    assert_equal(sparse.dumps(sparse.load[Integer, ndim=2](Source(sink.text.copy(), 0))), sink.text)
    var info = roots.diagnostics([-2, 0, 1], Rational(3, 2))
    var corrupted = info.evidence().copy()
    corrupted.backward_error_high = QQ(Rational(0))
    var bad_info = info.copy()
    bad_info._proof = ArcPointer(corrupted^)
    assert_true(not roots.check(bad_info))
    assert_true(roots.check(info))
    var cover = dense.real_root_cover([1, -2, 1], options=NumericOptions(absolute_bits=2), complete_only=True)
    var multiplicities = roots.multiplicities(cover)
    var proof = multiplicities.evidence().copy()
    for i in range(len(proof.multiplicities)):
        if proof.multiplicities[i]:
            proof.multiplicities[i] += 1
    var bad_multiplicities = multiplicities.copy()
    bad_multiplicities._proof = ArcPointer(proof^)
    assert_true(not roots.check(bad_multiplicities))
    assert_true(roots.check(multiplicities))
    var evaluated = dense.evaluate([1, 2], algebraic.asreal(3))
    var horner = evaluated._horner.take()
    horner[].proof.sums[0].operation = 2
    evaluated._horner = horner^
    assert_true(not algebraic.check(evaluated))
    var specialized = sparse.specialize([((1, 0), 1)], [(0, algebraic.asreal(2))])
    specialized._data[].proof.monomials.clear()
    assert_true(not algebraic.check(specialized))
    print("array gap boundaries checks passed")
