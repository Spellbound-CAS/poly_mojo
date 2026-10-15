"""Reusable Sturm chains count distinct roots and reject corrupted recurrence evidence."""
from poly import dense, sparse, roots, io
from apn_mojo import Rational
from polynomials import Budget, PolyError, PolyErrorKind, QQ, PolynomialTextSink, PolynomialTextSource
from std.testing import assert_equal, assert_true, assert_false


@fieldwise_init
struct Sink(PolynomialTextSink):
    var text: String

    def write(mut self, var chunk: String, mut budget: Budget) raises PolyError:
        budget.charge(UInt64(chunk.byte_length()) + 1)
        self.text.write(chunk)


@fieldwise_init
struct Source(PolynomialTextSource):
    var text: String
    var position: Int
    var empty: Bool

    def read(mut self, mut budget: Budget) raises PolyError -> Optional[String]:
        budget.charge(1)
        self.empty = not self.empty
        if self.empty: return String()
        if self.position == self.text.byte_length(): return None
        var chunk = String(self.text[byte=self.position:self.position + 1])
        self.position += 1
        return chunk^


def main() raises:
    var plan = dense.sturm_sequence([-2, 0, 1])
    assert_true(roots.check(plan))
    assert_equal(len(plan), 3)
    var chain = roots.polynomials(plan)
    assert_equal(chain[0].to_list(), List[Rational]([Rational(-2), Rational(0), Rational(1)]))
    assert_equal(chain[1].to_list(), List[Rational]([Rational(0), Rational(2)]))
    assert_equal(chain[2].to_list(), List[Rational]([Rational(2)]))
    assert_equal(roots.sign_changes(plan, -2), 2)
    assert_equal(roots.sign_changes(plan, 0), 1)
    assert_equal(roots.sign_changes(plan, 2), 0)
    var count = roots.count_roots(plan, -2, 2)
    assert_equal(count.count(), 2)
    assert_true(roots.check(count))
    assert_equal(roots.count_roots(plan, 0, Rational(3, 2)).count(), 1)
    assert_equal(roots.count_roots(plan, 2, 3).count(), 0)
    assert_true(roots.check(count.plan()))
    var repeated = dense.sturm_sequence([1, -2, 1])
    assert_equal(roots.count_roots(repeated, 0, 2).count(), 1)
    assert_equal(roots.sign_changes(repeated, 1), 0)
    var constant = dense.sturm_sequence([7])
    assert_true(roots.check(constant))
    assert_equal(roots.count_roots(constant, -10, 10).count(), 0)
    var sp = sparse.sturm_sequence([((0, 0), -2), ((0, 2), 1)], axis=1)
    assert_equal(roots.count_roots(sp, 1, 2).count(), 1)
    var saved = io.loads[roots.SturmPlan](io.dumps(plan))
    assert_true(roots.check(saved))
    var saved_count = io.loads[roots.SturmCount](io.dumps(count))
    assert_true(roots.check(saved_count))
    var sink = Sink(String())
    io.dump(count, sink)
    var source = Source(sink.text.copy(), 0, False)
    var streamed = io.load[roots.SturmCount](source^)
    assert_equal(streamed.count(), 2)
    assert_true(roots.check(streamed))
    var broken_source = Source(sink.text.copy(), 0, False)
    var field_limited = False
    try: _ = io.load[roots.SturmCount](broken_source^, max_field_bytes=1)
    except e: field_limited = e.kind == PolyErrorKind.CAPACITY
    assert_true(field_limited)
    saved_count._count = 1
    assert_false(roots.check(saved_count))
    saved._data[].sequence.quotients[0].coefficients[1] = QQ(Rational(9))
    assert_false(roots.check(saved))
    var forged_source = io.loads[roots.SturmPlan](io.dumps(plan))
    forged_source._data[].source = dense.to_core([Rational(-3), Rational(0), Rational(1)])
    assert_false(roots.check(forged_source))
    var forged_chain = io.loads[roots.SturmPlan](io.dumps(plan))
    forged_chain._data[].sequence.polynomials[2] = forged_chain._data[].sequence.polynomials[1].copy()
    assert_false(roots.check(forged_chain))
    for mode in range(6):
        var rejected = False
        try:
            if mode == 0: _ = dense.sturm_sequence([0])
            elif mode == 1: _ = roots.count_roots(plan, 2, 1)
            elif mode == 2: _ = roots.count_roots(repeated, 1, 2)
            elif mode == 3: _ = sparse.sturm_sequence([((1, 1), 1)], axis=0)
            elif mode == 4: _ = sparse.sturm_sequence([(1, 1)], axis=1)
            else: _ = io.dumps(saved)
        except e: rejected = e.kind == PolyErrorKind.DOMAIN or e.kind == PolyErrorKind.ENCODING
        assert_true(rejected)
    var meter = Budget()
    _ = roots.count_roots(plan, -2, 2, meter)
    var used = UInt64.MAX - meter.remaining()
    var exact = Budget(used)
    _ = roots.count_roots(plan, -2, 2, exact)
    assert_equal(exact.remaining(), UInt64(0))
    var limited = False
    try: _ = roots.count_roots(plan, -2, 2, max_steps=used - 1)
    except e: limited = e.kind == PolyErrorKind.LIMIT
    assert_true(limited)
    assert_true(roots.check(plan))
    for mode in range(2):
        var meter = Budget()
        if mode == 0: _ = dense.sturm_sequence([-2, 0, 1], meter)
        else: _ = roots.check(plan, meter)
        var work = UInt64.MAX - meter.remaining()
        var exact_budget = Budget(work)
        if mode == 0: _ = dense.sturm_sequence([-2, 0, 1], exact_budget)
        else: assert_true(roots.check(plan, exact_budget))
        assert_equal(exact_budget.remaining(), UInt64(0))
        var rejected = False
        try:
            if mode == 0: _ = dense.sturm_sequence([-2, 0, 1], max_steps=work - 1)
            else: _ = roots.check(plan, max_steps=work - 1)
        except e: rejected = e.kind == PolyErrorKind.LIMIT
        assert_true(rejected)
    print("access reusable Sturm checks passed")
