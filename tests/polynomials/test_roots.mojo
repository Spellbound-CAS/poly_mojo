"""Root isolation against independent corpus intervals and planted rational roots."""
from apn_mojo import Rational
from polynomials import Budget, PolyErrorKind, QQ, NoContext, DensePoly, MonomialOrder, RootInterval, RootIsolation, isolate_real_roots, refine, count_roots, check_root_count, sign_on_interval, evaluate_dense
from support.check import Checker
from support.corpus import read_cases, parse_polynomial


def validate(f: DensePoly[QQ], result: RootIsolation, mut c: Checker, mut work: Budget) raises:
    var width = QQ(Rational("1/1099511627776"))
    for i in range(len(result.exact_roots)):
        c.check(evaluate_dense(f, result.exact_roots[i], work).is_zero(), "exact root evaluates to zero")
        if i:
            c.check(result.exact_roots[i - 1].value < result.exact_roots[i].value, "exact roots strictly sorted")
    for i in range(len(result.intervals)):
        ref interval = result.intervals[i]
        var proof = count_roots(f, interval.low, interval.high, work)
        c.check(proof.count == 1 and check_root_count(f, proof, work), "one certified root per interval")
        c.check(interval.high.sub(interval.low, work).value < width.value, "isolation width below 2^-40")
        if i:
            c.check(result.intervals[i - 1].high.value <= interval.low.value, "intervals sorted and disjoint")
        for point in result.exact_roots:
            c.check(point.value < interval.low.value or point.value > interval.high.value, "interval excludes every exact root and endpoint")


def main() raises:
    var work = Budget()
    var c = Checker("roots")
    var count = 0
    for record in read_cases():
        if record.field("op ") != "isolate":
            continue
        var text = record.field("in f ")
        var f = DensePoly[QQ].from_sparse(parse_polynomial(text, MonomialOrder(), work), work)
        var got = isolate_real_roots(f, work)
        c.check(len(got.exact_roots) + len(got.intervals) == Int(record.field("out count ")), record.name)
        var exact = List[QQ]()
        var intervals = List[RootInterval]()
        for line in record.lines:
            if line.startswith("out root exact "):
                var number = String(line[byte=15:])
                exact.append(QQ(Rational(number)))
            elif line.startswith("out root interval "):
                var parts = line.split(" ")
                intervals.append(RootInterval(QQ(Rational(String(parts[3]))), QQ(Rational(String(parts[4])))))
        c.check(got.exact_roots == exact, "same exact rational roots")
        c.check(len(got.intervals) == len(intervals), "same number of irrational intervals")
        for i in range(len(intervals)):
            c.check(got.intervals[i].low.value < intervals[i].high.value and intervals[i].low.value < got.intervals[i].high.value, "corpus interval overlap")
        validate(f, got, c, work)
        count += 1
    c.check(count == 14, "all fourteen root isolation cases")
    var one = DensePoly[QQ].from_coefficients(4, NoContext(), [QQ.one(NoContext())], work)
    var no_real = DensePoly[QQ].from_coefficients(4, NoContext(), [QQ.one(NoContext()), QQ.zero(NoContext()), QQ.one(NoContext())], work)
    for n in range(1, 8):
        var expected = List[QQ]()
        var f = no_real.pow(2, work)
        for i in range(-n, n + 1):
            var root = QQ(Rational(i, 2 * n + 3))
            expected.append(root.copy())
            var linear = DensePoly[QQ].from_coefficients(4, NoContext(), [root.neg(work), QQ.one(NoContext())], work)
            f = f.mul(linear.pow(UInt64(abs(i) % 3 + 1), work), work)
        var got = isolate_real_roots(f, work)
        c.check(got.exact_roots == expected and len(got.intervals) == 0, "all planted rational roots including negative non-dyadics")
        validate(f, got, c, work)
    var rational = DensePoly[QQ].from_coefficients(4, NoContext(), [QQ(Rational(-1)), QQ(Rational(3))], work)
    var interval = RootInterval(QQ.zero(NoContext()), QQ.one(NoContext()))
    var narrow = QQ(Rational("1/1000000000000000000000000"))
    var refined = refine(rational, interval, narrow, work)
    c.check(refined.high.sub(refined.low, work).value < narrow.value, "strict refinement width")
    c.check(count_roots(rational, refined.low, refined.high, work).count == 1, "refined interval retains its root")
    var x = DensePoly[QQ].from_coefficients(4, NoContext(), [QQ.zero(NoContext()), QQ.one(NoContext())], work)
    var around_zero = RootInterval(QQ(Rational(-1)), QQ(Rational(1)))
    var centered = refine(x.pow(4, work), around_zero, narrow, work)
    c.check(centered.low.value.sign() < 0 and centered.high.value.sign() > 0 and centered.high.sub(centered.low, work).value < narrow.value, "exact midpoint remains strictly inside refinement")
    var minus_two = QQ(Rational(-2))
    var two = QQ(Rational(2))
    c.check(sign_on_interval(no_real, minus_two, two, work).value() == 1, "positive root-free interval")
    c.check(sign_on_interval(no_real.neg(work), minus_two, two, work).value() == -1, "negative root-free interval")
    c.check(not sign_on_interval(x, minus_two, two, work), "interior root prevents strict sign")
    c.check(not sign_on_interval(x, QQ.zero(NoContext()), two, work), "endpoint root prevents strict sign")
    c.check(sign_on_interval(one, two, two, work).value() == 1, "singleton interval sign")
    var zero = DensePoly[QQ](4, NoContext())
    c.check(not sign_on_interval(zero, minus_two, two, work), "zero has no strict sign")
    var constants = isolate_real_roots(one, work)
    c.check(len(constants.exact_roots) == 0 and len(constants.intervals) == 0, "constant has no roots")
    var budget_input = DensePoly[QQ].from_coefficients(4, NoContext(), [QQ(Rational(-2)), QQ.zero(NoContext()), QQ.zero(NoContext()), QQ.one(NoContext())], work)
    var before = budget_input.copy()
    var budgets: List[UInt64] = [0, 100, 10000]
    for limit in budgets:
        var tiny = Budget(limit)
        try:
            _ = isolate_real_roots(budget_input, tiny)
            c.check(False, "isolation budget ignored")
        except e:
            c.check(e.kind == PolyErrorKind.LIMIT and budget_input == before, "isolation fails closed and preserves input")
    for error in range(6):
        try:
            if error == 0:
                _ = isolate_real_roots(zero, work)
            elif error == 1:
                _ = refine(x, interval, narrow, work)
            elif error == 2:
                _ = refine(no_real, around_zero, narrow, work)
            elif error == 3:
                _ = refine(rational, interval, QQ.zero(NoContext()), work)
            elif error == 4:
                _ = refine(x.pow(2, work).sub(one, work), RootInterval(minus_two.copy(), two.copy()), narrow, work)
            else:
                _ = sign_on_interval(x, two, minus_two, work)
            c.check(False, "invalid root operation accepted")
        except e:
            c.check(e.kind == PolyErrorKind.DOMAIN, "invalid root operation DOMAIN")
    c.finish()
