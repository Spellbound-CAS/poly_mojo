"""Exact Sturm corpus sequences, distinct counts, and complete witness rejection."""
from apn_mojo import Rational
from polynomials import Budget, PolyErrorKind, QQ, NoContext, DensePoly, MonomialOrder, sturm_sequence, sign_changes, count_roots, check_root_count, SturmSequence
from support.corpus import read_cases, parse_polynomial
from support.check import Checker


def main() raises:
    var work = Budget()
    var c = Checker("sturm")
    var sequences = 0
    var counts = 0
    for record in read_cases():
        var op = record.field("op ")
        if op != "sturm" and op != "count_roots":
            continue
        var text = record.field("in f ")
        var f = DensePoly[QQ].from_sparse(parse_polynomial(text, MonomialOrder(), work), work)
        if op == "sturm":
            var expected = List[DensePoly[QQ]]()
            for line in record.lines:
                if line.startswith("out sturm "):
                    var term = String(line[byte=10:])
                    expected.append(DensePoly[QQ].from_sparse(parse_polynomial(term, MonomialOrder(), work), work))
            var got = sturm_sequence(f, work)
            c.check(got.polynomials == expected, record.name)
            c.check(len(got.quotients) == len(got.polynomials) - 1, "all quotients including terminal relation")
            sequences += 1
        else:
            var at = record.field("in a ")
            var bt = record.field("in b ")
            var a = QQ(Rational(at))
            var b = QQ(Rational(bt))
            var got = count_roots(f, a, b, work)
            c.check(got.count == Int(record.field("out count ")), record.name)
            c.check(check_root_count(f, got, work), "corpus root-count certificate")
            counts += 1
    c.check(sequences == 4 and counts == 6, "all four sequences and six counts")
    var one = DensePoly[QQ].from_coefficients(0, NoContext(), [QQ.one(NoContext())], work)
    var zero = DensePoly[QQ](0, NoContext())
    var low = QQ(Rational(-20))
    var high = QQ(Rational(20))
    for n in range(1, 7):
        var f = one.copy()
        for i in range(n):
            var linear = DensePoly[QQ].from_coefficients(0, NoContext(), [QQ(Rational(-i)), QQ.one(NoContext())], work)
            f = f.mul(linear.pow(UInt64(i % 3 + 1), work), work)
        var quadratic = DensePoly[QQ].from_coefficients(0, NoContext(), [QQ.one(NoContext()), QQ.zero(NoContext()), QQ.one(NoContext())], work)
        f = f.mul(quadratic.pow(2, work), work)
        var got = count_roots(f, low, high, work)
        c.check(got.count == n and check_root_count(f, got, work), "planted distinct roots and repeated factors")
        var corrupt = got.copy()
        corrupt.count += 1
        c.check(not check_root_count(f, corrupt, work), "wrong count")
        corrupt = got.copy()
        corrupt.sequence.polynomials[0] = one.copy()
        c.check(not check_root_count(f, corrupt, work), "wrong initial polynomial")
        corrupt = got.copy()
        corrupt.sequence.polynomials[1] = one.copy()
        c.check(not check_root_count(f, corrupt, work), "wrong derivative")
        corrupt = got.copy()
        corrupt.sequence.quotients[0] = corrupt.sequence.quotients[0].add(one, work)
        c.check(not check_root_count(f, corrupt, work), "wrong recurrence quotient")
        corrupt = got.copy()
        _ = corrupt.sequence.quotients.pop()
        c.check(not check_root_count(f, corrupt, work), "missing terminal relation")
        corrupt = got.copy()
        corrupt.sequence.polynomials[len(corrupt.sequence.polynomials) - 1] = zero.copy()
        c.check(not check_root_count(f, corrupt, work), "zero terminal polynomial")
        corrupt = got.copy()
        corrupt.a = high.copy()
        c.check(not check_root_count(f, corrupt, work), "unordered certificate endpoints")
        corrupt = got.copy()
        corrupt.a = QQ.zero(NoContext())
        c.check(not check_root_count(f, corrupt, work), "endpoint root in certificate")
        var tiny = Budget(50)
        var before = f.copy()
        try:
            _ = count_roots(f, low, high, tiny)
            c.check(False, "count budget ignored")
        except e:
            c.check(e.kind == PolyErrorKind.LIMIT and f == before, "count budget preserves input")
    var constant = count_roots(one, low, high, work)
    c.check(constant.count == 0 and check_root_count(one, constant, work), "nonzero constant count")
    var signs = SturmSequence([one.copy(), zero.copy(), one.neg(work), zero.copy(), one.copy()], List[DensePoly[QQ]]())
    c.check(sign_changes(signs, low, work) == 2, "zero signs are omitted")
    for i in range(3):
        try:
            if i == 0:
                _ = sturm_sequence(zero, work)
            elif i == 1:
                _ = count_roots(one, high, low, work)
            else:
                _ = count_roots(one, low, low, work)
            c.check(False, "invalid count accepted")
        except e:
            c.check(e.kind == PolyErrorKind.DOMAIN, "invalid count DOMAIN")
    c.finish()
