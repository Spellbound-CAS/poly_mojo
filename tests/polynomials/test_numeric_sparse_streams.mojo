"""Independent coefficient schedules for sparse merges and streamed products."""
from apn_mojo import Rational, Float, FloatFormat, ArithmeticContext, ComplexContext, BallContext
from polynomials import *
from polynomials.numeric_domains import _numeric_apply, _numeric_zero, _numeric_integer
from support.check import Checker


def same[D: NumericDomain](a: NumericSparsePoly[D], b: NumericSparsePoly[D], mut work: Budget) raises -> Bool:
    if a.order != b.order or len(a.terms) != len(b.terms):
        return False
    for i in range(len(a.terms)):
        if a.terms[i].monomial != b.terms[i].monomial or numeric_bounds[D](a.terms[i].coefficient, work) != numeric_bounds[D](b.terms[i].coefficient, work):
            return False
    return True


def grid_monomial(x: Int, y: Int, mut work: Budget) raises -> Monomial:
    var pairs = List[Tuple[Int64, UInt32]]()
    if x:
        pairs.append((Int64(-3), UInt32(x)))
    if y:
        pairs.append((Int64(8), UInt32(y)))
    return Monomial(pairs^, work)


def bivariate[D: NumericDomain](context: D.Context, mut check: Checker) raises:
    var work = Budget()
    var side = 8 if D == ComplexBallDomain else 4
    for kind in [MonomialOrder.LEX, MonomialOrder.GRLEX, MonomialOrder.GREVLEX, MonomialOrder.KERNEL]:
        var order = MonomialOrder(kind, [Int64(8), Int64(-3)], work) if kind != MonomialOrder.KERNEL else MonomialOrder()
        var terms = List[Term[QQ]]()
        var shifted = List[Term[QQ]]()
        for x in range(side):
            for y in range(side):
                terms.append(Term[QQ](grid_monomial(x, y, work), QQ(Rational(1))))
                shifted.append(Term[QQ](grid_monomial(x + 2, y + 2, work), QQ(Rational(1))))
        var a = numeric_from_sparse[D, QQ](SparsePoly[QQ].from_terms(terms^, order, NoContext(), work), context, work)
        var b = numeric_from_sparse[D, QQ](SparsePoly[QQ].from_terms(shifted^, order, NoContext(), work), context, work)
        terms = List[Term[QQ]]()
        # The squared geometric sum has triangular coefficients; multiply by the
        # corresponding y coefficients without using polynomial multiplication.
        for x in range(2 * side - 1):
            for y in range(2 * side - 1):
                terms.append(Term[QQ](grid_monomial(x, y, work), QQ(Rational(min(x + 1, 2 * side - 1 - x) * min(y + 1, 2 * side - 1 - y)))))
        var oracle = numeric_from_sparse[D, QQ](SparsePoly[QQ].from_terms(terms^, order, NoContext(), work), context, work)
        check.check(same(a.mul(a, context, work), oracle, work), "bivariate heap respects each monomial order")
        for subtract in [False, True]:
            terms = List[Term[QQ]]()
            for x in range(side + 2):
                for y in range(side + 2):
                    var coefficient = Int(x < side and y < side) + (-1 if subtract else 1) * Int(x >= 2 and y >= 2)
                    if coefficient:
                        terms.append(Term[QQ](grid_monomial(x, y, work), QQ(Rational(coefficient))))
            oracle = numeric_from_sparse[D, QQ](SparsePoly[QQ].from_terms(terms^, order, NoContext(), work), context, work)
            check.check(same(a.sub(b, context, work) if subtract else a.add(b, context, work), oracle, work), "interleaved bivariate merge retains unique and shared terms")


def domain[D: NumericDomain](context: D.Context, mut check: Checker) raises:
    bivariate[D](context, check)
    var work = Budget()
    for kind in [MonomialOrder.LEX, MonomialOrder.GRLEX, MonomialOrder.GREVLEX, MonomialOrder.KERNEL]:
        var order = MonomialOrder(kind, [Int64(-3)], work) if kind != MonomialOrder.KERNEL else MonomialOrder()
        for n in [1, 7, 8, 9, 17, 65]:
            var ac = List[QQ]()
            var bc = List[QQ]()
            for i in range(n):
                ac.append(QQ(Rational((1 if i % 3 else -1) * (2 * i + 1), 7)))
                bc.append(QQ(Rational((1 if i % 2 else -1) * (i + 3), 11)))
            var exact_a = DensePoly[QQ].from_coefficients(-3, NoContext(), ac^, work)
            var exact_b = DensePoly[QQ].from_coefficients(-3, NoContext(), bc^, work)
            var a = NumericSparsePoly[D].from_dense(numeric_from_dense[D, QQ](exact_a, context, work), order, context, work)
            var b = NumericSparsePoly[D].from_dense(numeric_from_dense[D, QQ](exact_b, context, work), order, context, work)
            var expected = List[NumericTerm[D]]()
            for k in range(2 * n - 1):
                var coefficient = _numeric_integer[D](0, context, work)
                var seen = False
                # Formula for coefficient k; iterate i in increasing order, as
                # mandated by left-major products. No polynomial multiplication.
                for i in range(n):
                    var j = k - i
                    if j < 0 or j >= n:
                        continue
                    var value = _numeric_apply[D](a.terms[i].coefficient, b.terms[j].coefficient, 2, context, work)
                    if _numeric_zero[D](value):
                        continue
                    coefficient = _numeric_apply[D](coefficient, value, 0, context, work) if seen else value^
                    seen = True
                if not _numeric_zero[D](coefficient):
                    var pairs = List[Tuple[Int64, UInt32]]()
                    if k:
                        pairs.append((Int64(-3), UInt32(k)))
                    expected.append(NumericTerm[D](Monomial(pairs^, work), coefficient^))
            var oracle = NumericSparsePoly[D].from_terms(expected^, order, context, work)
            check.check(same(a.mul(b, context, work), oracle, work), "streamed convolution retains every rounded coefficient")
            for subtract in [False, True]:
                expected = List[NumericTerm[D]]()
                for i in range(n):
                    var right = _numeric_apply[D](_numeric_integer[D](0, context, work), b.terms[i].coefficient, 1, context, work) if subtract else b.terms[i].coefficient.copy()
                    var value = _numeric_apply[D](a.terms[i].coefficient, right, 0, context, work)
                    expected.append(NumericTerm[D](a.terms[i].monomial.copy(), value^))
                oracle = NumericSparsePoly[D].from_terms(expected^, order, context, work)
                var got = a.sub(b, context, work) if subtract else a.add(b, context, work)
                check.check(same(got, oracle, work), "linear merge matches scalar sum/difference")
            var saved = a.copy()
            var tiny = Budget(30)
            try:
                _ = a.mul(b, context, tiny)
                check.check(False, "sparse multiplication ignored budget")
            except e:
                check.check(e.kind == PolyErrorKind.LIMIT and same(a, saved, work), "stream failure preserves borrowed input")
    # Noncanonical publicly edited storage still takes the normalizing path.
    var order = MonomialOrder(MonomialOrder.LEX, [Int64(0)], work)
    var one = _numeric_integer[D](1, context, work)
    var p = NumericSparsePoly[D](order.copy())
    p.terms.append(NumericTerm[D](Monomial([(Int64(0), UInt32(1))], work), one.copy()))
    p.terms.append(NumericTerm[D](Monomial(), one.copy()))
    var canonical = NumericSparsePoly[D].from_terms(p.terms.copy(), order, context, work)
    check.check(same(p.add(p, context, work), canonical.add(canonical, context, work), work), "edited input retains normalization fallback")


def main() raises:
    var check = Checker("numerical sparse streams")
    domain[FloatDomain](ArithmeticContext(format=FloatFormat(8)), check)
    domain[ComplexDomain](ComplexContext(ArithmeticContext(format=FloatFormat(8))), check)
    domain[BallDomain](BallContext(8), check)
    domain[ComplexBallDomain](BallContext(8), check)
    var work = Budget()
    var order = MonomialOrder(MonomialOrder.LEX, [Int64(0)], work)
    var high = ArithmeticContext(format=FloatFormat(64))
    var tiny = ArithmeticContext(format=FloatFormat(2, emin=-2, emax=2))
    var a = NumericSparsePoly[FloatDomain].from_terms([NumericTerm[FloatDomain](Monomial(), Float(1))], order, high, work)
    var b = NumericSparsePoly[FloatDomain].from_terms([NumericTerm[FloatDomain](Monomial(), Float(Rational(1, 1024)))], order, high, work)
    check.check(a.sub(b, tiny, work).terms[0].coefficient == Rational(1), "individually underflowed subtraction term is skipped")
    check.finish()
