"""Shared squares preserve every rounded product and original witness step."""
from apn_mojo import Rational, FloatFormat, ArithmeticContext, ComplexContext, BallContext
from polynomials import *
from polynomials._numeric_powers import _NumericPowerCache
from polynomials.numeric_domains import _numeric_power
from polynomials.numeric_evaluate import _exact_sparse_program, _evaluate_program, _evaluate_scheduled, _check_program, _enclosure_context
from polynomials._numeric_storage import _storage_plan
from support.check import Checker


def powers[D: NumericDomain](context: D.Context, mut check: Checker, mut work: Budget) raises:
    var cache = _NumericPowerCache[D]()
    var value = numeric_value[D](NumericBox.point(QQ(Rational(5, 7)), QQ(Rational(0)), work), context, work)
    for exponent in [0, 1, 2, 3, 4, 17, 255, 256, 257, 65535]:
        var original = _numeric_power[D](value, UInt64(exponent), context, work)
        var shared = cache.power(0, value, UInt64(exponent), context, work)
        check.check(numeric_bounds[D](original, work) == numeric_bounds[D](shared, work), "cached binary squares retain the original rounded exponentiation")
    var other = numeric_value[D](NumericBox.point(QQ(Rational(-3, 5)), QQ(Rational(0)), work), context, work)
    for exponent in [2, 3, 17, 255]:
        check.check(numeric_bounds[D](cache.power(1, other, UInt64(exponent), context, work), work) == numeric_bounds[D](_numeric_power[D](other, UInt64(exponent), context, work), work), "different assignment slots never share a ladder")
    var terms = List[Term[QQ]]()
    for i in range(512):
        terms.append(Term[QQ](Monomial([(Int64(0), UInt32(2 * i + 1))], work), QQ(Rational(1 if i % 2 else -1, 3))))
    var exact = SparsePoly[QQ].from_terms(terms^, MonomialOrder(MonomialOrder.LEX, [Int64(0)], work), NoContext(), work)
    var numeric = numeric_from_sparse[D, QQ](exact, context, work)
    var assignments: List[Tuple[Int64, D.Scalar]] = [(Int64(0), value.copy())]
    var original = numeric._evaluate[False](assignments, context, work)
    var shared = numeric.evaluate(assignments, context, work)
    check.check(numeric_bounds[D](original, work) == numeric_bounds[D](shared, work), "public sparse dispatch preserves the complete accumulation schedule")


def witness[D: NumericDomain](mut check: Checker, mut work: Budget) raises:
    var context = _enclosure_context[D](48, 8192)
    var point = NumericBox.real_interval(QQ(Rational(1, 2)), QQ(Rational(3, 4)), work)
    comptime if D == ComplexBallDomain:
        point.imag_low = QQ(Rational(-1, 8))
        point.imag_high = QQ(Rational(1, 8))
    var value = numeric_value[D](point, context, work)
    var terms = List[Term[QQ]]()
    for i in range(256):
        terms.append(Term[QQ](Monomial([(Int64(0), UInt32(2 * i + 1))], work), QQ(Rational(i + 1, 3))))
    var exact = SparsePoly[QQ].from_terms(terms^, MonomialOrder(MonomialOrder.LEX, [Int64(0)], work), NoContext(), work)
    var program = _exact_sparse_program[D](exact, [(Int64(0), value.copy())], work)
    check.check(len(program.aliases) == len(program.instructions), "large sparse graph enables internal reuse metadata")
    var count = 0
    for i in range(len(program.aliases)):
        if program.aliases[i] >= 0:
            count += 1
            check.check(program.aliases[i] < i and program.aliases[program.aliases[i]] == -1, "aliases point directly to an earlier canonical square")
    check.check(count > 1000, "odd powers share substantial squaring work")
    var original = _evaluate_scheduled[D, False](program, context, work)
    var shared = _evaluate_program[D](program, context, work)
    check.check(len(shared.steps) == len(original.steps), "reuse does not shorten existing certificates")
    for i in range(len(shared.steps)):
        check.check(shared.steps[i] == original.steps[i], "every reused witness is identical to independent execution")
    check.check(_check_program[D](program, shared, work), "original rational checker validates reused squares")
    check.check(numeric_bounds[D](shared.value(work), work) == numeric_bounds[D](original.value(work), work), "reused final value is identical")
    var changed = program.copy()
    changed.aliases[1] = len(changed.instructions)
    try:
        _ = _storage_plan[True](changed, work)
        check.check(False, "forward alias accepted")
    except e:
        check.check(e.kind == PolyErrorKind.DOMAIN, "invalid alias rejected before slot access")


def main() raises:
    var check = Checker("shared numerical powers")
    var work = Budget()
    powers[FloatDomain](ArithmeticContext(format=FloatFormat(24)), check, work)
    powers[ComplexDomain](ComplexContext(ArithmeticContext(format=FloatFormat(24))), check, work)
    powers[BallDomain](BallContext(24), check, work)
    powers[ComplexBallDomain](BallContext(24), check, work)
    witness[BallDomain](check, work)
    witness[ComplexBallDomain](check, work)
    check.finish()
