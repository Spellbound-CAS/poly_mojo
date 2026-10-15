"""Formal exponential/logarithm identities across exact fields and characteristics."""
from apn_mojo import Rational, ExactComplex
from polynomials import *
from polynomials._toolkit_exact import _ExactOps
from support.check import Checker


def field[C: CoefficientField](ctx: C.Context, mut check: Checker) raises:
    var work = Budget()
    for n in [1, 2, 15, 16, 63, 64, 65, 127, 128, 129, 255, 256, 257, 1024, 1025, 2049]:
        comptime if C != ModInt:
            if n > 257:
                continue
        var values = List[C](length=n, fill=C.one(ctx))
        var geometric = TruncatedSeries[C].from_polynomial(DensePoly[C].from_coefficients(0, ctx, values^, work), n, work)
        var logarithm = geometric.log(work)
        var valid = logarithm.polynomial.coefficients[0].is_zero() if n > 1 else logarithm.polynomial.is_zero()
        for k in range(1, n):
            valid = valid and logarithm.polynomial.coefficients[k] == _ExactOps[C].integer(k, ctx, work).inverse(work)
        check.check(valid, "log(1/(1-x)) has coefficients 1/k")
        var exponential = logarithm.exp(work)
        check.check(exponential.polynomial == geometric.polynomial, "exp(sum x^k/k) is the geometric series")
        values = List[C]()
        values.append(C.zero(ctx))
        if n > 1:
            values.append(C.one(ctx))
        var x = TruncatedSeries[C].from_polynomial(DensePoly[C].from_coefficients(0, ctx, values^, work), n, work)
        var exp_x = x.exp(work)
        var factorial = C.one(ctx)
        valid = True
        for k in range(n):
            if k:
                factorial = factorial.mul(_ExactOps[C].integer(k, ctx, work), work)
            valid = valid and exp_x.polynomial.coefficients[k] == factorial.inverse(work)
        check.check(valid, "exp(x) has inverse factorial coefficients")
        var tiny = Budget(10)
        try:
            _ = logarithm.exp(tiny)
            check.check(False, "exact exponential ignored budget")
        except e:
            check.check(e.kind == PolyErrorKind.LIMIT and geometric.order == n, "exhaustion preserves exact input")


def main() raises:
    var check = Checker("series transcendental")
    var work = Budget()
    field[QQ](NoContext(), check)
    field[QQi](NoContext(), check)
    field[ModInt](Modulus(65521, work), check)
    for prime in [2, 7, 257]:
        var ctx = Modulus(UInt64(prime), work)
        for order in [prime, prime + 1, 2 * prime + 1]:
            var unit = TruncatedSeries[ModInt].from_polynomial(DensePoly[ModInt].from_coefficients(0, ctx, [ModInt.one(ctx)], work), order, work)
            var zero = TruncatedSeries[ModInt].from_polynomial(DensePoly[ModInt](0, ctx), order, work)
            for operation in [0, 1]:
                try:
                    _ = unit.log(work) if operation == 0 else zero.exp(work)
                    check.check(order <= prime, "all required integers remain invertible at characteristic boundary")
                except e:
                    check.check(order > prime and e.kind == PolyErrorKind.DOMAIN, "characteristic obstruction rejected even for zero numerator")
    # Dense inputs cross the accelerated dispatch thresholds, so these cases
    # also exercise characteristic failures inside Newton iteration/division.
    for operation in [0, 1]:
        var ctx = Modulus(257, work)
        var values = List[ModInt](length=16, fill=ModInt.one(ctx))
        if operation:
            values[0] = ModInt.zero(ctx)
        var dense = TruncatedSeries[ModInt].from_polynomial(DensePoly[ModInt].from_coefficients(0, ctx, values^, work), 2049 if operation == 0 else 1024, work)
        try:
            _ = dense.log(work) if operation == 0 else dense.exp(work)
            check.check(False, "accelerated series ignored characteristic obstruction")
        except e:
            check.check(e.kind == PolyErrorKind.DOMAIN and len(dense.polynomial.coefficients) == 16, "accelerated characteristic failure preserves input")
    var ctx = NoContext()
    var n = 257
    var values = List[QQi]()
    values.append(QQi.zero(ctx))
    var power = QQi.one(ctx)
    var imaginary = QQi(ExactComplex(0, 1))
    for k in range(1, n):
        power = power.mul(imaginary, work)
        values.append(power.mul(_ExactOps[QQi].integer(k, ctx, work).inverse(work), work))
    var f = TruncatedSeries[QQi].from_polynomial(DensePoly[QQi].from_coefficients(0, ctx, values^, work), n, work)
    var exp_f = f.exp(work)
    power = QQi.one(ctx)
    var valid = True
    for k in range(n):
        valid = valid and exp_f.polynomial.coefficients[k] == power
        power = power.mul(imaginary, work)
    check.check(valid, "non-real exp(sum (i*x)^k/k) equals 1/(1-i*x)")
    check.finish()
