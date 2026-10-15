"""Truncated products across split boundaries and non-real exact series identities."""
from apn_mojo import Rational, ExactComplex
from polynomials import *
from support.check import Checker


def main() raises:
    var check = Checker("toolkit_products")
    var work = Budget()
    var ctx = Modulus(65521, work)
    for size in [63, 64, 65, 127, 128, 129, 257]:
        for short in [size, size - 3, size // 2]:
            var a = List[ModInt]()
            var b = List[ModInt]()
            for k in range(size):
                a.append(ModInt(UInt64((k * k + 3 * k + 1) % 17), ctx))
            for k in range(short):
                b.append(ModInt(UInt64((k * k + 5 * k + 2) % 19), ctx))
            var f = TruncatedSeries[ModInt].from_polynomial(DensePoly[ModInt].from_coefficients(-2, ctx, a.copy(), work), size, work)
            var g = TruncatedSeries[ModInt].from_polynomial(DensePoly[ModInt].from_coefficients(-2, ctx, b.copy(), work), size, work)
            var got = f.mul(g, work)
            # Independent native-residue oracle, with sums reduced each term.
            var valid = True
            for k in range(size):
                var expected = UInt64(0)
                for i in range(k + 1):
                    if i < len(a) and k - i < len(b):
                        expected = (expected + a[i].value * b[k - i].value) % ctx.value
                var actual = got.polynomial.coefficients[k].value if k < len(got.polynomial.coefficients) else UInt64(0)
                valid = valid and actual == expected
            check.check(valid, "truncated Karatsuba split agrees with modular oracle")
            var saved_f = f.polynomial.copy()
            var saved_g = g.polynomial.copy()
            var tiny = Budget(71)
            try:
                _ = f.mul(g, tiny)
                check.check(False, "truncated product ignored budget")
            except e:
                check.check(e.kind == PolyErrorKind.LIMIT and f.polynomial == saved_f and g.polynomial == saved_g, "product failure preserves operands")
    var imaginary = QQi(ExactComplex(0, 1))
    var long_values = List[ModInt](length=2049, fill=ModInt.one(ctx))
    var long_series = TruncatedSeries[ModInt].from_polynomial(DensePoly[ModInt].from_coefficients(0, ctx, long_values^, work), 2049, work)
    var long_inverse = long_series.inverse(work)
    var long_expected = DensePoly[ModInt].from_coefficients(0, ctx, [ModInt.one(ctx), ModInt.one(ctx).neg(work)], work)
    check.check(long_inverse.polynomial == long_expected, "large modular Newton inverse retains truncation boundary")
    var coefficients = List[QQi](length=65, fill=imaginary)
    var f = TruncatedSeries[QQi].from_polynomial(DensePoly[QQi].from_coefficients(0, NoContext(), coefficients^, work), 65, work)
    var inverse = f.inverse(work)
    var expected = DensePoly[QQi].from_coefficients(0, NoContext(), [QQi(ExactComplex(0, -1)), QQi(ExactComplex(0, 1))], work)
    check.check(inverse.polynomial == expected, "Gaussian rational inverse constant is not real")
    coefficients = List[QQi]()
    coefficients.append(QQi.zero(NoContext()))
    var power = QQi.one(NoContext())
    for _ in range(1, 65):
        coefficients.append(power.copy())
        power = power.mul(imaginary, work)
    f = TruncatedSeries[QQi].from_polynomial(DensePoly[QQi].from_coefficients(0, NoContext(), coefficients^, work), 65, work)
    var reverted = f.revert(work)
    power = QQi.one(NoContext())
    var valid = reverted.polynomial.coefficients[0].is_zero()
    for k in range(1, 65):
        valid = valid and reverted.polynomial.coefficients[k] == power
        power = power.mul(imaginary.neg(work), work)
    check.check(valid, "x/(1-i*x) reverts to x/(1+i*x)")
    check.finish()
