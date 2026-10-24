"""Dyadic intervals enclose the exact values of their operations, at every precision."""
from apn_mojo import Integer, Rational
from polynomials import Budget, QQ, DensePoly, NoContext
from polynomials._dyadic import _Dyadic, _bounds, _point, _add, _mul, _dense_value, _zero
from support.check import Checker


struct Lcg(Movable):
    var state: UInt64

    def __init__(out self, seed: UInt64):
        self.state = seed

    def next(mut self, bound: Int) -> Int:
        self.state = self.state * 6364136223846793005 + 1442695040888963407
        return Int((self.state >> 33) % UInt64(bound))


def value(x: Integer, e: Int) raises -> Rational:
    if e >= 0:
        return Rational(x << e)
    return Rational(x) / Rational(Integer(1) << (-e))


def contains(d: _Dyadic, q: Rational) raises -> Bool:
    return value(d.low, d.exponent) <= q and q <= value(d.high, d.exponent)


def random_q(mut rng: Lcg) raises -> Rational:
    var n = rng.next(2000001) - 1000000
    var den = 1 + rng.next(99999)
    var scale = rng.next(40)
    var num = Integer(n) << scale if rng.next(2) else Integer(n)
    return Rational(num) / Rational(Integer(den))


def main() raises:
    var work = Budget()
    var c = Checker("dyadic intervals")
    var rng = Lcg(42)
    for precision in [2, 8, 64, 200]:
        for _ in range(60):
            var a = random_q(rng)
            var b = random_q(rng)
            var da = _point(QQ(a), precision, work)
            var db = _point(QQ(b), precision, work)
            c.check(contains(da, a) and contains(db, b), "conversion encloses")
            c.check(contains(_add(da, db, precision, work), a + b), "sum encloses")
            c.check(contains(_mul(da, db, precision, work), a * b), "product encloses")
            var lo = a if a <= b else b
            var hi = b if a <= b else a
            var box = _bounds(QQ(lo), QQ(hi), precision, work)
            c.check(contains(box, lo) and contains(box, hi), "interval conversion encloses both ends")
            var coefficients = List[QQ]()
            for _ in range(4):
                coefficients.append(QQ(random_q(rng)))
            var p = DensePoly[QQ].from_coefficients(0, NoContext(), coefficients^, work)
            var exact = Rational(0)
            for i in range(len(p.coefficients) - 1, -1, -1):
                exact = exact * a + p.coefficients[i].value
            c.check(contains(_dense_value(p, da, precision, work), exact), "Horner encloses")
    var tiny = _point(QQ(Rational(Integer(1)) / Rational(Integer(3))), 64, work)
    c.check(tiny.sign() == 1 and _zero().sign() == 0, "signs of decided and undecided intervals")
    var neg = _point(QQ(Rational(Integer(-7)) / Rational(Integer(5))), 4, work)
    c.check(neg.sign() == -1, "negative interval at low precision")
    c.finish()
