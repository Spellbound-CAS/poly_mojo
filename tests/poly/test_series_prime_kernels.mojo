"""Prime-field series kernels: word-level square roots, prepared binary operations and the product kernels behind them."""
from poly import series
from polynomials import Budget, ModInt, Modulus
from std.testing import assert_equal, assert_true


def residues(n: Int, p: UInt64, seed: UInt64, ctx: Modulus) -> List[ModInt]:
    var out = List[ModInt](capacity=n)
    var x = seed
    for i in range(n):
        x = x * 6364136223846793005 + 1442695040888963407
        # Runs of p - 1 load every dot product with its largest terms.
        out.append(ModInt(p - 1 if (i // 5) % 4 == 1 else (x >> 1) % p, ctx))
    return out^


def naive(a: List[ModInt], b: List[ModInt], n: Int, p: UInt64) -> List[UInt64]:
    var out = List[UInt64](length=n, fill=0)
    for i in range(min(n, len(a))):
        for j in range(min(n - i, len(b))):
            out[i + j] = UInt64((UInt128(out[i + j]) + UInt128(a[i].value) * UInt128(b[j].value)) % UInt128(p))
    return out^


def values(s: series.Array[3]) -> List[UInt64]:
    var out = List[UInt64](capacity=len(s))
    for i in range(len(s)):
        out.append(s[i].value)
    return out^


def prefix(v: List[UInt64], n: Int) -> List[UInt64]:
    var out = List[UInt64](capacity=n)
    for i in range(n):
        out.append(v[i])
    return out^


def main() raises:
    var budget = Budget()
    var primes = [UInt64(3), UInt64(7), UInt64(65521), UInt64(1099511627791), UInt64(2305843009213693951), UInt64(4611686018427387847)]
    # Square roots: the root has constant term one and squares back to f,
    # including orders past a small prime, where only 2 must be invertible.
    for p in primes:
        var ctx = Modulus(p, budget)
        for n in [1, 2, 3, 17, 100, 333, 1000]:
            var f = residues(n, p, p + UInt64(n), ctx)
            f[0] = ModInt(1, ctx)
            var a = series.asarray(f.copy(), order=n)
            var root = series.sqrt(a)
            assert_equal(root[0].value, UInt64(1))
            var square = values(series.multiply(root, root))
            var expected = List[UInt64]()
            for c in f:
                expected.append(c.value)
            assert_true(square == expected, String("square root squares back over GF(", p, ") at order ", n))
    # Products across the schoolbook, Montgomery and transform kernels, with
    # balanced and lopsided operands: prepared pairs, the copying path at a
    # lower order, and the schoolbook reference agree.
    for p in primes:
        var ctx = Modulus(p, budget)
        for shape in [(63, 63), (64, 64), (100, 100), (200, 200), (350, 350), (400, 400), (700, 700), (10, 600), (600, 40), (330, 900)]:
            var n = max(shape[0], shape[1])
            var x = residues(shape[0], p, 11, ctx)
            var y = residues(shape[1], p, 29, ctx)
            var expected = naive(x, y, n, p)
            var a = series.asarray(x.copy(), order=n)
            var b = series.asarray(y.copy(), order=n)
            var product = values(series.multiply(a, b))
            assert_true(product == expected, String("prepared product over GF(", p, ") ", shape[0], "x", shape[1]))
            var lower = values(series.multiply(a, b, order=n - 1))
            assert_true(lower == prefix(expected, n - 1), String("copied product over GF(", p, ") ", shape[0], "x", shape[1]))
            var copied = values(series.multiply(x.copy(), y.copy(), order=n))
            assert_true(copied == expected, String("coefficient-list product over GF(", p, ") ", shape[0], "x", shape[1]))
    # The other binary operations read prepared pairs as the copying path does.
    var p = UInt64(4611686018427387847)
    var ctx = Modulus(p, budget)
    for n in [5, 120, 500]:
        var f = residues(n, p, 3, ctx)
        var g = residues(n, p, 5, ctx)
        g[0] = ModInt(0, ctx)
        var h = residues(n, p, 7, ctx)
        h[0] = ModInt(3, ctx)
        var a = series.asarray(f.copy(), order=n)
        var inner = series.asarray(g.copy(), order=n)
        var unit = series.asarray(h.copy(), order=n)
        assert_true(prefix(values(series.add(a, unit)), n - 1) == values(series.add(a, unit, order=n - 1)), "prepared sum")
        assert_true(prefix(values(series.subtract(a, unit)), n - 1) == values(series.subtract(a, unit, order=n - 1)), "prepared difference")
        assert_true(prefix(values(series.divide(a, unit)), n - 1) == values(series.divide(a, unit, order=n - 1)), "prepared quotient")
        assert_true(prefix(values(series.compose(a, inner)), n - 1) == values(series.compose(a, inner, order=n - 1)), "prepared composition")
    # Prepared operands over different primes still fail as before.
    var rejected = False
    try:
        _ = series.multiply(series.asarray([ModInt(1, Modulus(7, budget))], order=1), series.asarray([ModInt(1, Modulus(11, budget))], order=1))
    except:
        rejected = True
    assert_true(rejected, "prepared operands over different primes")
    print("prime-field series kernel checks passed")
