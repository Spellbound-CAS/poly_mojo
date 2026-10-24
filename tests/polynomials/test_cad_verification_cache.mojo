"""Complete verification must flush cached root jobs and reject corrupt data."""
from apn_mojo import Integer, Rational
from polynomials import Budget, QQ, DensePoly, NoContext, SparsePoly, MonomialOrder, PolyErrorKind
from polynomials._algebraic_filter import _constant
from polynomials._cad_final import _FinalRoot
from polynomials._cad_root_cache import _RootCountCache, _check_rational_roots
from polynomials._cad_tree_final import _construct_final_cad, _check_final_cad, _check_final_cad_with_cache
from polynomials._cad_final_signs import _SignStats, _prepare_signs, _pack_endpoints, _prepared_signs_at
from polynomials._root_ops import _sign_at
from support.check import Checker
from support.corpus import parse_polynomial


def member(values: List[Int], mut work: Budget) raises -> List[DensePoly[QQ]]:
    var result = List[DensePoly[QQ]]()
    for v in values:
        result.append(_constant(QQ(Rational(v)), work))
    return result^


def brackets(n: Int) raises -> List[_FinalRoot]:
    return [_FinalRoot(0, QQ(Rational(-n - 1)), QQ(Rational(-1, 2))),
            _FinalRoot(0, QQ(Rational(1, 2)), QQ(Rational(n + 1)))]


def main() raises:
    var c = Checker("CAD verification cache")
    var work = Budget()
    var basis: List[List[DensePoly[QQ]]] = [member([-2, 0, 1], work)]
    var roots = brackets(2)
    for vector in range(2):
        for capacity in [0, 1, 128]:
            for batch in [2, 8, 64]:
                var cache = _RootCountCache(capacity, batch)
                for _ in range(33):
                    c.check(_check_rational_roots(basis, roots, cache, work, Bool(vector)), "valid roots queue across repeated stacks")
                c.check(cache.flush(work, Bool(vector)), "all queued roots pass exact scalar/SIMD verification")
                c.check(cache.hits == (32 if capacity else 0), "only exact repeated polynomials reuse Sturm data")
                c.check(cache.misses == (1 if capacity else 33), "cache-disabled reference recomputes each sequence")
                if vector and capacity and batch >= 8:
                    c.check(cache.stats.vector_batches > 0, "compatible roots from multiple stacks fill SIMD batches")

        var cache = _RootCountCache(2, 64, 0)
        c.check(_check_rational_roots(basis, roots, cache, work, Bool(vector)), "first colliding key")
        var other: List[List[DensePoly[QQ]]] = [member([-3, 0, 1], work)]
        c.check(_check_rational_roots(other, brackets(3), cache, work, Bool(vector)), "different polynomial with same hash remains distinct")
        c.check(_check_rational_roots(basis, roots, cache, work, Bool(vector)), "exact comparison resolves a collision")
        c.check(cache.misses == 2 and cache.hits == 1 and cache.flush(work, Bool(vector)), "collisions cannot reuse incorrect root counts")

        cache = _RootCountCache(1, 64)
        var bad = roots.copy()
        bad[1].low = QQ(Rational(2))
        bad[1].high = QQ(Rational(3))
        c.check(_check_rational_roots(basis, bad, cache, work, Bool(vector)), "a queued bracket is not yet accepted")
        c.check(not cache.flush(work, Bool(vector)), "the final flush rejects a missed root")
        c.check(not cache.flush(work, Bool(vector)), "a rejected cache cannot later accept")
        cache = _RootCountCache(1, 64)
        _ = _check_rational_roots(basis, bad, cache, work, Bool(vector))
        c.check(not _check_rational_roots(other, brackets(3), cache, work, Bool(vector)), "eviction checks pending jobs before dropping them")

        cache = _RootCountCache()
        bad = roots.copy()
        _ = bad.pop()
        c.check(not _check_rational_roots(basis, bad, cache, work, Bool(vector)), "missing roots fail per-stack completeness before batching")
        cache = _RootCountCache()
        other = [member([1, 0, 1], work)]
        c.check(_check_rational_roots(other, List[_FinalRoot](), cache, work, Bool(vector)) and cache.flush(work, Bool(vector)), "a zero-real-root polynomial still has its total checked")
        other = [member([-1, 0, 1], work)]
        bad = [_FinalRoot(0, QQ(Rational(-1)), QQ(Rational(-1))), _FinalRoot(0, QQ(Rational(1)), QQ(Rational(1)))]
        c.check(_check_rational_roots(other, bad, cache, work, Bool(vector)) and cache.flush(work, Bool(vector)), "exact rational roots keep exact zero semantics")

    # Packing is reused for a whole sequence; mixed wide/native endpoints,
    # overflow, and every partial-vector length retain an exact oracle.
    var polynomial = DensePoly[QQ].from_coefficients(0, NoContext(), [QQ(Rational(-2, 7)), QQ(Rational(3, 11)), QQ(Rational(1))], work)
    var prepared = _prepare_signs(polynomial, work)
    for count in range(8, 25):
        var points = List[QQ]()
        for j in range(count):
            points.append(QQ(Rational(j - 12, j % 3 + 1)))
        points[count - 1] = QQ(Rational(Integer(1) << 90))
        var packed = _pack_endpoints(points, work)
        var stats = _SignStats()
        var signs = _prepared_signs_at(prepared, points, packed, stats, work)
        for j in range(count):
            c.check(signs[j] == _sign_at(polynomial, points[j], work), "cached packing and wide fallback match exact evaluation")

    var inputs: List[SparsePoly[QQ]] = [parse_polynomial("x0^4 - 10*x0^2 + 9", MonomialOrder(), work), parse_polynomial("x1^2 - 2", MonomialOrder(), work)]
    var variables: List[Int64] = [0, 1]
    var masks: List[Int] = [7, 2]
    var result = _construct_final_cad(inputs, variables, masks, work, interval=True)
    for vector in range(2):
        var cache = _RootCountCache()
        c.check(_check_final_cad_with_cache(inputs, variables, masks, result, cache, work, vectorized=Bool(vector)), "complete tree passes with a fresh verifier cache")
        c.check(cache.hits > 0, "full CAD reuses a repeated terminal polynomial")
        if vector:
            c.check(cache.stats.vector_batches > 0, "complete verification batches across stack boundaries")
        var bad = result.copy()
        ref last = bad.stacks[len(bad.stacks) - 1].stack
        last.roots[len(last.roots) - 1].low = QQ(Rational(1))
        last.roots[len(last.roots) - 1].high = QQ(Rational(1))
        c.check(not _check_final_cad(inputs, variables, masks, bad, work, vectorized=Bool(vector)), "complete gate rejects a late forged exact root")

    var cache = _RootCountCache()
    var limited = Budget(0)
    var exhausted = False
    try:
        _ = _check_rational_roots(basis, roots, cache, limited, True)
    except e:
        exhausted = e.kind == PolyErrorKind.LIMIT
    c.check(exhausted, "budget exhaustion cannot accept queued verification")
    c.finish()
