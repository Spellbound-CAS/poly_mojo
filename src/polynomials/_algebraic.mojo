"""Private producers for exact real algebraic arithmetic.

Source: Z3 2d08bfcad30cabb7b2819dd2f4c47b7bbf4c377f,
src/math/polynomial/algebraic_numbers.cpp, mk_binary and interval selection.
Pure Mojo adaptation using this package's rational polynomials and budgets.
"""
# Z3
# Copyright (c) Microsoft Corporation
# All rights reserved.
# MIT License
#
# Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED *AS IS*, WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

from apn_mojo import Rational, Integer
from std.collections import Dict
from .errors import Budget, PolyError, PolyErrorKind, _allocation
from .coefficients import QQ, NoContext
from .dense import DensePoly
from .certificates import Factorization, ResultantResult, _dense_equal
from ._dense_ops import _copy as _copy_poly, _gcd, _exact
from ._root_ops import _compare as _compare_q, _copy as _copy_q, _midpoint, _sign_at
from .roots import RootInterval, isolate_real_roots
from .descartes import _unique_root, _one_root, _narrow_to_one, _witness_intervals, AlgebraicDescartesCertificate, _algebraic_descartes
from .evaluate import _derivative_dense
from .coefficients import ZZ
from ._images import _integer_dense
from ._integer_ops import _primitive
from ._univariate_modular import _univariate_gcd
from .factor import factor
from .algebraic_types import RealAlgebraic as _Real, AlgebraicComparisonCertificate as _Comparison, _copy_real, _rational
from .algebraic_arithmetic import AlgebraicArithmeticCertificate


struct _Frame(Copyable, Movable):
    """Every coordinate of a point as a polynomial in one generator.

    `representations[j]` is reduced modulo the generator's polynomial and
    names `variables[j]`; without a generator every coordinate is rational
    and its representation constant. `valid` is false once a coordinate
    could not be expressed.
    """
    var valid: Bool
    var generator: Optional[_Real]
    var variables: List[Int64]
    var representations: List[DensePoly[QQ]]

    def __init__(out self):
        self.valid = True
        self.generator = None
        self.variables = List[Int64]()
        self.representations = List[DensePoly[QQ]]()


struct _NormalizationCache(Movable):
    var polynomials: List[DensePoly[QQ]]
    var factorizations: List[Factorization[QQ]]
    # Arithmetic results are a function of the operation and the exact operand
    # witnesses, so a repeated step reuses its certificate unchanged.
    var operations: List[Int]
    var operands: List[Tuple[_Real, _Real]]
    var results: List[AlgebraicArithmeticCertificate]
    # The frame of the point being lifted, when its caller knows one.
    var frame: Optional[_Frame]
    # Real-root lists of exact norm polynomials: the complete isolation and
    # the per-factor candidates are functions of the polynomial alone, so a
    # norm met again reuses both. Keyed by _norm_key, chained on collisions.
    var norm_heads: Dict[Int, Int]
    var norm_next: List[Int]
    var norms: List[DensePoly[QQ]]
    var norm_isolations: List[AlgebraicDescartesCertificate]
    var norm_candidates: List[List[_Real]]
    # Tower resultants Res_x(F, m_e) are functions of the element's polynomial
    # and F's coefficient representations alone; keyed like the norms.
    var tower_heads: Dict[Int, Int]
    var tower_next: List[Int]
    var tower_elements: List[DensePoly[QQ]]
    var tower_representations: List[List[DensePoly[QQ]]]
    var tower_results: List[Optional[ResultantResult[QQ]]]

    def __init__(out self):
        self.frame = None
        self.norm_heads = Dict[Int, Int]()
        self.norm_next = List[Int]()
        self.norms = List[DensePoly[QQ]]()
        self.norm_isolations = List[AlgebraicDescartesCertificate]()
        self.norm_candidates = List[List[_Real]]()
        self.tower_heads = Dict[Int, Int]()
        self.tower_next = List[Int]()
        self.tower_elements = List[DensePoly[QQ]]()
        self.tower_representations = List[List[DensePoly[QQ]]]()
        self.tower_results = List[Optional[ResultantResult[QQ]]]()
        self.polynomials = List[DensePoly[QQ]]()
        self.factorizations = List[Factorization[QQ]]()
        self.operations = List[Int]()
        self.operands = List[Tuple[_Real, _Real]]()
        self.results = List[AlgebraicArithmeticCertificate]()


def _same_witness(a: _Real, b: _Real, mut budget: Budget) raises PolyError -> Bool:
    budget.charge(1)
    if a.rational != b.rational or a.polynomial.degree() != b.polynomial.degree():
        return False
    return a.low == b.low and a.high == b.high and _dense_equal(a.polynomial, b.polynomial, budget)


def _isolated(p: DensePoly[QQ], low: QQ, high: QQ, mut budget: Budget) raises PolyError -> _Real:
    # A squarefree p's single root in (low, high), narrowed until Descartes'
    # bound certifies it.
    var a = _copy_q(low, budget)
    var b = _copy_q(high, budget)
    if not _one_root(p, a, b, budget):
        _narrow_to_one(p, a, b, budget)
    return _Real(_copy_poly(p, budget), a^, b^, False)


def _monic(var f: DensePoly[QQ], mut budget: Budget) raises PolyError -> DensePoly[QQ]:
    if f.degree() < 0 or _compare_q(f.coefficients[f.degree()], QQ.one(NoContext()), budget) == 0:
        return f^
    var inverse = f.coefficients[f.degree()].inverse(budget)
    return f.mul(DensePoly[QQ].from_coefficients(f.variable, NoContext(), [inverse^], budget), budget)


def _qq_gcd(f: DensePoly[QQ], g: DensePoly[QQ], mut budget: Budget) raises PolyError -> DensePoly[QQ]:
    # The monic gcd, from the modular gcd of the primitive integer parts:
    # no coefficient swell, and the same canonical result as rational Euclid.
    if f.degree() == g.degree() and _dense_equal(f, g, budget):
        return _monic(f.copy(), budget)
    var h = _univariate_gcd(_primitive(_integer_dense(f, budget), budget), _primitive(_integer_dense(g, budget), budget), budget)
    _allocation[QQ](len(h.coefficients), budget)
    var values = List[QQ](capacity=len(h.coefficients))
    for c in h.coefficients:
        values.append(QQ(Rational(c.value)))
    return _monic(DensePoly[QQ].from_coefficients(f.variable, NoContext(), values^, budget), budget)


def _squarefree_part(var f: DensePoly[QQ], mut budget: Budget) raises PolyError -> DensePoly[QQ]:
    if f.degree() < 2:
        return f^
    var repeated = _qq_gcd(f, _derivative_dense(f, budget), budget)
    if repeated.degree() < 1:
        return f^
    return _exact(f, repeated, budget)


def _root(f: DensePoly[QQ], low: QQ, high: QQ, mut budget: Budget) raises PolyError -> _Real:
    var cache = _NormalizationCache()
    return _root_cached(f, low, high, cache, budget)


def _root_cached(f: DensePoly[QQ], low: QQ, high: QQ, mut cache: _NormalizationCache, mut budget: Budget) raises PolyError -> _Real:
    f.check_invariants(budget)
    if f.degree() < 1:
        raise PolyError(PolyErrorKind.DOMAIN, "an algebraic root needs a nonconstant polynomial")
    var ordering = _compare_q(low, high, budget)
    if ordering == 0:
        if not _sign_at(f, low, budget) == 0:
            raise PolyError(PolyErrorKind.DOMAIN, "the rational point is not a root")
        return _rational(low, budget)
    if ordering > 0 or not _unique_root(f, low, high, budget):
        raise PolyError(PolyErrorKind.DOMAIN, "an algebraic interval must isolate exactly one root")
    var normalized = _copy_poly(f, budget)
    normalized.variable = 0
    var index = -1
    for i in range(len(cache.polynomials)):
        budget.charge(1)
        if _dense_equal(cache.polynomials[i], normalized, budget):
            index = i
            break
    if index < 0:
        var factors = factor(normalized, budget)
        index = len(cache.polynomials)
        _allocation[DensePoly[QQ]](1, budget)
        _allocation[Factorization[QQ]](1, budget)
        cache.polynomials.append(normalized^)
        cache.factorizations.append(factors^)
    for entry in cache.factorizations[index].factors:
        budget.charge(1)
        if _unique_root(entry[0], low, high, budget):
            if entry[0].degree() == 1:
                return _rational(entry[0].coefficients[0].neg(budget).mul(entry[0].coefficients[1].inverse(budget), budget), budget)
            return _isolated(entry[0], low, high, budget)
    raise PolyError(PolyErrorKind.DOMAIN, "no factor contains the selected root")


def _refine(mut a: _Real, mut budget: Budget) raises PolyError:
    budget.charge(1)
    if a.rational:
        return
    var middle = _midpoint(a.low, a.high, budget)
    var at_middle = _sign_at(a.polynomial, middle, budget)
    if at_middle == 0:
        # Keep the defining polynomial: refinement witnesses must remain
        # bound to their input, even for a nonminimal caller-built value.
        a.low = _midpoint(a.low, middle, budget)
        a.high = _midpoint(middle, a.high, budget)
        return
    # A valid witness changes sign across its interval at its single simple
    # root, so the half whose endpoint signs differ holds it.
    var at_low = _sign_at(a.polynomial, a.low, budget)
    if at_low == 0 or at_low == _sign_at(a.polynomial, a.high, budget):
        raise PolyError(PolyErrorKind.DOMAIN, "an algebraic witness must change sign across its interval")
    if at_middle == at_low:
        a.low = middle^
    else:
        a.high = middle^


def _interval_relation(a: _Real, b: _Real, mut budget: Budget) raises PolyError -> Optional[Int]:
    if a.rational and b.rational:
        return _compare_q(a.low, b.low, budget)
    if _compare_q(a.high, b.low, budget) <= 0:
        return -1
    if _compare_q(b.high, a.low, budget) <= 0:
        return 1
    return None


def _resolve_overlap(mut left: _Real, mut right: _Real, mut common: DensePoly[QQ], mut shared: Optional[RootInterval], mut budget: Budget) raises PolyError -> Int:
    if left.rational or right.rational:
        if _sign_at(right.polynomial if left.rational else left.polynomial, left.low if left.rational else right.low, budget) == 0:
            return 0
    else:
        common = _squarefree_part(_qq_gcd(left.polynomial, right.polynomial, budget), budget)
        if common.degree() > 0:
            # Both values are the common factor's root when it has one in the
            # overlap, whose ends are nonroots of one definition each.
            var low = _copy_q(left.low if _compare_q(left.low, right.low, budget) >= 0 else right.low, budget)
            var high = _copy_q(left.high if _compare_q(left.high, right.high, budget) <= 0 else right.high, budget)
            if _compare_q(low, high, budget) < 0 and _unique_root(common, low, high, budget):
                _narrow_to_one(common, low, high, budget)
                shared = RootInterval(low^, high^)
                return 0
    # Equality has been excluded on the original overlap. Refine once and keep
    # these endpoints for the ordering certificate instead of replaying the work.
    while True:
        budget.charge(1)
        _refine(left, budget)
        _refine(right, budget)
        var relation = _interval_relation(left, right, budget)
        if relation:
            return relation.value()


def _compare(a: _Real, b: _Real, mut budget: Budget) raises PolyError -> Int:
    var relation = _interval_relation(a, b, budget)
    if relation:
        return relation.value()
    var left = _copy_real(a, budget)
    var right = _copy_real(b, budget)
    var common = DensePoly[QQ](0, NoContext())
    var shared: Optional[RootInterval] = None
    return _resolve_overlap(left, right, common, shared, budget)


def _comparison(a: _Real, b: _Real, mut budget: Budget) raises PolyError -> _Comparison:
    var left = _copy_real(a, budget)
    var right = _copy_real(b, budget)
    var common = DensePoly[QQ](0, NoContext())
    var shared: Optional[RootInterval] = None
    var separated = _interval_relation(left, right, budget)
    var relation = separated.value() if separated else _resolve_overlap(left, right, common, shared, budget)
    if relation == 0 and not a.rational and not b.rational:
        var left_cofactor = _exact(a.polynomial, common, budget)
        var right_cofactor = _exact(b.polynomial, common, budget)
        return _Comparison(0, left^, right^, common^, left_cofactor^, right_cofactor^, shared^)
    return _Comparison(relation, left^, right^, DensePoly[QQ](0, NoContext()), DensePoly[QQ](0, NoContext()), DensePoly[QQ](0, NoContext()), None)


def _sign(a: _Real, mut budget: Budget) raises PolyError -> Int:
    if a.rational:
        budget.charge(1)
        return a.low.value.sign()
    var value = _copy_real(a, budget)
    while True:
        budget.charge(1)
        if value.low.value.sign() >= 0:
            return 1
        if value.high.value.sign() <= 0:
            return -1
        if value.polynomial.coefficients[0].is_zero():
            return 0
        _refine(value, budget)
        if value.rational:
            return value.low.value.sign()


def _roots(f: DensePoly[QQ], mut budget: Budget) raises PolyError -> List[_Real]:
    var cache = _NormalizationCache()
    return _roots_cached(f, cache, budget)


def _roots_cached(f: DensePoly[QQ], mut cache: _NormalizationCache, mut budget: Budget) raises PolyError -> List[_Real]:
    var isolation = isolate_real_roots(f, budget)
    var roots = List[_Real]()
    for q in isolation.exact_roots:
        budget.charge(1)
        _allocation[_Real](1, budget)
        roots.append(_rational(q, budget))
    for interval in isolation.intervals:
        budget.charge(1)
        _allocation[_Real](1, budget)
        roots.append(_root_cached(f, interval.low, interval.high, cache, budget))
    for i in range(1, len(roots)):
        budget.charge(1)
        var j = i
        while j > 0 and _compare(roots[j], roots[j - 1], budget) < 0:
            budget.charge(1)
            var moved = _copy_real(roots[j], budget)
            roots[j] = _copy_real(roots[j - 1], budget)
            roots[j - 1] = moved^
            j -= 1
    return roots^


def _next_group(mut cache: _NormalizationCache):
    """Start a new sibling group: drop everything but the norm root lists and tower resultants, which depend on their polynomials alone."""
    cache.polynomials = List[DensePoly[QQ]]()
    cache.factorizations = List[Factorization[QQ]]()
    cache.operations = List[Int]()
    cache.operands = List[Tuple[_Real, _Real]]()
    cache.results = List[AlgebraicArithmeticCertificate]()
    cache.frame = None


def _mix_key(key: UInt64, h: DensePoly[QQ]) raises -> UInt64:
    # Degree and every coefficient's numerator and denominator modulo 2^61 - 1.
    var prime = Integer((1 << 61) - 1)
    var mixed = key * 1000003 ^ UInt64(h.degree() + 1)
    for c in h.coefficients:
        mixed = mixed * 1000003 ^ UInt64(Int(c.value.numerator() % prime))
        mixed = mixed * 1000003 ^ UInt64(Int(c.value.denominator() % prime))
    return mixed


def _norm_key(h: DensePoly[QQ]) raises -> Int:
    return Int(_mix_key(0, h) >> 1)


def _tower_find(m: DensePoly[QQ], representations: List[DensePoly[QQ]], cache: _NormalizationCache, mut budget: Budget) raises PolyError -> Tuple[Int, Int]:
    """The key of a tower and the index of its stored resultant, or -1."""
    var key: Int
    try:
        var mixed = _mix_key(UInt64(len(representations)), m)
        for r in representations:
            mixed = _mix_key(mixed, r)
        key = Int(mixed >> 1)
    except e:
        raise PolyError(PolyErrorKind.CAPACITY, String(e))
    budget.charge(UInt64(len(m.coefficients) + len(representations)))
    var index = cache.tower_heads.get(key, -1)
    while index >= 0:
        budget.charge(1)
        if len(cache.tower_representations[index]) == len(representations) and _dense_equal(cache.tower_elements[index], m, budget):
            var same = True
            for i in range(len(representations)):
                if not _dense_equal(cache.tower_representations[index][i], representations[i], budget):
                    same = False
                    break
            if same:
                return (key, index)
        index = cache.tower_next[index]
    return (key, -1)


def _tower_store(key: Int, m: DensePoly[QQ], representations: List[DensePoly[QQ]], result: Optional[ResultantResult[QQ]], mut cache: _NormalizationCache, mut budget: Budget) raises PolyError:
    _allocation[DensePoly[QQ]](1 + len(representations), budget)
    _allocation[ResultantResult[QQ]](1, budget)
    cache.tower_next.append(cache.tower_heads.get(key, -1))
    cache.tower_heads[key] = len(cache.tower_elements)
    cache.tower_elements.append(_copy_poly(m, budget))
    cache.tower_representations.append(representations.copy())
    cache.tower_results.append(result.copy())


def _norm_roots(h: DensePoly[QQ], mut cache: _NormalizationCache, mut budget: Budget) raises PolyError -> Tuple[AlgebraicDescartesCertificate, List[_Real]]:
    """The complete isolation of h's real roots and its per-factor candidates, computed once per cache."""
    var key: Int
    try:
        key = _norm_key(h)
    except e:
        raise PolyError(PolyErrorKind.CAPACITY, String(e))
    budget.charge(UInt64(len(h.coefficients)))
    var index = cache.norm_heads.get(key, -1)
    while index >= 0:
        budget.charge(1)
        if _dense_equal(cache.norms[index], h, budget):
            _allocation[_Real](len(cache.norm_candidates[index]), budget)
            return (cache.norm_isolations[index].copy(), cache.norm_candidates[index].copy())
        index = cache.norm_next[index]
    var isolation = _algebraic_descartes(h, budget)
    var candidates = _roots_by_factors(h, cache, budget)
    _allocation[DensePoly[QQ]](1, budget)
    _allocation[AlgebraicDescartesCertificate](1, budget)
    _allocation[_Real](len(candidates), budget)
    cache.norm_next.append(cache.norm_heads.get(key, -1))
    cache.norm_heads[key] = len(cache.norms)
    cache.norms.append(_copy_poly(h, budget))
    cache.norm_isolations.append(isolation.copy())
    cache.norm_candidates.append(candidates.copy())
    return (isolation^, candidates^)


def _roots_by_factors(f: DensePoly[QQ], mut cache: _NormalizationCache, mut budget: Budget) raises PolyError -> List[_Real]:
    # Every distinct real root of f, isolated on the irreducible factor that
    # defines it: smaller polynomials than isolating f, the same roots.
    var normalized = _copy_poly(f, budget)
    normalized.variable = 0
    var index = -1
    for i in range(len(cache.polynomials)):
        budget.charge(1)
        if _dense_equal(cache.polynomials[i], normalized, budget):
            index = i
            break
    if index < 0:
        var factors = factor(normalized, budget)
        index = len(cache.polynomials)
        _allocation[DensePoly[QQ]](1, budget)
        _allocation[Factorization[QQ]](1, budget)
        cache.polynomials.append(normalized^)
        cache.factorizations.append(factors^)
    var roots = List[_Real]()
    for entry in cache.factorizations[index].factors:
        budget.charge(1)
        ref p = entry[0]
        if p.degree() == 1:
            _allocation[_Real](1, budget)
            roots.append(_rational(p.coefficients[0].neg(budget).mul(p.coefficients[1].inverse(budget), budget), budget))
            continue
        # An irreducible factor of degree two or more has no rational root;
        # its Descartes brackets already certify each root.
        for interval in _witness_intervals(p, budget):
            budget.charge(1)
            _allocation[_Real](1, budget)
            roots.append(_Real(_copy_poly(p, budget), interval.low.copy(), interval.high.copy(), False))
    for i in range(1, len(roots)):
        budget.charge(1)
        var j = i
        while j > 0 and _compare(roots[j], roots[j - 1], budget) < 0:
            budget.charge(1)
            var moved = _copy_real(roots[j], budget)
            roots[j] = _copy_real(roots[j - 1], budget)
            roots[j - 1] = moved^
            j -= 1
    # Witnesses come from separate factors; bisect neighbours until each
    # interval ends before the next begins, so all are disjoint. Two intervals
    # may share an endpoint, which is a root of neither, but none may reach a
    # rational root.
    for i in range(1, len(roots)):
        var strict = 0 if roots[i - 1].rational or roots[i].rational else 1
        while _compare_q(roots[i - 1].high, roots[i].low, budget) >= strict:
            budget.charge(1)
            _refine(roots[i - 1], budget)
            _refine(roots[i], budget)
    return roots^
