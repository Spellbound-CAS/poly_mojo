"""Verifier-owned rational Sturm data and bounded batches across CAD stacks.

Callers first establish squarefree coprime bases and disjoint ordered root
brackets. Checking every owned bracket and the total root count for each basis
member then proves completeness; roots of another member already have their
own disjoint brackets. A queued check is not acceptance: every caller must
flush before returning an accepted result. No producer data populates this cache.
"""
from .coefficients import QQ, NoContext
from .dense import DensePoly
from .errors import Budget, PolyError, PolyErrorKind, _allocation, _sum, _product
from .certificates import _dense_equal
from ._algebraic import _norm_key
from ._cad_final import _FinalRoot, _sequence, _total
from ._cad_final_signs import _SignStats, _PreparedSigns, _EndpointBatch, _prepare_signs, _pack_endpoints, _prepared_signs_at
from ._root_ops import _compare
from std.collections import Dict


@fieldwise_init
struct _RootCountEntry(Movable):
    var polynomial: DensePoly[QQ]
    var sequence: List[_PreparedSigns]
    var total: Int
    var endpoints: List[QQ]


struct _RootCountCache(Movable):
    """Bounded verified computations, exact collision checks, no cached answers."""
    var heads: Dict[Int, Int]
    var next: List[Int]
    var entries: List[_RootCountEntry]
    var capacity: Int
    var batch_points: Int
    var key_mask: Int
    var hits: Int
    var misses: Int
    var linear_checks: Int
    var flushes: Int
    var stats: _SignStats
    var valid: Bool

    def __init__(out self, capacity: Int = 128, batch_points: Int = 64, key_mask: Int = Int.MAX):
        self.heads = Dict[Int, Int]()
        self.next = List[Int]()
        self.entries = List[_RootCountEntry]()
        self.capacity = max(0, capacity)
        self.batch_points = max(2, batch_points)
        self.key_mask = key_mask
        self.hits = 0
        self.misses = 0
        self.linear_checks = 0
        self.flushes = 0
        self.stats = _SignStats()
        self.valid = True

    def flush_entry(mut self, index: Int, mut budget: Budget, vectorized: Bool) raises PolyError -> Bool:
        if not self.valid:
            return False
        ref entry = self.entries[index]
        ref points = entry.endpoints
        if len(points) == 0:
            return True
        self.flushes = _sum(self.flushes, 1)
        var batches = List[_EndpointBatch]()
        if vectorized and len(points) >= 8:
            batches = _pack_endpoints(points, budget)
        _allocation[Int](_product(3, len(points)), budget)
        var variations = List[Int](length=len(points), fill=0)
        var previous = List[Int](length=len(points), fill=0)
        var boundary = List[Int]()
        for i in range(len(entry.sequence)):
            var signs = _prepared_signs_at(entry.sequence[i], points, batches, self.stats, budget, vectorized=vectorized)
            if i == 0:
                boundary = signs.copy()
            for j in range(len(points)):
                budget.charge(1)
                if signs[j]:
                    variations[j] += Int(previous[j] != 0 and previous[j] != signs[j])
                    previous[j] = signs[j]
        for j in range(0, len(points), 2):
            budget.charge(1)
            if points[j] == points[j + 1]:
                if boundary[j] != 0:
                    self.valid = False
                    return False
            elif boundary[j] == 0 or boundary[j + 1] == 0 or variations[j] - variations[j + 1] != 1:
                self.valid = False
                return False
        entry.endpoints.clear()
        return True

    def flush(mut self, mut budget: Budget, vectorized: Bool) raises PolyError -> Bool:
        for i in range(len(self.entries)):
            if not self.flush_entry(i, budget, vectorized):
                return False
        return self.valid

    def find(mut self, p: DensePoly[QQ], mut budget: Budget, vectorized: Bool) raises PolyError -> Int:
        var key: Int
        try:
            key = _norm_key(p) & self.key_mask
        except e:
            raise PolyError(PolyErrorKind.CAPACITY, String(e))
        var index = self.heads.get(key, -1) if self.capacity else -1
        while index >= 0:
            budget.charge(1)
            if _dense_equal(p, self.entries[index].polynomial, budget):
                self.hits = _sum(self.hits, 1)
                return index
            index = self.next[index]
        if len(self.entries) >= max(1, self.capacity):
            if not self.flush(budget, vectorized):
                return -1
            self.heads.clear()
            self.next.clear()
            self.entries.clear()
        var sequence = _sequence(p, budget)
        var total = _total(sequence, budget)
        _allocation[_PreparedSigns](len(sequence), budget)
        var prepared = List[_PreparedSigns]()
        for f in sequence:
            prepared.append(_prepare_signs(f, budget))
        _allocation[_RootCountEntry](1, budget)
        _allocation[Int](1, budget)
        index = len(self.entries)
        self.next.append(self.heads.get(key, -1))
        self.heads[key] = index
        self.entries.append(_RootCountEntry(p.copy(), prepared^, total, List[QQ]()))
        self.misses = _sum(self.misses, 1)
        return index


def _check_rational_roots(basis: List[List[DensePoly[QQ]]], roots: List[_FinalRoot], mut cache: _RootCountCache, mut budget: Budget, vectorized: Bool) raises PolyError -> Bool:
    """Queue owned-root checks after the caller establishes algebra and geometry."""
    for k in range(len(basis)):
        ref member = basis[k]
        var count = 0
        var first = -1
        for i in range(len(roots)):
            budget.charge(1)
            if roots[i].basis == k:
                count += 1
                first = i
        if len(member) == 2:
            # The caller proved monicity. A linear polynomial has one real
            # root, so one exact bracket comparison establishes completeness.
            if count != 1:
                return False
            var zero = QQ.zero(NoContext())
            if not member[0].is_zero():
                zero = member[0].coefficients[0].neg(budget)
            ref root = roots[first]
            if root.low == root.high:
                if zero != root.low:
                    return False
            elif _compare(root.low, zero, budget) >= 0 or _compare(zero, root.high, budget) >= 0:
                return False
            cache.linear_checks = _sum(cache.linear_checks, 1)
            continue
        _allocation[QQ](len(member), budget)
        var coefficients = List[QQ]()
        for c in member:
            coefficients.append(c.coefficients[0].copy() if not c.is_zero() else QQ.zero(NoContext()))
        var p = DensePoly[QQ].from_coefficients(0, NoContext(), coefficients^, budget)
        var index = cache.find(p, budget, vectorized)
        if index < 0 or count != cache.entries[index].total:
            return False
        for root in roots:
            if root.basis != k:
                continue
            _allocation[QQ](2, budget)
            cache.entries[index].endpoints.append(root.low.copy())
            cache.entries[index].endpoints.append(root.high.copy())
            if len(cache.entries[index].endpoints) >= cache.batch_points:
                if not cache.flush_entry(index, budget, vectorized):
                    return False
        if cache.capacity == 0 and not cache.flush_entry(index, budget, vectorized):
            return False
    return cache.valid
