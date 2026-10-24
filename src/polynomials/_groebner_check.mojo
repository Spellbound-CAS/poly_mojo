"""Checking a proposed Gröbner basis on packed terms: the inputs and the S-pairs reduce to zero.

`sparse.check_groebner` proves both ideal containments and Buchberger's
criterion. Here the basis is packed once and moved into a packed engine,
which installs its pairs under the Gebauer–Möller criteria, as the default
mode's certificate does (`_rational_groebner._certify`):

- every input reduces to zero modulo the basis, one at a time, stopping at
  the first term no leading monomial divides;
- under a degree order, the S-pairs of each LCM degree share one F4 matrix:
  both halves of every pair and, from symbolic preprocessing, a reducer for
  every monomial some leading monomial divides. Over a prime field each
  S-polynomial is reduced alone in a dense row; over QQ the basis is scaled
  to integer polynomials and `_integer_f4._pairs_reduce_to_zero` reduces
  them fraction-free. Under LEX each S-polynomial is reduced by division.

Reducing to zero by some choice of reducers gives a standard representation,
which is all Buchberger's criterion needs. A basis with a zero element, or
with a leading monomial that divides another's, is left to the generic
check, as are fields and orders the packed engine does not take.

A basis over QQ known to lie in the ideal (`_proves_reduced`, for bases a
trace replay computed) is first offered to `_etale`: for n generators in n
variables whose image modulo a word prime is a reduced algebra of the
basis's dimension, Hensel's lemma bounds the ideal's colength from below,
and that with the basis inside the ideal makes it a Gröbner basis with no
S-pair reduced. Any unmet condition leaves the S-pairs to decide.
`check_groebner_basis` always uses the S-pairs.
"""
from std.bit import count_trailing_zeros, pop_count
from std.collections import Dict
from std.math import iota
from .coefficients import CoefficientField, CoefficientRing, ModInt, Modulus, NoContext, QQ, ZZ
from ._coefficient_ops import _exact, _gcd_coefficient
from .errors import Budget, PolyError, _allocation, _product, _sum
from .monomial import MonomialOrder
from .sparse import SparsePoly
from ._packed_groebner import _Buckets, _Engine, _Packed, _Pair, _Term, _compare, _div, _divides, _pack, _packed, _packed_supported
from ._f4 import _F4_PRIME_LIMIT, _Columns, _Matrix, _Rows, _Store, _f4_packed, _raw_hash, _reduce_row
from ._integer_f4 import _pairs_reduce_to_zero
from ._rational_groebner import _integral, _next_prime
from ._qq_trace import _admissible, _images


def _check_basis[C: CoefficientField](generators: List[SparsePoly[C]], basis: List[SparsePoly[C]], order: MonomialOrder, mut budget: Budget) raises PolyError -> Optional[Bool]:
    """Whether every generator and every S-pair reduces to zero modulo the basis, or None where the packed check does not apply.

    The basis must use `order`; the generators are reordered to it.
    """
    if not _packed_supported(order) or len(basis) == 0:
        return None
    for g in basis:
        if g.is_zero() or g.ring != basis[0].ring or g.order != order:
            return None
    for g in generators:
        if g.ring != basis[0].ring:
            return None
    comptime if C == ModInt:
        ref modular = rebind[List[SparsePoly[ModInt]]](basis)
        ref inputs = rebind[List[SparsePoly[ModInt]]](generators)
        if len(order.variables) <= 8:
            return _check_modular[8](inputs, modular, order, budget)
        return _check_modular[16](inputs, modular, order, budget)
    comptime if C == QQ:
        ref rational = rebind[List[SparsePoly[QQ]]](basis)
        ref inputs = rebind[List[SparsePoly[QQ]]](generators)
        if len(order.variables) <= 8:
            return _check_rational[8](inputs, rational, order, budget)
        return _check_rational[16](inputs, rational, order, budget)
    return None


def _check_modular[W: Int](generators: List[SparsePoly[ModInt]], basis: List[SparsePoly[ModInt]], order: MonomialOrder, mut budget: Budget) raises PolyError -> Optional[Bool]:
    var ring = basis[0].ring.copy()
    var engine = _Engine[ModInt, W, False](order.kind, ring)
    for g in basis:
        budget.charge(1)
        engine.insert(_pack[ModInt, W](g, order, budget), budget, heap=False)
    if len(engine.active) != len(basis):
        return None
    engine.heapify_pairs(budget)
    engine.sort_active(budget)
    for f in generators:
        budget.charge(1)
        if f.is_zero():
            continue
        var pending = _Buckets[ModInt, W, False](order.kind, ring)
        pending.add(_pack[ModInt, W](f.reorder(order, budget), order, budget), budget)
        if not engine.reduces_to_zero(pending^, engine.active, budget):
            return False
    if order.kind == MonomialOrder.LEX or ring.value >= _F4_PRIME_LIMIT:
        return _pairs_by_division(engine, budget)
    return _pairs_modular(engine, ring.value, budget)


def _check_rational[W: Int](generators: List[SparsePoly[QQ]], basis: List[SparsePoly[QQ]], order: MonomialOrder, mut budget: Budget) raises PolyError -> Optional[Bool]:
    var engine = _Engine[ZZ, W, True](order.kind, NoContext())
    for g in basis:
        budget.charge(1)
        engine.insert(_integral[W](_pack[QQ, W](g, order, budget), budget), budget, heap=False)
    if len(engine.active) != len(basis):
        return None
    engine.heapify_pairs(budget)
    engine.sort_active(budget)
    for f in generators:
        budget.charge(1)
        if f.is_zero():
            continue
        var pending = _Buckets[ZZ, W, True](order.kind, NoContext())
        pending.add(_integral[W](_pack[QQ, W](f.reorder(order, budget), order, budget), budget), budget)
        if not engine.reduces_to_zero(pending^, engine.active, budget):
            return False
    if order.kind == MonomialOrder.LEX:
        return _pairs_by_division(engine, budget)
    return _pairs_reduce_to_zero(engine, budget)


def _pairs_by_division[C: CoefficientRing, W: Int, integral: Bool](mut engine: _Engine[C, W, integral], mut budget: Budget) raises PolyError -> Bool:
    """Whether every pending pair's S-polynomial reduces to zero, each by division."""
    while len(engine.pairs):
        var pair = engine.select(budget)
        if not engine.reduces_to_zero(_s_polynomial(engine, pair, budget), engine.active, budget):
            return False
    return True


def _s_polynomial[C: CoefficientRing, W: Int, integral: Bool](engine: _Engine[C, W, integral], pair: _Pair[W], mut budget: Budget) raises PolyError -> _Buckets[C, W, integral]:
    """The S-polynomial of a pair, from its two translated tails, as `_Engine.reduce_pair` builds it."""
    var left: C
    var right: C
    comptime if integral:
        ref a = engine.polys[pair.first][len(engine.polys[pair.first]) - 1].coefficient
        ref b = engine.polys[pair.second][len(engine.polys[pair.second]) - 1].coefficient
        var common = _gcd_coefficient(a, b, budget)
        left = _exact(b, common, budget)
        right = _exact(a, common, budget).neg(budget)
    else:
        left = C.one(engine.ring)
        right = left.neg(budget)
    var pending = _Buckets[C, W, integral](engine.kind, engine.ring)
    comptime if C == ModInt:
        pending.add_multiple(engine.polys[pair.first], _div(pair.lcm, engine.leads[pair.first]), rebind[ModInt](left), budget)
        pending.add_multiple(engine.polys[pair.second], _div(pair.lcm, engine.leads[pair.second]), rebind[ModInt](right), budget)
    else:
        pending.add(engine._multiple(pair.first, _div(pair.lcm, engine.leads[pair.first]), left, budget), budget)
        pending.add(engine._multiple(pair.second, _div(pair.lcm, engine.leads[pair.second]), right, budget), budget)
    return pending^


def _pairs_modular[W: Int](mut engine: _Engine[ModInt, W, False], prime: UInt64, mut budget: Budget) raises PolyError -> Bool:
    """Whether every pending pair's S-polynomial reduces to zero, by one F4 matrix per LCM degree.

    The active elements are monic and sorted. Each S-polynomial is the
    difference of its two halves, rows that lead in one column; it is reduced
    by a row leading in each later column it reaches, and any entry left in a
    column no row leads in is a nonzero remainder, since symbolic
    preprocessing gave a row to every column a leading monomial divides.
    """
    var store = _Store[W]()
    store.extend(engine, budget)
    var found = _Rows()
    while len(engine.pairs):
        budget.charge(1)
        var degree = engine.pairs[0].lcm.degree
        var matrix = _Matrix[W]()
        var halves = List[Tuple[Int, Int]]()
        var seen = Dict[Int, Int]()
        var stride = len(engine.polys)
        while len(engine.pairs) and engine.pairs[0].lcm.degree == degree:
            var pair = engine.select(budget)
            matrix.columns.reserve(1, budget)
            var lead = Int(matrix.columns.find(pair.lcm, _raw_hash(pair.lcm)))
            var rows = List[Int](capacity=2)
            for source in [pair.first, pair.second]:
                budget.charge(1)
                var key = lead * stride + source
                var known = seen.get(key)
                if known:
                    rows.append(known.value())
                    continue
                seen[key] = len(matrix.sources)
                rows.append(len(matrix.sources))
                matrix.add_row(engine, store, source, _div(pair.lcm, engine.leads[source]), budget)
            halves.append((rows[0], rows[1]))
        matrix.preprocess(engine, store, budget)
        _ = matrix.number(engine.kind, budget)
        var ncols = len(matrix.columns.monomials)
        var nrows = len(matrix.sources)
        _allocation[Int](ncols, budget)
        var pivot = List[Int](length=ncols, fill=-1)
        for r in range(nrows - 1, -1, -1):
            pivot[Int(matrix.entries[matrix.starts[r]])] = r
        _allocation[UInt64](ncols, budget)
        var dense = List[UInt64](length=ncols, fill=0)
        for half in halves:
            var a = half[0]
            var b = half[1]
            # Both halves are monic and lead in the same column: their
            # difference starts after it.
            var sa = matrix.starts[a]
            var sb = matrix.starts[b]
            ref va = store.values[matrix.sources[a]]
            ref vb = store.values[matrix.sources[b]]
            budget.charge(UInt64(matrix.starts[a + 1] - sa + matrix.starts[b + 1] - sb))
            for k in range(1, matrix.starts[a + 1] - sa):
                dense[Int(matrix.entries[sa + k])] = UInt64(va[k])
            for k in range(1, matrix.starts[b + 1] - sb):
                var j = Int(matrix.entries[sb + k])
                dense[j] += prime - UInt64(vb[k])
            _reduce_row(dense, Int(matrix.entries[sa]) + 1, pivot, matrix, store.values, found, prime, budget)
            if len(found.columns):
                return False
    return True



def _proves_reduced(generators: List[SparsePoly[QQ]], basis: List[SparsePoly[QQ]], order: MonomialOrder, mut budget: Budget, etale_start: UInt64 = _ETALE_START) raises PolyError -> Bool:
    """Whether a basis inside the generators' ideal is its reduced Gröbner basis, by the checks of `_check_basis`.

    The basis must lie in the ideal, which the caller has proved. Its
    elements must be monic, with no term but the leading one divisible by a
    leading monomial, in decreasing leading-monomial order, as `groebner`
    returns them. Then either `_etale` shows the basis is a Gröbner basis
    (see there), or every generator reducing to zero gives equal ideals and
    every S-pair reducing to zero a Gröbner basis; reduced, it is the
    unique reduced basis. `etale_start` is `_etale`'s, for tests.
    """
    if len(basis) == 0:
        for g in generators:
            if not g.is_zero():
                return False
        return True
    if not _packed_supported(order):
        return False
    if len(order.variables) <= 8:
        return _proves_reduced_packed[8](generators, basis, order, budget, etale_start)
    return _proves_reduced_packed[16](generators, basis, order, budget, etale_start)


def _proves_reduced_packed[W: Int](generators: List[SparsePoly[QQ]], basis: List[SparsePoly[QQ]], order: MonomialOrder, mut budget: Budget, etale_start: UInt64) raises PolyError -> Bool:
    _allocation[_Packed[W]](len(basis), budget)
    var leads = List[_Packed[W]](capacity=len(basis))
    var packed = List[List[_Term[QQ, W]]](capacity=len(basis))
    for g in basis:
        budget.charge(1)
        if g.is_zero() or g.order != order:
            return False
        var p = _pack[QQ, W](g, order, budget)
        if not (p[len(p) - 1].coefficient == QQ.one(NoContext())):
            return False
        leads.append(p[len(p) - 1].monomial)
        packed.append(p^)
    for i in range(1, len(leads)):
        if _compare(order.kind, leads[i - 1], leads[i]) <= 0:
            return False
    # No tail term is divisible by a leading monomial: each distinct tail
    # monomial is tested once, as a reduced basis's tails repeat a few
    # standard monomials many times.
    var seen = _Columns[W]()
    for p in packed:
        seen.reserve(len(p), budget)
        for k in range(len(p) - 1):
            var m = p[k].monomial
            var before = len(seen.monomials)
            _ = seen.find(m, _raw_hash(m))
            if len(seen.monomials) == before:
                continue
            budget.charge(UInt64(len(leads)))
            for l in leads:
                if _divides(l, m):
                    return False
    var terms = 0
    for q in packed:
        terms += len(q)
    if _etale[W](generators, leads, terms, order, budget, etale_start):
        return True
    var checked = _check_basis(generators, basis, order, budget)
    return checked.value() if checked else False


comptime _ETALE_LIMIT = 1 << 14
"""Standard monomials beyond which `_etale` declines and the S-pairs decide."""

comptime _ETALE_ENTRIES = 1 << 24
"""Residues the dense normal forms of `_etale` may hold in all."""

comptime _ETALE_PRIMES = 2
"""Primes `_etale` tries for an image whose leading monomials are the basis's."""

comptime _ETALE_RATIO = 16
"""Work units `_etale` may spend per basis element and basis term, plus `_ETALE_BASE`, before it declines: well above what a quotient that is reduced needs, and far below the S-pair proof where it is not."""

comptime _ETALE_BASE = 1 << 20

comptime _ETALE_RUN = 64
"""Terms a short recurrence must have generated before `_separable` declines early."""

comptime _ETALE_START = (UInt64(1) << 31) + 1
"""Odd, so that `_next_prime` starts `_etale` from 2^31 - 1."""


def _etale[W: Int](generators: List[SparsePoly[QQ]], leads: List[_Packed[W]], terms: Int, order: MonomialOrder, mut budget: Budget, start: UInt64 = _ETALE_START) raises PolyError -> Bool:
    """Whether dim_Q Q[x]/I >= N is shown modulo a prime, for I the ideal of the generators and N the number of monomials no lead divides.

    With G a set of polynomials in I whose leading monomials are `leads`
    (the caller's G is reduced, and lies in I by its exact replay), this
    makes G the reduced basis of I without S-pairs (`_proves_reduced`):

    1. G in I gives <LT(G)> in LT(I), so the monomials no element of LT(I)
       divides, which form a basis of Q[x]/I (Macaulay), are among the N
       standard monomials of LT(G): dim Q[x]/I <= N.
    2. The generators f_1, ..., f_n are n nonzero polynomials in the n
       variables. For a prime p dividing no numerator or denominator of
       their coefficients, F4 gives the reduced basis G_p of their image
       I_p, and A = F_p[x]/I_p. When LT(G_p) = LT(G), the standard monomials
       of LT(G) form a basis of A: dim A = N.
    3. If A is reduced, A tensor the algebraic closure K of F_p is K^N, so
       I_p has exactly N zeros over K, each with local ring K. There the
       images of f_1, ..., f_n generate the maximal ideal, so their
       differentials span the n-dimensional cotangent space: the Jacobian
       of f_1, ..., f_n is invertible at each zero.
    4. Each zero a lies in F_q^n for a power q of p. By Hensel's lemma in
       W(F_q), the unramified extension of the p-adic integers with residue
       field F_q, a lifts to a unique zero of f_1, ..., f_n congruent to a;
       distinct zeros lift to distinct zeros. So I has at least N zeros
       over an algebraically closed field of characteristic 0, and
       dim Q[x]/I >= N, as dimension does not change under field extension.
    5. Then dim Q[x]/I = N = colength of <LT(G)>, so LT(I) = <LT(G)>: G is
       a Groebner basis of I, and, being reduced, its reduced basis. No step
       assumes anything of p beyond 2.

    A is reduced when some l in A has a minimal polynomial m of degree N
    with gcd(m, m') = 1: then F_p[l] has dimension N, so A = F_p[l] =
    F_p[t]/(m), a product of separable field extensions (`_reduced_algebra`
    finds m by Berlekamp-Massey and checks both conditions exactly).

    Every precondition is checked, and any failure, or work past
    `_ETALE_RATIO` units per basis element and basis term (`terms`, the
    basis's total), returns False: the caller then proves G with S-pairs.
    The prime search starts below `start`, for tests.
    """
    var n = len(order.variables)
    var count = 0
    for g in generators:
        if not g.is_zero():
            count += 1
    if count != n or n == 0 or n > W or order.kind == MonomialOrder.LEX:
        return False
    var allowance = UInt64(_sum(_product(_ETALE_RATIO, _product(len(leads), terms)), _ETALE_BASE))
    # Pure powers make Q[x]/<LT(G)> finite, and bound the exponents of
    # standard and border monomials.
    var bound = List[Int](length=n, fill=0)
    for m in leads:
        if Int(pop_count(m.support)) == 1:
            var lane = Int(count_trailing_zeros(m.support))
            if lane < n:
                bound[lane] = Int(m.exponents()[lane])
    var radix = List[Int](length=n, fill=0)
    var span = 1
    for k in range(n):
        if bound[k] == 0:
            return False
        radix[k] = span
        if span > (1 << 50) // (bound[k] + 1):
            return False
        span *= bound[k] + 1
    # The standard monomials, from 1 by multiplying by each variable.
    var standard = List[_Packed[W]]()
    var index = Dict[Int, Int]()
    standard.append(_Packed[W](SIMD[DType.uint32, W](0), 0, 0))
    index[0] = 0
    var next = 0
    while next < len(standard):
        var s = standard[next]
        next += 1
        for k in range(n):
            budget.charge(UInt64(len(leads)))
            var e = s.exponents()
            e[k] += 1
            var m = _packed(e)
            var key = _key(e, radix, n)
            if key in index:
                continue
            var divided = False
            for l in leads:
                if _divides(l, m):
                    divided = True
                    break
            if divided:
                continue
            index[key] = len(standard)
            standard.append(m)
            if len(standard) > _ETALE_LIMIT:
                return False
    var lead_keys = Dict[Int, Bool]()
    for l in leads:
        lead_keys[_key(l.exponents(), radix, n)] = True
    # The generators' images modulo a prime whose basis has these leads.
    _allocation[List[_Term[QQ, W]]](n, budget)
    var packed = List[List[_Term[QQ, W]]](capacity=n)
    for g in generators:
        if not g.is_zero():
            packed.append(_pack[QQ, W](g.reorder(order, budget), order, budget))
    var prime = start
    var tries = 0
    while tries < _ETALE_PRIMES:
        tries += 1
        prime = _next_prime(prime, budget)
        while not _admissible(packed, prime, budget):
            prime = _next_prime(prime, budget)
        var ring = Modulus(prime, budget)
        var image = _f4_packed[W](_images(packed, ring, budget), order.kind, ring, budget)
        # The allowance covers the work after the image, which the S-pairs
        # would not do.
        var stop = budget.remaining() - allowance if budget.remaining() > allowance else UInt64(0)
        # LT(G_p) = LT(G): as many elements, each leading in a lead of G,
        # and the leads of a reduced basis are distinct.
        if len(image) != len(leads):
            continue
        var by_lead = Dict[Int, Int]()
        var same = True
        for i in range(len(image)):
            var key = _key(image[i][len(image[i]) - 1].monomial.exponents(), radix, n)
            if key not in lead_keys or key in by_lead:
                same = False
                break
            by_lead[key] = i
        if not same:
            continue
        return _reduced_algebra[W](image, by_lead^, standard^, index^, radix^, n, prime, stop, budget)
    return False


@always_inline
def _key[W: Int](exponents: SIMD[DType.uint32, W], radix: List[Int], n: Int) -> Int:
    var key = 0
    for k in range(n):
        key += Int(exponents[k]) * radix[k]
    return key


comptime _AXPY = 8
"""Residues per vector in `_axpy`."""


def _axpy(mut dst: List[UInt64], src: List[UInt64], c: UInt64, square: UInt64):
    """dst += c src, for residues below p in src, c below p and entries of dst below `square` = 2 p^2 (or slightly above), kept below it again."""
    var n = len(dst)
    var d = dst.unsafe_ptr()
    var x = src.unsafe_ptr()
    var factor = SIMD[DType.uint64, _AXPY](c)
    var bound = SIMD[DType.uint64, _AXPY](square)
    var i = 0
    while i + _AXPY <= n:
        var y = d.unsafe_offset(i).unsafe_load[width=_AXPY]() + factor * x.unsafe_offset(i).unsafe_load[width=_AXPY]()
        d.unsafe_offset(i).unsafe_store[width=_AXPY](y.ge(bound).select(y - bound, y))
        i += _AXPY
    while i < n:
        var y = d.unsafe_offset(i).unsafe_load() + c * x.unsafe_offset(i).unsafe_load()
        d.unsafe_offset(i).unsafe_store(y - square if y >= square else y)
        i += 1


struct _Border[W: Int](Movable):
    """Normal forms modulo a reduced monic basis over F_p, as dense vectors over the standard monomials, for monomials whose exponents stay within the pure powers."""
    var standard: List[_Packed[Self.W]]
    var index: Dict[Int, Int]
    var by_lead: Dict[Int, Int]
    var tails: List[List[Tuple[Int, UInt64]]]
    var radix: List[Int]
    var n: Int
    var p: UInt64
    var memo: Dict[Int, Int]
    var vectors: List[List[UInt64]]
    var stop: UInt64
    """The budget's remaining units below which forms are no longer computed."""

    def __init__(out self, var standard: List[_Packed[Self.W]], var index: Dict[Int, Int], var by_lead: Dict[Int, Int], var tails: List[List[Tuple[Int, UInt64]]], var radix: List[Int], n: Int, p: UInt64, stop: UInt64):
        self.standard = standard^
        self.index = index^
        self.by_lead = by_lead^
        self.tails = tails^
        self.radix = radix^
        self.n = n
        self.p = p
        self.memo = Dict[Int, Int]()
        self.vectors = List[List[UInt64]]()
        self.stop = stop

    def form(mut self, m: _Packed[Self.W], mut budget: Budget) raises PolyError -> Int:
        """The slot of the normal form of a nonstandard monomial m, or -1.

        A leading monomial's form is minus its element's tail, which is
        standard as the basis is reduced. Any other m has a variable x_j
        with m / x_j nonstandard: then its form is x_j times that one's,
        each x_j t, t standard and below m / x_j, being standard or a
        smaller nonstandard monomial (FGLM's recursion).
        """
        var size = len(self.standard)
        var key = _key(m.exponents(), self.radix, self.n)
        var known = self.memo.get(key)
        if known:
            return known.value()
        if (len(self.vectors) + 1) * size > _ETALE_ENTRIES or budget.remaining() < self.stop:
            return -1
        _allocation[UInt64](size, budget)
        var vector = List[UInt64](length=size, fill=0)
        var g = self.by_lead.get(key)
        if g:
            for t in self.tails[g.value()]:
                vector[t[0]] = t[1]
        else:
            var e = m.exponents()
            var j = -1
            for k in range(self.n):
                if e[k] == 0:
                    continue
                var smaller = e
                smaller[k] -= 1
                if _key(smaller, self.radix, self.n) not in self.index:
                    j = k
                    break
            if j < 0:
                return -1
            var below = e
            below[j] -= 1
            var v = self.form(_packed(below), budget)
            if v < 0:
                return -1
            for t in range(size):
                var c = self.vectors[v][t]
                if c == 0:
                    continue
                var w = self.standard[t].exponents()
                w[j] += 1
                var s = self.index.get(_key(w, self.radix, self.n))
                if s:
                    vector[s.value()] += c
                    continue
                var u = self.form(_packed(w), budget)
                if u < 0:
                    return -1
                budget.charge(UInt64(size))
                _axpy(vector, self.vectors[u], c, 2 * self.p * self.p)
        for i in range(size):
            vector[i] %= self.p
        self.memo[key] = len(self.vectors)
        self.vectors.append(vector^)
        return len(self.vectors) - 1


@fieldwise_init
struct _Multiplication(Movable):
    """The multiplication matrix of one variable on the quotient, by columns: a standard monomial's multiple is standard (`unit`), a leading monomial whose negated tail is its normal form (`element`), or another border monomial with a dense normal form (`dense`)."""
    var unit: List[Int]
    var element: List[Int]
    var dense: List[Int]


def _columns[W: Int](mut border: _Border[W], k: Int, mut budget: Budget) raises PolyError -> Optional[_Multiplication]:
    var size = len(border.standard)
    var unit = List[Int](length=size, fill=-1)
    var element = List[Int](length=size, fill=-1)
    var dense = List[Int](length=size, fill=-1)
    for s in range(size):
        budget.charge(1)
        var e = border.standard[s].exponents()
        e[k] += 1
        var key = _key(e, border.radix, border.n)
        var found = border.index.get(key)
        if found:
            unit[s] = found.value()
            continue
        var g = border.by_lead.get(key)
        if g and _AXPY * len(border.tails[g.value()]) < size:
            element[s] = g.value()
            continue
        var slot = border.form(_packed(e), budget)
        if slot < 0:
            return None
        # A mostly zero form is applied as a list of its nonzero entries.
        ref vector = border.vectors[slot]
        var nonzero = 0
        for r in range(size):
            if vector[r] != 0:
                nonzero += 1
        if _AXPY * nonzero >= size:
            dense[s] = slot
            continue
        var tail = List[Tuple[Int, UInt64]](capacity=nonzero)
        for r in range(size):
            if vector[r] != 0:
                tail.append((r, vector[r]))
        element[s] = len(border.tails)
        border.tails.append(tail^)
    return _Multiplication(unit^, element^, dense^)


def _apply[W: Int](border: _Border[W], columns: _Multiplication, scale: UInt64, v: List[UInt64], mut w: List[UInt64], mut budget: Budget) raises PolyError:
    """w += scale M v modulo p, M the matrix of `columns`, with w's entries kept below 2 p^2 between steps."""
    var p = border.p
    var square = 2 * p * p
    var size = len(v)
    for s in range(size):
        var x = v[s] * scale % p
        if x == 0:
            continue
        if columns.unit[s] >= 0:
            var r = columns.unit[s]
            var y = w[r] + x
            w[r] = y - square if y >= square else y
        elif columns.element[s] >= 0:
            ref tail = border.tails[columns.element[s]]
            budget.charge(UInt64(len(tail)))
            for t in tail:
                var y = w[t[0]] + x * t[1]
                w[t[0]] = y - square if y >= square else y
        else:
            budget.charge(UInt64(size))
            _axpy(w, border.vectors[columns.dense[s]], x, square)


def _separable[W: Int](border: _Border[W], columns: List[_Multiplication], scales: List[UInt64], mut budget: Budget) raises PolyError -> Bool:
    """Whether l = sum scales[k] x_k has a squarefree minimal polynomial of degree N on the quotient (see `_etale`)."""
    var p = border.p
    var size = len(border.standard)
    # The scalar sequence u L^i v for i < 2N, u and v fixed pseudorandom.
    var state = UInt64(0x9E3779B97F4A7C15)
    var u = List[UInt64](length=size, fill=0)
    var v = List[UInt64](length=size, fill=0)
    for i in range(size):
        state = state * 6364136223846793005 + 1442695040888963407
        u[i] = (state >> 33) % p
        state = state * 6364136223846793005 + 1442695040888963407
        v[i] = (state >> 33) % p
    var sequence = List[UInt64](length=2 * size, fill=0)
    var w = List[UInt64](length=size, fill=0)
    var checkpoint = 2 * _ETALE_RUN
    for i in range(2 * size):
        var dot = UInt64(0)
        for j in range(size):
            dot = (dot + u[j] * v[j]) % p
        sequence[i] = dot
        if i + 1 == 2 * size:
            break
        # A sequence of minimal polynomial degree N has linear complexity
        # about half of every prefix; one that a far shorter recurrence has
        # generated for `_ETALE_RUN` terms almost surely stays short, and the
        # proof goes to the S-pairs without the remaining products (only a
        # degree of N, after all 2N terms, ever accepts).
        if i + 1 == checkpoint:
            checkpoint += checkpoint // 4
            var prefix = List[UInt64](capacity=i + 1)
            for j in range(i + 1):
                prefix.append(sequence[j])
            if 2 * (len(_berlekamp_massey(prefix, p, budget)) - 1) + _ETALE_RUN <= i + 1:
                return False
        if budget.remaining() < border.stop:
            return False
        for j in range(size):
            w[j] = 0
        for c in range(len(columns)):
            _apply(border, columns[c], scales[c], v, w, budget)
        for j in range(size):
            v[j] = w[j] % p
    var minimal = _berlekamp_massey(sequence, p, budget)
    if len(minimal) - 1 != size:
        return False
    # m is squarefree when gcd(m, m') = 1; N < p, so m' != 0.
    var derivative = List[UInt64](length=size, fill=0)
    for i in range(1, size + 1):
        derivative[i - 1] = minimal[i] * UInt64(i) % p
    return _coprime(minimal^, derivative^, p, budget)


def _reduced_algebra[W: Int](image: List[List[_Term[ModInt, W]]], var by_lead: Dict[Int, Int], var standard: List[_Packed[W]], var index: Dict[Int, Int], var radix: List[Int], n: Int, p: UInt64, stop: UInt64, mut budget: Budget) raises PolyError -> Bool:
    """Whether the quotient by the reduced monic basis `image` is reduced, by an element whose minimal polynomial has degree its dimension and no repeated factor (see `_etale`).

    The element is the last variable, whose matrix is the cheapest. A
    combination of all variables would separate more points, but where the
    last variable does not, the quotient has so far not been reduced either,
    and the other matrices cost more than the S-pairs.
    """
    # Tails, negated, as (row, value) pairs over the standard monomials.
    var tails = List[List[Tuple[Int, UInt64]]](capacity=len(image))
    for i in range(len(image)):
        ref terms = image[i]
        var tail = List[Tuple[Int, UInt64]](capacity=len(terms) - 1)
        for k in range(len(terms) - 1):
            budget.charge(1)
            var found = index.get(_key(terms[k].monomial.exponents(), radix, n))
            if not found:
                return False
            var c = terms[k].coefficient.value
            tail.append((found.value(), (p - c) % p))
        tails.append(tail^)
    var border = _Border[W](standard^, index^, by_lead^, tails^, radix^, n, p, stop)
    var columns = List[_Multiplication]()
    var scales = List[UInt64]()
    var last = _columns[W](border, n - 1, budget)
    if not last:
        return False
    columns.append(last.take())
    scales.append(1)
    return _separable(border, columns, scales, budget)


def _power_mod(var b: UInt64, var e: UInt64, p: UInt64) -> UInt64:
    var r = UInt64(1)
    b %= p
    while e > 0:
        if e & 1:
            r = r * b % p
        b = b * b % p
        e >>= 1
    return r


def _berlekamp_massey(sequence: List[UInt64], p: UInt64, mut budget: Budget) raises PolyError -> List[UInt64]:
    """The minimal polynomial of a linearly recurrent sequence over F_p, monic, coefficients from the constant term up."""
    var c: List[UInt64] = [UInt64(1)]
    var b: List[UInt64] = [UInt64(1)]
    var length = 0
    var shift = 1
    var last = UInt64(1)
    for i in range(len(sequence)):
        budget.charge(UInt64(length + 1))
        var d = sequence[i]
        for j in range(1, length + 1):
            if j < len(c):
                d = (d + c[j] * sequence[i - j]) % p
        if d == 0:
            shift += 1
            continue
        var coefficient = d * _power_mod(last, p - 2, p) % p
        var previous = c.copy()
        while len(c) < len(b) + shift:
            c.append(0)
        for j in range(len(b)):
            c[j + shift] = (c[j + shift] + p - coefficient * b[j] % p) % p
        if 2 * length <= i:
            length = i + 1 - length
            b = previous^
            last = d
            shift = 1
        else:
            shift += 1
    # c is the connection polynomial 1 + c_1 x + ... + c_L x^L; the minimal
    # polynomial is its reversal x^L + c_1 x^(L-1) + ... + c_L.
    var minimal = List[UInt64](length=length + 1, fill=0)
    for j in range(length + 1):
        minimal[length - j] = c[j] if j < len(c) else 0
    return minimal^


def _coprime(var a: List[UInt64], var b: List[UInt64], p: UInt64, mut budget: Budget) raises PolyError -> Bool:
    """Whether two polynomials over F_p, coefficients from the constant term up, have no common factor of positive degree."""
    while len(a) > 0 and a[len(a) - 1] == 0:
        _ = a.pop()
    while len(b) > 0 and b[len(b) - 1] == 0:
        _ = b.pop()
    while len(b) > 0:
        if len(b) == 1:
            return True
        # a mod b.
        var inverse = _power_mod(b[len(b) - 1], p - 2, p)
        while len(a) >= len(b):
            budget.charge(UInt64(len(b)))
            var q = a[len(a) - 1] * inverse % p
            var offset = len(a) - len(b)
            for j in range(len(b)):
                a[offset + j] = (a[offset + j] + p - q * b[j] % p) % p
            while len(a) > 0 and a[len(a) - 1] == 0:
                _ = a.pop()
        swap(a, b)
    return len(a) == 1
