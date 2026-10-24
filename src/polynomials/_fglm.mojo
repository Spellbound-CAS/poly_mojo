"""FGLM: the reduced LEX basis of a zero-dimensional ideal from its GREVLEX basis, over a prime field and over QQ.

A zero-dimensional ideal I leaves finitely many standard monomials, the D
monomials no GREVLEX leading monomial divides, and they are a basis of the
quotient R/I. Every polynomial has a normal form there, a vector of D
residues. Faugère, Gianni, Lazard and Mora walk the monomials in increasing
LEX order: each one's normal form is either independent of those of the LEX
standard monomials found so far, and the monomial is one of them, or a
combination of them, which gives an element of the LEX basis. The walk stops
once no monomial is left that no LEX leading monomial divides.

A normal form is found from an earlier one by a multiplication matrix: the
normal form of x_k * m is M_k applied to the normal form of m, whose columns
are the normal forms of x_k * d for the standard monomials d. Those are d
itself or border monomials, taken in increasing GREVLEX order: a border
monomial is a GREVLEX leading monomial, whose normal form is minus its
element's tail, or x_k times a smaller non-standard monomial w, whose normal
form is M_k applied to w's, all of whose columns come earlier. Only a
monomial all of whose divisors by one variable are standard can be neither,
and that is a leading monomial of the reduced basis.

The walk reduces each new normal form by the earlier rows in 64-bit sums
kept below p^2 and records the multiple of each row it takes out; only for
the few monomials that turn out dependent is their combination of the LEX
standard monomials rebuilt from those multiples. The linear algebra is dense
on D entries, so the cost grows as n D^3 for n variables. `_fglm` declines
(None) for an ideal that is not zero-dimensional or has more than
`_FGLM_LIMIT` standard monomials; the caller then computes the LEX basis
directly.

Over QQ (`_fglm_rational`) the GREVLEX basis G of I comes from the certified
multimodular computation. Its primitive integer elements modulo primes above
2^31, made monic, are converted by FGLM, and the LEX images are combined
and reconstructed as for GREVLEX bases (`_Images`), each candidate checked
against the next prime's image. With certification, a candidate C is proved
to be the reduced LEX basis of I in three steps, by fraction-free reduction
in the packed engine, which never reduces C's long elements by G:

1. Every element of G reduces to zero modulo C, so I lies in the ideal J
   that C generates.
2. The S-polynomial of every two elements of C whose leading monomials share
   a variable reduces to zero modulo C; those of coprime ones always do. So
   C is a Groebner basis of J, and dim R/J is the number of monomials that
   C's leading monomials leave.
3. That number is D = dim R/I, counted from G's leading monomials. So R/I,
   which maps onto R/J since I lies in J, has the same finite dimension,
   the map is an isomorphism, and J = I.

C is reduced and monic by construction, so it is then the reduced LEX basis.
A candidate that fails any step is not returned.
"""
from std.math import max, min
from apn_mojo import Integer, Rational, lcm
from .coefficients import ModInt, Modulus, NoContext, QQ, ZZ
from .division import divide_by_list
from ._nmod import _field, _to, _mulm, _addm
from ._word_simd import _column_axpy
from ._interpolate_rational import _word_primes
from ._integer_f4 import _CERT_LANES, _CertLanes, _CertPrimes, _cert_mod
from .errors import Budget, PolyError, PolyErrorKind, _allocation, _charge_size, _product
from .monomial import MonomialOrder
from .sparse import SparsePoly
from .primes import _residue as _integer_residue
from ._rational_groebner import _Candidate, _Images, _MISMATCH_LIMIT, _PRIME_LIMIT, _agrees, _copy, _grevlex_groebner_rational, _integral, _monic, _next_prime
from ._f4 import _Columns, _f4_packed, _raw_hash, _reciprocal, _residue
from ._packed_groebner import _Engine, _Packed, _Pair, _Term, _compare, _divides, _lcm, _pack, _packed, _priorities, _unpack

comptime _FGLM_LIMIT = 4096
"""The most standard monomials `_fglm` works with; its memory grows as n D^2."""


struct _Forms[W: Int](Movable):
    """Normal forms on the GREVLEX standard monomials, numbered first in increasing GREVLEX order, then the border monomials."""
    var index: _Columns[Self.W]
    var size: Int
    """The number of standard monomials."""
    var forms: List[List[UInt32]]
    """The normal form of each border monomial, by its number minus `size`."""
    var targets: List[Int32]
    """For variable k and standard monomial t, the number of x_k * t at k * size + t."""
    var prime: UInt64
    var reciprocal: UInt64

    def __init__(out self, prime: UInt64):
        self.index = _Columns[Self.W]()
        self.size = 0
        self.forms = List[List[UInt32]]()
        self.targets = List[Int32]()
        self.prime = prime
        self.reciprocal = _reciprocal(prime)

    def multiply(self, k: Int, v: List[UInt32], mut budget: Budget) raises PolyError -> List[UInt32]:
        """M_k v: the normal form of x_k times the polynomial whose normal form is v."""
        var d = self.size
        budget.charge(UInt64(d))
        _allocation[UInt64](d, budget)
        var acc = List[UInt64](length=d, fill=0)
        var out = acc.unsafe_ptr()
        var square = self.prime * self.prime
        for t in range(d):
            var c = UInt64(v.unsafe_get(t))
            if c == 0:
                continue
            var target = Int(self.targets.unsafe_get(k * d + t))
            if target < d:
                var x = out.unsafe_offset(target).unsafe_load() + c
                out.unsafe_offset(target).unsafe_store(x - square if x >= square else x)
                continue
            budget.charge(UInt64(d))
            var column = self.forms[target - d].unsafe_ptr()
            _column_axpy(out, column, d, c, square)
        _allocation[UInt32](d, budget)
        var result = List[UInt32](capacity=d)
        for s in range(d):
            result.append(UInt32(_residue(acc.unsafe_get(s), self.prime, self.reciprocal)))
        return result^

    def vector(self, number: Int, mut budget: Budget) raises PolyError -> List[UInt32]:
        """The normal form of the monomial numbered `number`, standard or border."""
        if number >= self.size:
            return self.forms[number - self.size].copy()
        _allocation[UInt32](self.size, budget)
        var v = List[UInt32](length=self.size, fill=0)
        v[number] = 1
        return v^


@always_inline
def _unit[W: Int](k: Int) -> _Packed[W]:
    var e = SIMD[DType.uint32, W](0)
    e[k] = 1
    return _packed(e)


def _standard[W: Int](leads: List[_Packed[W]], n: Int, mut budget: Budget) raises PolyError -> Optional[List[_Packed[W]]]:
    """The monomials no leading monomial divides, or None when there are more than `_FGLM_LIMIT`."""
    var seen = _Columns[W]()
    var stack = List[_Packed[W]]()
    var found = List[_Packed[W]]()
    var one = _packed(SIMD[DType.uint32, W](0))
    seen.reserve(1, budget)
    _ = seen.find(one, _raw_hash(one))
    stack.append(one)
    while len(stack):
        var m = stack.pop()
        var divided = False
        for lead in leads:
            budget.charge(1)
            if _divides(lead, m):
                divided = True
                break
        if divided:
            continue
        found.append(m)
        if len(found) > _FGLM_LIMIT:
            return None
        seen.reserve(n, budget)
        for k in range(n):
            var u = _packed(m.exponents() + _unit[W](k).exponents())
            var before = len(seen.monomials)
            _ = seen.find(u, _raw_hash(u))
            if len(seen.monomials) > before:
                stack.append(u)
    return found^


def _ascending[W: Int](kind: UInt8, mut items: List[_Packed[W]], mut budget: Budget) raises PolyError:
    """Sort monomials into increasing order by merge sort."""
    var n = len(items)
    _allocation[_Packed[W]](n, budget)
    var scratch = items.copy()
    var width = 1
    while width < n:
        budget.charge(UInt64(n))
        var start = 0
        while start < n:
            var middle = min(start + width, n)
            var end = min(start + 2 * width, n)
            var i = start
            var j = middle
            var k = start
            while i < middle and j < end:
                if _compare(kind, items[j], items[i]) < 0:
                    scratch[k] = items[j]
                    j += 1
                else:
                    scratch[k] = items[i]
                    i += 1
                k += 1
            while i < middle:
                scratch[k] = items[i]
                i += 1
                k += 1
            while j < end:
                scratch[k] = items[j]
                j += 1
                k += 1
            start = end
        swap(items, scratch)
        width *= 2


def _normal_forms[W: Int](basis: List[List[_Term[ModInt, W]]], standard: List[_Packed[W]], n: Int, prime: UInt64, mut budget: Budget) raises PolyError -> _Forms[W]:
    """The multiplication matrices' columns: normal forms of every x_k * d for standard d."""
    var forms = _Forms[W](prime)
    var d = len(standard)
    forms.size = d
    forms.index.reserve(_product(d, n + 1), budget)
    for m in standard:
        _ = forms.index.find(m, _raw_hash(m))
    _allocation[Int32](_product(n, d), budget)
    forms.targets = List[Int32](length=n * d, fill=0)
    for k in range(n):
        for t in range(d):
            budget.charge(1)
            var u = _packed(standard[t].exponents() + _unit[W](k).exponents())
            forms.targets[k * d + t] = forms.index.find(u, _raw_hash(u))
    # A reduced basis's leading monomials are border monomials: every divisor
    # by one variable is standard. Each basis element by its lead's number.
    _allocation[Int](len(forms.index.monomials), budget)
    var element_of = List[Int](length=len(forms.index.monomials), fill=-1)
    for e in range(len(basis)):
        ref lead = basis[e][len(basis[e]) - 1].monomial
        var count = len(forms.index.monomials)
        var number = Int(forms.index.find(lead, _raw_hash(lead)))
        if len(forms.index.monomials) > count or number < d:
            raise PolyError(PolyErrorKind.DOMAIN, "FGLM needs the reduced basis")
        element_of[number] = e
    # The border monomials in increasing GREVLEX order.
    var border = List[_Packed[W]]()
    for i in range(d, len(forms.index.monomials)):
        border.append(forms.index.monomials[i])
    _ascending(MonomialOrder.GREVLEX, border, budget)
    _allocation[List[UInt32]](len(border), budget)
    for _ in range(len(border)):
        forms.forms.append(List[UInt32]())
    for u in border:
        var number = Int(forms.index.find(u, _raw_hash(u)))
        var element = element_of[number]
        if element >= 0:
            # A leading monomial: its normal form is minus its element's tail.
            ref p = basis[element]
            budget.charge(UInt64(len(p)))
            _allocation[UInt32](d, budget)
            var v = List[UInt32](length=d, fill=0)
            for k in range(len(p) - 1):
                var t = Int(forms.index.find(p[k].monomial, _raw_hash(p[k].monomial)))
                if t >= d:
                    raise PolyError(PolyErrorKind.DOMAIN, "FGLM needs the reduced basis")
                var c = p[k].coefficient.value
                v[t] = UInt32(prime - c if c != 0 else 0)
            forms.forms[number - d] = v^
            continue
        var done = False
        for k in range(n):
            if u.exponents()[k] == 0:
                continue
            var w = _packed(u.exponents() - _unit[W](k).exponents())
            var count = len(forms.index.monomials)
            var wn = Int(forms.index.find(w, _raw_hash(w)))
            if len(forms.index.monomials) > count:
                raise PolyError(PolyErrorKind.DOMAIN, "FGLM border monomial without a predecessor")
            if wn < d:
                continue
            forms.forms[number - d] = forms.multiply(k, forms.forms[wn - d], budget)
            done = True
            break
        if not done:
            raise PolyError(PolyErrorKind.DOMAIN, "FGLM needs the reduced basis")
    return forms^


struct _Row(Movable):
    """An echelon row: the normal form of the r-th LEX standard monomial less the multiples `factors` of the rows before it, scaled by `inverse` to be 1 at `pivot`, its first nonzero position."""
    var pivot: Int
    var values: List[UInt32]
    var factors: List[UInt32]
    var inverse: UInt64

    def __init__(out self, pivot: Int, var values: List[UInt32], var factors: List[UInt32], inverse: UInt64):
        self.pivot = pivot
        self.values = values^
        self.factors = factors^
        self.inverse = inverse


def _fglm[W: Int](basis: List[List[_Term[ModInt, W]]], n: Int, ring: Modulus, mut budget: Budget) raises PolyError -> Optional[List[List[_Term[ModInt, W]]]]:
    """The reduced LEX basis, in decreasing leading-monomial order, from the reduced GREVLEX basis on the same lanes; None when the ideal is not zero-dimensional or too large.

    Elements are increasing term lists, leading term last, as the packed
    engine keeps them; lanes 0 to n - 1 are the variables in priority order.
    """
    var prime = ring.value
    var leads = List[_Packed[W]]()
    var pure = 0
    for p in basis:
        if len(p) == 0:
            continue
        ref lead = p[len(p) - 1].monomial
        if lead.degree == 0:
            # The unit ideal: its LEX basis is the same.
            return basis.copy()
        leads.append(lead)
        if lead.support & (lead.support - 1) == 0:
            pure |= Int(lead.support)
    if pure != (1 << n) - 1:
        return None
    var standard = _standard[W](leads, n, budget)
    if not standard:
        return None
    var ordered = standard.take()
    _ascending(MonomialOrder.GREVLEX, ordered, budget)
    var forms = _normal_forms[W](basis, ordered, n, prime, budget)
    var d = forms.size
    var reciprocal = _reciprocal(prime)
    var square = prime * prime
    # The walk: candidates by increasing LEX order in a binary heap, each
    # with the accepted monomial and variable it came from.
    var candidates = _Columns[W]()
    var heap = List[Int]()
    var parents = List[Int]()
    var variables = List[Int]()
    var accepted = List[_Packed[W]]()
    var normal = List[List[UInt32]]()
    var rows = List[_Row]()
    var found = List[List[_Term[ModInt, W]]]()
    var found_leads = List[_Packed[W]]()
    var one = _packed(SIMD[DType.uint32, W](0))
    candidates.reserve(1, budget)
    _ = candidates.find(one, _raw_hash(one))
    parents.append(-1)
    variables.append(-1)
    heap.append(0)
    while len(heap):
        budget.charge(1)
        var c = _pop_least(heap, candidates.monomials, budget)
        var t = candidates.monomials[c]
        var skip = False
        for lead in found_leads:
            budget.charge(1)
            if _divides(lead, t):
                skip = True
                break
        if skip:
            continue
        var v: List[UInt32]
        if parents[c] < 0:
            v = forms.vector(0, budget)
        else:
            v = forms.multiply(variables[c], normal[parents[c]], budget)
        var own = v.copy()
        # Reduce by the rows in 64-bit words, each kept below p^2, recording
        # the multiple of each row taken out.
        var count = len(rows)
        _allocation[UInt64](d, budget)
        var acc = List[UInt64](length=d, fill=0)
        var a = acc.unsafe_ptr()
        for s in range(d):
            a.unsafe_offset(s).unsafe_store(UInt64(v.unsafe_get(s)))
        _allocation[UInt32](count, budget)
        var factors = List[UInt32](length=count, fill=0)
        for r in range(count):
            ref row = rows[r]
            var f = _residue(a.unsafe_offset(row.pivot).unsafe_load(), prime, reciprocal)
            if f == 0:
                continue
            factors[r] = UInt32(f)
            budget.charge(UInt64(d - row.pivot))
            var g = prime - f
            var source = row.values.unsafe_ptr()
            for s in range(row.pivot, d):
                var x = a.unsafe_offset(s).unsafe_load() + g * UInt64(source.unsafe_offset(s).unsafe_load())
                a.unsafe_offset(s).unsafe_store(x - square if x >= square else x)
        var pivot = -1
        for s in range(d):
            var x = _residue(a.unsafe_offset(s).unsafe_load(), prime, reciprocal)
            v.unsafe_get(s) = UInt32(x)
            if x != 0 and pivot < 0:
                pivot = s
        if pivot < 0:
            # t's normal form is the sum of factors[r] times row r: rewrite
            # it in the LEX standard monomials' normal forms, from the last
            # row back, and t less that combination is a basis element.
            _allocation[UInt64](2 * count, budget)
            var x = List[UInt64](capacity=count)
            for r in range(count):
                x.append(UInt64(factors[r]))
            var coefficients = List[UInt64](length=count, fill=0)
            for r in range(count - 1, -1, -1):
                if x[r] == 0:
                    continue
                var c = (x[r] * rows[r].inverse) % prime
                coefficients[r] = c
                ref earlier = rows[r].factors
                budget.charge(UInt64(r + 1))
                for q in range(r):
                    var e = UInt64(earlier.unsafe_get(q))
                    if e != 0:
                        x[q] = (x[q] + (prime - c) * e) % prime
            var terms = List[_Term[ModInt, W]]()
            for r in range(count):
                if coefficients[r] != 0:
                    terms.append(_Term[ModInt, W](accepted[r], ModInt(prime - coefficients[r], ring)))
            terms.append(_Term[ModInt, W](t, ModInt(1, ring)))
            _sort_terms(terms, budget)
            found.append(terms^)
            found_leads.append(t)
            continue
        var inverse = ModInt(UInt64(v[pivot]), ring).inverse(budget).value
        for s in range(pivot, d):
            v.unsafe_get(s) = UInt32(_residue(UInt64(v.unsafe_get(s)) * inverse, prime, reciprocal))
        rows.append(_Row(pivot, v^, factors^, inverse))
        var position = len(accepted)
        accepted.append(t)
        normal.append(own^)
        if len(accepted) > d:
            raise PolyError(PolyErrorKind.DOMAIN, "FGLM found more standard monomials than the GREVLEX basis")
        candidates.reserve(n, budget)
        for k in range(n):
            var u = _packed(t.exponents() + _unit[W](k).exponents())
            var before = len(candidates.monomials)
            var number = Int(candidates.find(u, _raw_hash(u)))
            if len(candidates.monomials) > before:
                parents.append(position)
                variables.append(k)
                _push(heap, number, candidates.monomials, budget)
    if len(accepted) != d:
        raise PolyError(PolyErrorKind.DOMAIN, "FGLM found fewer standard monomials than the GREVLEX basis")
    # Decreasing leading-monomial order, as the packed engine returns bases.
    var result = List[List[_Term[ModInt, W]]](capacity=len(found))
    for i in range(len(found) - 1, -1, -1):
        var p = List[_Term[ModInt, W]]()
        swap(p, found[i])
        result.append(p^)
    return result^


def _sort_terms[W: Int](mut terms: List[_Term[ModInt, W]], mut budget: Budget) raises PolyError:
    """Sort terms into increasing LEX order by insertion; the lists are short."""
    for i in range(1, len(terms)):
        var j = i
        while j > 0 and _compare(MonomialOrder.LEX, terms[j].monomial, terms[j - 1].monomial) < 0:
            budget.charge(1)
            terms.swap_elements(j, j - 1)
            j -= 1


def _push[W: Int](mut heap: List[Int], item: Int, monomials: List[_Packed[W]], mut budget: Budget) raises PolyError:
    heap.append(item)
    var i = len(heap) - 1
    while i > 0:
        budget.charge(1)
        var parent = (i - 1) // 2
        if _compare(MonomialOrder.LEX, monomials[heap[i]], monomials[heap[parent]]) >= 0:
            break
        heap.swap_elements(i, parent)
        i = parent


def _pop_least[W: Int](mut heap: List[Int], monomials: List[_Packed[W]], mut budget: Budget) raises PolyError -> Int:
    var least = heap[0]
    var last = heap.pop()
    if len(heap) == 0:
        return least
    heap[0] = last
    var i = 0
    while True:
        budget.charge(1)
        var smallest = i
        var left = 2 * i + 1
        var right = left + 1
        if left < len(heap) and _compare(MonomialOrder.LEX, monomials[heap[left]], monomials[heap[smallest]]) < 0:
            smallest = left
        if right < len(heap) and _compare(MonomialOrder.LEX, monomials[heap[right]], monomials[heap[smallest]]) < 0:
            smallest = right
        if smallest == i:
            break
        heap.swap_elements(i, smallest)
        i = smallest
    return least


def _fglm_basis[W: Int](generators: List[SparsePoly[ModInt]], order: MonomialOrder, mut budget: Budget) raises PolyError -> Optional[List[SparsePoly[ModInt]]]:
    var ring = generators[0].ring.copy()
    var grevlex = MonomialOrder(MonomialOrder.GREVLEX, order.variables.copy(), budget)
    _allocation[List[_Term[ModInt, W]]](len(generators), budget)
    var packed = List[List[_Term[ModInt, W]]](capacity=len(generators))
    for g in generators:
        budget.charge(1)
        packed.append(_pack[ModInt, W](g.reorder(grevlex, budget), grevlex, budget))
    var reduced = _f4_packed[W](packed^, MonomialOrder.GREVLEX, ring, budget)
    var lex = _fglm[W](reduced, len(order.variables), ring, budget)
    if not lex:
        return None
    var basis = lex.take()
    var by_index = _priorities(order, budget)
    _allocation[SparsePoly[ModInt]](len(basis), budget)
    var result = List[SparsePoly[ModInt]](capacity=len(basis))
    for i in range(len(basis)):
        budget.charge(1)
        var p = List[_Term[ModInt, W]]()
        swap(p, basis[i])
        result.append(_unpack(p, order, ring, by_index, budget))
    return result^


def _fglm_supported(order: MonomialOrder, prime: UInt64) -> Bool:
    """LEX over at most sixteen variables and a prime F4 takes."""
    return order.kind == MonomialOrder.LEX and len(order.variables) > 0 and len(order.variables) <= 16 and prime < (UInt64(1) << 31)


def _fglm_groebner(generators: List[SparsePoly[ModInt]], order: MonomialOrder, mut budget: Budget) raises PolyError -> Optional[List[SparsePoly[ModInt]]]:
    """The reduced LEX basis by F4 under GREVLEX and FGLM, or None when the ideal is not zero-dimensional or too large; generators are nonempty and share a field."""
    if len(order.variables) <= 8:
        return _fglm_basis[8](generators, order, budget)
    return _fglm_basis[16](generators, order, budget)


@always_inline
def _lane_add(a: _CertLanes, b: _CertLanes, lanes: _CertPrimes) -> _CertLanes:
    """a + b for lanes below p^2 each, kept below p^2."""
    var x = a + b
    return x.ge(lanes.squares).select(x - lanes.squares, x)


def _lane_multiply(forms: List[List[_CertLanes]], targets: List[Int32], d: Int, k: Int, v: List[_CertLanes], lanes: _CertPrimes, mut budget: Budget) raises PolyError -> List[_CertLanes]:
    """M_k v in every lane: as `_Forms.multiply`, with border normal forms by lanes."""
    budget.charge(UInt64(d))
    _allocation[_CertLanes](d, budget)
    var acc = List[_CertLanes](length=d, fill=_CertLanes(0))
    for t in range(d):
        var c = v.unsafe_get(t)
        if c.reduce_or() == 0:
            continue
        var target = Int(targets.unsafe_get(k * d + t))
        if target < d:
            acc.unsafe_get(target) = _lane_add(acc.unsafe_get(target), c, lanes)
            continue
        budget.charge(UInt64(d))
        ref column = forms[target - d]
        for s in range(d):
            acc.unsafe_get(s) = _lane_add(acc.unsafe_get(s), c * column.unsafe_get(s), lanes)
    for s in range(d):
        acc.unsafe_get(s) = _cert_mod(acc.unsafe_get(s), lanes)
    return acc^


def _fglm_lanes[W: Int](bases: List[List[List[_Term[ModInt, W]]]], n: Int, rings: List[Modulus], mut budget: Budget) raises PolyError -> Optional[List[List[List[_Term[ModInt, W]]]]]:
    """`_fglm` for `_CERT_LANES` primes at once, one per lane; None when the lanes' leading monomials or pivots differ, or for a case `_fglm` itself declines.

    The GREVLEX images share their leading monomials and so the standard and
    border monomials and every multiplication target; the LEX walk shares
    its decisions while each lane's elimination puts its pivot in the same
    column. That is checked at every step: a lane that differs voids the
    batch, and the caller runs `_fglm` prime by prime.
    """
    var count = len(bases)
    if count != _CERT_LANES:
        return None
    var leads = List[_Packed[W]]()
    var pure = 0
    ref first = bases[0]
    for p in first:
        if len(p) == 0:
            return None
        ref lead = p[len(p) - 1].monomial
        if lead.degree == 0:
            return None
        leads.append(lead)
        if lead.support & (lead.support - 1) == 0:
            pure |= Int(lead.support)
    if pure != (1 << n) - 1:
        return None
    for l in range(1, count):
        if len(bases[l]) != len(first):
            return None
        for e in range(len(first)):
            ref other = bases[l][e]
            if len(other) == 0 or _compare(MonomialOrder.GREVLEX, other[len(other) - 1].monomial, leads[e]) != 0:
                return None
    var standard = _standard[W](leads, n, budget)
    if not standard:
        return None
    var ordered = standard.take()
    _ascending(MonomialOrder.GREVLEX, ordered, budget)
    var primes = List[UInt64](capacity=count)
    for ring in rings:
        primes.append(ring.value)
    var lanes = _CertPrimes(primes, 0)
    # Numbering and targets as `_normal_forms`.
    var index = _Columns[W]()
    var d = len(ordered)
    index.reserve(_product(d, n + 1), budget)
    for m in ordered:
        _ = index.find(m, _raw_hash(m))
    _allocation[Int32](_product(n, d), budget)
    var targets = List[Int32](length=n * d, fill=0)
    for k in range(n):
        for t in range(d):
            budget.charge(1)
            var u = _packed(ordered[t].exponents() + _unit[W](k).exponents())
            targets[k * d + t] = index.find(u, _raw_hash(u))
    _allocation[Int](len(index.monomials), budget)
    var element_of = List[Int](length=len(index.monomials), fill=-1)
    for e in range(len(first)):
        var before = len(index.monomials)
        var number = Int(index.find(leads[e], _raw_hash(leads[e])))
        if len(index.monomials) > before or number < d:
            return None
        element_of[number] = e
    var border = List[_Packed[W]]()
    for i in range(d, len(index.monomials)):
        border.append(index.monomials[i])
    _ascending(MonomialOrder.GREVLEX, border, budget)
    var forms = List[List[_CertLanes]](length=len(border), fill=List[_CertLanes]())
    for u in border:
        var number = Int(index.find(u, _raw_hash(u)))
        var element = element_of[number]
        if element >= 0:
            # Minus the element's tail, lane by lane.
            _allocation[_CertLanes](d, budget)
            var v = List[_CertLanes](length=d, fill=_CertLanes(0))
            for l in range(count):
                ref p = bases[l][element]
                budget.charge(UInt64(len(p)))
                for k in range(len(p) - 1):
                    var t = Int(index.find(p[k].monomial, _raw_hash(p[k].monomial)))
                    if t >= d:
                        return None
                    var c = p[k].coefficient.value
                    v[t][l] = primes[l] - c if c != 0 else 0
            forms[number - d] = v^
            continue
        var done = False
        for k in range(n):
            if u.exponents()[k] == 0:
                continue
            var w = _packed(u.exponents() - _unit[W](k).exponents())
            var before = len(index.monomials)
            var wn = Int(index.find(w, _raw_hash(w)))
            if len(index.monomials) > before:
                return None
            if wn < d:
                continue
            forms[number - d] = _lane_multiply(forms, targets, d, k, forms[wn - d], lanes, budget)
            done = True
            break
        if not done:
            return None
    # The walk, shared by the lanes.
    var candidates = _Columns[W]()
    var heap = List[Int]()
    var parents = List[Int]()
    var variables = List[Int]()
    var accepted = List[_Packed[W]]()
    var normal = List[List[_CertLanes]]()
    var pivots = List[Int]()
    var rows = List[List[_CertLanes]]()
    var row_factors = List[List[_CertLanes]]()
    var inverses = List[_CertLanes]()
    var found = List[List[List[_Term[ModInt, W]]]](length=count, fill=List[List[_Term[ModInt, W]]]())
    var found_leads = List[_Packed[W]]()
    var one = _packed(SIMD[DType.uint32, W](0))
    candidates.reserve(1, budget)
    _ = candidates.find(one, _raw_hash(one))
    parents.append(-1)
    variables.append(-1)
    heap.append(0)
    while len(heap):
        budget.charge(1)
        var c = _pop_least(heap, candidates.monomials, budget)
        var t = candidates.monomials[c]
        var skip = False
        for lead in found_leads:
            budget.charge(1)
            if _divides(lead, t):
                skip = True
                break
        if skip:
            continue
        var v: List[_CertLanes]
        if parents[c] < 0:
            v = List[_CertLanes](length=d, fill=_CertLanes(0))
            v[0] = _CertLanes(1)
        else:
            v = _lane_multiply(forms, targets, d, variables[c], normal[parents[c]], lanes, budget)
        var own = v.copy()
        var rcount = len(rows)
        _allocation[_CertLanes](rcount, budget)
        var factors = List[_CertLanes](length=rcount, fill=_CertLanes(0))
        for r in range(rcount):
            var pivot = pivots[r]
            var f = _cert_mod(v.unsafe_get(pivot), lanes)
            if f.reduce_or() == 0:
                continue
            factors[r] = f
            budget.charge(UInt64(d - pivot))
            var g = f.eq(0).select(f, lanes.primes - f)
            ref source = rows[r]
            for s in range(pivot, d):
                v.unsafe_get(s) = _lane_add(v.unsafe_get(s), g * source.unsafe_get(s), lanes)
        var pivot = -1
        var mixed = False
        for s in range(d):
            var x = _cert_mod(v.unsafe_get(s), lanes)
            v.unsafe_get(s) = x
            if pivot < 0:
                var nonzero = x.ne(0)
                if nonzero.reduce_or():
                    if not nonzero.reduce_and():
                        mixed = True
                    pivot = s
        if mixed:
            return None
        if pivot < 0:
            # Back-substitution in every lane, as `_fglm`.
            var x = factors.copy()
            var coefficients = List[_CertLanes](length=rcount, fill=_CertLanes(0))
            for r in range(rcount - 1, -1, -1):
                if x[r].reduce_or() == 0:
                    continue
                var cc = _cert_mod(_cert_mod(x[r], lanes) * inverses[r], lanes)
                coefficients[r] = cc
                ref earlier = row_factors[r]
                budget.charge(UInt64(r + 1))
                var neg = cc.eq(0).select(cc, lanes.primes - cc)
                for q in range(r):
                    var e = earlier.unsafe_get(q)
                    if e.reduce_or() != 0:
                        x[q] = _cert_mod(x[q] + _cert_mod(neg * e, lanes), lanes)
            for l in range(count):
                var terms = List[_Term[ModInt, W]]()
                for r in range(rcount):
                    var value = coefficients[r][l]
                    if value != 0:
                        terms.append(_Term[ModInt, W](accepted[r], ModInt(primes[l] - value, rings[l])))
                terms.append(_Term[ModInt, W](t, ModInt(1, rings[l])))
                _sort_terms(terms, budget)
                found[l].append(terms^)
            found_leads.append(t)
            continue
        var inverse = _CertLanes(0)
        var lead = v[pivot]
        for l in range(count):
            inverse[l] = ModInt(lead[l], rings[l]).inverse(budget).value
        for s in range(pivot, d):
            v.unsafe_get(s) = _cert_mod(v.unsafe_get(s) * inverse, lanes)
        pivots.append(pivot)
        rows.append(v^)
        row_factors.append(factors^)
        inverses.append(inverse)
        var position = len(accepted)
        accepted.append(t)
        normal.append(own^)
        if len(accepted) > d:
            return None
        candidates.reserve(n, budget)
        for k in range(n):
            var u = _packed(t.exponents() + _unit[W](k).exponents())
            var before = len(candidates.monomials)
            var number = Int(candidates.find(u, _raw_hash(u)))
            if len(candidates.monomials) > before:
                parents.append(position)
                variables.append(k)
                _push(heap, number, candidates.monomials, budget)
    if len(accepted) != d:
        return None
    var results = List[List[List[_Term[ModInt, W]]]](capacity=count)
    for l in range(count):
        var result = List[List[_Term[ModInt, W]]](capacity=len(found[l]))
        for i in range(len(found[l]) - 1, -1, -1):
            var p = List[_Term[ModInt, W]]()
            swap(p, found[l][i])
            result.append(p^)
        results.append(result^)
    return results^


def _image[W: Int](integral: List[List[_Term[ZZ, W]]], p: UInt64, ring: Modulus, mut budget: Budget) raises PolyError -> Optional[List[List[_Term[ModInt, W]]]]:
    """Each primitive integer element modulo p, made monic; None when p divides a leading coefficient."""
    var image = List[List[_Term[ModInt, W]]](capacity=len(integral))
    for q in integral:
        var lead = _integer_residue(q[len(q) - 1].coefficient.value, p, budget)
        if lead == 0:
            return None
        var inverse = ModInt(lead, ring).inverse(budget).value
        _allocation[_Term[ModInt, W]](len(q), budget)
        var terms = List[_Term[ModInt, W]](capacity=len(q))
        for t in q:
            var r = _integer_residue(t.coefficient.value, p, budget)
            if r != 0:
                terms.append(_Term[ModInt, W](t.monomial, ModInt(UInt64((UInt128(r) * UInt128(inverse)) % UInt128(p)), ring)))
        image.append(terms^)
    return image^


def _rational_lex[W: Int](reduced: List[SparsePoly[QQ]], grevlex: MonomialOrder, order: MonomialOrder, certify: Bool, mut budget: Budget) raises PolyError -> Optional[List[SparsePoly[QQ]]]:
    """The reduced LEX basis over QQ from the reduced GREVLEX one, by FGLM modulo primes and reconstruction; None when the ideal is not zero-dimensional or reconstruction gives up.

    With `certify`, the result is proved: every element reduces to zero
    modulo the GREVLEX basis, so the ideal J it generates lies in I; and its
    leading monomials leave as many standard monomials as I's, D, so R/J has
    dimension at most D = dim R/I, which forces J = I with the result's
    leading monomials generating LT(I). It is reduced and monic by
    construction, hence the reduced LEX basis.
    """
    var n = len(order.variables)
    var integral = List[List[_Term[ZZ, W]]](capacity=len(reduced))
    var leads = List[_Packed[W]]()
    var pure = 0
    for g in reduced:
        budget.charge(1)
        var q = _integral(_pack[QQ, W](g, grevlex, budget), budget)
        ref lead = q[len(q) - 1].monomial
        if lead.degree == 0:
            return None
        leads.append(lead)
        if lead.support & (lead.support - 1) == 0:
            pure |= Int(lead.support)
        integral.append(q^)
    if pure != (1 << n) - 1:
        return None
    var standard = _standard[W](leads, n, budget)
    if not standard:
        return None
    var dimension = len(standard.value())
    var images = _Images[W](MonomialOrder.LEX)
    var candidate = Optional[_Candidate[W]](None)
    var mismatches = 0
    var attempt = 1
    var next = (UInt64(1) << 31) + 1
    # Primes come in batches of `_CERT_LANES` images, one FGLM for all of
    # them; each image's basis is then taken in turn as before.
    var queued = List[UInt64]()
    var bases = List[Optional[List[List[_Term[ModInt, W]]]]]()
    var head = 0
    var tried = 0
    while True:
        if head == len(queued):
            if tried >= _PRIME_LIMIT:
                break
            queued = List[UInt64]()
            bases = List[Optional[List[List[_Term[ModInt, W]]]]]()
            head = 0
            var batch = List[List[List[_Term[ModInt, W]]]]()
            var rings = List[Modulus]()
            while len(batch) < _CERT_LANES and tried < _PRIME_LIMIT:
                tried += 1
                budget.charge(1)
                next = _next_prime(next, budget)
                var image = _image[W](integral, next, Modulus(_prime_value=next), budget)
                if image:
                    batch.append(image.take())
                    rings.append(Modulus(_prime_value=next))
                    queued.append(next)
            var lanes = _fglm_lanes[W](batch, n, rings, budget)
            if lanes:
                var results = lanes.take()
                for i in range(len(results)):
                    var one = List[List[_Term[ModInt, W]]]()
                    swap(one, results[i])
                    bases.append(Optional(one^))
            else:
                for i in range(len(batch)):
                    bases.append(_fglm[W](batch[i], n, rings[i], budget))
            if len(queued) == 0:
                break
        var p = queued[head]
        var lex = bases[head].take() if bases[head] else Optional[List[List[_Term[ModInt, W]]]](None)
        head += 1
        if not lex:
            continue
        var basis = lex.take()
        if images.count > 0 and not images.matches(basis):
            mismatches += 1
            if mismatches < _MISMATCH_LIMIT:
                continue
            images.reset()
            candidate = None
            attempt = 1
        mismatches = 0
        if candidate:
            if _agrees(MonomialOrder.LEX, candidate.value(), images.table, basis, p, budget):
                var result = _lex_result[W](images.written(candidate.take(), budget), reduced, order, grevlex, dimension, n, certify, budget)
                if result or not certify:
                    return result^
                return None
            candidate = None
        images.add(basis, p, budget)
        if images.count >= attempt:
            candidate = images.reconstruct(False, budget)
            attempt = 2 * images.count if images.count < 8 else images.count + max(8, images.count // 8)
    return None


def _exact_form[W: Int](m: _Packed[W], reduced: List[SparsePoly[QQ]], grevlex: MonomialOrder, by_index: List[Int], mut columns: _Columns[W], dimension: Int, mut budget: Budget) raises PolyError -> Optional[List[Tuple[Int, Rational]]]:
    """The normal form of m modulo the GREVLEX basis, as (standard index, coefficient) pairs; None if a term is not standard."""
    columns.reserve(1, budget)
    var index = Int(columns.find(m, _raw_hash(m)))
    var result = List[Tuple[Int, Rational]]()
    if index < dimension:
        result.append((index, Rational(1)))
        return result^
    var terms: List[_Term[QQ, W]] = [_Term[QQ, W](m, QQ(Rational(1)))]
    var poly = _unpack[QQ, W](terms, grevlex, NoContext(), by_index, budget)
    var division = divide_by_list(poly, reduced, budget)
    for t in _pack[QQ, W](division.remainder, grevlex, budget):
        columns.reserve(1, budget)
        var row = Int(columns.find(t.monomial, _raw_hash(t.monomial)))
        if row >= dimension:
            return None
        result.append((row, t.coefficient.value))
    return result^


def _scaled_form(form: List[Tuple[Int, Rational]], scale: Integer) raises -> List[Tuple[Int, Integer]]:
    var result = List[Tuple[Int, Integer]](capacity=len(form))
    for entry in form:
        ref v = entry[1]
        result.append((entry[0], v.numerator() * (scale // v.denominator())))
    return result^


def _shape_certified[W: Int](candidate: List[List[_Term[ZZ, W]]], leads: List[_Packed[W]], reduced: List[SparsePoly[QQ]], grevlex: MonomialOrder, n: Int, dimension: Int, mut budget: Budget) raises PolyError -> Optional[Bool]:
    """Whether a LEX candidate in shape position lies in I; None when it is not in shape position.

    Shape position: one element leads in x_n^D, with x_n the last variable and
    D = dim R/I, and one leads in each other variable to the first power; all
    other terms are powers of x_n below D. With M the matrix of
    multiplication by x_n on the GREVLEX standard monomials and e the normal
    form of 1, the normal form of h(x_n) is h(M) e, so f(x_n) lies in I
    exactly when f(M) e = 0, and a x_i + h(x_n) exactly when
    a NF(x_i) + h(M) e = 0. With delta the common denominator of M and of
    the NF(x_i), multiplying by delta^D makes each identity integral, with
    entries bounded by the infinity-norm sum over its terms; vanishing modulo
    word primes whose product exceeds twice that bound proves it. Every
    element in I gives <LM(C)> inside LT(I); both leave D standard
    monomials, so they are equal and the candidate is a Groebner basis of I.
    """
    var last = n - 1
    if len(candidate) != n or dimension < 1:
        return None
    var power = UInt32(1) << UInt32(last)
    var owner = List[Int](length=n, fill=-1)
    for e in range(len(candidate)):
        budget.charge(1)
        var lead = leads[e]
        var v = -1
        if lead.support == power:
            if Int(lead.degree) != dimension:
                return None
            v = last
        elif lead.degree == 1:
            var ex = lead.exponents()
            for k in range(n):
                if ex[k] == 1:
                    v = k
            if v == last:
                return None
        else:
            return None
        if owner[v] >= 0:
            return None
        owner[v] = e
        for t in range(len(candidate[e]) - 1):
            ref m = candidate[e][t].monomial
            if m.support & ~power != 0 or Int(m.degree) >= dimension:
                return None
    # The GREVLEX standard monomials, numbered in a column table.
    var glead = List[_Packed[W]](capacity=len(reduced))
    for g in reduced:
        var q = _pack[QQ, W](g, grevlex, budget)
        glead.append(q[len(q) - 1].monomial)
    var found = _standard[W](glead, n, budget)
    if not found or len(found.value()) != dimension:
        return None
    var standard = found.take()
    var columns = _Columns[W]()
    columns.reserve(dimension, budget)
    for s in standard:
        _ = columns.find(s, _raw_hash(s))
    var by_index = _priorities(grevlex, budget)
    var one = _packed(SIMD[DType.uint32, W](0))
    columns.reserve(1, budget)
    var origin = Int(columns.find(one, _raw_hash(one)))
    if origin >= dimension:
        return None
    # M by columns, and NF(x_i), exactly.
    _allocation[List[Tuple[Int, Rational]]](dimension + n, budget)
    var exact = List[List[Tuple[Int, Rational]]](capacity=dimension + n)
    for j in range(dimension):
        var m = _packed(standard[j].exponents() + _unit[W](last).exponents())
        var form = _exact_form[W](m, reduced, grevlex, by_index, columns, dimension, budget)
        if not form:
            return None
        exact.append(form.take())
    for i in range(last):
        var form = _exact_form[W](_unit[W](i), reduced, grevlex, by_index, columns, dimension, budget)
        if not form:
            return None
        exact.append(form.take())
    var delta = Integer(1)
    var matrix = List[List[Tuple[Int, Integer]]](capacity=dimension + n)
    var top = 0
    try:
        for form in exact:
            for entry in form:
                ref den = entry[1].denominator()
                if den != Integer(1):
                    delta = lcm(delta, den)
        for form in exact:
            matrix.append(_scaled_form(form, delta))
        # Row sums of |delta M| bound its infinity norm by 2^norm.
        var sums = List[Integer](length=dimension, fill=Integer(0))
        for j in range(dimension):
            for entry in matrix[j]:
                sums[entry[0]] = sums[entry[0]] + abs(entry[1])
        var norm = 0
        for x in sums:
            norm = max(norm, x.magnitude_bit_length())
        var dbits = delta.magnitude_bit_length()
        # Bits of each element's identity, scaled by delta^D.
        for v in range(n):
            ref element = candidate[owner[v]]
            var b = 0
            var count = len(element) + 1
            for t in range(len(element)):
                ref term = element[t]
                if v != last and t == len(element) - 1:
                    var u = 0
                    for entry in matrix[dimension + v]:
                        u = max(u, entry[1].magnitude_bit_length())
                    b = max(b, term.coefficient.value.magnitude_bit_length() + (dimension - 1) * dbits + u)
                else:
                    var k = Int(term.monomial.degree)
                    b = max(b, term.coefficient.value.magnitude_bit_length() + (dimension - k) * dbits + k * norm)
            # The sum of `count` terms adds at most log2(count) bits.
            while count > 1:
                b += 1
                count = (count + 1) // 2
            top = max(top, b)
    except e:
        raise PolyError(PolyErrorKind.CAPACITY, String(e))
    var primes = _word_primes(top // 61 + 2, UInt64(1) << 61, budget)
    _allocation[UInt64](_product(dimension + 1, dimension), budget)
    for p in primes:
        budget.charge(1)
        var f = _field(p)
        var reduced_matrix = List[List[Tuple[Int, UInt64]]](capacity=len(matrix))
        for column in matrix:
            var c = List[Tuple[Int, UInt64]](capacity=len(column))
            for entry in column:
                c.append((entry[0], _to(_integer_residue(entry[1], p, budget), f)))
            reduced_matrix.append(c^)
        # Krylov vectors M^k e, k = 0 ... D.
        var krylov = List[List[UInt64]](capacity=dimension + 1)
        var start = List[UInt64](length=dimension, fill=0)
        start[origin] = f.one
        krylov.append(start^)
        _charge_size(_product(dimension, dimension), budget)
        for k in range(1, dimension + 1):
            var next = List[UInt64](length=dimension, fill=0)
            ref previous = krylov[k - 1]
            for j in range(dimension):
                var x = previous[j]
                if x == 0:
                    continue
                for entry in reduced_matrix[j]:
                    next[entry[0]] = _addm(next[entry[0]], _mulm(entry[1], x, f), f)
            krylov.append(next^)
        var powers = List[UInt64](length=dimension + 1, fill=f.one)
        var d = _to(_integer_residue(delta, p, budget), f)
        for t in range(1, dimension + 1):
            powers[t] = _mulm(powers[t - 1], d, f)
        for v in range(n):
            ref element = candidate[owner[v]]
            var residual = List[UInt64](length=dimension, fill=0)
            for t in range(len(element)):
                ref term = element[t]
                var c = _to(_integer_residue(term.coefficient.value, p, budget), f)
                if v != last and t == len(element) - 1:
                    var scale = _mulm(c, powers[dimension - 1], f)
                    for entry in reduced_matrix[dimension + v]:
                        residual[entry[0]] = _addm(residual[entry[0]], _mulm(scale, entry[1], f), f)
                else:
                    var k = Int(term.monomial.degree)
                    var scale = _mulm(c, powers[dimension - k], f)
                    ref vector = krylov[k]
                    for r in range(dimension):
                        if vector[r] != 0:
                            residual[r] = _addm(residual[r], _mulm(scale, vector[r], f), f)
            for x in residual:
                if x != 0:
                    return False
    return True


def _lex_result[W: Int](var candidate: List[List[_Term[ZZ, W]]], reduced: List[SparsePoly[QQ]], order: MonomialOrder, grevlex: MonomialOrder, dimension: Int, n: Int, certify: Bool, mut budget: Budget) raises PolyError -> Optional[List[SparsePoly[QQ]]]:
    """The candidate as monic LEX polynomials, certified when asked; None when certification fails.

    The proof reduces only the small GREVLEX elements and S-polynomials,
    never the candidate's large elements by the GREVLEX basis. Let J be the
    ideal the candidate generates. Every GREVLEX element reducing to zero
    modulo the candidate gives I in J. Every S-polynomial of two elements
    whose leading monomials share a variable reducing to zero (those of
    coprime ones always do) makes the candidate a Groebner basis of J, so
    dim R/J is the number of monomials its leading monomials leave. That
    number being D = dim R/I, I in J forces J = I, and the reduced, monic
    candidate is the reduced LEX basis of I.
    """
    var leads = List[_Packed[W]]()
    for q in candidate:
        leads.append(q[len(q) - 1].monomial)
    var proved = False
    if certify:
        var standard = _standard[W](leads, n, budget)
        if not standard or len(standard.value()) != dimension:
            return None
        var shaped = _shape_certified[W](candidate, leads, reduced, grevlex, n, dimension, budget)
        if shaped:
            if not shaped.value():
                return None
            proved = True
    if certify and not proved:
        # Fraction-free reduction in the packed engine, over the candidate.
        var engine = _Engine[ZZ, W, True](MonomialOrder.LEX, NoContext())
        for e in range(len(candidate)):
            engine.leads.append(leads[e])
            engine.active.append(e)
            engine.polys.append(_copy(candidate[e], budget))
        engine.sort_active(budget)
        for g in reduced:
            budget.charge(1)
            var q = _integral(_pack[QQ, W](g.reorder(order, budget), order, budget), budget)
            if len(engine.reduce_base(q^, engine.active, budget)) != 0:
                return None
        for a in range(len(candidate)):
            for b in range(a + 1, len(candidate)):
                budget.charge(1)
                if leads[a].support & leads[b].support == 0:
                    continue
                if len(engine.reduce_pair(_Pair[W](a, b, _lcm(leads[a], leads[b])), budget)) != 0:
                    return None
    var by_index = _priorities(order, budget)
    var result = List[SparsePoly[QQ]](capacity=len(candidate))
    for e in range(len(candidate)):
        budget.charge(1)
        var q = List[_Term[ZZ, W]]()
        swap(q, candidate[e])
        result.append(_monic(q, order, by_index, budget))
    return result^


def _fglm_rational_supported(order: MonomialOrder) -> Bool:
    """LEX over at most fifteen variables, which the GREVLEX basis over QQ takes."""
    return order.kind == MonomialOrder.LEX and len(order.variables) > 0 and len(order.variables) + 1 <= 16


def _fglm_rational(generators: List[SparsePoly[QQ]], order: MonomialOrder, certify: Bool, mut budget: Budget) raises PolyError -> Optional[List[SparsePoly[QQ]]]:
    """The reduced LEX basis over QQ from the GREVLEX one by multimodular FGLM, or None when the ideal is not zero-dimensional or reconstruction gives up."""
    var grevlex = MonomialOrder(MonomialOrder.GREVLEX, order.variables.copy(), budget)
    var reduced = _grevlex_groebner_rational(generators, grevlex, budget, certify=certify)
    if len(reduced) == 0:
        return None
    if len(reduced) == 1 and len(reduced[0].terms) == 1 and len(reduced[0].terms[0].monomial.pairs) == 0:
        # The unit ideal: its LEX basis is the same.
        var result = List[SparsePoly[QQ]]()
        result.append(reduced[0].reorder(order, budget))
        return result^
    if len(order.variables) <= 8:
        return _rational_lex[8](reduced, grevlex, order, certify, budget)
    return _rational_lex[16](reduced, grevlex, order, certify, budget)
