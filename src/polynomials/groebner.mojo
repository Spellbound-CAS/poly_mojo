"""Deterministic Buchberger bases with Gebauer–Möller pair installation.

Pair installation follows SymPy groebnertools.py _buchberger.update.
Set iteration is replaced by stable indices; cofactor propagation is added.

Copyright (c) 2006-2023 SymPy Development Team

All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

  a. Redistributions of source code must retain the above copyright notice,
     this list of conditions and the following disclaimer.
  b. Redistributions in binary form must reproduce the above copyright
     notice, this list of conditions and the following disclaimer in the
     documentation and/or other materials provided with the distribution.
  c. Neither the name of SymPy nor the names of its contributors
     may be used to endorse or promote products derived from this software
     without specific prior written permission.


THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
ARE DISCLAIMED. IN NO EVENT SHALL THE REGENTS OR CONTRIBUTORS BE LIABLE FOR
ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY
OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
DAMAGE.
"""
from std.bit import count_leading_zeros
from apn_mojo import Integer, Rational, div_rem_euclid, div_rem_trunc, gcd as integer_gcd
from .coefficients import CoefficientField, ModInt, Modulus, NoContext, QQ, ZZ, _work
from ._images import _denominator, _scaled_integer
from ._product_stream import _subtract_column
from .errors import Budget, PolyError, PolyErrorKind, _allocation, _charge_size, _product, _sum
from .monomial import Monomial, MonomialOrder
from .sparse import SparsePoly, Term
from .division import _divide_prepared, _ReducerData, divide_by_list
from .certificates import ListDivision, MembershipCertificate
from ._sparse_ops import _subtract_multiples
from ._packed_groebner import _packed_groebner, _packed_groebner_rational, _packed_supported
from ._f4 import _f4_groebner, _f4_supported
from ._f4_trace import _f4_traced, _f4_recorded, _lanes_trace, _Replay, _TracedBasis
from ._f4 import _LANES, _Columns, _raw_hash, _sort_descending, _times
from ._groebner_check import _check_basis, _proves_reduced
from ._packed_reduce import _packed_divide, _packed_divide_rational, _add_multiple, _drain, _pack_monomial
from ._packed_groebner import _Buckets, _Packed, _equal as _packed_equal, _mul as _packed_mul, _packed
from ._packed_groebner import _monomial, _pack as _pack_terms, _unpack, _priorities, _Term as _PackedTerm
from ._fglm import _fglm_groebner, _fglm_rational, _fglm_rational_supported, _fglm_supported
from ._rational_groebner import _grevlex_groebner_rational, _grevlex_groebner_selected, _grevlex_supported, _integer, _next_prime
from ._integer_f4 import _accumulate, _magnitude
from ._qq_trace import _Lifted, _admissible as _qq_admissible, _images as _qq_images


@fieldwise_init
struct TraceTerm[C: CoefficientField](Copyable, Movable, Writable):
    """One multiple in a trace step: a term times an earlier trace polynomial.

    Parameters:
        C: The coefficient field.
    """
    var source: Int
    """The trace polynomial multiplied: the input generators first, then earlier steps."""
    var multiplier: Term[Self.C]
    """The coefficient and monomial it is multiplied by."""


@fieldwise_init
struct GroebnerTrace[C: CoefficientField](Copyable, Movable, Writable):
    """How every basis element arises from the input generators, as a straight-line combination.

    Trace polynomial i, for i below `generators`, is input generator i in
    the basis's monomial order; polynomial `generators` + k is the sum of
    the multiples in step k, which refer only to earlier polynomials. Each
    basis element equals the trace polynomial `basis_sources` names, so it
    lies in the ideal of the generators, without its combination of them
    ever being expanded.

    Parameters:
        C: The coefficient field.
    """
    var generators: Int
    """The number of input generators, the first trace polynomials."""
    var steps: List[List[TraceTerm[Self.C]]]
    """The multiples summed into each later trace polynomial."""
    var basis_sources: List[Int]
    """The trace polynomial equal to each basis element, in basis order."""

    def replay(self, generators: List[SparsePoly[Self.C]], basis: List[SparsePoly[Self.C]], mut budget: Budget) raises PolyError -> Bool:
        """Recompute every trace polynomial and compare the basis with its sources.

        Args:
            generators: The input generators, in trace order.
            basis: The basis the trace describes.
            budget: The shared work allowance.

        Returns:
            Whether every basis element equals its source trace polynomial.

        Raises:
            PolyError: DOMAIN for mixed fields or orders, LIMIT or CAPACITY on work limits.
        """
        if len(generators) != self.generators or len(basis) != len(self.basis_sources):
            return False
        if len(basis) == 0:
            return True
        ref order = basis[0].order
        if _packed_supported(order) and _packable(generators, basis):
            comptime if Self.C == QQ:
                ref rational = rebind[GroebnerTrace[QQ]](self)
                ref inputs = rebind[List[SparsePoly[QQ]]](generators)
                ref elements = rebind[List[SparsePoly[QQ]]](basis)
                if len(order.variables) <= 8:
                    return _rational_replay[8](rational, inputs, elements, budget)
                return _rational_replay[16](rational, inputs, elements, budget)
            comptime if Self.C == ModInt:
                ref modular = rebind[GroebnerTrace[ModInt]](self)
                ref inputs = rebind[List[SparsePoly[ModInt]]](generators)
                ref elements = rebind[List[SparsePoly[ModInt]]](basis)
                if elements[0].ring.value < (UInt64(1) << 31):
                    if len(order.variables) <= 8:
                        return _modular_replay[8](modular, inputs, elements, budget)
                    return _modular_replay[16](modular, inputs, elements, budget)
            if len(order.variables) <= 8:
                return _packed_replay[Self.C, 8](self, generators, basis, budget)
            return _packed_replay[Self.C, 16](self, generators, basis, budget)
        var polys = _trace_polynomials(self, generators, order, len(self.steps), budget)
        for i in range(len(basis)):
            budget.charge(1)
            var source = self.basis_sources[i]
            if source < 0 or source >= len(polys) or not (polys[source] == basis[i]):
                return False
        return True

    def cofactors(self, index: Int, generators: List[SparsePoly[Self.C]], mut budget: Budget) raises PolyError -> List[SparsePoly[Self.C]]:
        """Expand one basis element's combination of the input generators.

        Args:
            index: The basis element.
            generators: The input generators, in trace order; their order is that of the cofactors.
            budget: The shared work allowance.

        Returns:
            One multiplier per input generator, summing with them to the basis element.

        Raises:
            PolyError: DOMAIN for an index outside the basis, LIMIT or CAPACITY on work limits.
        """
        if index < 0 or index >= len(self.basis_sources) or len(generators) != self.generators:
            raise PolyError(PolyErrorKind.DOMAIN, "no such basis element in the trace")
        var target = self.basis_sources[index]
        var s = self.generators
        if len(generators) == 0:
            return List[SparsePoly[Self.C]]()
        ref order = generators[0].order
        ref ring = generators[0].ring
        # Only the polynomials the target depends on are expanded.
        _allocation[Bool](target + 1, budget)
        var needed = List[Bool](length=target + 1, fill=False)
        needed[target] = True
        for k in range(target, s - 1, -1):
            budget.charge(1)
            if not needed[k]:
                continue
            for t in self.steps[k - s]:
                needed[t.source] = True
        var rows = List[List[SparsePoly[Self.C]]](capacity=target + 1)
        for k in range(target + 1):
            budget.charge(1)
            var row = List[SparsePoly[Self.C]]()
            if needed[k]:
                _allocation[SparsePoly[Self.C]](s, budget)
                for g in range(s):
                    if k < s and g == k:
                        row.append(SparsePoly[Self.C].constant(Self.C.one(ring), order, budget))
                    else:
                        row.append(SparsePoly[Self.C].zero(order, ring, budget))
                if k >= s:
                    for t in self.steps[k - s]:
                        for g in range(s):
                            if not rows[t.source][g].is_zero():
                                row[g] = row[g].add(rows[t.source][g]._mul_term(t.multiplier, budget), budget)
            rows.append(row^)
        var result = List[SparsePoly[Self.C]]()
        swap(result, rows[target])
        return result^


def _packable[C: CoefficientField](generators: List[SparsePoly[C]], basis: List[SparsePoly[C]]) -> Bool:
    """Whether every polynomial shares the first basis element's ring, and the basis its order."""
    for g in generators:
        if g.ring != basis[0].ring:
            return False
    for b in basis:
        if b.ring != basis[0].ring or b.order != basis[0].order:
            return False
    return True


def _packed_replay[C: CoefficientField, W: Int](trace: GroebnerTrace[C], generators: List[SparsePoly[C]], basis: List[SparsePoly[C]], mut budget: Budget) raises PolyError -> Bool:
    """`trace.replay` on packed terms, with the same verdicts.

    Each basis element is compared with its source as soon as that is made,
    and each trace polynomial is freed after the last step that reads it.
    """
    ref order = basis[0].order
    var s = trace.generators
    var count = _sum(s, len(trace.steps))
    if s == 0 and len(trace.steps):
        raise PolyError(PolyErrorKind.DOMAIN, "a trace step without generators")
    # The last step reading each polynomial, and the basis elements it makes.
    _allocation[Int](count, budget)
    var last = List[Int](length=count, fill=-1)
    for k in range(len(trace.steps)):
        for t in trace.steps[k]:
            budget.charge(1)
            if t.source < 0 or t.source >= s + k:
                raise PolyError(PolyErrorKind.DOMAIN, "a trace step refers to a later polynomial")
            last[t.source] = k
    _allocation[List[Int]](count, budget)
    var wanted = List[List[Int]](length=count, fill=List[Int]())
    for i in range(len(trace.basis_sources)):
        budget.charge(1)
        var source = trace.basis_sources[i]
        if source < 0 or source >= count:
            return False
        wanted[source].append(i)
    _allocation[List[_PackedTerm[C, W]]](count, budget)
    var polys = List[List[_PackedTerm[C, W]]](capacity=count)
    for g in generators:
        budget.charge(1)
        polys.append(_pack_terms[C, W](g.reorder(order, budget), order, budget))
    for k in range(count):
        if k >= s:
            ref step = trace.steps[k - s]
            var pending = _Buckets[C, W, False](order.kind, basis[0].ring)
            for t in step:
                budget.charge(1)
                _add_multiple(pending, polys[t.source], _pack_monomial[W](t.multiplier.monomial, order, budget), t.multiplier.coefficient, budget)
            polys.append(_drain(pending, budget))
        for i in wanted[k]:
            var expected = _pack_terms[C, W](basis[i], order, budget)
            ref made = polys[k]
            if len(expected) != len(made):
                return False
            for j in range(len(made)):
                budget.charge(1)
                if not _packed_equal(made[j].monomial, expected[j].monomial) or not (made[j].coefficient == expected[j].coefficient):
                    return False
        if k >= s:
            for t in trace.steps[k - s]:
                if last[t.source] == k - s:
                    polys[t.source] = List[_PackedTerm[C, W]]()
    return True


struct _Hashed[W: Int](Movable):
    """A trace polynomial over a word prime for `_modular_replay`: increasing terms, each one's raw hash and exponent bytes, and the greatest exponent of each variable."""
    var terms: List[_PackedTerm[ModInt, Self.W]]
    var hashes: List[UInt64]
    var keys: List[SIMD[DType.uint8, Self.W]]
    var top: _Packed[Self.W]

    def __init__(out self, var terms: List[_PackedTerm[ModInt, Self.W]]):
        self.hashes = List[UInt64](capacity=len(terms))
        self.keys = List[SIMD[DType.uint8, Self.W]](capacity=len(terms))
        var top = SIMD[DType.uint32, Self.W](0)
        for t in terms:
            self.hashes.append(_raw_hash(t.monomial))
            self.keys.append(t.monomial.exponents().cast[DType.uint8]())
            top = max(top, t.monomial.exponents())
        self.top = _packed(top)
        self.terms = terms^


def _modular_replay[W: Int](trace: GroebnerTrace[ModInt], generators: List[SparsePoly[ModInt]], basis: List[SparsePoly[ModInt]], mut budget: Budget) raises PolyError -> Bool:
    """`trace.replay` over a prime below 2^31, each step summed in a hashed dense row, with the same verdicts.

    A step's products are numbered by a hash table of their monomials, as
    an F4 matrix numbers its columns (`_f4._Matrix.add_row`): a product's
    hash is the sum of its factors', and when every exponent stays below 256
    its exponent bytes are the sum of theirs. Residues are summed below p^2
    in one accumulator per monomial; the nonzero sums, sorted once, make the
    step's polynomial. As in `_packed_replay`, each basis element is
    compared with its source as soon as that is made, and each trace
    polynomial is freed after the last step that reads it.
    """
    ref order = basis[0].order
    ref ring = basis[0].ring
    var prime = ring.value
    var square = prime * prime
    var s = trace.generators
    var count = _sum(s, len(trace.steps))
    if s == 0 and len(trace.steps):
        raise PolyError(PolyErrorKind.DOMAIN, "a trace step without generators")
    _allocation[Int](count, budget)
    var last = List[Int](length=count, fill=-1)
    for k in range(len(trace.steps)):
        for t in trace.steps[k]:
            budget.charge(1)
            if t.source < 0 or t.source >= s + k:
                raise PolyError(PolyErrorKind.DOMAIN, "a trace step refers to a later polynomial")
            last[t.source] = k
    _allocation[List[Int]](count, budget)
    var wanted = List[List[Int]](length=count, fill=List[Int]())
    for i in range(len(trace.basis_sources)):
        budget.charge(1)
        var source = trace.basis_sources[i]
        if source < 0 or source >= count:
            return False
        wanted[source].append(i)
    _allocation[_Hashed[W]](count, budget)
    var polys = List[_Hashed[W]](capacity=count)
    for g in generators:
        budget.charge(1)
        polys.append(_Hashed[W](_pack_terms[ModInt, W](g.reorder(order, budget), order, budget)))
    var sums = List[UInt64]()
    for k in range(count):
        if k >= s:
            ref step = trace.steps[k - s]
            var table = _Columns[W]()
            sums.clear()
            for t in step:
                budget.charge(1)
                if t.multiplier.coefficient.modulus != prime:
                    raise PolyError(PolyErrorKind.DOMAIN, "trace coefficients in another field")
                ref source = polys[t.source]
                var n = len(source.terms)
                var c = t.multiplier.coefficient.value
                if n == 0 or c == 0:
                    continue
                var shift = _pack_monomial[W](t.multiplier.monomial, order, budget)
                # Some product's exponent overflows exactly when the greatest exponents' does.
                var top = _packed_mul(source.top, shift)
                budget.charge(UInt64(n))
                table.reserve(n, budget)
                var base = _raw_hash(shift)
                var narrow = top.exponents().reduce_max() < 256
                var offset = shift.exponents().cast[DType.uint8]()
                for u in range(n):
                    var j: Int
                    if narrow:
                        j = Int(table.find_narrow(source.keys.unsafe_get(u) + offset, base + source.hashes.unsafe_get(u)))
                    else:
                        j = Int(table.find(_times(source.terms.unsafe_get(u).monomial, shift), base + source.hashes.unsafe_get(u)))
                    if j == len(sums):
                        sums.append(0)
                    var x = sums.unsafe_get(j) + c * source.terms.unsafe_get(u).coefficient.value
                    sums.unsafe_get(j) = x - square if x >= square else x
            var descending = _sort_descending(order.kind, table.monomials, budget)
            _allocation[_PackedTerm[ModInt, W]](len(descending), budget)
            var terms = List[_PackedTerm[ModInt, W]](capacity=len(descending))
            for i in range(len(descending) - 1, -1, -1):
                var j = descending[i]
                var r = sums[j] % prime
                if r:
                    terms.append(_PackedTerm[ModInt, W](table.monomials[j], ModInt(r, ring)))
            polys.append(_Hashed[W](terms^))
        for i in wanted[k]:
            var expected = _pack_terms[ModInt, W](basis[i], order, budget)
            ref made = polys[k].terms
            if len(expected) != len(made):
                return False
            for j in range(len(made)):
                budget.charge(1)
                if not _packed_equal(made[j].monomial, expected[j].monomial) or not (made[j].coefficient == expected[j].coefficient):
                    return False
        if k >= s:
            for t in trace.steps[k - s]:
                if last[t.source] == k - s:
                    polys[t.source] = _Hashed[W](List[_PackedTerm[ModInt, W]]())
    return True


struct _Scaled[W: Int](Movable):
    """A trace polynomial over QQ as a rational scale times a primitive integer polynomial."""
    var terms: List[_PackedTerm[ZZ, Self.W]]
    var scale: QQ

    def __init__(out self, var terms: List[_PackedTerm[ZZ, Self.W]], var scale: QQ):
        self.terms = terms^
        self.scale = scale^


def _integer_div_rem(a: Integer, b: Integer) raises PolyError -> Tuple[Integer, Integer]:
    """a divided by b, truncated: the quotient and the remainder."""
    try:
        var qr = div_rem_trunc(a, b)
        return (qr[0].copy(), qr[1].copy())
    except e:
        raise PolyError(PolyErrorKind.CAPACITY, String(e))


def _integer_gcd(a: Integer, b: Integer) raises PolyError -> Integer:
    """The nonnegative gcd of a and b."""
    try:
        var g = integer_gcd(a, b)
        return -g if g.sign() < 0 else g^
    except e:
        raise PolyError(PolyErrorKind.CAPACITY, String(e))


def _primitive[W: Int](var terms: List[_PackedTerm[ZZ, W]], var denominator: Integer, mut budget: Budget) raises PolyError -> _Scaled[W]:
    """terms / denominator, as its positive content over the denominator times the primitive part.

    Each coefficient is divided by the content of the ones before it: where
    that divides, the quotient is already its primitive coefficient, and
    only a remainder takes a gcd (with the remainder, the smaller one),
    after which the quotients so far are multiplied by the content's drop.
    The content settles within a few terms, so most coefficients cost one
    division instead of a gcd and a division.
    """
    var content = Integer(0)
    for i in range(len(terms)):
        budget.charge(1)
        ref c = terms[i].coefficient.value
        if content.sign() == 0:
            if c.sign() != 0:
                content = -c if c.sign() < 0 else c.copy()
                terms[i].coefficient = ZZ(Integer(c.sign()))
            continue
        if content == Integer(1):
            break
        _work(terms[i].coefficient.size(), content.magnitude_bit_length() // 8 + 1, True, budget)
        var qr = _integer_div_rem(c, content)
        if qr[1].sign() == 0:
            terms[i].coefficient = ZZ(qr[0].copy())
            continue
        _work(content.magnitude_bit_length() // 8 + 1, qr[1].magnitude_bit_length() // 8 + 1, True, budget)
        var g = _integer_gcd(content, qr[1])
        var quotient = _exact_quotient(c, g, budget)
        var drop = _exact_quotient(content, g, budget)
        for k in range(i):
            budget.charge(1)
            if terms[k].coefficient.value.sign() != 0:
                terms[k].coefficient = ZZ(_integer_product(terms[k].coefficient.value, drop, budget))
        terms[i].coefficient = ZZ(quotient^)
        content = g^
    if content.sign() == 0:
        return _Scaled[W](List[_PackedTerm[ZZ, W]](), QQ.zero(NoContext()))
    var scale: QQ
    try:
        scale = QQ(Rational(content^, denominator^))
    except e:
        raise PolyError(PolyErrorKind.CAPACITY, String(e))
    return _Scaled[W](terms^, scale^)


@always_inline
def _bits_of(n: Int) -> Int:
    """The bit length of a positive n."""
    return 64 - Int(count_leading_zeros(UInt64(n)))


struct _Limbs(Movable):
    """A trace polynomial's integer coefficients as magnitudes in limbs, for `_rational_made`'s sums: term i's are limbs[starts[i] : starts[i + 1]], negated when `negative[i]`; `bits` bounds every magnitude."""
    var limbs: List[UInt64]
    var starts: List[Int]
    var negative: List[Bool]
    var bits: Int
    var sizes: Int
    """The sum over terms of `ZZ.size`, and below of its square, for charging a multiple as `_work` would term by term."""
    var squares: Int

    def __init__(out self):
        self.limbs = List[UInt64]()
        self.starts = [0]
        self.negative = List[Bool]()
        self.bits = 0
        self.sizes = 0
        self.squares = 0

    def __init__[W: Int](out self, terms: List[_PackedTerm[ZZ, W]], mut budget: Budget) raises PolyError:
        self.limbs = List[UInt64]()
        self.starts = List[Int](capacity=len(terms) + 1)
        self.starts.append(0)
        self.negative = List[Bool](capacity=len(terms))
        self.bits = 0
        self.sizes = 0
        self.squares = 0
        for t in terms:
            budget.charge(1)
            _magnitude(t.coefficient.value, self.limbs)
            self.starts.append(len(self.limbs))
            self.negative.append(t.coefficient.value.sign() < 0)
            var bits = t.coefficient.value.magnitude_bit_length()
            self.bits = max(self.bits, bits)
            var size = t.coefficient.size()
            self.sizes = _sum(self.sizes, size)
            self.squares = _sum(self.squares, _product(size, size))


def _signed_integer(data: Pointer[mut=False, UInt64, _], width: Int) -> Integer:
    """The two's complement value of `width` limbs."""
    if data.unsafe_offset(width - 1).unsafe_load() >> 63 == 0:
        return _integer(data, width)
    var magnitude = List[UInt64](length=width, fill=0)
    var carry = UInt64(1)
    for l in range(width):
        var x = ~data.unsafe_offset(l).unsafe_load() + carry
        carry = UInt64(1) if carry == 1 and x == 0 else UInt64(0)
        magnitude[l] = x
    return -_integer(magnitude.unsafe_ptr(), width)


def _exact_quotient(common: Integer, d: Integer, mut budget: Budget) raises PolyError -> Integer:
    """common / d, for d dividing common."""
    _work(common.magnitude_bit_length() // 8 + 1, d.magnitude_bit_length() // 8 + 1, True, budget)
    try:
        return div_rem_euclid(common, d)[0].copy()
    except e:
        raise PolyError(PolyErrorKind.CAPACITY, String(e))


def _integer_product(a: Integer, b: Integer, mut budget: Budget) raises PolyError -> Integer:
    """a * b, charged as `_work` charges a product."""
    _work(a.magnitude_bit_length() // 8 + 1, b.magnitude_bit_length() // 8 + 1, True, budget)
    try:
        return a * b
    except e:
        raise PolyError(PolyErrorKind.CAPACITY, String(e))


def _rational_made[W: Int](trace: GroebnerTrace[QQ], generators: List[SparsePoly[QQ]], order: MonomialOrder, mut budget: Budget) raises PolyError -> Optional[List[_Scaled[W]]]:
    """The trace polynomials `trace.basis_sources` names, in `order`, recomputed over QQ fraction-free; None for a source outside the trace.

    Each trace polynomial is kept as a rational scale times a primitive
    integer polynomial. A step brings its multiples' rational factors, each
    term's coefficient times its source's scale, over one common
    denominator, sums the integer multiples, and divides out the content
    once: no fraction is reduced per term. Each trace polynomial is freed
    after the last step that reads it.
    """
    var s = trace.generators
    var count = _sum(s, len(trace.steps))
    if s == 0 and len(trace.steps):
        raise PolyError(PolyErrorKind.DOMAIN, "a trace step without generators")
    _allocation[Int](count, budget)
    var last = List[Int](length=count, fill=-1)
    for k in range(len(trace.steps)):
        for t in trace.steps[k]:
            budget.charge(1)
            if t.source < 0 or t.source >= s + k:
                raise PolyError(PolyErrorKind.DOMAIN, "a trace step refers to a later polynomial")
            last[t.source] = k
    _allocation[List[Int]](count, budget)
    var wanted = List[List[Int]](length=count, fill=List[Int]())
    for i in range(len(trace.basis_sources)):
        budget.charge(1)
        var source = trace.basis_sources[i]
        if source < 0 or source >= count:
            return None
        wanted[source].append(i)
    _allocation[_Scaled[W]](_sum(count, len(trace.basis_sources)), budget)
    var made = List[_Scaled[W]](capacity=len(trace.basis_sources))
    for _ in range(len(trace.basis_sources)):
        made.append(_Scaled[W](List[_PackedTerm[ZZ, W]](), QQ.zero(NoContext())))
    var polys = List[_Scaled[W]](capacity=count)
    var limbed = List[_Limbs](capacity=count)
    var sums = List[UInt64]()
    for g in generators:
        budget.charge(1)
        var packed = _pack_terms[QQ, W](g.reorder(order, budget), order, budget)
        var common = Integer(1)
        for t in packed:
            budget.charge(1)
            common = _denominator(common, t.coefficient, budget)
        _allocation[_PackedTerm[ZZ, W]](len(packed), budget)
        var terms = List[_PackedTerm[ZZ, W]](capacity=len(packed))
        for t in packed:
            budget.charge(1)
            terms.append(_PackedTerm[ZZ, W](t.monomial, _scaled_integer(t.coefficient, common, budget)))
        polys.append(_primitive[W](terms^, common^, budget))
        limbed.append(_Limbs(polys[len(polys) - 1].terms, budget))
    for k in range(count):
        if k >= s:
            ref step = trace.steps[k - s]
            _allocation[QQ](len(step), budget)
            var factors = List[QQ](capacity=len(step))
            # Each factor's denominator d joins the common one c as it
            # comes: where d divides c, the quotient c / d is kept; otherwise
            # c grows to c / g * d, g = gcd(c, d), and c / g is kept. A
            # factor's quotient into the final common one is then its kept
            # quotient times the final one over c as it was: one division
            # each time c grows, not one per factor.
            _allocation[Integer](len(step), budget)
            var ratios = List[Integer](capacity=len(step))
            var segments = List[Int](capacity=len(step))
            var commons: List[Integer] = [Integer(1)]
            for j in range(len(step)):
                ref t = step[j]
                budget.charge(1)
                var factor = t.multiplier.coefficient.mul(polys[t.source].scale, budget)
                ref d = factor.value.denominator()
                var last = len(commons) - 1
                if j > 0 and d == factors[j - 1].value.denominator():
                    # Consecutive factors often share a denominator.
                    ratios.append(ratios[j - 1].copy())
                elif d == Integer(1):
                    ratios.append(commons[last].copy())
                else:
                    _work(commons[last].magnitude_bit_length() // 8 + 1, d.magnitude_bit_length() // 8 + 1, True, budget)
                    var qr = _integer_div_rem(commons[last], d)
                    if qr[1].sign() == 0:
                        ratios.append(qr[0].copy())
                    else:
                        _work(d.magnitude_bit_length() // 8 + 1, qr[1].magnitude_bit_length() // 8 + 1, True, budget)
                        var g = _integer_gcd(d, qr[1])
                        var ratio = _exact_quotient(commons[last], g, budget)
                        var grown = _integer_product(ratio, d, budget)
                        ratios.append(ratio^)
                        commons.append(grown^)
                segments.append(len(commons) - 1)
                factors.append(factor^)
            var final = len(commons) - 1
            _allocation[Integer](len(commons), budget)
            var scales = List[Integer](capacity=len(commons))
            for i in range(final):
                scales.append(_exact_quotient(commons[final], commons[i], budget))
            scales.append(Integer(1))
            var common = commons[final].copy()
            # The integer multiples summed in a hashed dense row, one
            # two's complement accumulator of `width` limbs per monomial,
            # wide enough for every product and their sum.
            _allocation[List[UInt64]](len(step), budget)
            var multiples = List[List[UInt64]](capacity=len(step))
            var negative = List[Bool](capacity=len(step))
            var bytes = List[Int](capacity=len(step))
            var top = 0
            var products = 1
            for j in range(len(step)):
                ref t = step[j]
                var limbs = List[UInt64]()
                var minus = False
                var size = 0
                if len(polys[t.source].terms) and not factors[j].is_zero():
                    var a = _integer_product(factors[j].value.numerator(), ratios[j], budget)
                    if segments[j] != final:
                        a = _integer_product(a, scales[segments[j]], budget)
                    _magnitude(a, limbs)
                    minus = a.sign() < 0
                    size = a.magnitude_bit_length() // 8 + 1
                    top = max(top, a.magnitude_bit_length() + limbed[t.source].bits)
                    products += len(polys[t.source].terms)
                multiples.append(limbs^)
                negative.append(minus)
                bytes.append(size)
            var width = (top + _bits_of(products) + 1) // 64 + 1
            var table = _Columns[W]()
            sums.clear()
            for j in range(len(step)):
                ref a = multiples[j]
                if len(a) == 0:
                    continue
                ref t = step[j]
                ref source = polys[t.source].terms
                ref limbs = limbed[t.source]
                var shift = _pack_monomial[W](t.multiplier.monomial, order, budget)
                budget.charge(UInt64(len(source)))
                # As `_work` charges each product, (|a| + |c|)^2, summed.
                var size = bytes[j]
                _charge_size(_sum(_sum(_product(len(source), _product(size, size)), _product(2 * size, limbs.sizes)), limbs.squares), budget)
                table.reserve(len(source), budget)
                for i in range(len(source)):
                    var m = _packed_mul(source[i].monomial, shift)
                    var c = Int(table.find(m, _raw_hash(m)))
                    if width * c == len(sums):
                        _allocation[UInt64](_sum(len(sums), width), budget)
                        sums.resize(len(sums) + width, 0)
                    var start = limbs.starts[i]
                    var at = sums.unsafe_ptr().unsafe_offset(width * c)
                    if negative[j] != limbs.negative[i]:
                        _accumulate[True](at, width, a.unsafe_ptr(), len(a), limbs.limbs.unsafe_ptr().unsafe_offset(start), limbs.starts[i + 1] - start)
                    else:
                        _accumulate[False](at, width, a.unsafe_ptr(), len(a), limbs.limbs.unsafe_ptr().unsafe_offset(start), limbs.starts[i + 1] - start)
            var descending = _sort_descending(order.kind, table.monomials, budget)
            _allocation[_PackedTerm[ZZ, W]](len(descending), budget)
            var sum = List[_PackedTerm[ZZ, W]](capacity=len(descending))
            for i in range(len(descending) - 1, -1, -1):
                var c = descending[i]
                var at = sums.unsafe_ptr().unsafe_offset(width * c)
                var zero = True
                for l in range(width):
                    if at.unsafe_offset(l).unsafe_load() != 0:
                        zero = False
                        break
                if not zero:
                    sum.append(_PackedTerm[ZZ, W](table.monomials[c], ZZ(_signed_integer(at, width))))
            polys.append(_primitive[W](sum^, common^, budget))
            limbed.append(_Limbs(polys[len(polys) - 1].terms, budget))
        for i in wanted[k]:
            ref source = polys[k]
            _allocation[_PackedTerm[ZZ, W]](len(source.terms), budget)
            var terms = List[_PackedTerm[ZZ, W]](capacity=len(source.terms))
            for t in source.terms:
                budget.charge(1)
                terms.append(_PackedTerm[ZZ, W](t.monomial, t.coefficient.copy()))
            made[i] = _Scaled[W](terms^, source.scale.copy())
        if k >= s:
            for t in trace.steps[k - s]:
                if last[t.source] == k - s:
                    polys[t.source].terms = List[_PackedTerm[ZZ, W]]()
                    limbed[t.source] = _Limbs()
    return made^


def _rational_replay[W: Int](trace: GroebnerTrace[QQ], generators: List[SparsePoly[QQ]], basis: List[SparsePoly[QQ]], mut budget: Budget) raises PolyError -> Bool:
    """`trace.replay` over QQ, fraction-free (`_rational_made`), with the same verdicts."""
    ref order = basis[0].order
    var made = _rational_made[W](trace, generators, order, budget)
    if not made:
        return False
    for i in range(len(basis)):
        var expected = _pack_terms[QQ, W](basis[i], order, budget)
        ref element = made.value()[i]
        if len(expected) != len(element.terms):
            return False
        for j in range(len(expected)):
            budget.charge(1)
            if not _packed_equal(element.terms[j].monomial, expected[j].monomial):
                return False
            var value = QQ(Rational(element.terms[j].coefficient.value.copy()))
            if not (value.mul(element.scale, budget) == expected[j].coefficient):
                return False
    return True


def _trace_polynomials[C: CoefficientField](trace: GroebnerTrace[C], generators: List[SparsePoly[C]], order: MonomialOrder, count: Int, mut budget: Budget) raises PolyError -> List[SparsePoly[C]]:
    """The input generators in `order`, then the first `count` steps' sums."""
    _allocation[SparsePoly[C]](_sum(len(generators), count), budget)
    var polys = List[SparsePoly[C]](capacity=_sum(len(generators), count))
    for g in generators:
        polys.append(g.reorder(order, budget))
    for k in range(count):
        if len(polys) == 0:
            raise PolyError(PolyErrorKind.DOMAIN, "a trace step without generators")
        var total = SparsePoly[C].zero(order, polys[0].ring, budget)
        for t in trace.steps[k]:
            budget.charge(1)
            if t.source < 0 or t.source >= len(polys):
                raise PolyError(PolyErrorKind.DOMAIN, "a trace step refers to a later polynomial")
            total = total.add(polys[t.source]._mul_term(t.multiplier, budget), budget)
        polys.append(total^)
    return polys^


struct GroebnerBasis[C: CoefficientField](Copyable, Movable, Writable):
    """A reduced monic basis and optional combinations of its input generators.

    Parameters:
        C: The coefficient field.
    """
    var order: MonomialOrder
    """The admissible monomial order used by the basis."""
    var basis: List[SparsePoly[Self.C]]
    """The reduced basis in decreasing leading-monomial order."""
    var cofactors: Optional[List[List[SparsePoly[Self.C]]]]
    """One row per basis element, with columns in input generator order."""
    var trace: Optional[GroebnerTrace[Self.C]]
    """How each basis element arises from the input generators, when traced."""

    def __init__(out self, var order: MonomialOrder, var basis: List[SparsePoly[Self.C]], var cofactors: Optional[List[List[SparsePoly[Self.C]]]]):
        """Hold a basis and optional cofactor rows, without a trace.

        Args:
            order: The admissible monomial order used by the basis.
            basis: The reduced basis in decreasing leading-monomial order.
            cofactors: One row per basis element, with columns in input generator order.
        """
        self.order = order^
        self.basis = basis^
        self.cofactors = cofactors^
        self.trace = None

    def __init__(out self, var order: MonomialOrder, var basis: List[SparsePoly[Self.C]], var cofactors: Optional[List[List[SparsePoly[Self.C]]]], var trace: Optional[GroebnerTrace[Self.C]]):
        """Hold a basis with optional cofactor rows and an optional trace.

        Args:
            order: The admissible monomial order used by the basis.
            basis: The reduced basis in decreasing leading-monomial order.
            cofactors: One row per basis element, with columns in input generator order.
            trace: How each basis element arises from the input generators.
        """
        self.order = order^
        self.basis = basis^
        self.cofactors = cofactors^
        self.trace = trace^


@fieldwise_init
struct _Pair(Movable):
    var first: Int
    var second: Int
    var lcm: Monomial


def _monomial_equal(a: Monomial, b: Monomial, mut budget: Budget) raises PolyError -> Bool:
    _charge_size(_sum(len(a.pairs), len(b.pairs)), budget)
    return a == b


def _sort_indices(mut indices: List[Int], leads: List[Monomial], order: MonomialOrder, mut budget: Budget) raises PolyError:
    for i in range(1, len(indices)):
        budget.charge(1)
        var value = indices[i]
        var j = i
        while j > 0 and order._compare(leads[value], leads[indices[j - 1]], budget) < 0:
            indices[j] = indices[j - 1]
            j -= 1
        indices[j] = value


def _install(leads: List[Monomial], index: Int, mut active: List[Int], mut pairs: List[_Pair], mut budget: Budget) raises PolyError:
    # Gebauer–Möller installation, Becker–Weispfenning p.230 / SymPy update.
    # Stable active order replaces unordered sets. Keep product pairs during
    # divisibility filtering; discard them only after they remove redundant
    # new pairs. Old pairs survive when either replacement has the same LCM.
    ref mh = leads[index]
    _allocation[Monomial](len(active), budget)
    var lcms = List[Monomial](capacity=len(active))
    for old in active:
        budget.charge(1)
        lcms.append(mh.lcm(leads[old], budget))
    _allocation[Int](len(active), budget)
    var kept = List[Int](capacity=len(active))
    for i in range(len(active)):
        budget.charge(1)
        var product = _monomial_equal(mh.mul(leads[active[i]], budget), lcms[i], budget)
        var redundant = False
        if not product:
            for j in range(i + 1, len(active)):
                budget.charge(1)
                if lcms[j].divides(lcms[i], budget):
                    redundant = True
                    break
            if not redundant:
                for j in kept:
                    budget.charge(1)
                    if lcms[j].divides(lcms[i], budget):
                        redundant = True
                        break
        if not redundant:
            kept.append(i)
    _allocation[_Pair](_sum(len(pairs), len(kept)), budget)
    var next_pairs = List[_Pair](capacity=_sum(len(pairs), len(kept)))
    for i in range(len(pairs)):
        budget.charge(1)
        ref pair = pairs[i]
        if not mh.divides(pair.lcm, budget) or _monomial_equal(mh.lcm(leads[pair.first], budget), pair.lcm, budget) or _monomial_equal(mh.lcm(leads[pair.second], budget), pair.lcm, budget):
            var moved = _Pair(0, 0, Monomial())
            swap(moved, pairs[i])
            next_pairs.append(moved^)
    for i in kept:
        budget.charge(1)
        if not _monomial_equal(mh.mul(leads[active[i]], budget), lcms[i], budget):
            var lcm = Monomial()
            swap(lcm, lcms[i])
            next_pairs.append(_Pair(min(index, active[i]), max(index, active[i]), lcm^))
    pairs = next_pairs^
    _allocation[Int](_sum(len(active), 1), budget)
    var next_active = List[Int](capacity=_sum(len(active), 1))
    for old in active:
        budget.charge(1)
        if not mh.divides(leads[old], budget):
            next_active.append(old)
    next_active.append(index)
    active = next_active^


def _pair_before(a: _Pair, b: _Pair, order: MonomialOrder, mut budget: Budget) raises PolyError -> Bool:
    var comparison = order._compare(a.lcm, b.lcm, budget)
    return comparison < 0 or (comparison == 0 and (a.first < b.first or (a.first == b.first and a.second < b.second)))


def _sift_pairs(mut pairs: List[_Pair], start: Int, order: MonomialOrder, mut budget: Budget) raises PolyError:
    var value = _Pair(0, 0, Monomial())
    swap(value, pairs[start])
    var position = start
    while position < len(pairs) // 2:
        budget.charge(1)
        var child = position * 2 + 1
        if child + 1 < len(pairs) and _pair_before(pairs[child + 1], pairs[child], order, budget):
            child += 1
        if not _pair_before(pairs[child], value, order, budget):
            break
        var moved = _Pair(0, 0, Monomial())
        swap(moved, pairs[child])
        swap(moved, pairs[position])
        position = child
    swap(value, pairs[position])


def _heapify_pairs(mut pairs: List[_Pair], order: MonomialOrder, mut budget: Budget) raises PolyError:
    for i in range(len(pairs) // 2 - 1, -1, -1):
        budget.charge(1)
        _sift_pairs(pairs, i, order, budget)


def _select(mut pairs: List[_Pair], order: MonomialOrder, mut budget: Budget) raises PolyError -> _Pair:
    var selected = pairs.pop()
    if len(pairs):
        swap(selected, pairs[0])
        _sift_pairs(pairs, 0, order, budget)
    return selected^


def _zero_row[C: CoefficientField](p: SparsePoly[C], size: Int, mut budget: Budget) raises PolyError -> List[SparsePoly[C]]:
    _allocation[SparsePoly[C]](size, budget)
    var row = List[SparsePoly[C]](capacity=size)
    for _ in range(size):
        budget.charge(1)
        row.append(SparsePoly[C].zero(p.order, p.ring, budget))
    return row^


def _reduce_row[C: CoefficientField](mut row: List[SparsePoly[C]], quotients: List[SparsePoly[C]], indices: List[Int], rows: List[List[SparsePoly[C]]], mut budget: Budget) raises PolyError:
    var terms = 0
    for quotient in quotients:
        budget.charge(1)
        terms = _sum(terms, len(quotient.terms))
    if terms >= 512:
        for k in range(len(row)):
            row[k] = _subtract_column(row[k], quotients, indices, rows, k, budget)
        return
    for j in range(len(indices)):
        budget.charge(1)
        if quotients[j].is_zero():
            continue
        for k in range(len(row)):
            budget.charge(1)
            if not rows[indices[j]][k].is_zero():
                row[k] = row[k].sub(quotients[j].mul(rows[indices[j]][k], budget), budget)


def _normalize[C: CoefficientField](mut p: SparsePoly[C], mut row: List[SparsePoly[C]], mut budget: Budget) raises PolyError:
    var inverse = p.terms[len(p.terms) - 1].coefficient.inverse(budget)
    p = p.scale(inverse, budget)
    for i in range(len(row)):
        budget.charge(1)
        row[i] = row[i].scale(inverse, budget)


def _insert[C: CoefficientField](var p: SparsePoly[C], var row: List[SparsePoly[C]], mut arena: List[SparsePoly[C]], mut rows: List[List[SparsePoly[C]]], mut data: _ReducerData[C], mut active: List[Int], mut pairs: List[_Pair], mut budget: Budget) raises PolyError:
    _normalize(p, row, budget)
    _allocation[SparsePoly[C]](1, budget)
    data.append(p, budget)
    arena.append(p^)
    if len(row):
        _allocation[List[SparsePoly[C]]](1, budget)
        rows.append(row^)
    _install(data.leads, len(arena) - 1, active, pairs, budget)
    _heapify_pairs(pairs, arena[len(arena) - 1].order, budget)


def _record[C: CoefficientField](var step: List[TraceTerm[C]], quotients: List[SparsePoly[C]], indices: List[Int], ids: List[Int], remainder: SparsePoly[C], mut steps: List[List[TraceTerm[C]]], mut budget: Budget) raises PolyError:
    """Finish a trace step: subtract the quotients' multiples, scale by the remainder's monic factor, and append it.

    The step then sums to the remainder made monic, the polynomial `_insert`
    and `_normalize` keep.
    """
    for j in range(len(indices)):
        for t in quotients[j].terms:
            budget.charge(1)
            step.append(TraceTerm[C](ids[indices[j]], Term[C](t.monomial.copy(), t.coefficient.neg(budget))))
    var inverse = remainder.terms[len(remainder.terms) - 1].coefficient.inverse(budget)
    if not (inverse == C.one(remainder.ring)):
        for k in range(len(step)):
            budget.charge(1)
            step[k].multiplier.coefficient = step[k].multiplier.coefficient.mul(inverse, budget)
    steps.append(step^)


def _f4_trace_basis[W: Int](generators: List[SparsePoly[ModInt]], order: MonomialOrder, mut budget: Budget) raises PolyError -> GroebnerBasis[ModInt]:
    """The reduced basis by F4 with the trace of every polynomial the run makes."""
    var ring = generators[0].ring.copy()
    _allocation[List[_PackedTerm[ModInt, W]]](len(generators), budget)
    var packed = List[List[_PackedTerm[ModInt, W]]](capacity=len(generators))
    for g in generators:
        budget.charge(1)
        packed.append(_pack_terms[ModInt, W](g.reorder(order, budget), order, budget))
    var run = _f4_traced[W](packed, order.kind, ring, budget)
    var by_index = _priorities(order, budget)
    _allocation[SparsePoly[ModInt]](len(run.basis), budget)
    var basis = List[SparsePoly[ModInt]](capacity=len(run.basis))
    for p in run.basis:
        budget.charge(1)
        basis.append(_unpack(p, order, ring, by_index, budget))
    ref traced = run.traced
    var count = len(traced.starts) - 1
    _allocation[List[TraceTerm[ModInt]]](count, budget)
    var steps = List[List[TraceTerm[ModInt]]](capacity=count)
    for k in range(count):
        var start = traced.starts[k]
        var end = traced.starts[k + 1]
        _allocation[TraceTerm[ModInt]](end - start, budget)
        var step = List[TraceTerm[ModInt]](capacity=end - start)
        for t in range(start, end):
            budget.charge(1)
            step.append(TraceTerm[ModInt](traced.sources[t], Term[ModInt](_monomial(traced.monomials[t], order, by_index, budget), ModInt(traced.coefficients[t], ring))))
        steps.append(step^)
    var trace = GroebnerTrace[ModInt](len(generators), steps^, run.sources.copy())
    return GroebnerBasis[ModInt](order.copy(), basis^, Optional[List[List[SparsePoly[ModInt]]]](None), Optional[GroebnerTrace[ModInt]](trace^))


comptime _TRACE_PRIMES = 4096
"""Word primes tried for a trace over QQ before the generic traced computation."""
comptime _TRACE_MISMATCHES = 3
"""Consecutive traces of another structure before they replace the first."""


def _qq_f4_trace[W: Int](generators: List[SparsePoly[QQ]], order: MonomialOrder, mut budget: Budget) raises PolyError -> Optional[GroebnerBasis[QQ]]:
    """The reduced basis of the generators over QQ with a trace of it, from traced F4 over word primes, or None.

    The lifted trace's exact replay (`_rational_made`) computes the basis
    elements it names, which therefore lie in the ideal; they are the
    reduced basis once `_proves_reduced` accepts them.

    The first prime's run is recorded, and later primes follow its matrices
    `_LANES` at a time (`_lanes_trace`), each lane giving that prime's trace
    coefficients over the recorded structure; a run too large to record
    leaves each prime to a traced run of its own.
    """
    _allocation[List[_PackedTerm[QQ, W]]](len(generators), budget)
    var packed = List[List[_PackedTerm[QQ, W]]](capacity=len(generators))
    for g in generators:
        budget.charge(1)
        packed.append(_pack_terms[QQ, W](g.reorder(order, budget), order, budget))
    var by_index = _priorities(order, budget)
    # Odd, so that _next_prime starts from 2^31 - 1, the largest prime F4 takes.
    var prime = (UInt64(1) << 31) + 1
    var lifted = Optional[_Lifted[W]](None)
    var replay = _Replay[W]()
    var recorded = Optional[_TracedBasis[W]](None)
    var mismatches = 0
    var failures = 0
    var used = 0
    while used < _TRACE_PRIMES:
        if not lifted or not replay.complete:
            prime = _next_prime(prime, budget)
            while not _qq_admissible(packed, prime, budget):
                prime = _next_prime(prime, budget)
            used += 1
            var ring = Modulus(prime, budget)
            if not lifted:
                var run = _f4_recorded[W](_qq_images(packed, ring, budget), order.kind, ring, replay, budget)
                lifted = _Lifted[W](run, prime)
                recorded = run^
                mismatches = 0
                continue
            var run = _f4_traced[W](_qq_images(packed, ring, budget), order.kind, ring, budget)
            if not lifted.value().agrees(run):
                mismatches += 1
                if mismatches >= _TRACE_MISMATCHES:
                    lifted = _Lifted[W](run, prime)
                    mismatches = 0
                continue
            mismatches = 0
            lifted.value().combine(run, prime, budget)
        else:
            var primes = List[UInt64](capacity=_LANES)
            var rings = List[Modulus](capacity=_LANES)
            var images = List[List[List[_PackedTerm[ModInt, W]]]](capacity=_LANES)
            for _ in range(_LANES):
                prime = _next_prime(prime, budget)
                while not _qq_admissible(packed, prime, budget):
                    prime = _next_prime(prime, budget)
                var ring = Modulus(prime, budget)
                primes.append(prime)
                images.append(_qq_images(packed, ring, budget))
                rings.append(ring^)
            used += _LANES
            var lanes = _lanes_trace[W](replay, recorded.value(), images, order.kind, rings, budget)
            if lanes.live == 0:
                # The recorded run may be the unlucky one: after a few
                # batches that all leave it, the next prime starts over.
                mismatches += 1
                if mismatches >= _TRACE_MISMATCHES:
                    lifted = None
                continue
            mismatches = 0
            lifted.value().combine_lanes(lanes.values, lanes.live, primes, budget)
        if not lifted.value().confirmed():
            continue
        # Every coefficient is a fraction two primes agree on. Replay decides.
        ref structure = lifted.value()
        var count = len(structure.starts) - 1
        _allocation[List[TraceTerm[QQ]]](count, budget)
        var steps = List[List[TraceTerm[QQ]]](capacity=count)
        for k in range(count):
            var start = structure.starts[k]
            var end = structure.starts[k + 1]
            _allocation[TraceTerm[QQ]](end - start, budget)
            var step = List[TraceTerm[QQ]](capacity=end - start)
            for t in range(start, end):
                budget.charge(1)
                step.append(TraceTerm[QQ](structure.sources[t], Term[QQ](_monomial(structure.monomials[t], order, by_index, budget), structure.fractions[t].copy())))
            steps.append(step^)
        var candidate = GroebnerTrace[QQ](len(generators), steps^, structure.basis_sources.copy())
        var made = _rational_made[W](candidate, generators, order, budget)
        if made:
            _allocation[SparsePoly[QQ]](len(made.value()), budget)
            var basis = List[SparsePoly[QQ]](capacity=len(made.value()))
            for e in range(len(made.value())):
                ref element = made.value()[e]
                _allocation[_PackedTerm[QQ, W]](len(element.terms), budget)
                var terms = List[_PackedTerm[QQ, W]](capacity=len(element.terms))
                for t in element.terms:
                    budget.charge(1)
                    terms.append(_PackedTerm[QQ, W](t.monomial, QQ(Rational(t.coefficient.value.copy())).mul(element.scale, budget)))
                basis.append(_unpack(terms, order, NoContext(), by_index, budget))
            if _proves_reduced(generators, basis, order, budget):
                return GroebnerBasis[QQ](order.copy(), basis^, Optional[List[List[SparsePoly[QQ]]]](None), Optional[GroebnerTrace[QQ]](candidate^))
        failures += 1
        if failures >= _TRACE_MISMATCHES:
            return None
        lifted = None
    return None


def groebner[C: CoefficientField](generators: List[SparsePoly[C]], order: MonomialOrder, mut budget: Budget, track_cofactors: Bool = False, certify: Bool = True, trace: Bool = False) raises PolyError -> GroebnerBasis[C]:
    """Compute the reduced monic basis with deterministic Gebauer–Möller pruning.

    Parameters:
        C: The coefficient field.

    Args:
        generators: Polynomials in one field; zero generators are allowed.
        order: LEX, GRLEX, or GREVLEX with every input variable listed.
        budget: The shared work allowance.
        track_cofactors: Whether to retain combinations of the original generators.
        certify: Whether a basis over QQ computed from modular images is certified exactly; without it the result is the reduced basis with high probability, and is not proved.
        trace: Whether to record how each basis element arises from the generators, without expanding cofactor rows.

    Returns:
        The reduced basis, optional cofactor rows in generator order, and an optional trace.

    Raises:
        PolyError: DOMAIN for KERNEL order, mixed fields or malformed input, LIMIT or CAPACITY on work limits.
    """
    if order.kind == MonomialOrder.KERNEL:
        raise PolyError(PolyErrorKind.DOMAIN, "Groebner bases require a monomial order")
    order.validate(Monomial(), budget)
    comptime if C == ModInt:
        if trace and not track_cofactors and len(generators) and _packed_supported(order):
            ref modular = rebind[List[SparsePoly[ModInt]]](generators)
            if _f4_supported(order, modular[0].ring.value):
                for i in range(len(modular)):
                    budget.charge(1)
                    modular[i].check_invariants(budget)
                    if modular[i].ring != modular[0].ring:
                        raise PolyError(PolyErrorKind.DOMAIN, "Groebner generators in different fields")
                var traced = _f4_trace_basis[8](modular, order, budget) if len(order.variables) <= 8 else _f4_trace_basis[16](modular, order, budget)
                return rebind_var[GroebnerBasis[C]](traced^)
    comptime if C == QQ:
        if trace and not track_cofactors and len(generators) and _grevlex_supported(order):
            # The basis comes from a trace, found by traced F4 over word primes
            # and lifted to QQ (`_qq_f4_trace`), and is proved, not trusted.
            # The exact replay over QQ applies the recorded operations to the
            # actual inputs and takes no modular value as given, so every
            # element G it computes lies in the input ideal I. Then (a) every
            # input reduces to zero modulo G, so I is in <G> and <G> = I;
            # (b) every S-pair the Gebauer-Moller criteria keep reduces to
            # zero, so G is a Gröbner basis of I; (d) G is monic and no term
            # but a leading one is divisible by a leading monomial, so G is
            # the reduced basis of I, which is unique: the result is the
            # basis the default mode returns. An unlucky first prime only
            # fails (a), (b) or (d), and the search starts over.
            ref rational = rebind[List[SparsePoly[QQ]]](generators)
            for i in range(len(rational)):
                budget.charge(1)
                rational[i].check_invariants(budget)
            var found = _qq_f4_trace[8](rational, order, budget) if len(order.variables) <= 8 else _qq_f4_trace[16](rational, order, budget)
            if found:
                return rebind_var[GroebnerBasis[C]](found.take())
    if not track_cofactors and not trace and _packed_supported(order):
        # The same result on packed exponent vectors: F4 over small prime
        # fields under degree orders, Buchberger's algorithm otherwise.
        for i in range(len(generators)):
            budget.charge(1)
            generators[i].check_invariants(budget)
            if generators[i].ring != generators[0].ring:
                raise PolyError(PolyErrorKind.DOMAIN, "Groebner generators in different fields")
        comptime if C == QQ:
            # Over QQ under GREVLEX, the homogenized computation by certified
            # reconstruction; otherwise the packed engine fraction-free on
            # integer coefficients.
            ref inputs = rebind[List[SparsePoly[QQ]]](generators)
            if _grevlex_supported(order):
                if certify:
                    # The homogenized proof, unless the homogenized ideal's
                    # first image takes far more work than the input ideal's
                    # (`_grevlex_groebner_selected`); then the basis proved
                    # from a trace, as in certified mode (see the trace branch
                    # above), its witnesses dropped, and the homogenized
                    # proof only if no trace is found.
                    var selected = _grevlex_groebner_selected(inputs, order, budget)
                    if selected:
                        return GroebnerBasis[C](order.copy(), rebind_var[List[SparsePoly[C]]](selected.take()), Optional[List[List[SparsePoly[C]]]](None))
                    var found = _qq_f4_trace[8](inputs, order, budget) if len(order.variables) <= 8 else _qq_f4_trace[16](inputs, order, budget)
                    if found:
                        var proved = found.take()
                        var elements = List[SparsePoly[QQ]]()
                        swap(elements, proved.basis)
                        return GroebnerBasis[C](order.copy(), rebind_var[List[SparsePoly[C]]](elements^), Optional[List[List[SparsePoly[C]]]](None))
                var homogenized = _grevlex_groebner_rational(inputs, order, budget, certify=certify)
                return GroebnerBasis[C](order.copy(), rebind_var[List[SparsePoly[C]]](homogenized^), Optional[List[List[SparsePoly[C]]]](None))
            if _fglm_rational_supported(order):
                # A zero-dimensional ideal's LEX basis from its GREVLEX one.
                var lex = _fglm_rational(inputs, order, certify, budget)
                if lex:
                    return GroebnerBasis[C](order.copy(), rebind_var[List[SparsePoly[C]]](lex.take()), Optional[List[List[SparsePoly[C]]]](None))
            var rational = _packed_groebner_rational(inputs, order, budget)
            return GroebnerBasis[C](order.copy(), rebind_var[List[SparsePoly[C]]](rational^), Optional[List[List[SparsePoly[C]]]](None))
        elif C == ModInt:
            ref modular = rebind[List[SparsePoly[ModInt]]](generators)
            if len(modular) and _fglm_supported(order, modular[0].ring.value):
                # A zero-dimensional ideal's LEX basis from its GREVLEX one.
                var lex = _fglm_groebner(modular, order, budget)
                if lex:
                    return GroebnerBasis[C](order.copy(), rebind_var[List[SparsePoly[C]]](lex.take()), Optional[List[List[SparsePoly[C]]]](None))
            if len(modular) and _f4_supported(order, modular[0].ring.value):
                var f4 = _f4_groebner(modular, order, budget)
                return GroebnerBasis[C](order.copy(), rebind_var[List[SparsePoly[C]]](f4^), Optional[List[List[SparsePoly[C]]]](None))
            return GroebnerBasis[C](order.copy(), _packed_groebner(generators, order, budget), Optional[List[List[SparsePoly[C]]]](None))
        else:
            return GroebnerBasis[C](order.copy(), _packed_groebner(generators, order, budget), Optional[List[List[SparsePoly[C]]]](None))
    var arena = List[SparsePoly[C]]()
    var rows = List[List[SparsePoly[C]]]()
    # With `trace`, each arena polynomial's trace number and the steps so far.
    var carry = track_cofactors or trace
    var ids = List[Int]()
    var steps = List[List[TraceTerm[C]]]()
    var sources = List[Int]()
    var data = _ReducerData[C]()
    var active = List[Int]()
    var pairs = List[_Pair]()
    for i in range(len(generators)):
        budget.charge(1)
        generators[i].check_invariants(budget)
        if generators[i].ring != generators[0].ring:
            raise PolyError(PolyErrorKind.DOMAIN, "Groebner generators in different fields")
        var p = generators[i].reorder(order, budget)
        var row = List[SparsePoly[C]]()
        if track_cofactors:
            row = _zero_row(p, len(generators), budget)
            row[i] = SparsePoly[C].constant(C.one(p.ring), order, budget)
        var reduced = _divide_prepared(p, arena, active, carry, data, budget)
        if reduced.remainder.is_zero():
            continue
        if track_cofactors:
            _reduce_row(row, reduced.quotients, active, rows, budget)
        swap(p, reduced.remainder)
        if trace:
            var step: List[TraceTerm[C]] = [TraceTerm[C](i, Term[C](Monomial(), C.one(p.ring)))]
            _record(step^, reduced.quotients, active, ids, p, steps, budget)
            ids.append(len(generators) + len(steps) - 1)
        _insert(p^, row^, arena, rows, data, active, pairs, budget)
    var dirty_order = True
    while len(pairs):
        budget.charge(1)
        var pair = _select(pairs, order, budget)
        var left = Term[C](pair.lcm.div(data.leads[pair.first], budget), C.one(arena[pair.first].ring))
        var right = Term[C](pair.lcm.div(data.leads[pair.second], budget), C.one(arena[pair.second].ring))
        var s = _subtract_multiples(arena[pair.first], left, arena[pair.second], right, budget)
        if dirty_order:
            _sort_indices(active, data.leads, order, budget)
            dirty_order = False
        var reduced = _divide_prepared(s, arena, active, carry, data, budget)
        if reduced.remainder.is_zero():
            continue
        var row = List[SparsePoly[C]]()
        if track_cofactors:
            row = _zero_row(s, len(generators), budget)
            for k in range(len(generators)):
                budget.charge(1)
                row[k] = _subtract_multiples(rows[pair.first][k], left, rows[pair.second][k], right, budget)
            _reduce_row(row, reduced.quotients, active, rows, budget)
        swap(s, reduced.remainder)
        if trace:
            var step: List[TraceTerm[C]] = [TraceTerm[C](ids[pair.first], left.copy()), TraceTerm[C](ids[pair.second], Term[C](right.monomial.copy(), right.coefficient.neg(budget)))]
            _record(step^, reduced.quotients, active, ids, s, steps, budget)
            ids.append(len(generators) + len(steps) - 1)
        _insert(s^, row^, arena, rows, data, active, pairs, budget)
        dirty_order = True
    _sort_indices(active, data.leads, order, budget)
    _allocation[SparsePoly[C]](len(active), budget)
    var basis = List[SparsePoly[C]](capacity=len(active))
    _allocation[List[SparsePoly[C]]](len(active) if track_cofactors else 0, budget)
    var cofactors = List[List[SparsePoly[C]]]()
    # Active leading monomials are mutually nondividing after installation.
    # Reducing each tail against all the others therefore preserves every
    # leading term, including the equal-leading-monomial regression case.
    for position in range(len(active) - 1, -1, -1):
        budget.charge(1)
        var index = active[position]
        _allocation[Int](len(active) - 1, budget)
        var others = List[Int](capacity=len(active) - 1)
        for j in active:
            budget.charge(1)
            if j != index:
                others.append(j)
        var reduced = _divide_prepared(arena[index], arena, others, carry, data, budget)
        var row = List[SparsePoly[C]]()
        if track_cofactors:
            row = _zero_row(arena[index], len(generators), budget)
            for k in range(len(generators)):
                budget.charge(1)
                row[k] = rows[index][k]._copy(budget)
            _reduce_row(row, reduced.quotients, others, rows, budget)
        var p = SparsePoly[C].zero(order, arena[index].ring, budget)
        swap(p, reduced.remainder)
        if trace:
            var step: List[TraceTerm[C]] = [TraceTerm[C](ids[index], Term[C](Monomial(), C.one(p.ring)))]
            _record(step^, reduced.quotients, others, ids, p, steps, budget)
            sources.append(len(generators) + len(steps) - 1)
        _normalize(p, row, budget)
        basis.append(p^)
        if track_cofactors:
            cofactors.append(row^)
    _charge_size(len(order.variables), budget)
    var evidence = Optional[List[List[SparsePoly[C]]]](None)
    if track_cofactors:
        evidence = cofactors^
    var recorded = Optional[GroebnerTrace[C]](None)
    if trace:
        recorded = GroebnerTrace[C](len(generators), steps^, sources^)
    return GroebnerBasis[C](order.copy(), basis^, evidence^, recorded^)


def reduce[C: CoefficientField](f: SparsePoly[C], basis: GroebnerBasis[C], mut budget: Budget) raises PolyError -> ListDivision[C]:
    """Reduce a polynomial by a Gröbner basis in its chosen monomial order.

    Parameters:
        C: The coefficient field.

    Args:
        f: The polynomial to reduce.
        basis: The basis and its monomial order.
        budget: The shared work allowance.

    Returns:
        Quotients in basis order and the unique normal form.

    Raises:
        PolyError: DOMAIN for incompatible fields or invalid orders, LIMIT or CAPACITY on work limits.
    """
    var dividend = f.reorder(basis.order, budget)
    if _packed_supported(basis.order) and len(basis.basis):
        # The same reductions on packed terms, once the divisors are known
        # to be nonzero and in the dividend's ring and order.
        var packable = True
        for g in basis.basis:
            budget.charge(1)
            if g.is_zero():
                packable = False
                break
            dividend._same(g, budget)
        if packable:
            comptime if C == QQ:
                ref rational = rebind[SparsePoly[QQ]](dividend)
                ref divisors = rebind[List[SparsePoly[QQ]]](basis.basis)
                if len(basis.order.variables) <= 8:
                    return rebind_var[ListDivision[C]](_packed_divide_rational[8](rational, divisors, budget))
                return rebind_var[ListDivision[C]](_packed_divide_rational[16](rational, divisors, budget))
            if len(basis.order.variables) <= 8:
                return _packed_divide[C, 8](dividend, basis.basis, budget)
            return _packed_divide[C, 16](dividend, basis.basis, budget)
    return divide_by_list(dividend, basis.basis, budget)


def membership[C: CoefficientField](f: SparsePoly[C], generators: List[SparsePoly[C]], order: MonomialOrder, mut budget: Budget) raises PolyError -> Optional[MembershipCertificate[C]]:
    """Find a certificate expressing f as a combination of the given generators.

    Parameters:
        C: The coefficient field.

    Args:
        f: The proposed ideal member.
        generators: The ordered original generators.
        order: The admissible monomial order used for the computation.
        budget: The shared work allowance.

    Returns:
        Cofactors in generator order when f belongs to the ideal, otherwise None.

    Raises:
        PolyError: DOMAIN for mixed fields or invalid orders, LIMIT or CAPACITY on work limits.
    """
    f.check_invariants(budget)
    if len(generators) and f.ring != generators[0].ring:
        raise PolyError(PolyErrorKind.DOMAIN, "ideal member and generators in different fields")
    var basis = groebner(generators, order, budget, track_cofactors=True)
    var division = reduce(f, basis, budget)
    if not division.remainder.is_zero():
        return None
    var row = _zero_row(division.remainder, len(generators), budget)
    for j in range(len(basis.basis)):
        budget.charge(1)
        if division.quotients[j].is_zero():
            continue
        for k in range(len(generators)):
            budget.charge(1)
            row[k] = row[k].add(division.quotients[j].mul(basis.cofactors.value()[j][k], budget), budget)
    return MembershipCertificate[C](row^)


def check_groebner_basis[C: CoefficientField](generators: List[SparsePoly[C]], basis: GroebnerBasis[C], mut budget: Budget) raises PolyError -> Bool:
    """Check that the generators reduce to zero modulo a basis and that the basis satisfies Buchberger's criterion.

    Together these give <generators> ⊆ <basis> and make the basis a Gröbner
    basis. The check does not show that the basis lies in the generators'
    ideal (a trace or cofactor rows do that) or that it is reduced.

    Parameters:
        C: The coefficient field.

    Args:
        generators: The polynomials that must reduce to zero.
        basis: The basis and its monomial order.
        budget: The shared work allowance.

    Returns:
        Whether every generator and every S-pair reduces to zero.

    Raises:
        PolyError: DOMAIN for mixed fields or orders, LIMIT or CAPACITY on work limits.
    """
    # On packed terms, with the basis packed once, where that applies.
    var packed = _check_basis(generators, basis.basis, basis.order, budget)
    if packed:
        return packed.value()
    for p in generators:
        if not reduce(p, basis, budget).remainder.is_zero():
            return False
    # Buchberger's criterion, with the two pair criteria that are theorems:
    # coprime leading monomials give an S-polynomial reducing to zero, and a
    # pair (i, j) has a standard representation when some leading monomial
    # m_k divides lcm(m_i, m_j) and both lcm(m_i, m_k) and lcm(m_j, m_k) are
    # proper divisors of it, given that every pair with a smaller lcm does;
    # induction on the lcm then covers every skipped pair.
    ref elements = basis.basis
    _allocation[Monomial](len(elements), budget)
    var leads = List[Monomial](capacity=len(elements))
    for q in elements:
        leads.append(q.leading(budget).monomial.copy())
    for i in range(len(elements)):
        for j in range(i):
            budget.charge(1)
            var a = elements[i].leading(budget)
            var b = elements[j].leading(budget)
            var common = a.monomial.lcm(b.monomial, budget)
            if common.degree(budget) == a.monomial.degree(budget) + b.monomial.degree(budget):
                continue
            var chained = False
            for k in range(len(leads)):
                budget.charge(1)
                if k == i or k == j or not leads[k].divides(common, budget):
                    continue
                if leads[i].lcm(leads[k], budget) != common and leads[j].lcm(leads[k], budget) != common:
                    chained = True
                    break
            if chained:
                continue
            _allocation[Term[C]](2, budget)
            var left = SparsePoly[C].from_terms([Term[C](common.div(a.monomial, budget), a.coefficient.inverse(budget))], basis.order, elements[i].ring, budget)
            var right = SparsePoly[C].from_terms([Term[C](common.div(b.monomial, budget), b.coefficient.inverse(budget))], basis.order, elements[j].ring, budget)
            var s_polynomial = left.mul(elements[i], budget).sub(right.mul(elements[j], budget), budget)
            if not reduce(s_polynomial, basis, budget).remainder.is_zero():
                return False
    return True
