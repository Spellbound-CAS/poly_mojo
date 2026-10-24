"""Budgeted Horner evaluation, polynomial substitution, and formal derivatives."""
from .errors import Budget, PolyError, PolyErrorKind, _charge_size, _allocation, _sum
from .coefficients import CoefficientRing, ZZ, QQ, _pow
from apn_mojo import Integer, Rational, gcd as _gcd_integer
from .monomial import Monomial, MonomialOrder
from .sparse import Term, SparsePoly
from .dense import DensePoly
from ._grouping import _coefficient_groups, _CoefficientGroups
from ._integer_product import _zz_mul
from ._ntt import _integer

def _natural[C: CoefficientRing](n: UInt64, ctx: C.Context, mut budget: Budget) raises PolyError -> C:
    var result = C.zero(ctx)
    var bit = C.one(ctx)
    var value = n
    while value:
        budget.charge(1)
        if value & 1:
            result = result.add(bit, budget)
        value >>= 1
        if value:
            bit = bit.add(bit, budget)
    return result^

@fieldwise_init
struct _HornerTask(Movable):
    var lo: Int
    var hi: Int
    var level: Int
    var exponents: List[UInt32]

def evaluate[C: CoefficientRing](p: SparsePoly[C], point: List[Tuple[Int64, C]], mut budget: Budget) raises PolyError -> C:
    """Evaluate by nested sparse Horner schemes using an explicit work stack.

    Parameters:
        C: The coefficient ring type.

    Args:
        p: The polynomial.
        point: Distinct variable indices and values in the same ring.
        budget: The shared work allowance.

    Returns:
        The exact value in the coefficient ring.

    Raises:
        PolyError: DOMAIN for missing or repeated variables or mixed rings, LIMIT or CAPACITY on work limits.
    """
    var variables = List[Int64]()
    for entry in point:
        budget.charge(1)
        if entry[1].context() != p.ring:
            raise PolyError(PolyErrorKind.DOMAIN, "evaluation point in another ring")
        variables.append(entry[0])
    var order = MonomialOrder(MonomialOrder.LEX, variables^, budget)
    var sorted_poly = p.reorder(order, budget)
    var tasks = List[_HornerTask]()
    var values = List[C]()
    tasks.append(_HornerTask(0, len(sorted_poly.terms), 0, List[UInt32]()))
    while len(tasks):
        budget.charge(1)
        var task = tasks.pop()
        if len(task.exponents):
            var result = values.pop()
            var last = task.exponents[len(task.exponents) - 1]
            for i in range(len(task.exponents) - 1):
                budget.charge(1)
                var e = task.exponents[len(task.exponents) - 2 - i]
                result = result.mul(_pow(point[task.level][1], UInt64(last - e), budget), budget)
                var next_value = values.pop()
                result = result.add(next_value, budget)
                last = e
            if last:
                result = result.mul(_pow(point[task.level][1], UInt64(last), budget), budget)
            values.append(result^)
        elif task.lo == task.hi:
            values.append(C.zero(p.ring))
        elif task.level == len(point):
            _charge_size(sorted_poly.terms[task.lo].coefficient.size(), budget)
            values.append(sorted_poly.terms[task.lo].coefficient.copy())
        else:
            var groups = List[Tuple[Int, Int]]()
            var exponents = List[UInt32]()
            var start = task.lo
            while start < task.hi:
                budget.charge(1)
                var e = sorted_poly.terms[start].monomial.exponent(point[task.level][0], budget)
                var end = start + 1
                while end < task.hi:
                    budget.charge(1)
                    if sorted_poly.terms[end].monomial.exponent(point[task.level][0], budget) != e:
                        break
                    end += 1
                groups.append((start, end))
                exponents.append(e)
                start = end
            tasks.append(_HornerTask(0, 0, task.level, exponents^))
            for i in range(len(groups)):
                budget.charge(1)
                var group = groups[len(groups) - 1 - i]
                tasks.append(_HornerTask(group[0], group[1], task.level + 1, List[UInt32]()))
    return values.pop()

comptime _FRACTION_FREE = 32
"""Integer polynomial length from which evaluation at an integer or a fraction runs by halves."""


def _cancel(numerator: Integer, b: Integer) raises -> Integer:
    """The common factor of a long numerator and a short denominator base."""
    return _gcd_integer(numerator % b, b)


comptime _WORD_BLOCK = 32
"""Coefficients summed per block by word-level Horner before the balanced halves (measured)."""


def _word_blocks(coefficients: List[ZZ], a: ZZ, mut budget: Budget) raises PolyError -> Optional[List[ZZ]]:
    """The sums sum_(i<w) e_(kw+i) a^i over consecutive blocks of w = _WORD_BLOCK coefficients, or None unless a and every e_i are signed words.

    With e_i (-1)^i folded in for negative a, each block runs Horner's rule by
    |a| on a growing run of 64-bit limbs in two's complement, so the bottom
    levels of the halves allocate one integer per block instead of one per
    product and sum.
    """
    var e = List[Int64](capacity=len(coefficients))
    var m: UInt64
    var flip: Bool
    try:
        if a.value.magnitude_bit_length() > 63:
            return None
        var x = a.value.to_native_exact[DType.int64]()
        flip = x < 0
        m = UInt64(-x) if flip else UInt64(x)
        for i in range(len(coefficients)):
            ref c = coefficients[i]
            if c.value.magnitude_bit_length() > 63:
                return None
            var v = c.value.to_native_exact[DType.int64]()
            # e_i a^i = (e_i (-1)^i) |a|^i.
            e.append(-v if flip and i % 2 == 1 else v)
    except err:
        raise PolyError(PolyErrorKind.CAPACITY, String(err))
    var a_bits = a.value.magnitude_bit_length()
    # |block| < 2^63 sum |a|^i < 2^(63 + (w - 1) a_bits + 5), plus a sign bit.
    var limbs = (64 + (_WORD_BLOCK - 1) * a_bits + 5) // 64 + 1
    var count = (len(e) + _WORD_BLOCK - 1) // _WORD_BLOCK
    _allocation[ZZ](count, budget)
    _charge_size(len(e) * limbs, budget)
    var result = List[ZZ](capacity=count)
    var acc = List[UInt64](length=limbs, fill=0)
    var bytes = List[UInt8](length=8 * limbs, fill=0)
    var pa = acc.unsafe_ptr()
    for k in range(count):
        var top = min(len(e), (k + 1) * _WORD_BLOCK)
        # Only the low `used` limbs can differ from the sign fill above them.
        var used = 1
        var fill = UInt64(0)
        pa.unsafe_store(UInt64(0))
        for i in range(top - 1, k * _WORD_BLOCK - 1, -1):
            # acc = acc |a| modulo 2^(64 limbs), on the used limbs and one more.
            var carry = UInt64(0)
            for w in range(used):
                var t = UInt128(pa.unsafe_offset(w).unsafe_load()) * UInt128(m) + UInt128(carry)
                pa.unsafe_offset(w).unsafe_store(UInt64(t))
                carry = UInt64(t >> 64)
            if used < limbs:
                # The sign fill times |a| continues the product above `used`.
                var t = UInt128(fill) * UInt128(m) + UInt128(carry)
                pa.unsafe_offset(used).unsafe_store(UInt64(t))
                used += 1
            # acc += e_i, borrowing or carrying through the used limbs.
            var add = e[i]
            if add >= 0:
                var c = UInt64(add)
                var w = 0
                while c and w < used:
                    var s = pa.unsafe_offset(w).unsafe_load() + c
                    c = UInt64(1) if s < c else UInt64(0)
                    pa.unsafe_offset(w).unsafe_store(s)
                    w += 1
            else:
                var b = UInt64(0) - UInt64(add)
                var w = 0
                while b and w < used:
                    var v = pa.unsafe_offset(w).unsafe_load()
                    pa.unsafe_offset(w).unsafe_store(v - b)
                    b = UInt64(1) if v < b else UInt64(0)
                    w += 1
            fill = UInt64(0) - (pa.unsafe_offset(used - 1).unsafe_load() >> 63)
        var negative = fill != 0
        if negative:
            var c = UInt64(1)
            for w in range(used):
                var v = ~pa.unsafe_offset(w).unsafe_load() + c
                c = UInt64(1) if (c == 1 and v == 0) else UInt64(0)
                pa.unsafe_offset(w).unsafe_store(v)
        while used and pa.unsafe_offset(used - 1).unsafe_load() == 0:
            used -= 1
        result.append(ZZ(_integer(acc, used, negative, bytes)) if used else ZZ(Integer(0)))
    return result^


def _halves(var level: List[ZZ], a: ZZ, b: ZZ, mut budget: Budget) raises PolyError -> ZZ:
    """sum e_i a^i b^(L-1-i) for the L integers e_i in `level`, by balanced halves.

    A block of width w holds the same sum over its coefficients; a full left
    block and its right neighbour of width v combine as
    S_left * b^v + a^w * S_right. Powers of b are skipped when b is one.
    """
    var whole = b.value == Integer(1)
    var power_a = a.copy()
    var power_b = b.copy()
    # Every block has width `width` except the last, whose width is `last`.
    var width = 1
    var last = 1
    while len(level) > 1:
        var next = List[ZZ](capacity=(len(level) + 1) // 2)
        for j in range(0, len(level) - 1, 2):
            budget.charge(1)
            var right = _zz_mul(power_a, level[j + 1], budget)
            if whole:
                next.append(level[j].add(right, budget))
            else:
                var shift = power_b.copy() if j + 2 < len(level) or last == width else _pow(b, UInt64(last), budget)
                next.append(_zz_mul(level[j], shift, budget).add(right, budget))
        if len(level) % 2:
            next.append(level[len(level) - 1].copy())
        else:
            last += width
        level = next^
        width *= 2
        if len(level) > 1:
            power_a = _zz_mul(power_a, power_a, budget)
            if not whole:
                power_b = _zz_mul(power_b, power_b, budget)
    return level[0].copy()


def _evaluate_fraction(p: DensePoly[QQ], x: QQ, mut budget: Budget) raises PolyError -> Optional[QQ]:
    # For integer coefficients e_i and x = a/b in lowest terms, the value is
    # S / b^(L-1) with S = sum e_i a^i b^(L-1-i), summed by _halves. A prime that
    # divides b and S divides the leading coefficient only (a and b are coprime),
    # so cancelling word gcds with b reduces the fraction without a long gcd.
    var count = len(p.coefficients)
    _allocation[ZZ](count, budget)
    var level = List[ZZ](capacity=count)
    for i in range(count):
        budget.charge(1)
        ref c = p.coefficients[i].value
        if c.denominator() != Integer(1):
            return None
        level.append(ZZ(c.numerator()))
    var b = ZZ(x.value.denominator())
    var numerator = _halves(level^, ZZ(x.value.numerator()), b, budget)
    var denominator = _pow(b, UInt64(count - 1), budget)
    if numerator.value.sign() == 0:
        return QQ(Rational(0))
    if b.value != Integer(1):
        while True:
            budget.charge(1)
            _charge_size(_sum(numerator.size(), denominator.size()), budget)
            var common: Integer
            try:
                common = _cancel(numerator.value, b.value)
            except e:
                raise PolyError(PolyErrorKind.CAPACITY, String(e))
            if common == Integer(1):
                break
            try:
                numerator = ZZ(numerator.value // common)
                denominator = ZZ(denominator.value // common)
            except e:
                raise PolyError(PolyErrorKind.CAPACITY, String(e))
    return QQ(Rational(_numerator=numerator.value, _denominator=denominator.value))


def evaluate_dense[C: CoefficientRing](p: DensePoly[C], x: C, mut budget: Budget) raises PolyError -> C:
    """Evaluate a dense polynomial by Horner's rule.

    Parameters:
        C: The coefficient ring type.

    Args:
        p: The polynomial.
        x: The value in the coefficient ring.
        budget: The shared work allowance.

    Returns:
        The exact value.

    Raises:
        PolyError: DOMAIN for another ring, LIMIT or CAPACITY on work limits.
    """
    if x.context() != p.ring:
        raise PolyError(PolyErrorKind.DOMAIN, "evaluation point in another ring")
    comptime if C == ZZ:
        if len(p.coefficients) >= _FRACTION_FREE:
            _allocation[ZZ](len(p.coefficients), budget)
            var level = List[ZZ](capacity=len(p.coefficients))
            for i in range(len(p.coefficients)):
                budget.charge(1)
                level.append(rebind[ZZ](p.coefficients[i]).copy())
            var blocks = _word_blocks(level, rebind[ZZ](x), budget)
            if blocks:
                # Block k weighs (a^w)^k: the halves continue from width w.
                var step = _pow(rebind[ZZ](x), UInt64(_WORD_BLOCK), budget)
                return rebind_var[C](_halves(blocks.take(), step, ZZ(Integer(1)), budget))
            return rebind_var[C](_halves(level^, rebind[ZZ](x), ZZ(Integer(1)), budget))
    comptime if C == QQ:
        if len(p.coefficients) >= _FRACTION_FREE and rebind[QQ](x).value.denominator() != Integer(1):
            var value = _evaluate_fraction(rebind[DensePoly[QQ]](p), rebind[QQ](x), budget)
            if value:
                return rebind_var[C](value.take())
    var result = C.zero(p.ring)
    for i in range(len(p.coefficients)):
        budget.charge(1)
        result = result.mul(x, budget).add(p.coefficients[len(p.coefficients) - 1 - i], budget)
    return result^

def substitute[C: CoefficientRing](p: SparsePoly[C], index: Int64, q: SparsePoly[C], mut budget: Budget) raises PolyError -> SparsePoly[C]:
    """Substitute a polynomial for one variable.

    Parameters:
        C: The coefficient ring type.

    Args:
        p: The source polynomial.
        index: The variable to replace.
        q: The replacement in the same ring and order.
        budget: The shared work allowance.

    Returns:
        The canonical substituted polynomial.

    Raises:
        PolyError: DOMAIN for incompatible operands, LIMIT or CAPACITY on work limits.
    """
    p._same(q, budget)
    if len(p.terms) >= 64:
        var groups = _CoefficientGroups(p, index, budget)
        var count = len(groups.starts) - 1
        var result = groups.coefficient(p, count - 1, budget)
        var last = groups.exponent(count - 1)
        for i in range(count - 2, -1, -1):
            var exponent = groups.exponent(i)
            result = result.mul(q.pow(UInt64(last - exponent), budget), budget).add(groups.coefficient(p, i, budget), budget)
            last = exponent
        if last:
            result = result.mul(q.pow(UInt64(last), budget), budget)
        return result^
    var grouped = _coefficient_groups(p, index, False, budget, first=False)
    var result = SparsePoly[C].zero(p.order, p.ring, budget)
    if len(grouped[0]) == 0:
        return result^
    var last = grouped[0][len(grouped[0]) - 1].exponent(index, budget)
    swap(result, grouped[1][len(grouped[1]) - 1])
    # Sparse Horner steps depend on occupied exponents, never the degree box.
    for i in range(len(grouped[0]) - 2, -1, -1):
        budget.charge(1)
        var exponent = grouped[0][i].exponent(index, budget)
        result = result.mul(q.pow(UInt64(last - exponent), budget), budget).add(grouped[1][i], budget)
        last = exponent
    if last:
        result = result.mul(q.pow(UInt64(last), budget), budget)
    return result^


def derivative[C: CoefficientRing](p: SparsePoly[C], index: Int64, mut budget: Budget) raises PolyError -> SparsePoly[C]:
    """Differentiate formally in one variable in any coefficient ring.

    Parameters:
        C: The coefficient ring type.

    Args:
        p: The polynomial.
        index: The variable to differentiate.
        budget: The shared work allowance.

    Returns:
        The canonical formal derivative.

    Raises:
        PolyError: LIMIT or CAPACITY on work limits.
    """
    _allocation[Term[C]](len(p.terms), budget)
    var terms = List[Term[C]](capacity=len(p.terms))
    for t in p.terms:
        budget.charge(1)
        var exponent = t.monomial.exponent(index, budget)
        if exponent == 0:
            continue
        var pairs = List[Tuple[Int64, UInt32]](capacity=len(t.monomial.pairs))
        for pair in t.monomial.pairs:
            budget.charge(1)
            if pair[0] != index:
                pairs.append(pair)
            elif pair[1] > 1:
                pairs.append((index, pair[1] - 1))
        var c = t.coefficient.mul(_natural[C](UInt64(exponent), p.ring, budget), budget)
        if not c.is_zero():
            terms.append(Term[C](Monomial(pairs^, budget), c^))
    if p.order.kind == MonomialOrder.KERNEL:
        # Storage order compares pair lists, which dropping a variable can reorder.
        return SparsePoly[C].from_terms(terms^, p.order, p.ring, budget)
    # Dividing every monomial by the variable keeps distinct monomials
    # distinct and, monomial orders being multiplicative, in order; dropped
    # terms and zero products leave the rest canonical.
    var result = SparsePoly[C].zero(p.order, p.ring, budget)
    result.terms = terms^
    return result^


def _derivative_dense[C: CoefficientRing](f: DensePoly[C], mut budget: Budget) raises PolyError -> DensePoly[C]:
    _allocation[C](max(f.degree(), 0), budget)
    var values = List[C](capacity=max(f.degree(), 0))
    for i in range(1, len(f.coefficients)):
        budget.charge(1)
        values.append(f.coefficients[i].mul(_natural[C](UInt64(i), f.ring, budget), budget))
    return DensePoly[C].from_coefficients(f.variable, f.ring, values^, budget)
