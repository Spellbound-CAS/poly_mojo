"""The kernel's canonical polynomial byte encoding, with bounded parsing.

Magnitude conversion uses explicit-stack binary splitting over apn_mojo's
public API.
"""
from apn_mojo import Integer, Rational
from .errors import Budget, PolyError, PolyErrorKind, _sum, _product, _charge_size
from .coefficients import QQ, NoContext, _work
from .monomial import Monomial, MonomialOrder
from .sparse import Term, SparsePoly

comptime MAX_KERNEL_BYTES = 1048576
"""The kernel's one-mebibyte literal size ceiling."""

def _bad(message: String) -> PolyError:
    return PolyError(PolyErrorKind.ENCODING, message)

def _length(value: Integer) -> Int:
    return value.magnitude_bit_length() // 8 + Int(value.magnitude_bit_length() % 8 != 0)

def _write_word(mut output: List[UInt8], value: UInt64, width: Int, mut budget: Budget) raises PolyError:
    _charge_size(width, budget)
    for i in range(width):
        output.append(UInt8((value >> UInt64(8 * i)) & 255))

def _read_word(data: Span[UInt8, _], mut position: Int, width: Int, mut budget: Budget) raises PolyError -> UInt64:
    if width > len(data) - position:
        raise _bad("truncated field")
    _charge_size(width, budget)
    var result = UInt64(0)
    for i in range(width):
        result |= UInt64(data[position + i]) << UInt64(8 * i)
    position += width
    return result

def _magnitude(value: Integer, mut budget: Budget) raises PolyError -> List[UInt8]:
    var n = _length(value)
    _charge_size(n, budget)
    var output = List[UInt8](length=n, fill=0)
    var values = List[Integer]()
    var offsets = List[Int]()
    var counts = List[Int]()
    values.append(-value if value.sign() < 0 else value.copy())
    offsets.append(0)
    counts.append(n)
    while len(values):
        budget.charge(1)
        var v = values.pop()
        var offset = offsets.pop()
        var count = counts.pop()
        _charge_size(count, budget)
        if count <= 7:
            try:
                var x = UInt64(Int(v))
                for i in range(count):
                    output[offset + i] = UInt8((x >> UInt64(8 * i)) & 255)
            except e:
                raise PolyError(PolyErrorKind.CAPACITY, String(e))
        else:
            var half = count // 2
            var high: Integer
            var low: Integer
            try:
                high = v >> (8 * half)
                low = v - (high << (8 * half))
            except e:
                raise PolyError(PolyErrorKind.CAPACITY, String(e))
            values.append(high^)
            offsets.append(offset + half)
            counts.append(count - half)
            values.append(low^)
            offsets.append(offset)
            counts.append(half)
    return output^

def _from_magnitude(data: Span[UInt8, _], mut budget: Budget) raises PolyError -> Integer:
    if len(data) == 0:
        return Integer(0)
    if data[len(data) - 1] == 0:
        raise _bad("padded integer magnitude")
    var level = List[Integer]()
    var position = 0
    while position < len(data):
        budget.charge(1)
        var count = min(7, len(data) - position)
        var value = _read_word(data, position, count, budget)
        level.append(Integer(value))
    var width = 7
    while len(level) > 1:
        var next_level = List[Integer]()
        var i = 0
        while i + 1 < len(level):
            _charge_size(_product(width, 2), budget)
            try:
                next_level.append(level[i] + (level[i + 1] << (8 * width)))
            except e:
                raise PolyError(PolyErrorKind.CAPACITY, String(e))
            i += 2
        if i < len(level):
            _charge_size(_length(level[i]), budget)
            next_level.append(level[i].copy())
        level = next_level^
        width = _product(width, 2)
    return level.pop()

def _rational_bytes(q: Rational, mut budget: Budget) raises PolyError -> List[UInt8]:
    var numerator = q.numerator()
    var num = _magnitude(numerator, budget)
    var den = _magnitude(q.denominator(), budget)
    var output = List[UInt8]()
    _write_word(output, UInt64(len(num) + 1), 4, budget)
    budget.charge(1)
    output.append(UInt8(numerator.sign() < 0))
    for byte in num:
        budget.charge(1)
        output.append(byte)
    for byte in den:
        budget.charge(1)
        output.append(byte)
    return output^

def _read_rational(data: Span[UInt8, _], mut budget: Budget) raises PolyError -> QQ:
    var position = 0
    var length = Int(_read_word(data, position, 4, budget))
    if length < 1 or length >= len(data) - position:
        raise _bad("invalid rational numerator length")
    budget.charge(1)
    var sign = data[position]
    if sign > 1 or (sign == 1 and length == 1):
        raise _bad("invalid sign or negative zero")
    var numerator = _from_magnitude(data[position + 1:position + length], budget)
    if sign:
        numerator = -numerator
    var denominator = _from_magnitude(data[position + length:], budget)
    if denominator.sign() <= 0:
        raise _bad("zero denominator")
    _work(_length(numerator) + 1, _length(denominator) + 1, True, budget)
    var rational: Rational
    try:
        rational = Rational(numerator, denominator)
    except e:
        raise PolyError(PolyErrorKind.CAPACITY, String(e))
    _charge_size(len(data), budget)
    if rational.numerator() != numerator or rational.denominator() != denominator:
        raise _bad("unreduced rational")
    return QQ(rational^)

def to_kernel_bytes(p: SparsePoly[QQ], mut budget: Budget) raises PolyError -> List[UInt8]:
    """Encode a polynomial as a canonical kernel literal.

    Args:
        p: The rational polynomial, in any monomial order.
        budget: The shared work allowance.

    Returns:
        The canonical little-endian payload.

    Raises:
        PolyError: ENCODING for zero or one, DOMAIN for broken invariants, CAPACITY beyond the kernel ceiling, LIMIT on exhaustion.
    """
    p.check_invariants(budget)
    if p.is_zero() or (len(p.terms) == 1 and len(p.terms[0].monomial.pairs) == 0 and p.terms[0].coefficient.value == Rational(1)):
        raise _bad("zero and one are not kernel literal atoms")
    var size = 4
    for term in p.terms:
        budget.charge(1)
        var coefficient_size = _sum(5, _sum(_length(term.coefficient.value.numerator()), _length(term.coefficient.value.denominator())))
        size = _sum(size, _sum(8, _sum(_product(len(term.monomial.pairs), 12), coefficient_size)))
        if size > MAX_KERNEL_BYTES:
            raise PolyError(PolyErrorKind.CAPACITY, "kernel literal exceeds size ceiling")
    if p.order.kind == MonomialOrder.KERNEL:
        return _encode_canonical(p, size, budget)
    var canonical = p.reorder(MonomialOrder(), budget)
    return _encode_canonical(canonical, size, budget)

def _encode_canonical(canonical: SparsePoly[QQ], size: Int, mut budget: Budget) raises PolyError -> List[UInt8]:
    _charge_size(size, budget)
    var output = List[UInt8](capacity=size)
    _write_word(output, UInt64(len(canonical.terms)), 4, budget)
    for term in canonical.terms:
        budget.charge(1)
        _write_word(output, UInt64(len(term.monomial.pairs)), 4, budget)
        for pair in term.monomial.pairs:
            budget.charge(1)
            _write_word(output, UInt64(pair[0]), 8, budget)
            _write_word(output, UInt64(pair[1]), 4, budget)
        var coefficient = _rational_bytes(term.coefficient.value, budget)
        _write_word(output, UInt64(len(coefficient)), 4, budget)
        for byte in coefficient:
            budget.charge(1)
            output.append(byte)
    return output^

def from_kernel_bytes(data: Span[UInt8, _], mut budget: Budget) raises PolyError -> SparsePoly[QQ]:
    """Decode canonical kernel bytes without repairing malformed input.

    Args:
        data: The byte payload.
        budget: The shared work allowance.

    Returns:
        The canonical polynomial in KERNEL order.

    Raises:
        PolyError: ENCODING for malformed or noncanonical bytes, CAPACITY beyond the size ceiling, LIMIT on exhaustion.
    """
    if len(data) > MAX_KERNEL_BYTES:
        raise PolyError(PolyErrorKind.CAPACITY, "kernel literal exceeds size ceiling")
    var position = 0
    var count = Int(_read_word(data, position, 4, budget))
    if count == 0 or count > (len(data) - position) // 15:
        raise _bad("invalid term count")
    var result = SparsePoly[QQ](MonomialOrder(), NoContext())
    for _ in range(count):
        budget.charge(1)
        var pair_count = Int(_read_word(data, position, 4, budget))
        if pair_count > (len(data) - position) // 12:
            raise _bad("invalid monomial pair count")
        var monomial = Monomial()
        for i in range(pair_count):
            budget.charge(1)
            var index = Int64(_read_word(data, position, 8, budget))
            var exponent = UInt32(_read_word(data, position, 4, budget))
            if exponent == 0 or (i > 0 and index <= monomial.pairs[i - 1][0]):
                raise _bad("noncanonical monomial")
            monomial.pairs.append((index, exponent))
        var length = Int(_read_word(data, position, 4, budget))
        if length > len(data) - position:
            raise _bad("truncated coefficient")
        var c = _read_rational(data[position:position + length], budget)
        position += length
        if c.is_zero():
            raise _bad("zero coefficient")
        if len(result.terms) and result.order._compare(result.terms[len(result.terms) - 1].monomial, monomial, budget) >= 0:
            raise _bad("nonincreasing terms")
        result.terms.append(Term[QQ](monomial^, c^))
    if position != len(data):
        raise _bad("trailing bytes")
    if len(result.terms) == 1 and len(result.terms[0].monomial.pairs) == 0 and result.terms[0].coefficient.value == Rational(1):
        raise _bad("one is not a kernel literal atom")
    return result^
