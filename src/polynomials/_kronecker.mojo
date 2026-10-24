"""Products of short integer polynomials with large coefficients by one integer product.

Signed coefficients c_i of magnitude below 2^(32 w - 1) pack as N = P - Q,
with P holding the positive coefficients and Q the magnitudes of the
negative ones in slots of w words. One product of the packed integers
carries every coefficient of the polynomial product in its slot; signed
digits come back by borrowing from the slot above whenever a slot's top
bit is set.
"""
from apn_mojo import Integer
from .coefficients import ZZ, NoContext
from .dense import DensePoly
from .errors import Budget, PolyError, PolyErrorKind, _allocation, _charge_size, _product
from ._integer_product import _zz_mul


def _bits(a: DensePoly[ZZ]) -> Int:
    var top = 0
    for i in range(len(a.coefficients)):
        top = max(top, a.coefficients[i].value.magnitude_bit_length())
    return top


def _pack(a: DensePoly[ZZ], w: Int) -> Integer:
    var n = len(a.coefficients)
    var positive = List[UInt32](length=n * w, fill=0)
    var negative = List[UInt32](length=n * w, fill=0)
    var any_negative = False
    for i in range(n):
        ref c = a.coefficients[i].value
        var count = c._word_count()
        if c.sign() < 0:
            any_negative = True
            for j in range(count):
                negative[i * w + j] = c._word(j)
        else:
            for j in range(count):
                positive[i * w + j] = c._word(j)
    var p = Integer._from_words(positive^, False)
    if not any_negative:
        return p^
    try:
        return p - Integer._from_words(negative^, False)
    except:
        # Subtraction of canonical integers does not fail.
        return Integer(0)


def _kronecker_mul(a: DensePoly[ZZ], b: DensePoly[ZZ], mut budget: Budget) raises PolyError -> DensePoly[ZZ]:
    """a * b for nonzero integer polynomials, by packing both into integers."""
    var la = len(a.coefficients)
    var lb = len(b.coefficients)
    var length = la + lb - 1
    # Every product coefficient is below min(la, lb) 2^(bits a + bits b) in magnitude.
    var span = _bits(a) + _bits(b) + Integer(min(la, lb)).magnitude_bit_length() + 1
    var w = span // 32 + 1
    _charge_size(_product(length, w), budget)
    _allocation[UInt32](_product(length + 1, w), budget)
    var product = _zz_mul(ZZ(_pack(a, w)), ZZ(_pack(b, w)), budget).value
    var flip = product.sign() < 0
    var count = product._word_count()
    _allocation[ZZ](length, budget)
    var coefficients = List[ZZ](capacity=length)
    var slot = List[UInt32](length=w, fill=0)
    var carry = UInt64(0)
    for i in range(length):
        # The slot's digit plus the borrow carried from the slot below.
        var c = carry
        for j in range(w):
            var k = i * w + j
            c += UInt64(product._word(k)) if k < count else 0
            slot[j] = UInt32(c & 0xFFFFFFFF)
            c >>= 32
        # A borrow that overflows the slot moves on to the next one.
        carry = c
        var negative = (slot[w - 1] >> 31) != 0
        if negative:
            # The digit is slot - 2^(32 w): its magnitude is the two's complement.
            var t = UInt64(1)
            for j in range(w):
                t += UInt64(~slot[j])
                slot[j] = UInt32(t & 0xFFFFFFFF)
                t >>= 32
            carry += 1
        var words = List[UInt32](length=w, fill=0)
        for j in range(w):
            words[j] = slot[j]
        coefficients.append(ZZ(Integer._from_words(words^, negative != flip)))
    return DensePoly[ZZ].from_coefficients(a.variable, NoContext(), coefficients^, budget)
