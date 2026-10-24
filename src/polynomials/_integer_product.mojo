"""Products of large integers by one segmented convolution.

A product of two integers is the product of two one-coefficient
polynomials, so Kronecker segmentation over the float primes multiplies
them in time linear in their bits. Below the crossover the apn product
is faster and runs instead; both give the same integer.
"""
from .coefficient_buffer import _view
from .coefficients import ZZ, _work
from .errors import Budget, PolyError, PolyErrorKind
from ._ntt import _ks_range

comptime _SEGMENTED_WORDS = 250
"""The shorter operand's 64-bit words from which segmentation can win."""
comptime _SEGMENTED_AREA = 680 * 680
"""The product of both word counts from which it does win (measured crossovers)."""


def _zz_mul(a: ZZ, b: ZZ, mut budget: Budget) raises PolyError -> ZZ:
    """a * b, charged as ZZ.mul charges, by segmentation past the crossover and apn's product below it."""
    _work(a.size(), b.size(), True, budget)
    var wa = (a.value.magnitude_bit_length() + 63) // 64
    var wb = (b.value.magnitude_bit_length() + 63) // 64
    if min(wa, wb) >= _SEGMENTED_WORDS and wa * wb >= _SEGMENTED_AREA:
        var square = wa == wb and a == b
        var left = List[ZZ]()
        left.append(a.copy())
        var right = List[ZZ]()
        if not square:
            right.append(b.copy())
        var product = _ks_range[ZZ](_view(left), _view(left), 0, 1, True, budget) if square else _ks_range[ZZ](_view(left), _view(right), 0, 1, False, budget)
        if product:
            return product.value()[0].copy()
    try:
        return ZZ(a.value * b.value)
    except e:
        raise PolyError(PolyErrorKind.CAPACITY, String(e))
