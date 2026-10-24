"""Paired native-integer transforms for bounded CAD interval coefficients.

Each coefficient stores the independent lower and upper transforms in two
Int64 lanes. Callers establish the signed overflow bound before entering the
kernel and charge the complete coefficient work. No precision is discarded.
"""


def _shift_native[vectorized: Bool = True](mut coefficients: List[SIMD[DType.int64, 2]]):
    """Replace both coefficient vectors by their Taylor shifts by one.

    The inner recurrence must run from high to low degree. The two interval
    bounds are independent, so SIMD acts across bounds rather than across
    dependent coefficients. The scalar specialization uses the same layout
    and visits each bound separately as a measurement/reference path.
    """
    var n = len(coefficients) - 1
    comptime if vectorized:
        for i in range(n):
            for j in range(n - 1, i - 1, -1):
                coefficients[j] += coefficients[j + 1]
    else:
        for lane in range(2):
            for i in range(n):
                for j in range(n - 1, i - 1, -1):
                    coefficients[j][lane] += coefficients[j + 1][lane]
