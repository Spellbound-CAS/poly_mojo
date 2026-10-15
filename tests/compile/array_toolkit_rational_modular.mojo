# output: GF(7)
from poly import rational
from polynomials import ModInt, Modulus, Budget
from std.testing import assert_equal, assert_true

def main() raises:
    var budget = Budget()
    var field = Modulus(7, budget)
    var value = rational.asarray([ModInt(1, field)], [1])
    assert_equal(value.modulus(), UInt64(7))
    assert_true(rational.check(value))
    assert_equal(rational.evaluate(value, ModInt(0, field)).value, UInt64(1))
    print(value.domain())
