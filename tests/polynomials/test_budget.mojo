"""Budget boundaries and failure atomicity."""
from polynomials import Budget, PolyErrorKind
from support.check import Checker

def main() raises:
    var c = Checker("budget")
    var b = Budget(3)
    b.charge(0)
    b.charge(2)
    c.check(b.remaining() == 1, "successful charges")
    try:
        b.charge(2)
        c.check(False, "must exhaust")
    except e:
        c.check(e.kind == PolyErrorKind.LIMIT, "typed LIMIT")
    c.check(b.remaining() == 1, "failed charge unchanged")
    b.charge(1)
    c.check(b.remaining() == 0, "exact exhaustion")
    var full = Budget()
    full.charge(UInt64.MAX)
    c.check(full.remaining() == 0, "full unsigned range")
    c.finish()
