"""Independent prime table, CRT integers, and cache recovery after exhaustion."""
from apn_mojo import Integer, div_rem_euclid
from polynomials import Budget, PolyErrorKind
from polynomials.primes import prime, crt
from support.check import Checker
from support.corpus import read_cases

def main() raises:
    var c = Checker("primes")
    var work = Budget()
    var interrupted = Budget(80)
    try:
        _ = prime(0, interrupted)
        c.check(False, "prime ignored exhausted budget")
    except e:
        c.check(e.kind == PolyErrorKind.LIMIT, "prime budget exhaustion")
    var count = 0
    for record in read_cases():
        if record.field("op ") != "primes":
            continue
        for line in record.lines:
            if line.startswith("out prime "):
                var expected = UInt64(Int(String(line[byte=10:])))
                c.check(prime(count, work) == expected, String("prime ", count))
                count += 1
    c.check(count == 100, "all independent prime table entries")
    for i in range(count):
        var small = Budget(200)
        c.check(prime(i, small) == prime(i, work), "cached prime uses bounded work")
    var moduli: List[UInt64] = [prime(0, work), prime(1, work), prime(2, work)]
    var large = Integer(1) << 150
    for i in range(-8, 9):
        var value = large * Integer(i) + Integer(71)
        var residues = List[UInt64]()
        for p in moduli:
            residues.append(UInt64(Int(div_rem_euclid(value, Integer(Int(p)))[1])))
        c.check(crt(residues, moduli, work) == value, "CRT exact signed large integer")
    var invalid: List[List[UInt64]] = [List[UInt64](), [4], [17, 17], [0], [1], [UInt64(1) << 62]]
    for values in invalid:
        var residues = List[UInt64](length=len(values), fill=0)
        try:
            _ = crt(residues, values, work)
            c.check(False, "invalid CRT moduli accepted")
        except e:
            c.check(e.kind == PolyErrorKind.DOMAIN, "invalid CRT DOMAIN")
    try:
        _ = prime(-1, work)
        c.check(False, "negative prime index accepted")
    except e:
        c.check(e.kind == PolyErrorKind.DOMAIN, "negative prime index DOMAIN")
    c.check(crt([UInt64(100)], [UInt64(7)], work) == Integer(2), "CRT reduces residues")
    c.check(crt([UInt64(1)], [UInt64(2)], work) == Integer(1), "CRT positive midpoint")
    c.finish()
