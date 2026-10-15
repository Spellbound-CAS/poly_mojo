"""Independent multiplicity blocks, product certificates, and squarefreeness."""
from apn_mojo import Rational
from polynomials import Budget, PolyErrorKind, QQ, ModInt, Modulus, NoContext, DensePoly, MonomialOrder, squarefree, check_squarefree_product
from polynomials._dense_ops import _gcd, _derivative
from support.corpus import read_cases, parse_polynomial
from support.check import Checker
from support.gen import coefficient


def main() raises:
    var work = Budget()
    var c = Checker("squarefree")
    var count = 0
    for record in read_cases():
        if record.field("op ") != "squarefree":
            continue
        var input_text = record.field("in f ")
        var f = DensePoly[QQ].from_sparse(parse_polynomial(input_text, MonomialOrder(), work), work)
        var result = squarefree(f, work)
        var unit_text = record.field("out unit ")
        var unit = QQ(Rational(unit_text))
        c.check(result.unit == unit, record.name)
        var i = 0
        for line in record.lines:
            if not line.startswith("out factor "):
                continue
            var parts = line.split(" ")
            var exponent = Int(String(parts[2]))
            var text = String(line[byte=12 + String(exponent).byte_length():])
            var expected = DensePoly[QQ].from_sparse(parse_polynomial(text, MonomialOrder(), work), work)
            c.check(i < len(result.factors), "enough blocks")
            c.check(result.factors[i][1] == exponent and result.factors[i][0] == expected, "corpus block")
            i += 1
        c.check(i == len(result.factors), "no extra blocks")
        c.check(check_squarefree_product(f, result, work), "product certificate")
        for i in range(len(result.factors)):
            ref block = result.factors[i][0]
            c.check(_gcd(block, _derivative(block, work), work).degree() == 0, "squarefree block")
            for j in range(i):
                c.check(_gcd(block, result.factors[j][0], work).degree() == 0, "pairwise coprime blocks")
        result.unit = result.unit.add(QQ.one(NoContext()), work)
        c.check(not check_squarefree_product(f, result, work), "corrupt product")
        count += 1
    c.check(count == 9, "all nine squarefree cases")
    for degree in range(1, 9):
        var xminus = DensePoly[QQ].from_coefficients(7, NoContext(), [coefficient[QQ](NoContext(), -degree, work), QQ.one(NoContext())], work)
        var f = xminus.pow(UInt64(degree), work)
        var result = squarefree(f, work)
        c.check(len(result.factors) == 1 and result.factors[0][0] == xminus and result.factors[0][1] == degree, "planted multiplicity")
    var zero = DensePoly[QQ](0, NoContext())
    c.check(check_squarefree_product(zero, squarefree(zero, work), work), "zero decomposition")
    var modulus = Modulus(2, work)
    var finite = DensePoly[ModInt](0, modulus)
    c.check(check_squarefree_product(finite, squarefree(finite, work), work), "finite-field zero decomposition")
    var one = DensePoly[QQ].from_coefficients(0, NoContext(), [QQ.one(NoContext())], work)
    var tiny = Budget(0)
    try:
        _ = squarefree(one, tiny)
        c.check(False, "budget ignored")
    except e:
        c.check(e.kind == PolyErrorKind.LIMIT, "squarefree budget")
    c.finish()
