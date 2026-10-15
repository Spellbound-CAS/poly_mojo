"""P2 domain boundaries, degree checks, and borrowed-input preservation."""
from apn_mojo import Integer
from polynomials import Budget, PolyErrorKind, ZZ, QQ, ModInt, Modulus, NoContext, MonomialOrder, DensePoly, SparsePoly, DivisionResult, CancellationCertificate, divmod, pseudo_divmod, divide_exact, divide_by_list, gcd, extended_gcd, resultant, crt, check_division, check_cancellation
from support.check import Checker

def main() raises:
    var c = Checker("P2 boundaries")
    var work = Budget()
    var order = MonomialOrder(MonomialOrder.LEX, [0, 1], work)
    var x = SparsePoly[QQ].variable(0, order, NoContext(), work)
    var one = SparsePoly[QQ].constant(QQ.one(NoContext()), order, work)
    var zero = SparsePoly[QQ].zero(order, NoContext(), work)
    var d = DensePoly[QQ].from_sparse(x, work)
    var zd = DensePoly[QQ](0, NoContext())
    try:
        _ = divmod(d, zd, work)
        c.check(False, "zero field divisor accepted")
    except e:
        c.check(e.kind == PolyErrorKind.DOMAIN, "zero field divisor DOMAIN")
    try:
        _ = pseudo_divmod(d, zd, work)
        c.check(False, "zero pseudo-divisor accepted")
    except e:
        c.check(e.kind == PolyErrorKind.DOMAIN, "zero pseudo-divisor DOMAIN")
    var other = DensePoly[QQ].from_coefficients(1, NoContext(), [QQ.one(NoContext())], work)
    try:
        _ = extended_gcd(d, other, work)
        c.check(False, "mixed dense variables accepted")
    except e:
        c.check(e.kind == PolyErrorKind.DOMAIN, "mixed dense variables DOMAIN")
    var unreduced = DivisionResult[QQ](zd.copy(), d.copy())
    c.check(not check_division(d, d, unreduced, work), "identity alone does not establish remainder degree")
    try:
        _ = divide_by_list(x, [zero.copy()], work)
        c.check(False, "zero divisor in list accepted")
    except e:
        c.check(e.kind == PolyErrorKind.DOMAIN, "zero list divisor DOMAIN")
    var stored = x.reorder(MonomialOrder(), work)
    try:
        _ = divide_by_list(stored, List[SparsePoly[QQ]](), work)
        c.check(False, "KERNEL list order accepted")
    except e:
        c.check(e.kind == PolyErrorKind.DOMAIN, "KERNEL list order DOMAIN")
    c.check(divide_by_list(x, List[SparsePoly[QQ]](), work).remainder == x, "empty divisor list remainder")
    var different = MonomialOrder(MonomialOrder.GREVLEX, [0, 1], work)
    var reordered = x.reorder(different, work)
    try:
        _ = gcd(x, reordered, work)
        c.check(False, "mixed gcd orders accepted")
    except e:
        c.check(e.kind == PolyErrorKind.DOMAIN, "mixed gcd orders DOMAIN")
    var m5 = Modulus(5, work)
    var m7 = Modulus(7, work)
    var f5 = SparsePoly[ModInt].variable(0, order, m5, work)
    var f7 = SparsePoly[ModInt].variable(0, order, m7, work)
    try:
        _ = resultant(f5, f7, 0, work)
        c.check(False, "mixed resultant rings accepted")
    except e:
        c.check(e.kind == PolyErrorKind.DOMAIN, "mixed resultant rings DOMAIN")
    var two = SparsePoly[ZZ].constant(ZZ(Integer(2)), order, work)
    var three = SparsePoly[ZZ].constant(ZZ(Integer(3)), order, work)
    c.check(not divide_exact(two, three, work), "nonintegral integer quotient refused")
    c.check(gcd(two, three, work).gcd == SparsePoly[ZZ].constant(ZZ(Integer(1)), order, work), "integer constant gcd")
    var invalid_cancel = CancellationCertificate(one.copy(), (one.copy(), zero.copy()))
    c.check(not check_cancellation(one, zero, invalid_cancel, work), "zero-denominator cancellation rejected")
    try:
        _ = crt(List[UInt64](), [UInt64(7)], work)
        c.check(False, "mismatched CRT lengths accepted")
    except e:
        c.check(e.kind == PolyErrorKind.DOMAIN, "mismatched CRT lengths DOMAIN")
    var before = x.pow(6, work).sub(one, work)
    var divisor = x.pow(2, work).sub(one, work)
    var saved = before.copy()
    var saved_divisor = divisor.copy()
    var limits: List[UInt64] = [0, 1, 20, 100, 1000]
    for limit in limits:
        var allowance = Budget(limit)
        try:
            _ = gcd(before, divisor, allowance, certify_maximal=True)
        except e:
            c.check(e.kind == PolyErrorKind.LIMIT, "gcd finite budget has typed exhaustion")
        c.check(before == saved and divisor == saved_divisor, "failed gcd leaves borrowed inputs unchanged")
    var constant = SparsePoly[QQ].constant(QQ.one(NoContext()).add(QQ.one(NoContext()), work), order, work)
    c.check(resultant(x.pow(3, work), constant, 0, work).resultant == constant.pow(3, work), "constant resultant power")
    c.check(resultant(zero, constant, 0, work).resultant.is_zero(), "zero resultant convention")
    c.finish()
