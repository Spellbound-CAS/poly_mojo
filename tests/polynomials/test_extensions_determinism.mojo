"""Deterministic replay of the additive values and complete certificate structures."""
from apn_mojo import Rational
from polynomials import Budget, QQ, DensePoly, SparsePoly, RealAlgebraic, MonomialOrder, real_root, algebraic_add, sparse_gcd, factor_multivariate, cad_cell
from support.corpus import parse_polynomial
from support.check import Checker
from polynomials.cylindrical_covering import cylindrical_covering


def main() raises:
    var work = Budget()
    var c = Checker("extension replay")
    var order = MonomialOrder()
    var x = parse_polynomial("x0", order, work)
    var a = parse_polynomial("x0^2 - 2", order, work)
    var root = real_root(DensePoly[QQ].from_sparse(a, work), QQ(Rational(1)), QQ(Rational(2)), work)
    var first = algebraic_add(root, root, work)
    var second = algebraic_add(root, root, work)
    c.check(String(first) == String(second), "complete arithmetic proof replay")
    var common = parse_polynomial("x0*x1 + 1", order, work)
    var left = common.mul(x, work)
    var right = common.mul(parse_polynomial("x0 + 1", order, work), work)
    c.check(String(sparse_gcd(left, right, work)) == String(sparse_gcd(left, right, work)), "complete sparse gcd proof replay")
    var product = common.pow(2, work).mul(parse_polynomial("x1^2 + x0 + 1", order, work), work)
    c.check(String(factor_multivariate(product, work)) == String(factor_multivariate(product, work)), "multivariate factors and ordering replay")
    var inputs: List[SparsePoly[QQ]] = [a.copy()]
    var point: List[RealAlgebraic] = [root.copy()]
    var cell1 = cad_cell(inputs, [Int64(0)], point, work, lazard=True)
    var cell2 = cad_cell(inputs, [Int64(0)], point, work, lazard=True)
    c.check(String(cell1) == String(cell2), "complete projection lifting and cell certificate replay")
    c.check(String(cylindrical_covering(inputs, [UInt8(2)], [Int64(0)], work, lazard=True)) == String(cylindrical_covering(inputs, [UInt8(2)], [Int64(0)], work, lazard=True)), "selective covering certificate replay")
    c.finish()
