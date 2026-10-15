"""Compute and check a gcd and a rational-function cancellation."""
from polynomials import Budget, QQ, NoContext, MonomialOrder, SparsePoly, gcd, check_common_divisor, check_bezout, RationalFunction, check_cancellation

def main() raises:
    var budget = Budget()
    var order = MonomialOrder(MonomialOrder.LEX, [0], budget)
    var x = SparsePoly[QQ].variable(0, order, NoContext(), budget)
    var one = SparsePoly[QQ].constant(QQ.one(NoContext()), order, budget)
    var numerator = x.pow(2, budget).sub(one, budget)
    var denominator = x.sub(one, budget)
    var result = gcd(numerator, denominator, budget, certify_maximal=True)
    if result.gcd != denominator or not check_common_divisor(numerator, denominator, result, budget) or not check_bezout(numerator, denominator, result, budget):
        raise Error("gcd certificate failed")
    var reduced = RationalFunction.normalize(numerator, denominator, budget)
    if reduced[0].numerator != x.add(one, budget) or reduced[0].denominator != one or not check_cancellation(numerator, denominator, reduced[1], budget):
        raise Error("cancellation certificate failed")
    print("gcd and cancellation certificates checked")
