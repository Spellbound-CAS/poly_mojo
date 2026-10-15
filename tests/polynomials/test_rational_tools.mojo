"""Rational-function calculus, partial fractions and Laurent identities."""
from apn_mojo import Rational
from polynomials import Budget, PolyErrorKind, QQ, NoContext, DensePoly, MonomialOrder, RationalFunction, rational_evaluate, rational_derivative, rational_substitute, rational_compose, partial_fractions, check_partial_fractions, rational_poles, laurent_expansion, check_laurent_expansion, rational_residue
from support.check import Checker


def poly(values: List[Int], mut work: Budget) raises -> DensePoly[QQ]:
    var coefficients = List[QQ]()
    for value in values:
        coefficients.append(QQ(Rational(value)))
    return DensePoly[QQ].from_coefficients(0, NoContext(), coefficients^, work)


def fraction(a: List[Int], b: List[Int], mut work: Budget) raises -> RationalFunction:
    var order = MonomialOrder()
    return RationalFunction.normalize(poly(a, work).to_sparse(order, work), poly(b, work).to_sparse(order, work), work)[0].copy()


def main() raises:
    var c = Checker("rational_tools")
    var work = Budget()
    # (2x+3)/((x-1)^2*(x+1)) = -1/(4(x-1)) + 5/(2(x-1)^2) + 1/(4(x+1)).
    var f = fraction([3, 2], [1, -1, -1, 1], work)
    var assignments: List[Tuple[Int64, QQ]] = [(Int64(0), QQ(Rational(0)))]
    c.check(rational_evaluate(f, assignments, work) == QQ(Rational(3)), "rational evaluation")
    var derivative = rational_derivative(f, 0, work)
    c.check(rational_evaluate(derivative, assignments, work) == QQ(Rational(5)), "quotient rule at zero")
    var decomposition = partial_fractions(f, 0, work)
    c.check(len(decomposition.terms) == 3 and check_partial_fractions(f, 0, decomposition, work), "partial fraction recomposition identity")
    var corrupt = decomposition.copy()
    corrupt.terms[0].numerator = poly([99], work)
    c.check(not check_partial_fractions(f, 0, corrupt, work), "corrupt partial fraction rejected")
    var poles = rational_poles(f, 0, work)
    c.check(len(poles.factors) == 2, "both exact pole factors retained")
    var expansion = laurent_expansion(f, 0, QQ(Rational(1)), 4, work)
    c.check(expansion.valuation == -2 and expansion.coefficients[0] == QQ(Rational(5, 2)) and expansion.coefficients[1] == QQ(Rational(-1, 4)), "Laurent principal part")
    c.check(check_laurent_expansion(f, expansion, work), "Laurent truncated identity")
    var wrong = expansion.copy()
    wrong.coefficients[1] = QQ(Rational(7))
    c.check(not check_laurent_expansion(f, wrong, work), "corrupt Laurent residue rejected")
    c.check(rational_residue(f, 0, QQ(Rational(1)), work) == QQ(Rational(-1, 4)), "rational residue")
    var polynomial_fraction = fraction([0, 0, 1], [-1, 1], work)
    var nonproper = partial_fractions(polynomial_fraction, 0, work)
    c.check(nonproper.polynomial == poly([1, 1], work) and check_partial_fractions(polynomial_fraction, 0, nonproper, work), "polynomial part retained")
    var inner = fraction([1, 1], [1], work)
    var composed = rational_compose(f, inner, 0, work)
    var sample: List[Tuple[Int64, QQ]] = [(Int64(0), QQ(Rational(1)))]
    c.check(rational_evaluate(composed, sample, work) == QQ(Rational(7, 3)), "rational composition")
    var substitutions = List[Tuple[Int64, type_of(inner.numerator)]]()
    substitutions.append((Int64(0), inner.numerator.copy()))
    c.check(rational_substitute(f, substitutions, work) == composed, "simultaneous polynomial substitution")
    var pole: List[Tuple[Int64, QQ]] = [(Int64(0), QQ(Rational(1)))]
    try:
        _ = rational_evaluate(f, pole, work)
        c.check(False, "evaluation at pole accepted")
    except e:
        c.check(e.kind == PolyErrorKind.DOMAIN, "pole evaluation rejected")
    var tiny = Budget(0)
    try:
        _ = partial_fractions(f, 0, tiny)
        c.check(False, "partial fractions ignored budget")
    except e:
        c.check(e.kind == PolyErrorKind.LIMIT, "partial fractions budget")
    c.finish()
