"""Value adapters letting identity tests exercise the certified public producers."""
from polynomials import Budget, PolyError, QQ, DensePoly, SparsePoly, RealAlgebraic, real_rational, algebraic_add, algebraic_multiply, algebraic_negate, algebraic_inverse, algebraic_evaluate, algebraic_roots, algebraic_roots_at
from polynomials.algebraic_types import _copy_real


def _add(a: RealAlgebraic, b: RealAlgebraic, mut budget: Budget) raises PolyError -> RealAlgebraic:
    var proof = algebraic_add(a, b, budget)
    return _copy_real(proof.value, budget)


def _mul(a: RealAlgebraic, b: RealAlgebraic, mut budget: Budget) raises PolyError -> RealAlgebraic:
    var proof = algebraic_multiply(a, b, budget)
    return _copy_real(proof.value, budget)


def _neg(a: RealAlgebraic, mut budget: Budget) raises PolyError -> RealAlgebraic:
    var proof = algebraic_negate(a, budget)
    return _copy_real(proof.value, budget)


def _inverse(a: RealAlgebraic, mut budget: Budget) raises PolyError -> RealAlgebraic:
    var proof = algebraic_inverse(a, budget)
    return _copy_real(proof.value, budget)


def _evaluate(f: DensePoly[QQ], a: RealAlgebraic, mut budget: Budget) raises PolyError -> RealAlgebraic:
    var coefficients = List[RealAlgebraic]()
    for value in f.coefficients:
        coefficients.append(real_rational(value, budget))
    var proof = algebraic_evaluate(coefficients, a, budget)
    return _copy_real(proof.value, budget)


def _roots_algebraic(coefficients: List[RealAlgebraic], mut budget: Budget) raises PolyError -> List[RealAlgebraic]:
    var proof = algebraic_roots(coefficients, budget)
    var roots = List[RealAlgebraic]()
    swap(roots, proof.roots)
    return roots^


def _roots_at(f: SparsePoly[QQ], variable: Int64, point: List[Tuple[Int64, RealAlgebraic]], mut budget: Budget) raises PolyError -> List[RealAlgebraic]:
    var proof = algebraic_roots_at(f, variable, point, budget)
    var roots = List[RealAlgebraic]()
    swap(roots, proof.isolation.roots)
    return roots^
