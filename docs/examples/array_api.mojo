"""Work with exact dense and sparse polynomial arrays."""
import poly
from apn_mojo import Rational, ExactComplex
from polynomials import Budget, PolyErrorKind


def main() raises:
    # Coefficients run from the constant term to the highest power.
    var p = poly.dense.asarray([1, 2, 3])
    print(p.domain(), p.to_list())
    print(poly.dense.evaluate(p, 2))
    print(poly.dense.multiply([1, 1], [1, -1]).to_list())

    # Promote heterogeneous values before making a homogeneous collection.
    var mixed = poly.coefficients((1, Rational(1, 2), ExactComplex(0, 1)))
    print(mixed.domain(), poly.dense.evaluate(mixed, 2))

    # 2*x + 3*y^2: tuple positions identify the variables.
    var terms = [((1, 0), 2), ((0, 2), 3)]
    print(poly.sparse.evaluate(terms, (2, 3)))
    print(poly.sparse.derivative(terms, axis=1).to_list())

    var modular = poly.dense.asarray([-1, 2], modulus=7)
    print(modular.domain(), poly.dense.evaluate(modular, 3).value)

    # A tuple keeps sparse coefficient types intact until promotion.
    var mixed_terms = poly.sparse.asarray(((0, 3), (2, Rational(1, 2))))
    print(mixed_terms.domain(), poly.sparse.evaluate(mixed_terms, 2))

    # A shared budget covers a sequence of operations, including conversion.
    var budget = Budget(10000)
    var slope = poly.dense.derivative(p, budget)
    print(poly.dense.evaluate(slope, 2, budget))
    try:
        _ = poly.dense.multiply(p, p, max_steps=1)
    except error:
        if error.kind != PolyErrorKind.LIMIT:
            raise error
        print("work limit reached")
