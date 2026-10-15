"""Divide, integrate, evaluate, and change representation using coefficient arrays."""
import poly


def main() raises:
    # (x^2 + 1) / (2*x + 1) has rational quotient and remainder.
    var result = poly.dense.divmod([1, 0, 1], [1, 2])
    print("Quotient:", result[0].to_list())
    print("Remainder:", result[1].to_list())
    print("Integral:", poly.dense.integral([1, 1]).to_list())

    var p = poly.dense.from_roots([0, 1, 1])
    print("Coefficients:", p.to_list())
    print("Values:", poly.dense.evaluate_many(p, [0, 1, 2, 0]).to_list())
    var terms = poly.to_sparse(p)
    print("Terms:", terms.to_list())
    print("Dense again:", poly.to_dense(terms).to_list())

    # 2*x + 3*y^2, evaluated at three coordinate pairs.
    var surface = [((1, 0), 2), ((0, 2), 3)]
    print("Surface:", poly.sparse.evaluate_many(surface, [(0, 0), (2, 3), (0, 0)]).to_list())
    print("Total degree:", poly.sparse.degree(surface))
