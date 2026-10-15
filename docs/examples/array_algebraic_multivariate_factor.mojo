from poly import sparse, algebraic, exact


def main() raises:
    # Powers (a, b) represent x**a * y**b. Choose the coefficient field
    # with one algebraic value; integer entries promote automatically.
    var real = sparse.asarray((((2, 0), algebraic.asreal(1)), ((0, 2), 1)))
    var complex = sparse.asarray((((2, 0), algebraic.ascomplex(1)), ((0, 2), 1)))
    var real_result = sparse.factor_checked(real)
    var complex_result = sparse.factor_checked(complex)
    print(real.domain(), "factors:", len(exact.factors(real_result)))
    print(complex.domain(), "factors:", len(exact.factors(complex_result)))
    print("verified:", exact.check(real_result), exact.check(complex_result))
    var line = sparse.asarray((((1, 0), algebraic.asreal(1)), ((0, 1), 1)))
    var repeated = sparse.factor_checked(sparse.power(line, 3))
    print("multiplicity:", exact.factors(repeated)[0][1])
