"""Complex algebraic multivariate factors need not be linear."""
from poly import sparse, dense, algebraic, exact
from std.testing import assert_equal, assert_true


def main() raises:
    var i = algebraic.ascomplex(0, 1)
    var minus = sparse.asarray((((1, 0), 1), ((0, 1), algebraic.negative(i))))
    var plus = sparse.asarray((((1, 0), 1), ((0, 1), i.copy())))
    var p = sparse.asarray((((2, 0), algebraic.ascomplex(1)), ((0, 2), 1)))
    var r = sparse.factor_checked(p)
    assert_true(exact.check(r))
    var f = exact.factors(r)
    assert_equal(len(f), 2)
    assert_equal(f[0][0].to_list(), minus.to_list())
    assert_equal(f[1][0].to_list(), plus.to_list())
    # Coefficients can combine unrelated real extensions and an imaginary part.
    var root2 = dense.root_between([-2, 0, 1], 1, 2)
    var z = algebraic.ascomplex(root2, 1)
    var a = sparse.asarray((((1, 0), 1), ((0, 1), algebraic.negative(z))))
    var scaled = sparse.multiply(sparse.power(a, 2), sparse.asarray((((0, 0), z.copy()),)))
    var repeated = sparse.factor_checked(scaled)
    assert_true(exact.check(repeated))
    assert_true(exact.unit(repeated) == z)
    assert_equal(exact.factors(repeated)[0][0].to_list(), a.to_list())
    assert_equal(exact.factors(repeated)[0][1], 2)
    # y+x^2 is irreducible over QQbar: it is monic and linear in y.
    var curve = sparse.asarray((((2, 0), algebraic.ascomplex(1)), ((0, 1), 1)))
    var cr = sparse.factor_checked(curve)
    assert_true(exact.check(cr))
    assert_equal(len(exact.factors(cr)), 1)
    assert_equal(exact.factors(cr)[0][0].to_list(), curve.to_list())
    # xy+1 has no scalar leading coefficient in either variable. Its split
    # univariate image does not yield any proper sparse divisor.
    var hyperbola = sparse.asarray((((1, 1), algebraic.ascomplex(1)), ((0, 0), 1)))
    var hr = sparse.factor_checked(hyperbola)
    assert_true(exact.check(hr))
    assert_equal(len(exact.factors(hr)), 1)
    assert_equal(exact.factors(hr)[0][0].to_list(), hyperbola.to_list())
    # A nonhomogeneous product requires recombining several image roots.
    var line = sparse.asarray((((1, 0), algebraic.ascomplex(1)), ((0, 1), 1)))
    var nonhomogeneous = sparse.factor_checked(sparse.multiply(hyperbola, line))
    assert_true(exact.check(nonhomogeneous))
    var nf = exact.factors(nonhomogeneous)
    assert_equal(len(nf), 2)
    assert_equal(nf[0][0].to_list(), line.to_list())
    assert_equal(nf[1][0].to_list(), hyperbola.to_list())
    print("array multivariate complex algebraic factor checks passed")
