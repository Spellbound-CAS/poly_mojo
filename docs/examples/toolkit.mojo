"""A compact tour of the expanded polynomial toolkit."""
from apn_mojo import Rational, BallContext
from polynomials import *


def main() raises:
    var work = Budget()
    var ring = PolynomialRing[QQ].create(["x"], NoContext(), work)
    var p = DensePoly[QQ].from_sparse(ring.parse("(x+1)^3", work), work, variable=Int64(0))
    if derivative_dense(integrate_dense(p, work), work) != p:
        raise Error("calculus identity failed")
    var nodes: List[QQ] = [QQ(Rational(0)), QQ(Rational(1)), QQ(Rational(2)), QQ(Rational(3))]
    var plan = MultipointPlan[QQ].prepare(0, NoContext(), nodes, work)
    var values = plan.evaluate(p, work)
    if not check_interpolation(plan.interpolate(values, work), nodes, values, work):
        raise Error("interpolation witness failed")
    if deserialize_dense[QQ](serialize_dense(p, work), NoContext(), work) != p:
        raise Error("interchange failed")
    print("Calculus, interpolation and interchange checked")

    var basis = BasisPolynomial[QQ].bernstein(p, 5, work)
    var series = TruncatedSeries[QQ].from_polynomial(p, 6, work)
    if basis.to_polynomial(work) != p or series.log(work).exp(work).polynomial != p:
        raise Error("basis or formal-series identity failed")
    print("Basis conversion and formal series checked")

    var fraction = RationalFunction.normalize(ring.parse("x^2", work), ring.parse("x-1", work), work)
    var parts = partial_fractions(fraction[0], 0, work)
    var local = laurent_expansion(fraction[0], 0, QQ(Rational(1)), 4, work)
    if not check_partial_fractions(fraction[0], 0, parts, work) or not check_laurent_expansion(fraction[0], local, work):
        raise Error("rational-function witness failed")
    print("Partial fractions and Laurent expansion checked")

    var context = BallContext(128)
    var numerical = numeric_from_dense[BallDomain, QQ](p, context, work)
    var divisor = numeric_from_dense[BallDomain, QQ](DensePoly[QQ].from_sparse(ring.parse("x+1", work), work, variable=Int64(0)), context, work)
    var division = numeric_divmod(numerical, divisor, context, work)
    if not check_numeric_division(numerical, divisor, division, work):
        raise Error("numerical division witness failed")
    var coarse = NumericOptions(absolute_bits=2, max_subdivisions=32)
    var fine = NumericOptions(absolute_bits=5, max_subdivisions=64)
    var roots = real_roots_exact(p, coarse, work, complete_only=True)
    var refined = refine_real_roots_exact(p, roots, coarse, fine, work, complete_only=True)
    var multiplicity = real_root_multiplicities(p, refined, fine, work)
    if not check_real_root_multiplicities(p, refined, fine, multiplicity, work):
        raise Error("root multiplicity witness failed")
    print("Numerical division, root refinement and multiplicity checked")
