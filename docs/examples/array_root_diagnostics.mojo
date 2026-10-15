"""Array root diagnostics."""
from apn_mojo import Rational
from polynomials import NumericOptions
from poly import dense, roots


def main() raises:
    var cover = dense.real_root_cover([1, -2, 1], options=NumericOptions(absolute_bits=2), complete_only=True)
    var repeated = roots.multiplicities(cover)
    print("Multiplicity:", repeated.counts()[0])
    print("Multiplicity checked:", roots.check(repeated))
    var diagnostic = roots.diagnostics([-2, 0, 1], Rational(3, 2))
    print("Residual:", diagnostic.evidence().residual.real_low.value)
    print("Diagnostics checked:", roots.check(diagnostic))
