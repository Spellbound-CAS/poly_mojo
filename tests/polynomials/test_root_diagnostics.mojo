"""Root residual bounds and Newton candidates require independent validation."""
from apn_mojo import Rational, Ball, Complex, ComplexContext, ArithmeticContext, FloatFormat
from polynomials import Budget, QQ, NoContext, DensePoly, NumericBox, NumericOptions, NumericOutcome, root_diagnostics_exact, check_root_diagnostics_exact, polish_root_exact, check_polished_root_exact, NumericDensePoly, BallDomain, ComplexDomain, root_diagnostics_numeric, check_root_diagnostics_numeric, polish_root_numeric, check_polished_root_numeric
from support.check import Checker


def main() raises:
    var c = Checker("root_diagnostics")
    var work = Budget()
    var f = DensePoly[QQ].from_coefficients(0, NoContext(), [QQ(Rational(-2)), QQ(Rational(0)), QQ(Rational(1))], work)
    var point = NumericBox.point(QQ(Rational(3, 2)), QQ(Rational(0)), work)
    var diagnostic = root_diagnostics_exact(f, point, work)
    c.check(check_root_diagnostics_exact(f, point, diagnostic, work), "checked exact diagnostic")
    c.check(diagnostic.backward_error_low == QQ(Rational(1, 19)) and diagnostic.backward_error_high == QQ(Rational(1, 19)), "absolute coefficient backward error |1/4|/(1+3/2+9/4)")
    c.check(diagnostic.newton_correction_upper.value() == QQ(Rational(1, 12)), "Newton correction at 3/2")
    c.check(diagnostic.coefficient_sensitivity_upper.value() == QQ(Rational(19, 12)), "local absolute coefficient sensitivity")
    var bad = diagnostic.copy()
    bad.backward_error_high = QQ(Rational(0))
    c.check(not check_root_diagnostics_exact(f, point, bad, work), "understated residual diagnostic rejected")
    var singular = root_diagnostics_exact(f, NumericBox.point(QQ(Rational(0)), QQ(Rational(0)), work), work)
    c.check(not singular.newton_correction_upper and not singular.coefficient_sensitivity_upper, "zero derivative exposes unavailable bounds")
    var options = NumericOptions(absolute_bits=6, max_subdivisions=128)
    var context = ComplexContext(ArithmeticContext(format=FloatFormat(128)))
    var polished = polish_root_exact(f, Complex(1), QQ(Rational(1, 256)), 7, context, options, work)
    c.check(polished.outcome == NumericOutcome.SUCCESS and check_polished_root_exact(f, polished, options, work), "Newton candidate certified by local contour")
    var wrong = polished.copy()
    wrong.candidate = Complex(99)
    c.check(not check_polished_root_exact(f, wrong, options, work), "candidate outside root witness rejected")
    var stationary = polish_root_exact(f, Complex(0), QQ(Rational(1, 256)), 7, context, options, work)
    c.check(stationary.outcome != NumericOutcome.SUCCESS and not check_polished_root_exact(f, stationary, options, work), "stationary Newton seed cannot invent a root")
    var family = NumericDensePoly[BallDomain].from_coefficients(0, [Ball(-2, Rational(1, 4096)), Ball(0), Ball(1)], work)
    var family_diagnostic = root_diagnostics_numeric(family, point, work)
    c.check(check_root_diagnostics_numeric(family, point, family_diagnostic, work), "uniform family residual bounds")
    var family_polished = polish_root_numeric(family, Complex(1), QQ(Rational(1, 256)), 7, context, options, work)
    c.check(check_polished_root_numeric(family, family_polished, options, work), "Newton representative validated against original ball family")
    var complex_poly = NumericDensePoly[ComplexDomain].from_coefficients(0, [Complex(1), Complex(0), Complex(1)], work)
    var imaginary = polish_root_numeric(complex_poly, Complex(1, 1), QQ(Rational(1, 256)), 8, context, options, work)
    c.check(check_polished_root_numeric(complex_poly, imaginary, options, work), "imaginary root candidate retains complex components")
    c.finish()
