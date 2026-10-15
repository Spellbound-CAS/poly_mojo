"""Independent exact-value and corruption tests for ball evaluation traces."""
from apn_mojo import Integer, Rational, Ball, ComplexBall, BallContext
from polynomials import Budget, QQ, NoContext, DensePoly, MonomialOrder, PolyErrorKind
from polynomials.numeric_domains import BallDomain, ComplexBallDomain
from polynomials.numeric_options import NumericOptions, NumericOutcome
from polynomials.numeric_bounds import NumericBox
from polynomials.numeric_convert import numeric_from_dense, numeric_from_sparse
from polynomials.numeric_evaluate import enclose_dense, check_enclose_dense, enclose_numeric, check_enclose_numeric, enclose_sparse, check_enclose_sparse, enclose_numeric_sparse, check_enclose_numeric_sparse, adaptive_enclose_dense, adaptive_enclose_sparse, check_numeric_accuracy
from support.check import Checker
from support.corpus import parse_polynomial


def main() raises:
    var c = Checker("numerical enclosures")
    var work = Budget()
    var ctx = BallContext(64)
    var p = DensePoly[QQ].from_coefficients(0, NoContext(), [QQ(Rational(1, 3)), QQ(Rational(-2)), QQ(Rational(1))], work)
    var x = Ball(Rational(3, 2), Rational(1, 8))
    var proof = enclose_dense[BallDomain](p, x, ctx, work)
    c.check(check_enclose_dense[BallDomain](p, x, proof, work), "real ball trace checked independently")
    c.check(proof.value(work).lower_rational() <= Rational(-5, 12) and proof.value(work).upper_rational() >= Rational(-5, 12), "hand computed midpoint value enclosed")
    var corrupt = proof.copy()
    corrupt.steps[len(corrupt.steps) - 1] = NumericBox.real_interval(QQ(Rational(99)), QQ(Rational(100)), work)
    c.check(not check_enclose_dense[BallDomain](p, x, corrupt, work), "wrong intermediate bound rejected")
    corrupt = proof.copy()
    corrupt._value[0] = Ball(99, precision=64)
    c.check(not check_enclose_dense[BallDomain](p, x, corrupt, work), "wrong reported value rejected")
    corrupt = proof.copy()
    _ = corrupt.steps.pop()
    c.check(not check_enclose_dense[BallDomain](p, x, corrupt, work), "missing trace step rejected")
    var numerical = numeric_from_dense[BallDomain, QQ](p, ctx, work)
    var numeric_proof = enclose_numeric[BallDomain](numerical, x, ctx, work)
    c.check(check_enclose_numeric[BallDomain](numerical, x, numeric_proof, work), "coefficient-family trace checked")
    var z = ComplexBall(Ball(0, Rational(1, 16)), Ball(1, Rational(1, 16)))
    var complex_proof = enclose_dense[ComplexBallDomain](p, z, ctx, work)
    c.check(check_enclose_dense[ComplexBallDomain](p, z, complex_proof, work), "complex rectangle trace checked")
    c.check(complex_proof.value(work).real().lower_rational() <= Rational(-2, 3) and complex_proof.value(work).real().upper_rational() >= Rational(-2, 3), "complex exact real component enclosed")
    c.check(complex_proof.value(work).imag().lower_rational() <= Rational(-2) and complex_proof.value(work).imag().upper_rational() >= Rational(-2), "complex exact imaginary component enclosed")
    var sparse = parse_polynomial("x0^3 + 2*x0*x1 + 1/3", MonomialOrder(), work)
    var assignments: List[Tuple[Int64, Ball]] = [(Int64(1), Ball(2)), (Int64(0), Ball(3))]
    var sparse_proof = enclose_sparse[BallDomain](sparse, assignments, ctx, work)
    c.check(check_enclose_sparse[BallDomain](sparse, assignments, sparse_proof, work), "sparse shared power trace checked")
    c.check(sparse_proof.value(work).lower_rational() <= Rational(118, 3) and sparse_proof.value(work).upper_rational() >= Rational(118, 3), "hand computed sparse value")
    var sparse_numeric = numeric_from_sparse[BallDomain, QQ](sparse, ctx, work)
    var sparse_numeric_proof = enclose_numeric_sparse[BallDomain](sparse_numeric, assignments, ctx, work)
    c.check(check_enclose_numeric_sparse[BallDomain](sparse_numeric, assignments, sparse_numeric_proof, work), "sparse coefficient-family trace checked")
    var options = NumericOptions(initial_precision=8, max_precision=128, absolute_bits=40)
    var point = Ball(1)
    var adaptive = adaptive_enclose_dense[BallDomain](p, point, options, work)
    c.check(adaptive.outcome == NumericOutcome.SUCCESS and adaptive.evaluation.precision > 8, "adaptive precision reconstructs rational coefficients")
    c.check(check_enclose_dense[BallDomain](p, point, adaptive.evaluation, work) and check_numeric_accuracy(adaptive, options, work), "adaptive value and accuracy independently checked")
    var adaptive_sparse = adaptive_enclose_sparse[BallDomain](sparse, assignments, options, work)
    c.check(check_enclose_sparse[BallDomain](sparse, assignments, adaptive_sparse.evaluation, work) and check_numeric_accuracy(adaptive_sparse, options, work), "adaptive sparse enclosure")
    var uncertain = adaptive_enclose_dense[BallDomain](p, x, options, work)
    c.check(uncertain.outcome == NumericOutcome.PRECISION_CEILING and check_numeric_accuracy(uncertain, options, work), "input uncertainty survives higher precision")
    uncertain.outcome = NumericOutcome.SUCCESS
    c.check(not check_numeric_accuracy(uncertain, options, work), "unearned accuracy claim rejected")
    var tiny = Budget(0)
    try:
        _ = check_enclose_dense[BallDomain](p, x, proof, tiny)
        c.check(False, "checker ignored work exhaustion")
    except e:
        c.check(e.kind == PolyErrorKind.LIMIT, "checker charges work")
    var huge = Integer(1) << 300
    var cancellation = DensePoly[QQ].from_coefficients(0, NoContext(), [QQ(Rational(huge)), QQ(Rational(1, 3))], work)
    var negative = Ball(Rational(-3 * huge))
    var precision = NumericOptions(initial_precision=16, max_precision=1024, absolute_bits=40)
    var cancelled = adaptive_enclose_dense[BallDomain](cancellation, negative, precision, work)
    c.check(cancelled.outcome == NumericOutcome.SUCCESS and cancelled.evaluation.precision >= 512 and check_numeric_accuracy(cancelled, precision, work), "large cancellation refines from original exact coefficients")
    c.check(check_enclose_dense[BallDomain](cancellation, negative, cancelled.evaluation, work), "large cancellation certificate")
    var tiny_exact = DensePoly[QQ].from_coefficients(0, NoContext(), [QQ(Rational(Integer(1), Integer(1) << 500))], work)
    var tiny_proof = enclose_dense[BallDomain](tiny_exact, point, ctx, work)
    c.check(check_enclose_dense[BallDomain](tiny_exact, point, tiny_proof, work) and tiny_proof.value(work).lower_rational() == Rational(Integer(1), Integer(1) << 500), "very small exact value retained")
    var altered = p.copy()
    altered.coefficients[0] = QQ(Rational(50))
    c.check(not check_enclose_dense[BallDomain](altered, x, proof, work), "changed source coefficient rejected")
    var wrong_x = Ball(30)
    c.check(not check_enclose_dense[BallDomain](p, wrong_x, proof, work), "changed argument rejected")
    var wrong_assignments = assignments.copy()
    wrong_assignments[0] = (Int64(1), Ball(30))
    c.check(not check_enclose_sparse[BallDomain](sparse, wrong_assignments, sparse_proof, work), "changed sparse assignment rejected")
    corrupt = proof.copy()
    corrupt.precision = 128
    c.check(not check_enclose_dense[BallDomain](p, x, corrupt, work), "forged precision rejected")
    c.finish()
