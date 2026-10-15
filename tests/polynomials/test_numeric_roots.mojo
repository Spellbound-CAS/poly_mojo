"""Theorem-derived root counts, family guarantees and hostile certificate checks."""
from apn_mojo import Rational, Float, Complex, Ball, ComplexBall
from polynomials import Budget, QQ, ZZ, NoContext, DensePoly, PolyErrorKind
from polynomials.numeric_domains import FloatDomain, ComplexDomain, BallDomain, ComplexBallDomain
from polynomials.numeric_dense import NumericDensePoly
from polynomials.numeric_bounds import NumericBox
from polynomials.numeric_options import NumericOptions, NumericOutcome
from polynomials.numeric_real_roots import real_roots_exact, real_roots_numeric, check_real_roots_exact, check_real_roots_numeric
from polynomials.numeric_complex_roots import complex_roots_exact, complex_roots_numeric, check_complex_roots_exact, check_complex_roots_numeric
from polynomials.numeric_real_roots import _real_predicate
from support.check import Checker


def main() raises:
    var c = Checker("validated numerical roots")
    var b = Budget()
    var options = NumericOptions(absolute_bits=4, max_subdivisions=2048)
    var p = DensePoly[QQ].from_coefficients(0, NoContext(), [QQ(Rational(-2)), QQ(Rational(0)), QQ(Rational(1))], b)
    var real = real_roots_exact(p, options, b, complete_only=True)
    c.check(real.complete and real.accurate and real.isolated and check_real_roots_exact(p, real, options, b), "sqrt two complete real isolation")
    var count = 0
    for cell in real.cells:
        count += cell.count
    c.check(count == 2, "two distinct real square roots")
    var bad_real = real.copy()
    bad_real.cells[0].count += 1
    c.check(not check_real_roots_exact(p, bad_real, options, b), "forged Sturm count rejected")
    bad_real = real.copy()
    _ = bad_real.cells.pop()
    c.check(not check_real_roots_exact(p, bad_real, options, b), "missing real cover cell rejected")
    bad_real = real.copy()
    bad_real.cells[0].interval.real_low = QQ(Rational(-99))
    c.check(not check_real_roots_exact(p, bad_real, options, b), "incorrect global bound rejected")
    var repeated = DensePoly[QQ].from_coefficients(0, NoContext(), [QQ(Rational(0)), QQ(Rational(1)), QQ(Rational(-2)), QQ(Rational(1))], b)
    var repeated_real = real_roots_exact(repeated, options, b, complete_only=True)
    count = 0
    for cell in repeated_real.cells:
        count += cell.count
    c.check(count == 2 and check_real_roots_exact(repeated, repeated_real, options, b), "zero, boundary and repeated real roots counted distinctly")
    var fp = NumericDensePoly[FloatDomain].from_coefficients(0, [Float(-2), Float(0), Float(1)], b)
    var fr = real_roots_numeric(fp, options, b, complete_only=True)
    c.check(check_real_roots_numeric(fp, fr, options, b) and Bool(fr.sequence), "stored binary points take exact Sturm path")
    var family = NumericDensePoly[BallDomain].from_coefficients(0, [Ball(-2, Rational(1, 4096)), Ball(0), Ball(1)], b)
    var family_roots = real_roots_numeric(family, options, b, complete_only=True)
    c.check(check_real_roots_numeric(family, family_roots, options, b) and not family_roots.sequence, "uniform real family root counts")
    var bad_family = family_roots.copy()
    bad_family.cells[0].evaluations[0].steps[0].real_low = QQ(Rational(999))
    c.check(not check_real_roots_numeric(family, bad_family, options, b), "malformed family enclosure rejected")
    var unstable = NumericDensePoly[BallDomain].from_coefficients(0, [Ball(0, Rational(1, 8)), Ball(0), Ball(1)], b)
    var small = NumericOptions(initial_precision=32, max_precision=32, absolute_bits=4, max_subdivisions=12)
    var unresolved = real_roots_numeric(unstable, small, b)
    c.check(not unresolved.complete and unresolved.outcome == NumericOutcome.INPUT_UNCERTAINTY, "x squared minus uncertain a has no fixed real count")
    var degree_family = NumericDensePoly[BallDomain].from_coefficients(0, [Ball(1), Ball(0, 1)], b)
    c.check(real_roots_numeric(degree_family, small, b).outcome == NumericOutcome.UNCERTAIN_DEGREE, "uncertain real degree explicit")
    c.check(complex_roots_numeric(degree_family, small, b).outcome == NumericOutcome.UNCERTAIN_DEGREE, "uncertain complex degree explicit")
    try:
        _ = real_roots_numeric(degree_family, small, b, complete_only=True)
        c.check(False, "complete-only uncertain degree accepted")
    except e:
        c.check(e.kind == PolyErrorKind.UNSUPPORTED, "complete-only uncertain degree typed failure")
    var constant = DensePoly[ZZ].from_coefficients(0, NoContext(), [ZZ(3)], b)
    var no_real = real_roots_exact(constant, options, b)
    var no_complex = complex_roots_exact(constant, options, b)
    c.check(no_real.complete and check_real_roots_exact(constant, no_real, options, b), "constant complete real empty root cover")
    c.check(no_complex.complete and len(no_complex.regions) == 0 and check_complex_roots_exact(constant, no_complex, options, b), "constant complete complex empty cover")
    var zero = DensePoly[QQ](0, NoContext())
    try:
        _ = real_roots_exact(zero, options, b)
        c.check(False, "zero polynomial accepted for finite isolation")
    except e:
        c.check(e.kind == PolyErrorKind.DOMAIN, "zero polynomial rejected")
    var z2p1 = DensePoly[QQ].from_coefficients(0, NoContext(), [QQ(Rational(1)), QQ(Rational(0)), QQ(Rational(1))], b)
    var complex = complex_roots_exact(z2p1, options, b, complete_only=True)
    c.check(complex.complete and complex.accurate and complex.isolated and len(complex.regions) == 2 and check_complex_roots_exact(z2p1, complex, options, b), "conjugate simple roots with checked winding")
    var plus_i = NumericBox.point(QQ(Rational(0)), QQ(Rational(1)), b)
    var minus_i = NumericBox.point(QQ(Rational(0)), QQ(Rational(-1)), b)
    var hits = 0
    for region in complex.regions:
        hits += Int(region.rectangle.contains(plus_i, b)) + Int(region.rectangle.contains(minus_i, b))
    c.check(hits == 2, "both planted conjugate roots enclosed")
    var bad = complex.copy()
    bad.regions[0].count += 1
    c.check(not check_complex_roots_exact(z2p1, bad, options, b), "forged winding count rejected")
    bad = complex.copy()
    _ = bad.regions[0].boundary.pop()
    c.check(not check_complex_roots_exact(z2p1, bad, options, b), "missing boundary segment rejected")
    bad = complex.copy()
    bad.regions[0].boundary[0].end = plus_i.copy()
    c.check(not check_complex_roots_exact(z2p1, bad, options, b), "noncontiguous boundary rejected")
    bad = complex.copy()
    bad.regions[0].boundary[0].evaluation._value[0] = ComplexBall(Ball(0, 100), Ball(0, 100))
    c.check(not check_complex_roots_exact(z2p1, bad, options, b), "zero-containing image rejected even when enclosing true values")
    bad = complex.copy()
    bad.regions.append(complex.regions[0].copy())
    c.check(not check_complex_roots_exact(z2p1, bad, options, b), "overlapping root regions rejected")
    bad = complex.copy()
    _ = bad.regions.pop()
    c.check(not check_complex_roots_exact(z2p1, bad, options, b), "omitted cluster rejected by degree sum")
    var replay = complex_roots_exact(z2p1, options, b)
    c.check(replay.regions[0].rectangle == complex.regions[0].rectangle and replay.regions[1].rectangle == complex.regions[1].rectangle, "root region replay deterministic")
    var nonreal = NumericDensePoly[ComplexDomain].from_coefficients(0, [Complex(-1, -2), Complex(1)], b)
    var nr = complex_roots_numeric(nonreal, options, b, complete_only=True)
    c.check(len(nr.regions) == 1 and nr.regions[0].rectangle.contains(NumericBox.point(QQ(Rational(1)), QQ(Rational(2)), b), b) and check_complex_roots_numeric(nonreal, nr, options, b), "nonreal point coefficient root")
    var square = NumericDensePoly[ComplexDomain].from_coefficients(0, [Complex(-1), Complex(0, -2), Complex(1)], b)
    var clusters = NumericOptions(absolute_bits=4, max_subdivisions=2048, require_isolation=False)
    var cluster = complex_roots_numeric(square, clusters, b, complete_only=True)
    c.check(cluster.complete and cluster.accurate and not cluster.isolated and len(cluster.regions) == 1 and cluster.regions[0].count == 2 and cluster.regions[0].rectangle.contains(plus_i, b) and check_complex_roots_numeric(square, cluster, clusters, b), "double complex root certified as multiplicity-two cluster")
    var isolation = NumericOptions(absolute_bits=3, max_subdivisions=180)
    var unseparated = complex_roots_numeric(square, isolation, b)
    c.check(unseparated.outcome != NumericOutcome.SUCCESS and unseparated.complete and not unseparated.isolated and check_complex_roots_numeric(square, unseparated, isolation, b), "repeated root isolation inconclusive with valid cluster cover")
    var cbp = NumericDensePoly[ComplexBallDomain].from_coefficients(0, [ComplexBall(Ball(-1, Rational(1, 4096)), Ball(-2, Rational(1, 4096))), ComplexBall(1)], b)
    var cbr = complex_roots_numeric(cbp, options, b, complete_only=True)
    c.check(check_complex_roots_numeric(cbp, cbr, options, b) and cbr.regions[0].count == 1, "complex coefficient family uniform root cover")
    var limited = NumericOptions(absolute_bits=32, max_subdivisions=0)
    c.check(real_roots_exact(p, limited, b).outcome == NumericOutcome.SUBDIVISION_LIMIT, "real subdivision ceiling explicit")
    c.check(complex_roots_exact(z2p1, limited, b).outcome == NumericOutcome.UNRESOLVED_BOUNDARY, "complex unresolved initial boundary explicit")
    var tiny = Budget(0)
    try:
        _ = check_complex_roots_exact(z2p1, complex, options, tiny)
        c.check(False, "checker ignored budget")
    except e:
        c.check(e.kind == PolyErrorKind.LIMIT and z2p1.degree() == 2, "checker exhaustion preserves source")
    for i in range(len(family_roots.cells)):
        if family_roots.cells[i].count == 1:
            var invalid_derivative = family_roots.copy()
            invalid_derivative.cells[i].evaluations[1]._value[0] = Ball(0, 100)
            c.check(not check_real_roots_numeric(family, invalid_derivative, options, b), "zero-containing derivative cannot prove uniqueness")
            break
    var unit = NumericBox.real_interval(QQ(Rational(-1)), QQ(Rational(1)), b)
    var wide = NumericBox.real_interval(QQ(Rational(-3)), QQ(Rational(3)), b)
    var positive_derivative = NumericBox.real_interval(QQ(Rational(1)), QQ(Rational(2)), b)
    var origin = NumericBox.point(QQ(Rational(0)), QQ(Rational(0)), b)
    var images: List[NumericBox] = [wide.copy(), positive_derivative.copy(), origin.copy(), wide.copy(), wide.copy()]
    c.check(_real_predicate(unit, images, b)[1] == 3, "strict Newton inclusion establishes one root")
    images[2] = NumericBox.real_interval(QQ(Rational(3)), QQ(Rational(4)), b)
    c.check(_real_predicate(unit, images, b)[1] == 4, "disjoint Newton interval excludes roots")
    images[2] = NumericBox.point(QQ(Rational(-1)), QQ(Rational(0)), b)
    images[1] = NumericBox.point(QQ(Rational(1)), QQ(Rational(0)), b)
    c.check(_real_predicate(unit, images, b)[0] == -1, "Newton boundary contact is inconclusive")
    var relative_only = NumericOptions(absolute_bits=0, relative_bits=8, max_subdivisions=8)
    var at_zero = DensePoly[QQ].from_coefficients(0, NoContext(), [QQ(Rational(0)), QQ(Rational(1))], b)
    c.check(real_roots_exact(at_zero, relative_only, b).outcome == NumericOutcome.SUBDIVISION_LIMIT, "relative-only tolerance near zero is inconclusive")
    var float_complex = complex_roots_numeric(fp, options, b, complete_only=True)
    c.check(check_complex_roots_numeric(fp, float_complex, options, b) and len(float_complex.regions) == 2, "Float coefficients embedded as exact stored complex points")
    var ball_complex = complex_roots_numeric(family, options, b, complete_only=True)
    c.check(check_complex_roots_numeric(family, ball_complex, options, b) and len(ball_complex.regions) == 2, "real ball family complex roots")
    var invalid_options = NumericOptions(initial_precision=128, max_precision=64)
    try:
        _ = real_roots_exact(p, invalid_options, b)
        c.check(False, "invalid root options accepted")
    except e:
        c.check(e.kind == PolyErrorKind.DOMAIN, "invalid root options rejected")
    var no_work = Budget(0)
    try:
        _ = complex_roots_exact(z2p1, options, no_work)
        c.check(False, "root producer ignored work budget")
    except e:
        c.check(e.kind == PolyErrorKind.LIMIT and z2p1.degree() == 2, "root producer budget failure preserves source")
    c.finish()
