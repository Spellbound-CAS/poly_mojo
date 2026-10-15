"""Export actual certificates for the independent Python Fraction verifier."""
from apn_mojo import Rational, Ball, ComplexBall, BallContext
from polynomials import Budget, QQ, DensePoly, NoContext, BallDomain, ComplexBallDomain, NumericDomain, NumericBox, enclose_dense, check_enclose_dense
from polynomials.numeric_evaluate import _exact_dense_boxes
from polynomials.numeric_bounds import numeric_bounds
from polynomials._numeric_program import _dense_program


def row(box: NumericBox):
    print(box.real_low.value, box.real_high.value, box.imag_low.value, box.imag_high.value, sep=",")


def emit[D: NumericDomain](p: DensePoly[QQ], x: D.Scalar, mut budget: Budget) raises:
    var proof = enclose_dense[D](p, x, rebind_var[D.Context](BallContext(32)), budget)
    if not check_enclose_dense[D](p, x, proof, budget):
        raise Error("producer certificate rejected")
    var program = _dense_program(_exact_dense_boxes(p, budget), numeric_bounds[D](x, budget), budget)
    print("trace", len(program.inputs), len(program.instructions), program.result, sep=",")
    for input in program.inputs:
        row(input)
    for i in range(len(program.instructions)):
        ref instruction = program.instructions[i]
        print(instruction.operation, instruction.left, instruction.right, sep=",")
        row(proof.steps[i])
    row(numeric_bounds[D](proof.value(budget), budget))


def main() raises:
    var b = Budget()
    var p = DensePoly[QQ].from_coefficients(0, NoContext(), [QQ(Rational(1, 3)), QQ(Rational(-2)), QQ(Rational(1))], b)
    emit[BallDomain](p, Ball(Rational(3, 2), Rational(1, 8)), b)
    emit[ComplexBallDomain](p, ComplexBall(Ball(0, Rational(1, 16)), Ball(1, Rational(1, 16))), b)
    var q = DensePoly[QQ].from_coefficients(0, NoContext(), [QQ(Rational(1)), QQ(Rational(-4)), QQ(Rational(6)), QQ(Rational(-4)), QQ(Rational(1))], b)
    emit[BallDomain](q, Ball(1, Rational(1, 1024)), b)
