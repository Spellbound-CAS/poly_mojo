"""Storage lifetime, shared operands and unchanged certified instruction traces."""
from apn_mojo import Rational, BallContext
from polynomials import *
from polynomials._numeric_program import _NumericProgram, _dense_program, _sparse_program
from polynomials.numeric_evaluate import _evaluate_program, _evaluate_scheduled, _check_program
from polynomials.numeric_domains import _numeric_apply
from polynomials._numeric_storage import _storage_plan
from support.check import Checker


def exercise[D: NumericDomain](program: _NumericProgram, context: D.Context, mut check: Checker, mut work: Budget) raises:
    var proof = _evaluate_program[D](program, context, work)
    var compact = _evaluate_scheduled[D, True](program, context, work)
    var values = List[D.Scalar]()
    var valid = True
    # Direct scalar semantics with one owning slot per instruction. This is an
    # independent lifetime oracle for the compact-storage interpreter.
    for i in range(len(program.instructions)):
        ref instruction = program.instructions[i]
        var value = numeric_value[D](program.inputs[instruction.left], context, work) if instruction.operation == 0 else _numeric_apply[D](values[instruction.left], values[instruction.right], 0 if instruction.operation == 1 else 2, context, work)
        valid = valid and proof.steps[i] == numeric_bounds[D](value, work) and compact.steps[i] == proof.steps[i]
        values.append(value^)
    check.check(valid and numeric_bounds[D](proof.value(work), work) == numeric_bounds[D](values[program.result], work) and numeric_bounds[D](compact.value(work), work) == numeric_bounds[D](proof.value(work), work), "all trace bounds and final value retain scalar schedule")
    check.check(_check_program[D](program, proof, work), "compact evaluation passes independent rational checker")
    var corrupt = proof.copy()
    corrupt.steps[program.result] = NumericBox.point(QQ(Rational(99999)), QQ(Rational(0)), work)
    check.check(not _check_program[D](program, corrupt, work), "altered live result trace rejected")
    var tiny = Budget(20)
    try:
        _ = _evaluate_program[D](program, context, tiny)
        check.check(False, "storage analysis ignored budget")
    except e:
        check.check(e.kind == PolyErrorKind.LIMIT, "storage analysis and evaluation are budgeted")


def main() raises:
    var check = Checker("numerical storage")
    var work = Budget()
    var context = BallContext(12)
    for size in [0, 1, 2, 20, 21, 65, 257]:
        var coefficients = List[NumericBox]()
        for i in range(size):
            coefficients.append(NumericBox.point(QQ(Rational(i + 1, 3)), QQ(Rational(0)), work))
        var program = _dense_program(coefficients, NumericBox.point(QQ(Rational(1, 7)), QQ(Rational(0)), work), work)
        var storage = _storage_plan(program, work)
        check.check(storage.count <= size + 3, "dense temporary slots grow with live coefficients, not all instructions")
        exercise[BallDomain](program, context, check, work)
        exercise[ComplexBallDomain](program, context, check, work)
    var dag = _NumericProgram()
    var one = dag.input(NumericBox.point(QQ(Rational(1)), QQ(Rational(0)), work), work)
    _ = dag.input(NumericBox.point(QQ(Rational(42)), QQ(Rational(0)), work), work)
    var square = dag.operation(2, one, one, work)
    var double = dag.operation(1, square, square, work)
    _ = dag.operation(2, double, one, work)
    dag.result = square  # Result survives later consumers and dead instructions.
    exercise[BallDomain](dag, context, check, work)
    exercise[ComplexBallDomain](dag, context, check, work)
    var complex_dag = _NumericProgram()
    var z = complex_dag.input(NumericBox.point(QQ(Rational(2)), QQ(Rational(3)), work), work)
    _ = complex_dag.input(NumericBox.point(QQ(Rational(42)), QQ(Rational(1)), work), work)
    var z_squared = complex_dag.operation(2, z, z, work)
    _ = complex_dag.operation(1, z_squared, z, work)
    complex_dag.result = z_squared
    exercise[ComplexBallDomain](complex_dag, context, check, work)
    var coefficients: List[NumericBox] = [NumericBox.point(QQ(Rational(2)), QQ(Rational(0)), work), NumericBox.point(QQ(Rational(3)), QQ(Rational(0)), work), NumericBox.point(QQ(Rational(5)), QQ(Rational(0)), work)]
    var monomials: List[Monomial] = [Monomial([(Int64(-3), UInt32(2)), (Int64(8), UInt32(3))], work), Monomial([(Int64(-3), UInt32(2)), (Int64(99), UInt32(1))], work), Monomial([(Int64(8), UInt32(3))], work)]
    var assignments = List[Tuple[Int64, NumericBox]]()
    for variable in [99, 8, -3, 123]:
        assignments.append((Int64(variable), NumericBox.point(QQ(Rational(1)), QQ(Rational(0)), work)))
    var sparse = _sparse_program(coefficients, monomials, assignments, work)
    var expected: List[Tuple[Int, Int, Int]] = [(2, 5, 5), (2, 7, 9), (2, 0, 10), (2, 7, 4), (2, 4, 4), (2, 12, 13), (2, 11, 14), (1, 8, 15), (2, 1, 10), (2, 7, 3), (2, 17, 18), (1, 16, 19), (2, 2, 14), (1, 20, 21)]
    var unchanged = len(sparse.instructions) == 23 and sparse.result == 22
    for i in range(len(expected)):
        ref instruction = sparse.instructions[9 + i]
        unchanged = unchanged and Int(instruction.operation) == expected[i][0] and instruction.left == expected[i][1] and instruction.right == expected[i][2]
    check.check(unchanged, "indexed lookup preserves the specified sparse certificate instruction stream")
    exercise[BallDomain](sparse, context, check, work)
    exercise[ComplexBallDomain](sparse, context, check, work)
    assignments.append(assignments[0].copy())
    try:
        _ = _sparse_program(coefficients, monomials, assignments, work)
        check.check(False, "duplicate enclosure assignment accepted")
    except e:
        check.check(e.kind == PolyErrorKind.DOMAIN, "duplicate enclosure assignment rejected by shared index")
    coefficients = List[NumericBox]()
    monomials = List[Monomial]()
    assignments = List[Tuple[Int64, NumericBox]]()
    for k in range(20):
        coefficients.append(NumericBox.point(QQ(Rational(1)), QQ(Rational(0)), work))
        monomials.append(Monomial([(Int64(3 * k - 15), UInt32(2))], work))
        var j = 19 - k
        assignments.append((Int64(3 * j - 15), NumericBox.point(QQ(Rational(j + 1)), QQ(Rational(0)), work)))
    sparse = _sparse_program(coefficients, monomials, assignments, work)
    exercise[BallDomain](sparse, context, check, work)
    exercise[ComplexBallDomain](sparse, context, check, work)
    var sum_squares = _evaluate_program[BallDomain](sparse, context, work)
    check.check(sum_squares.value(work).lower_rational() == Rational(2870) and sum_squares.value(work).upper_rational() == Rational(2870), "lookup promotion preserves sum of twenty integer squares")
    assignments.append(assignments[0].copy())
    try:
        _ = _sparse_program(coefficients, monomials, assignments, work)
        check.check(False, "promoted duplicate assignment accepted")
    except e:
        check.check(e.kind == PolyErrorKind.DOMAIN, "promoted index rejects duplicate assignments")
    _ = assignments.pop()
    _ = assignments.pop()
    try:
        _ = _sparse_program(coefficients, monomials, assignments, work)
        check.check(False, "promoted missing assignment accepted")
    except e:
        check.check(e.kind == PolyErrorKind.DOMAIN, "promoted index rejects missing assignments")
    check.finish()
