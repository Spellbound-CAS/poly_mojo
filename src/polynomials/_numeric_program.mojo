"""Deterministic evaluation graphs shared by numerical producers and exact checkers."""
from std.collections import Dict
from .errors import Budget, PolyError, PolyErrorKind, _allocation
from .monomial import Monomial
from ._assignment import _AssignmentPlan
from .numeric_bounds import NumericBox, _point, _box_charge


@fieldwise_init
struct _NumericInstruction(ImplicitlyCopyable):
    var operation: UInt8
    var left: Int
    var right: Int


struct _NumericProgram(Copyable, Movable):
    var inputs: List[NumericBox]
    var instructions: List[_NumericInstruction]
    var result: Int
    var aliases: List[Int]
    var _share_squares: Bool

    def __init__(out self):
        self.inputs = List[NumericBox]()
        self.instructions = List[_NumericInstruction]()
        self.result = -1
        self.aliases = List[Int]()
        self._share_squares = False

    def input(mut self, value: NumericBox, mut budget: Budget) raises PolyError -> Int:
        _box_charge(value, budget)
        _allocation[NumericBox](1, budget)
        _allocation[_NumericInstruction](1, budget)
        var index = len(self.inputs)
        self.inputs.append(value.copy())
        self.instructions.append(_NumericInstruction(0, index, 0))
        if self._share_squares:
            _allocation[Int](1, budget)
            self.aliases.append(-1)
        return len(self.instructions) - 1

    def operation(mut self, operation: UInt8, left: Int, right: Int, mut budget: Budget) raises PolyError -> Int:
        _allocation[_NumericInstruction](1, budget)
        self.instructions.append(_NumericInstruction(operation, left, right))
        if self._share_squares:
            _allocation[Int](1, budget)
            self.aliases.append(-1)
        return len(self.instructions) - 1


def _dense_program(coefficients: List[NumericBox], argument: NumericBox, mut budget: Budget) raises PolyError -> _NumericProgram:
    var program = _NumericProgram()
    for c in coefficients:
        _ = program.input(c, budget)
    var x = program.input(argument, budget)
    var result = program.input(_point(0), budget)
    for i in range(len(coefficients) - 1, -1, -1):
        var product = program.operation(2, result, x, budget)
        result = program.operation(1, product, i, budget)
    program.result = result
    return program^


def _sparse_program(coefficients: List[NumericBox], monomials: List[Monomial], assignments: List[Tuple[Int64, NumericBox]], mut budget: Budget) raises PolyError -> _NumericProgram:
    var program = _NumericProgram()
    # Reuse repeated squaring only for substantial power-sharing workloads.
    # The original instructions, result index and witness order are untouched.
    program._share_squares = len(monomials) >= 256 and len(assignments) <= len(monomials) // 8
    var squares = Dict[Int, Int]()
    for c in coefficients:
        _ = program.input(c, budget)
    var plan = _AssignmentPlan()
    for i in range(len(assignments)):
        plan.assign(assignments[i][0], i, budget)
        _ = program.input(assignments[i][1], budget)
    var one = program.input(_point(1), budget)
    var result = program.input(_point(0), budget)
    var powers = List[Int]()
    for i in range(len(monomials)):
        var term = i
        for pair in monomials[i].pairs:
            var slot = plan.power(pair, budget)
            if slot == len(powers):
                var assigned = plan.assignment(pair[0], budget)
                if assigned < 0:
                    raise PolyError(PolyErrorKind.DOMAIN, "missing enclosure assignment")
                var base = len(coefficients) + assigned
                var exponent = pair[1]
                var power = one
                while exponent:
                    budget.charge(1)
                    if exponent & 1:
                        power = program.operation(2, power, base, budget)
                    exponent >>= 1
                    if exponent:
                        var key = base
                        if program._share_squares and program.aliases[base] >= 0:
                            key = program.aliases[base]
                        var step = program.operation(2, base, base, budget)
                        if program._share_squares:
                            budget.charge(1)
                            var previous = squares.get(key, -1)
                            if previous >= 0:
                                program.aliases[step] = previous
                            else:
                                _allocation[Tuple[Int, Int]](1, budget)
                                squares[key] = step
                        base = step
                _allocation[Int](1, budget)
                powers.append(power)
            term = program.operation(2, term, powers[slot], budget)
        result = program.operation(1, result, term, budget)
    program.result = result
    return program^
