"""Certified Lazard evaluation using exact Taylor coefficients.

Source: CArL 8afa65d, ran/interval/helper/LazardEvaluation.h.
Differentiating to the first nonzero specialization implements division by
vanishing linear factors, up to a positive factorial scalar.
"""
# Copyright (c) 2015-2022 Florian Corzilius, Gereon Kremer, Jasper Nalbach,
# Sebastian Junges, Stefan Schupp, Erika Abraham.
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
# THE SOFTWARE.
from ._algebraic import _NormalizationCache
from .errors import Budget, PolyError, PolyErrorKind, _allocation, _sum
from .coefficients import QQ
from .sparse import SparsePoly
from .evaluate import derivative
from .algebraic_types import RealAlgebraic, _copy_real
from .certificates import _dense_equal
from ._root_ops import _compare as _compare_q
from .algebraic_specialization import AlgebraicSpecializationCertificate, _specialize_cached, check_algebraic_specialization, _zero_specialization, _coefficients, _point


@fieldwise_init
struct LazardEvaluationCertificate(Copyable, Movable, Writable):
    """Successive vanishing orders and all zero/nonzero specialization witnesses."""
    var orders: List[Int]
    """The least nonvanishing derivative order for each assigned variable."""
    var evaluations: List[AlgebraicSpecializationCertificate]
    """All tested derivatives in prefix order; the last is the lifting polynomial."""


def lazard_evaluate(f: SparsePoly[QQ], variable: Int64, point: List[Tuple[Int64, RealAlgebraic]], mut budget: Budget) raises PolyError -> LazardEvaluationCertificate:
    """Remove vanishing factors before each algebraic substitution in cylindrical order.

    Args:
        f: The nonzero rational polynomial.
        variable: The final unassigned lifting variable.
        point: Earlier variable assignments in base-to-top cylindrical order.
        budget: The shared work allowance.

    Returns:
        Certified vanishing orders and a nonzero lifting specialization, up to positive factorials.

    Raises:
        PolyError: DOMAIN for zero input or invalid assignments, LIMIT or CAPACITY on work limits.
    """
    var cache = _NormalizationCache()
    return _lazard_cached(f, variable, point, cache, budget)


@fieldwise_init
struct _LazardSteps(Movable):
    """One polynomial's last Lazard evaluation, step by step along its point."""
    var f: SparsePoly[QQ]
    var variable: Int64
    var coordinates: List[Tuple[Int64, RealAlgebraic]]
    """The assigned coordinates, one per completed step."""
    var evaluations: List[List[AlgebraicSpecializationCertificate]]
    """The specializations each step tested."""
    var orders: List[Int]
    var currents: List[SparsePoly[QQ]]
    """The polynomial after each step's derivatives."""


struct _LazardMemo(Movable):
    """The last evaluation of each lifted polynomial. Breadth-first lifting
    visits sibling cells together, and their points share every coordinate but
    the last, so each sibling's evaluation repeats all steps but the last."""
    var steps: List[_LazardSteps]

    def __init__(out self):
        self.steps = List[_LazardSteps]()


def _same_coordinate(a: Tuple[Int64, RealAlgebraic], b: Tuple[Int64, RealAlgebraic], mut budget: Budget) raises PolyError -> Bool:
    budget.charge(1)
    ref x = a[1]
    ref y = b[1]
    return a[0] == b[0] and x.rational == y.rational and _compare_q(x.low, y.low, budget) == 0 and _compare_q(x.high, y.high, budget) == 0 and _dense_equal(x.polynomial, y.polynomial, budget)


def _lazard_memo(f: SparsePoly[QQ], variable: Int64, point: List[Tuple[Int64, RealAlgebraic]], mut cache: _NormalizationCache, mut memo: _LazardMemo, mut budget: Budget) raises PolyError -> LazardEvaluationCertificate:
    """As _lazard_cached, reusing the steps of the memo's evaluation of f on the shared leading coordinates.

    CAD lifting builds and checks the point, so its coordinates' witnesses
    are not checked again here or at each step's specialization.
    """
    f.check_invariants(budget)
    _point(point, budget, trusted=True)
    if f.is_zero() or len(point) == 0:
        return _lazard_cached(f, variable, point, cache, budget)
    var slot = -1
    for i in range(len(memo.steps)):
        budget.charge(1)
        if memo.steps[i].variable == variable and memo.steps[i].f == f:
            slot = i
            break
    if slot < 0:
        memo.steps.append(_LazardSteps(f._copy(budget), variable, List[Tuple[Int64, RealAlgebraic]](), List[List[AlgebraicSpecializationCertificate]](), List[Int](), List[SparsePoly[QQ]]()))
        slot = len(memo.steps) - 1
    ref steps = memo.steps[slot]
    var shared = 0
    while shared < len(steps.coordinates) and shared < len(point) and _same_coordinate(steps.coordinates[shared], point[shared], budget):
        shared += 1
    while len(steps.coordinates) > shared:
        _ = steps.coordinates.pop()
        _ = steps.evaluations.pop()
        _ = steps.orders.pop()
        _ = steps.currents.pop()
    var current = steps.currents[shared - 1]._copy(budget) if shared else f._copy(budget)
    var prefix = List[Tuple[Int64, RealAlgebraic]]()
    for i in range(shared):
        prefix.append((point[i][0], _copy_real(point[i][1], budget)))
    for i in range(shared, len(point)):
        budget.charge(1)
        ref entry = point[i]
        if entry[0] == variable:
            raise PolyError(PolyErrorKind.DOMAIN, "Lazard lifting variable was assigned")
        _allocation[Tuple[Int64, RealAlgebraic]](1, budget)
        prefix.append((entry[0], _copy_real(entry[1], budget)))
        var order = 0
        var tested = List[AlgebraicSpecializationCertificate]()
        while True:
            budget.charge(1)
            var evaluation = _specialize_cached(current, prefix, cache, budget, trusted=True)
            var zero = _zero_specialization(evaluation, budget)
            _allocation[AlgebraicSpecializationCertificate](1, budget)
            tested.append(evaluation^)
            if not zero:
                break
            current = derivative(current, entry[0], budget)
            if current.is_zero():
                raise PolyError(PolyErrorKind.DOMAIN, "Lazard derivative chain unexpectedly vanished")
            order = _sum(order, 1)
        steps.coordinates.append((entry[0], _copy_real(entry[1], budget)))
        steps.evaluations.append(tested^)
        steps.orders.append(order)
        steps.currents.append(current._copy(budget))
    var orders = steps.orders.copy()
    var evaluations = List[AlgebraicSpecializationCertificate]()
    for tested in steps.evaluations:
        for evaluation in tested:
            budget.charge(1)
            evaluations.append(evaluation.copy())
    _ = _coefficients(evaluations[len(evaluations) - 1], variable, budget)
    return LazardEvaluationCertificate(orders^, evaluations^)


def _lazard_cached(f: SparsePoly[QQ], variable: Int64, point: List[Tuple[Int64, RealAlgebraic]], mut cache: _NormalizationCache, mut budget: Budget) raises PolyError -> LazardEvaluationCertificate:
    f.check_invariants(budget)
    _point(point, budget)
    if f.is_zero():
        raise PolyError(PolyErrorKind.DOMAIN, "Lazard evaluation requires nonzero input")
    var current = f._copy(budget)
    var prefix = List[Tuple[Int64, RealAlgebraic]]()
    var orders = List[Int]()
    var evaluations = List[AlgebraicSpecializationCertificate]()
    for entry in point:
        budget.charge(1)
        if entry[0] == variable:
            raise PolyError(PolyErrorKind.DOMAIN, "Lazard lifting variable was assigned")
        _allocation[Tuple[Int64, RealAlgebraic]](1, budget)
        prefix.append((entry[0], _copy_real(entry[1], budget)))
        var order = 0
        while True:
            budget.charge(1)
            var evaluation = _specialize_cached(current, prefix, cache, budget)
            var zero = _zero_specialization(evaluation, budget)
            _allocation[AlgebraicSpecializationCertificate](1, budget)
            evaluations.append(evaluation^)
            if not zero:
                break
            current = derivative(current, entry[0], budget)
            if current.is_zero():
                raise PolyError(PolyErrorKind.DOMAIN, "Lazard derivative chain unexpectedly vanished")
            order = _sum(order, 1)
        _allocation[Int](1, budget)
        orders.append(order)
    if len(point) == 0:
        _allocation[AlgebraicSpecializationCertificate](1, budget)
        evaluations.append(_specialize_cached(current, prefix, cache, budget))
    _ = _coefficients(evaluations[len(evaluations) - 1], variable, budget)
    return LazardEvaluationCertificate(orders^, evaluations^)


def check_lazard_evaluation(f: SparsePoly[QQ], variable: Int64, point: List[Tuple[Int64, RealAlgebraic]], certificate: LazardEvaluationCertificate, mut budget: Budget) raises PolyError -> Bool:
    """Verify minimal vanishing orders and every Taylor-coefficient specialization.

    Args:
        f: The original nonzero polynomial.
        variable: The final unassigned variable.
        point: The original ordered assignments.
        certificate: The proposed orders and specialization witnesses.
        budget: The shared work allowance.

    Returns:
        Whether the final specialization is the nonzero Lazard evaluation up to positive factorials.

    Raises:
        PolyError: DOMAIN for malformed inputs, LIMIT or CAPACITY on work limits.
    """
    f.check_invariants(budget)
    _point(point, budget)
    if f.is_zero() or len(certificate.orders) != len(point):
        return False
    var current = f._copy(budget)
    var prefix = List[Tuple[Int64, RealAlgebraic]]()
    var cursor = 0
    for i in range(len(point)):
        budget.charge(1)
        if point[i][0] == variable or certificate.orders[i] < 0 or certificate.orders[i] > Int(current.degree_in(point[i][0], budget)):
            return False
        _allocation[Tuple[Int64, RealAlgebraic]](1, budget)
        prefix.append((point[i][0], _copy_real(point[i][1], budget)))
        for k in range(_sum(certificate.orders[i], 1)):
            budget.charge(1)
            if cursor >= len(certificate.evaluations):
                return False
            ref evaluation = certificate.evaluations[cursor]
            if not check_algebraic_specialization(current, prefix, evaluation, budget) or _zero_specialization(evaluation, budget) != (k < certificate.orders[i]):
                return False
            cursor += 1
            if k < certificate.orders[i]:
                current = derivative(current, point[i][0], budget)
    if len(point) == 0:
        if len(certificate.evaluations) != 1 or not check_algebraic_specialization(current, prefix, certificate.evaluations[0], budget) or _zero_specialization(certificate.evaluations[0], budget):
            return False
        cursor = 1
    if cursor != len(certificate.evaluations):
        return False
    _ = _coefficients(certificate.evaluations[cursor - 1], variable, budget)
    return True
