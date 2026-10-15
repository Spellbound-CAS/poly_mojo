"""Fraction-free determinant transcripts for certified CAD projection."""
from .errors import Budget, PolyError, _allocation, _product, _sum
from .coefficients import QQ
from .sparse import SparsePoly
from ._views import _split, _zero, _one
from ._tower import _exact_poly
from .certificates import _equal


@fieldwise_init
struct DeterminantStep(Copyable, Movable, Writable):
    """One exact fraction-free elimination step on the current trailing matrix."""
    var pivot_row: Int
    """The row swapped into position zero, or minus one for an all-zero first column."""
    var trailing: List[List[SparsePoly[QQ]]]
    """The next trailing matrix, empty for a singular termination."""


@fieldwise_init
struct ProjectionResultantCertificate(Copyable, Movable, Writable):
    """A Sylvester determinant value with every exact elimination identity."""
    var variable: Int64
    """The eliminated variable."""
    var resultant: SparsePoly[QQ]
    """The exact determinant, including its sign."""
    var steps: List[DeterminantStep]
    """The fraction-free elimination transcript in forward order."""


def _sylvester(f: SparsePoly[QQ], g: SparsePoly[QQ], variable: Int64, mut budget: Budget) raises PolyError -> List[List[SparsePoly[QQ]]]:
    f._same(g, budget)
    var a = _split(f, variable, budget)
    var b = _split(g, variable, budget)
    var m = len(a) - 1
    var n = len(b) - 1
    var size = _sum(m, n)
    _allocation[List[SparsePoly[QQ]]](size, budget)
    _allocation[SparsePoly[QQ]](_product(size, size), budget)
    var matrix = List[List[SparsePoly[QQ]]]()
    for i in range(size):
        budget.charge(1)
        var row = List[SparsePoly[QQ]]()
        for j in range(size):
            budget.charge(1)
            var k = j - i if i < n else j - (i - n)
            if i < n and 0 <= k <= m:
                row.append(a[m - k]._copy(budget))
            elif i >= n and 0 <= k <= n:
                row.append(b[n - k]._copy(budget))
            else:
                row.append(_zero(f, budget))
        matrix.append(row^)
    return matrix^


def _row(index: Int, pivot: Int) -> Int:
    return pivot if index == 0 else (0 if index == pivot else index)


def _next_step(matrix: List[List[SparsePoly[QQ]]], previous: SparsePoly[QQ], mut budget: Budget) raises PolyError -> DeterminantStep:
    var size = len(matrix)
    var pivot = 0
    while pivot < size and matrix[pivot][0].is_zero():
        budget.charge(1)
        pivot += 1
    if pivot == size:
        return DeterminantStep(-1, List[List[SparsePoly[QQ]]]())
    _allocation[List[SparsePoly[QQ]]](size - 1, budget)
    _allocation[SparsePoly[QQ]](_product(size - 1, size - 1), budget)
    var trailing = List[List[SparsePoly[QQ]]](capacity=size - 1)
    for i in range(1, size):
        budget.charge(1)
        var row = List[SparsePoly[QQ]](capacity=size - 1)
        var source = _row(i, pivot)
        for j in range(1, size):
            budget.charge(1)
            var numerator = matrix[pivot][0].mul(matrix[source][j], budget).sub(matrix[source][0].mul(matrix[pivot][j], budget), budget)
            row.append(_exact_poly(numerator, previous, budget))
        trailing.append(row^)
    return DeterminantStep(pivot, trailing^)


def _check_step(matrix: List[List[SparsePoly[QQ]]], previous: SparsePoly[QQ], step: DeterminantStep, mut budget: Budget) raises PolyError -> Bool:
    var size = len(matrix)
    if size <= 1:
        return False
    if step.pivot_row == -1:
        if len(step.trailing):
            return False
        for row in matrix:
            budget.charge(1)
            if not row[0].is_zero():
                return False
        return True
    var pivot = step.pivot_row
    if pivot < 0 or pivot >= size or len(step.trailing) != size - 1 or matrix[pivot][0].is_zero():
        return False
    for i in range(1, size):
        budget.charge(1)
        if len(step.trailing[i - 1]) != size - 1:
            return False
        var source = _row(i, pivot)
        for j in range(1, size):
            budget.charge(1)
            var numerator = matrix[pivot][0].mul(matrix[source][j], budget).sub(matrix[source][0].mul(matrix[pivot][j], budget), budget)
            if not _equal(numerator, previous.mul(step.trailing[i - 1][j - 1], budget), budget):
                return False
    return True


def projection_resultant(f: SparsePoly[QQ], g: SparsePoly[QQ], variable: Int64, mut budget: Budget) raises PolyError -> ProjectionResultantCertificate:
    """Compute a resultant with a complete Sylvester determinant transcript.

    Args:
        f: The first rational polynomial.
        g: The second polynomial in the same stored order.
        variable: The variable to eliminate.
        budget: The shared work allowance.

    Returns:
        The exact resultant and fraction-free elimination steps.

    Raises:
        PolyError: DOMAIN for incompatible inputs, LIMIT or CAPACITY on work limits.
    """
    f._same(g, budget)
    var steps = List[DeterminantStep]()
    if f.is_zero() or g.is_zero():
        return ProjectionResultantCertificate(variable, _zero(f, budget), steps^)
    var initial = _sylvester(f, g, variable, budget)
    var size = len(initial)
    var previous = _one(f, budget)
    var negative = False
    while size > 1:
        budget.charge(1)
        var step = _next_step(initial, previous, budget) if len(steps) == 0 else _next_step(steps[len(steps) - 1].trailing, previous, budget)
        _allocation[DeterminantStep](1, budget)
        if step.pivot_row < 0:
            steps.append(step^)
            return ProjectionResultantCertificate(variable, _zero(f, budget), steps^)
        var pivot = step.pivot_row
        previous = initial[pivot][0]._copy(budget) if len(steps) == 0 else steps[len(steps) - 1].trailing[pivot][0]._copy(budget)
        if pivot:
            negative = not negative
        # No borrow of steps survives its growth. The certificate owns the only
        # copy of each trailing matrix; pivoting is an immutable row mapping.
        steps.append(step^)
        initial = List[List[SparsePoly[QQ]]]()
        size -= 1
    var value = _one(f, budget)
    if size:
        value = initial[0][0]._copy(budget) if len(steps) == 0 else steps[len(steps) - 1].trailing[0][0]._copy(budget)
    if negative:
        value = value.neg(budget)
    return ProjectionResultantCertificate(variable, value^, steps^)


def check_projection_resultant(f: SparsePoly[QQ], g: SparsePoly[QQ], certificate: ProjectionResultantCertificate, mut budget: Budget) raises PolyError -> Bool:
    """Check a Sylvester determinant by multiplication identities without division.

    Args:
        f: The original first polynomial.
        g: The original second polynomial.
        certificate: The proposed determinant and elimination transcript.
        budget: The shared work allowance.

    Returns:
        Whether the complete transcript proves the exact resultant.

    Raises:
        PolyError: DOMAIN for malformed inputs, LIMIT or CAPACITY on work limits.
    """
    f._same(g, budget)
    if f.is_zero() or g.is_zero():
        return len(certificate.steps) == 0 and _equal(certificate.resultant, _zero(f, budget), budget)
    var initial = _sylvester(f, g, certificate.variable, budget)
    var size = len(initial)
    var previous = _one(f, budget)
    var negative = False
    for index in range(len(certificate.steps)):
        budget.charge(1)
        ref step = certificate.steps[index]
        var valid = _check_step(initial, previous, step, budget) if index == 0 else _check_step(certificate.steps[index - 1].trailing, previous, step, budget)
        if not valid:
            return False
        if step.pivot_row == -1:
            return index == len(certificate.steps) - 1 and _equal(certificate.resultant, _zero(f, budget), budget)
        var pivot = step.pivot_row
        if pivot:
            negative = not negative
        previous = initial[pivot][0]._copy(budget) if index == 0 else certificate.steps[index - 1].trailing[pivot][0]._copy(budget)
        initial = List[List[SparsePoly[QQ]]]()
        size -= 1
    if size > 1:
        return False
    var value = _one(f, budget)
    if size:
        value = initial[0][0]._copy(budget) if len(certificate.steps) == 0 else certificate.steps[len(certificate.steps) - 1].trailing[0][0]._copy(budget)
    if negative:
        value = value.neg(budget)
    return _equal(value, certificate.resultant, budget)
