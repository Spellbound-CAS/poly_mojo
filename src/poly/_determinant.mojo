"""Generic Sylvester transcripts using public exact polynomial operations."""
from polynomials import CoefficientRing, SparsePoly, Term, Monomial, Budget, PolyError, PolyErrorKind, divide_exact
from ._common import _allocation, _product, _sum, _charge_size


def _same[C: CoefficientRing](a: SparsePoly[C], b: SparsePoly[C], mut budget: Budget) raises PolyError:
    _ = a.add(b, budget)


def _copy[C: CoefficientRing](p: SparsePoly[C], mut budget: Budget) raises PolyError -> SparsePoly[C]:
    return p.reorder(p.order, budget)


def _zero[C: CoefficientRing](p: SparsePoly[C], mut budget: Budget) raises PolyError -> SparsePoly[C]:
    return SparsePoly[C].zero(p.order, p.ring, budget)


def _one[C: CoefficientRing](p: SparsePoly[C], mut budget: Budget) raises PolyError -> SparsePoly[C]:
    return SparsePoly[C].constant(C.one(p.ring), p.order, budget)


def _equal[C: CoefficientRing](a: SparsePoly[C], b: SparsePoly[C], mut budget: Budget) raises PolyError -> Bool:
    return a.sub(b, budget).is_zero()


def _exact_poly[C: CoefficientRing](a: SparsePoly[C], b: SparsePoly[C], mut budget: Budget) raises PolyError -> SparsePoly[C]:
    var q = divide_exact(a, b, budget)
    if not q:
        raise PolyError(PolyErrorKind.DOMAIN, "determinant division was not exact")
    return q.take()


def _split[C: CoefficientRing](p: SparsePoly[C], variable: Int64, mut budget: Budget) raises PolyError -> List[SparsePoly[C]]:
    var count = _sum(Int(p.degree_in(variable, budget)), 1)
    _allocation[List[Term[C]]](count, budget)
    var groups = List[List[Term[C]]](capacity=count)
    for _ in range(count):
        groups.append(List[Term[C]]())
    for t in p.terms:
        _allocation[Tuple[Int64, UInt32]](len(t.monomial.pairs), budget)
        var pairs = List[Tuple[Int64, UInt32]]()
        var exponent = 0
        for pair in t.monomial.pairs:
            budget.charge(1)
            if pair[0] == variable:
                exponent = Int(pair[1])
            else:
                pairs.append(pair)
        _allocation[Term[C]](1, budget)
        _charge_size(t.coefficient.size(), budget)
        groups[exponent].append(Term[C](Monomial(pairs^, budget), t.coefficient.copy()))
    _allocation[SparsePoly[C]](count, budget)
    var result = List[SparsePoly[C]](capacity=count)
    for i in range(count):
        var terms = List[Term[C]]()
        swap(terms, groups[i])
        result.append(SparsePoly[C].from_terms(terms^, p.order, p.ring, budget))
    return result^


@fieldwise_init
struct Step[C: CoefficientRing](Copyable, Movable, Writable):
    """One exact fraction-free elimination step on the current trailing matrix."""
    var pivot_row: Int
    """The row swapped into position zero, or minus one for an all-zero first column."""
    var trailing: List[List[SparsePoly[Self.C]]]
    """The next trailing matrix, empty for a singular termination."""


@fieldwise_init
struct Certificate[C: CoefficientRing](Copyable, Movable, Writable):
    """A Sylvester determinant value with every exact elimination identity."""
    var variable: Int64
    """The eliminated variable."""
    var resultant: SparsePoly[Self.C]
    """The exact determinant, including its sign."""
    var steps: List[Step[Self.C]]
    """The fraction-free elimination transcript in forward order."""


def _sylvester[C: CoefficientRing](f: SparsePoly[C], g: SparsePoly[C], variable: Int64, mut budget: Budget) raises PolyError -> List[List[SparsePoly[C]]]:
    _same(f, g, budget)
    var a = _split(f, variable, budget)
    var b = _split(g, variable, budget)
    var m = len(a) - 1
    var n = len(b) - 1
    var size = _sum(m, n)
    _allocation[List[SparsePoly[C]]](size, budget)
    _allocation[SparsePoly[C]](_product(size, size), budget)
    var matrix = List[List[SparsePoly[C]]]()
    for i in range(size):
        budget.charge(1)
        var row = List[SparsePoly[C]]()
        for j in range(size):
            budget.charge(1)
            var k = j - i if i < n else j - (i - n)
            if i < n and 0 <= k <= m:
                row.append(_copy(a[m - k], budget))
            elif i >= n and 0 <= k <= n:
                row.append(_copy(b[n - k], budget))
            else:
                row.append(_zero(f, budget))
        matrix.append(row^)
    return matrix^


def _row(index: Int, pivot: Int) -> Int:
    return pivot if index == 0 else (0 if index == pivot else index)


def _next_step[C: CoefficientRing](matrix: List[List[SparsePoly[C]]], previous: SparsePoly[C], mut budget: Budget) raises PolyError -> Step[C]:
    var size = len(matrix)
    var pivot = 0
    while pivot < size and matrix[pivot][0].is_zero():
        budget.charge(1)
        pivot += 1
    if pivot == size:
        return Step[C](-1, List[List[SparsePoly[C]]]())
    _allocation[List[SparsePoly[C]]](size - 1, budget)
    _allocation[SparsePoly[C]](_product(size - 1, size - 1), budget)
    var trailing = List[List[SparsePoly[C]]](capacity=size - 1)
    for i in range(1, size):
        budget.charge(1)
        var row = List[SparsePoly[C]](capacity=size - 1)
        var source = _row(i, pivot)
        for j in range(1, size):
            budget.charge(1)
            var numerator = matrix[pivot][0].mul(matrix[source][j], budget).sub(matrix[source][0].mul(matrix[pivot][j], budget), budget)
            row.append(_exact_poly(numerator, previous, budget))
        trailing.append(row^)
    return Step[C](pivot, trailing^)


def _check_step[C: CoefficientRing](matrix: List[List[SparsePoly[C]]], previous: SparsePoly[C], step: Step[C], mut budget: Budget) raises PolyError -> Bool:
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


def _produce[C: CoefficientRing](f: SparsePoly[C], g: SparsePoly[C], variable: Int64, mut budget: Budget) raises PolyError -> Certificate[C]:
    """Compute a resultant with a complete Sylvester determinant transcript.

    Args:
        f: The first exact polynomial.
        g: The second polynomial in the same stored order.
        variable: The variable to eliminate.
        budget: The shared work allowance.

    Returns:
        The exact resultant and fraction-free elimination steps.

    Raises:
        PolyError: DOMAIN for incompatible inputs, LIMIT or CAPACITY on work limits.
    """
    _same(f, g, budget)
    var steps = List[Step[C]]()
    if f.is_zero() or g.is_zero():
        return Certificate[C](variable, _zero(f, budget), steps^)
    var initial = _sylvester(f, g, variable, budget)
    var size = len(initial)
    var previous = _one(f, budget)
    var negative = False
    while size > 1:
        budget.charge(1)
        var step = _next_step(initial, previous, budget) if len(steps) == 0 else _next_step(steps[len(steps) - 1].trailing, previous, budget)
        _allocation[Step[C]](1, budget)
        if step.pivot_row < 0:
            steps.append(step^)
            return Certificate[C](variable, _zero(f, budget), steps^)
        var pivot = step.pivot_row
        previous = _copy(initial[pivot][0], budget) if len(steps) == 0 else _copy(steps[len(steps) - 1].trailing[pivot][0], budget)
        if pivot:
            negative = not negative
        # No borrow of steps survives its growth. The certificate owns the only
        # copy of each trailing matrix; pivoting is an immutable row mapping.
        steps.append(step^)
        initial = List[List[SparsePoly[C]]]()
        size -= 1
    var value = _one(f, budget)
    if size:
        value = _copy(initial[0][0], budget) if len(steps) == 0 else _copy(steps[len(steps) - 1].trailing[0][0], budget)
    if negative:
        value = value.neg(budget)
    return Certificate[C](variable, value^, steps^)


def _check[C: CoefficientRing](f: SparsePoly[C], g: SparsePoly[C], certificate: Certificate[C], mut budget: Budget) raises PolyError -> Bool:
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
    _same(f, g, budget)
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
        previous = _copy(initial[pivot][0], budget) if index == 0 else _copy(certificate.steps[index - 1].trailing[pivot][0], budget)
        initial = List[List[SparsePoly[C]]]()
        size -= 1
    if size > 1:
        return False
    var value = _one(f, budget)
    if size:
        value = _copy(initial[0][0], budget) if len(certificate.steps) == 0 else _copy(certificate.steps[len(certificate.steps) - 1].trailing[0][0], budget)
    if negative:
        value = value.neg(budget)
    return _equal(value, certificate.resultant, budget)
