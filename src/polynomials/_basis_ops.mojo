"""Scalar-vector recurrences for Newton, Chebyshev and Bernstein bases."""
from .coefficient_buffer import _CoefficientView, _view
from .errors import Budget, PolyError, PolyErrorKind, _allocation, _sum
from ._toolkit_ops import _ScalarOps, _values, _add_values, _mul_values, _newton_values


def _basis_shape(kind: UInt8, count: Int, nodes: Int, mut budget: Budget) raises PolyError:
    budget.charge(1)
    if kind > 2 or (kind == 0 and nodes != max(0, count - 1)) or (kind != 0 and nodes != 0):
        raise PolyError(PolyErrorKind.DOMAIN, "invalid polynomial basis or node count")
    if kind == 2 and count == 0:
        raise PolyError(PolyErrorKind.DOMAIN, "Bernstein coefficients must retain a nonnegative degree")


def _binomial_row[O: _ScalarOps](n: Int, context: O.Context, mut budget: Budget) raises PolyError -> List[O.Scalar]:
    _allocation[O.Scalar](_sum(n, 1), budget)
    var row = List[O.Scalar](length=n + 1, fill=O.integer(0, context, budget))
    row[0] = O.integer(1, context, budget)
    for k in range(1, n + 1):
        for j in range(k, 0, -1):
            row[j] = O.apply(row[j], row[j - 1], 0, context, budget)
    return row^


def _newton_coefficients[O: _ScalarOps](values: _CoefficientView[O.Scalar, _], nodes: _CoefficientView[O.Scalar, _], context: O.Context, mut budget: Budget) raises PolyError -> List[O.Scalar]:
    var count = _sum(len(nodes), 1)
    if len(values) > count:
        raise PolyError(PolyErrorKind.DOMAIN, "insufficient Newton nodes for polynomial degree")
    var remaining = _values[O](_view(values), count, context, budget)
    _allocation[O.Scalar](count, budget)
    var result = List[O.Scalar](capacity=count)
    for k in range(len(nodes)):
        var quotient = _values[O](_view(remaining), len(remaining) - 1, context, budget)
        var carry = O.copy(remaining[len(remaining) - 1], budget)
        for j in range(len(remaining) - 2, -1, -1):
            quotient[j] = O.copy(carry, budget)
            carry = O.apply(remaining[j], O.apply(nodes[k], carry, 2, context, budget), 0, context, budget)
        result.append(carry^)
        remaining = quotient^
    result.append(O.copy(remaining[0], budget))
    return result^


def _chebyshev_exact[O: _ScalarOps](values: _CoefficientView[O.Scalar, _], context: O.Context, mut budget: Budget) raises PolyError -> List[O.Scalar]:
    var n = len(values)
    if n <= 2:
        return _values[O](_view(values), n, context, budget)
    var previous: List[O.Scalar] = [O.integer(1, context, budget)]
    var current: List[O.Scalar] = [O.integer(0, context, budget), O.integer(1, context, budget)]
    for k in range(2, n):
        var next = _values[O](_view(List[O.Scalar]()), k + 1, context, budget)
        for j in range(k % 2, k + 1, 2):
            if j:
                next[j] = O.apply(current[j - 1], current[j - 1], 0, context, budget)
            if j < len(previous):
                next[j] = O.apply(next[j], previous[j], 1, context, budget)
        previous = current^
        current = next^
    var remaining = _values[O](_view(values), n, context, budget)
    var result = _values[O](_view(List[O.Scalar]()), n, context, budget)
    for k in range(n - 1, -1, -1):
        result[k] = O.apply(remaining[k], current[k], 3, context, budget)
        for j in range(k % 2, k + 1, 2):
            remaining[j] = O.apply(remaining[j], O.apply(result[k], current[j], 2, context, budget), 1, context, budget)
        if k >= 2:
            # Reverse the exact recurrence; numerical domains use the forward
            # Horner schedule to avoid cancellation-induced enclosure growth.
            var next = _values[O](_view(List[O.Scalar]()), k - 1, context, budget)
            for j in range(k % 2, k - 1, 2):
                if j:
                    next[j] = O.apply(previous[j - 1], previous[j - 1], 0, context, budget)
                next[j] = O.apply(next[j], current[j], 1, context, budget)
            current = previous^
            previous = next^
        elif k == 1:
            swap(current, previous)
    return result^


def _chebyshev_coefficients[O: _ScalarOps](values: _CoefficientView[O.Scalar, _], context: O.Context, mut budget: Budget) raises PolyError -> List[O.Scalar]:
    comptime if O.EXACT:
        return _chebyshev_exact[O](_view(values), context, budget)
    # Horner in the Chebyshev basis: x*T0=T1 and
    # x*Tk=(T(k-1)+T(k+1))/2 for k>0. Only the current row survives.
    var row = List[O.Scalar]()
    var two = O.integer(2, context, budget)
    for i in range(len(values) - 1, -1, -1):
        var next = _values[O](_view(List[O.Scalar]()), _sum(len(row), 1), context, budget)
        if len(row):
            next[1] = O.copy(row[0], budget)
        for k in range(1, len(row)):
            var half = O.apply(row[k], two, 3, context, budget)
            next[k - 1] = O.apply(next[k - 1], half, 0, context, budget)
            next[k + 1] = O.copy(half, budget)
        next[0] = O.apply(next[0], values[i], 0, context, budget)
        row = next^
    return row^


def _bernstein_coefficients[O: _ScalarOps](values: _CoefficientView[O.Scalar, _], degree: Int, context: O.Context, mut budget: Budget) raises PolyError -> List[O.Scalar]:
    if degree < 0 or degree < len(values) - 1:
        raise PolyError(PolyErrorKind.DOMAIN, "Bernstein degree below polynomial degree")
    var n = _sum(degree, 1)
    var denominators = _binomial_row[O](degree, context, budget)
    var scaled = _values[O](_view(values), len(values), context, budget)
    for j in range(len(values)):
        scaled[j] = O.apply(scaled[j], denominators[j], 3, context, budget)
    var row = _values[O](_view(List[O.Scalar]()), n, context, budget)
    row[0] = O.integer(1, context, budget)
    var result = _values[O](_view(List[O.Scalar]()), n, context, budget)
    for k in range(n):
        for j in range(min(k + 1, len(scaled))):
            result[k] = O.apply(result[k], O.apply(scaled[j], row[j], 2, context, budget), 0, context, budget)
        for j in range(min(k + 1, degree), 0, -1):
            row[j] = O.apply(row[j], row[j - 1], 0, context, budget)
    return result^


def _basis_values[O: _ScalarOps](kind: UInt8, coefficients: _CoefficientView[O.Scalar, _], nodes: _CoefficientView[O.Scalar, _], context: O.Context, mut budget: Budget) raises PolyError -> List[O.Scalar]:
    _basis_shape(kind, len(coefficients), len(nodes), budget)
    if kind == 0:
        return _newton_values[O](_view(coefficients), _view(nodes), context, budget)
    var n = len(coefficients)
    if kind == 2:
        var row = _values[O](_view(coefficients), n, context, budget)
        var choose = _binomial_row[O](n - 1, context, budget)
        _allocation[O.Scalar](n, budget)
        var result = List[O.Scalar](capacity=n)
        for k in range(n):
            result.append(O.apply(choose[k], row[0], 2, context, budget))
            for j in range(n - k - 1):
                row[j] = O.apply(row[j + 1], row[j], 1, context, budget)
        return result^
    var b1 = List[O.Scalar]()
    var b2 = List[O.Scalar]()
    var twice_x: List[O.Scalar] = [O.integer(0, context, budget), O.integer(2, context, budget)]
    for k in range(n - 1, 0, -1):
        var current = _add_values[O](_view(_mul_values[O](_view(b1), _view(twice_x), n, context, budget)), _view(b2), n, True, context, budget)
        if not len(current):
            _allocation[O.Scalar](1, budget)
            current.append(O.copy(coefficients[k], budget))
        else:
            current[0] = O.apply(current[0], coefficients[k], 0, context, budget)
        b2 = b1^
        b1 = current^
    var x: List[O.Scalar] = [O.integer(0, context, budget), O.integer(1, context, budget)]
    var result = _add_values[O](_view(_mul_values[O](_view(b1), _view(x), n, context, budget)), _view(b2), n, True, context, budget)
    if n:
        if not len(result):
            _allocation[O.Scalar](1, budget)
            result.append(O.copy(coefficients[0], budget))
        else:
            result[0] = O.apply(result[0], coefficients[0], 0, context, budget)
    return result^


def _bernstein_step[O: _ScalarOps](mut row: List[O.Scalar], count: Int, point: O.Scalar, complement: O.Scalar, context: O.Context, mut budget: Budget) raises PolyError:
    for j in range(count):
        row[j] = O.apply(O.apply(complement, row[j], 2, context, budget), O.apply(point, row[j + 1], 2, context, budget), 0, context, budget)


def _bernstein_split[O: _ScalarOps](coefficients: _CoefficientView[O.Scalar, _], point: O.Scalar, context: O.Context, mut budget: Budget) raises PolyError -> Tuple[List[O.Scalar], List[O.Scalar]]:
    if not len(coefficients):
        raise PolyError(PolyErrorKind.DOMAIN, "empty Bernstein representation")
    var n = len(coefficients)
    var row = _values[O](_view(coefficients), n, context, budget)
    var left = _values[O](_view(List[O.Scalar]()), n, context, budget)
    var right = _values[O](_view(List[O.Scalar]()), n, context, budget)
    var complement = O.apply(O.integer(1, context, budget), point, 1, context, budget)
    for k in range(n):
        left[k] = O.copy(row[0], budget)
        right[n - 1 - k] = O.copy(row[n - 1 - k], budget)
        _bernstein_step[O](row, n - 1 - k, point, complement, context, budget)
    return (left^, right^)


def _basis_evaluate[O: _ScalarOps](kind: UInt8, coefficients: _CoefficientView[O.Scalar, _], nodes: _CoefficientView[O.Scalar, _], point: O.Scalar, context: O.Context, mut budget: Budget) raises PolyError -> O.Scalar:
    _basis_shape(kind, len(coefficients), len(nodes), budget)
    if kind == 2:
        var row = _values[O](_view(coefficients), len(coefficients), context, budget)
        var complement = O.apply(O.integer(1, context, budget), point, 1, context, budget)
        for count in range(len(row) - 1, 0, -1):
            _bernstein_step[O](row, count, point, complement, context, budget)
        return O.copy(row[0], budget)
    var b1 = O.integer(0, context, budget)
    var b2 = O.integer(0, context, budget)
    if kind == 0:
        for k in range(len(coefficients) - 1, -1, -1):
            if k < len(nodes):
                b1 = O.apply(b1, O.apply(point, nodes[k], 1, context, budget), 2, context, budget)
            b1 = O.apply(b1, coefficients[k], 0, context, budget)
        return b1^
    var twice = O.apply(O.integer(2, context, budget), point, 2, context, budget)
    for k in range(len(coefficients) - 1, 0, -1):
        var current = O.apply(coefficients[k], O.apply(O.apply(twice, b1, 2, context, budget), b2, 1, context, budget), 0, context, budget)
        b2 = b1^
        b1 = current^
    var value = O.apply(O.apply(point, b1, 2, context, budget), b2, 1, context, budget)
    return O.apply(coefficients[0], value, 0, context, budget) if len(coefficients) else value^
