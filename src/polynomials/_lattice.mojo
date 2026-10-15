"""Exact rational Gram-Schmidt and integer LLL for van Hoeij recombination.

Fresh implementation of the lattice step in van Hoeij (2002), section 2.1.
All Gram-Schmidt calculations are exact; no floating-point decisions occur.
"""
from apn_mojo import Rational
from .coefficients import ZZ, QQ, NoContext, _work
from .errors import Budget, PolyError, PolyErrorKind, _allocation, _product, _charge_size
from ._root_ops import _compare
from ._coefficient_ops import _copy


@fieldwise_init
struct _Orthogonal(Movable):
    var mu: List[List[QQ]]
    var norms: List[QQ]


@fieldwise_init
struct _LLL(Movable):
    var basis: List[List[ZZ]]
    var transform: List[List[ZZ]]
    var mu: List[List[QQ]]
    var norms: List[QQ]


def _q(a: ZZ, mut budget: Budget) raises PolyError -> QQ:
    _charge_size(a.size(), budget)
    return QQ(Rational(a.value.copy()))


def _copy_matrix(a: List[List[ZZ]], mut budget: Budget) raises PolyError -> List[List[ZZ]]:
    _allocation[List[ZZ]](len(a), budget)
    var result = List[List[ZZ]]()
    for row in a:
        budget.charge(1)
        _allocation[ZZ](len(row), budget)
        var copied = List[ZZ]()
        for value in row:
            _charge_size(value.size(), budget)
            copied.append(value.copy())
        result.append(copied^)
    return result^


def _orthogonal(a: List[List[ZZ]], mut budget: Budget) raises PolyError -> _Orthogonal:
    var rows = len(a)
    var columns = len(a[0])
    _allocation[List[QQ]](rows, budget)
    _allocation[QQ](_product(rows, rows), budget)
    _allocation[QQ](_product(rows, columns), budget)
    _allocation[QQ](rows, budget)
    var mu = List[List[QQ]](length=rows, fill=List[QQ](length=rows, fill=QQ.zero(NoContext())))
    var orthogonal = List[List[QQ]]()
    var norms = List[QQ]()
    _allocation[QQ](rows, budget)
    var inverse_norms = List[QQ](capacity=rows)
    for i in range(rows):
        budget.charge(1)
        var vector = List[QQ]()
        for value in a[i]:
            budget.charge(1)
            vector.append(_q(value, budget))
        for j in range(i):
            budget.charge(1)
            var dot = QQ.zero(NoContext())
            for k in range(columns):
                budget.charge(1)
                dot = dot.add(_q(a[i][k], budget).mul(orthogonal[j][k], budget), budget)
            mu[i][j] = dot.mul(inverse_norms[j], budget)
            for k in range(columns):
                budget.charge(1)
                vector[k] = vector[k].sub(mu[i][j].mul(orthogonal[j][k], budget), budget)
        var norm = QQ.zero(NoContext())
        for value in vector:
            budget.charge(1)
            norm = norm.add(value.mul(value, budget), budget)
        if norm.is_zero():
            raise PolyError(PolyErrorKind.DOMAIN, "LLL requires independent nonzero rows")
        if i + 1 < rows:
            inverse_norms.append(norm.inverse(budget))
        norms.append(norm^)
        orthogonal.append(vector^)
    return _Orthogonal(mu^, norms^)


def _nearest(q: QQ, mut budget: Budget) raises PolyError -> ZZ:
    _work(q.size(), 1, True, budget)
    try:
        return ZZ(q.value.round())
    except e:
        raise PolyError(PolyErrorKind.CAPACITY, String(e))


def _swap_rows(mut a: List[List[ZZ]], i: Int, j: Int):
    var temporary = List[ZZ]()
    swap(temporary, a[i])
    swap(temporary, a[j])
    swap(temporary, a[i])


def _lll[carry: Bool = True](input: List[List[ZZ]], mut budget: Budget) raises PolyError -> _LLL:
    if len(input) == 0 or len(input) > len(input[0]):
        raise PolyError(PolyErrorKind.DOMAIN, "LLL requires a nonempty full-row-rank matrix")
    var rows = len(input)
    var columns = len(input[0])
    for row in input:
        budget.charge(1)
        if len(row) != columns:
            raise PolyError(PolyErrorKind.DOMAIN, "LLL matrix has unequal row lengths")
    var basis = _copy_matrix(input, budget)
    var transform = List[List[ZZ]]()
    comptime if carry:
        _allocation[List[ZZ]](rows, budget)
        _allocation[ZZ](_product(rows, rows), budget)
        transform = List[List[ZZ]](length=rows, fill=List[ZZ](length=rows, fill=ZZ.zero(NoContext())))
        for i in range(rows):
            budget.charge(1)
            transform[i][i] = ZZ.one(NoContext())
    var orthogonal = _orthogonal(basis, budget)
    var delta = QQ.one(NoContext()).sub(QQ(Rational(4)).inverse(budget), budget)
    var k = 1
    while k < rows:
        budget.charge(1)
        for j in range(k - 1, -1, -1):
            budget.charge(1)
            var nearest = _nearest(orthogonal.mu[k][j], budget)
            if nearest.is_zero():
                continue
            for column in range(columns):
                budget.charge(1)
                basis[k][column] = basis[k][column].sub(nearest.mul(basis[j][column], budget), budget)
            comptime if carry:
                for column in range(rows):
                    budget.charge(1)
                    transform[k][column] = transform[k][column].sub(nearest.mul(transform[j][column], budget), budget)
            # Subtracting an earlier row leaves every orthogonal vector
            # unchanged. Update only the affected Gram-Schmidt coefficients.
            var reduction = _q(nearest, budget)
            for l in range(j):
                budget.charge(1)
                orthogonal.mu[k][l] = orthogonal.mu[k][l].sub(reduction.mul(orthogonal.mu[j][l], budget), budget)
            orthogonal.mu[k][j] = orthogonal.mu[k][j].sub(reduction, budget)
        var rhs = delta.sub(orthogonal.mu[k][k - 1].mul(orthogonal.mu[k][k - 1], budget), budget).mul(orthogonal.norms[k - 1], budget)
        if _compare(orthogonal.norms[k], rhs, budget) >= 0:
            k += 1
        else:
            _swap_rows(basis, k, k - 1)
            comptime if carry:
                _swap_rows(transform, k, k - 1)
            # Exact adjacent-swap identities avoid a cubic recomputation.
            var old_mu = _copy(orthogonal.mu[k][k - 1], budget)
            var beta = orthogonal.norms[k].add(old_mu.mul(old_mu, budget).mul(orthogonal.norms[k - 1], budget), budget)
            var inverse_beta = beta.inverse(budget)
            var next_mu = old_mu.mul(orthogonal.norms[k - 1], budget).mul(inverse_beta, budget)
            orthogonal.norms[k] = orthogonal.norms[k - 1].mul(orthogonal.norms[k], budget).mul(inverse_beta, budget)
            orthogonal.norms[k - 1] = beta^
            for j in range(k - 1):
                budget.charge(1)
                var temporary = _copy(orthogonal.mu[k][j], budget)
                orthogonal.mu[k][j] = _copy(orthogonal.mu[k - 1][j], budget)
                orthogonal.mu[k - 1][j] = temporary^
            orthogonal.mu[k][k - 1] = _copy(next_mu, budget)
            for i in range(k + 1, rows):
                budget.charge(1)
                var old = _copy(orthogonal.mu[i][k], budget)
                orthogonal.mu[i][k] = orthogonal.mu[i][k - 1].sub(old_mu.mul(old, budget), budget)
                orthogonal.mu[i][k - 1] = old.add(next_mu.mul(orthogonal.mu[i][k], budget), budget)
            k = max(k - 1, 1)
    var mu = List[List[QQ]]()
    var norms = List[QQ]()
    swap(mu, orthogonal.mu)
    swap(norms, orthogonal.norms)
    return _LLL(basis^, transform^, mu^, norms^)
