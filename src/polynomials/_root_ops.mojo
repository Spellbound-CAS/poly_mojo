"""Exact scalar comparisons and sign-variation checks, independent of producers."""
from ._coefficient_ops import _copy
from .coefficients import QQ, NoContext
from .dense import DensePoly
from .errors import Budget, PolyError
from .evaluate import evaluate_dense, _derivative_dense


def _compare(a: QQ, b: QQ, mut budget: Budget) raises PolyError -> Int:
    return a.sub(b, budget).value.sign()


def _midpoint(a: QQ, b: QQ, mut budget: Budget) raises PolyError -> QQ:
    var two = QQ.one(NoContext()).add(QQ.one(NoContext()), budget)
    return a.add(b, budget).mul(two.inverse(budget), budget)


def _variations(polynomials: List[DensePoly[QQ]], a: QQ, mut budget: Budget, side: Int = 0) raises PolyError -> Int:
    var previous = 0
    var count = 0
    for p in polynomials:
        budget.charge(1)
        var sign = evaluate_dense(p, a, budget).value.sign()
        if sign == 0 and side != 0 and not p.is_zero():
            # At an endpoint root, the first nonzero Taylor coefficient
            # gives the exact one-sided sign, including repeated roots.
            var derivative = _derivative_dense(p, budget)
            var parity = 1
            while True:
                budget.charge(1)
                sign = evaluate_dense(derivative, a, budget).value.sign()
                if sign != 0:
                    if side < 0 and parity:
                        sign = -sign
                    break
                derivative = _derivative_dense(derivative, budget)
                parity = 1 - parity
        if sign != 0:
            if previous != 0 and previous != sign:
                count += 1
            previous = sign
    return count
