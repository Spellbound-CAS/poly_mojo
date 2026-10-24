"""Private prototype: construct a rational CAD stack, verify its final data.

This first slice covers one variable over QQ, including shared roots and
multiplicities. The producer stores no specialization, division, comparison,
or isolation histories. The independent verifier reconstructs identities and
Sturm counts from the input and final stack. No public CAD path uses it yet.
"""
from apn_mojo import Integer, Rational
from .coefficients import QQ, NoContext
from .dense import DensePoly
from .division import divmod
from .squarefree import squarefree
from .evaluate import _derivative_dense
from ._dense_ops import _gcd, _exact, _remainder
from ._root_ops import _compare, _midpoint, _sign_at
from .errors import Budget, PolyError, PolyErrorKind, _allocation, _sum, _product
from ._cad_final_signs import _SignStats, _signs_at


@fieldwise_init
struct _FinalRoot(Copyable, Movable):
    var basis: Int
    var low: QQ
    var high: QQ


@fieldwise_init
struct _FinalStack(Copyable, Movable):
    var basis: List[DensePoly[QQ]]
    var leading: List[QQ]
    var multiplicities: List[List[Int]]
    var roots: List[_FinalRoot]
    var sectors: List[QQ]
    var signs: List[List[Int]]


def _constant(variable: Int64, value: QQ, mut budget: Budget) raises PolyError -> DensePoly[QQ]:
    return DensePoly[QQ].from_coefficients(variable, NoContext(), [value.copy()], budget)


def _sequence(f: DensePoly[QQ], mut budget: Budget) raises PolyError -> List[DensePoly[QQ]]:
    """An ordinary Sturm computation; no quotient witnesses or replay trace."""
    _allocation[DensePoly[QQ]](1, budget)
    var sequence: List[DensePoly[QQ]] = [f.copy()]
    if f.degree() > 0:
        _allocation[DensePoly[QQ]](1, budget)
        sequence.append(_derivative_dense(f, budget))
        while True:
            budget.charge(1)
            var n = len(sequence)
            var remainder = _remainder(sequence[n - 2], sequence[n - 1], budget)
            if remainder.is_zero():
                break
            _allocation[DensePoly[QQ]](1, budget)
            sequence.append(remainder.neg(budget))
    return sequence^


def _variation(sequence: List[DensePoly[QQ]], x: QQ, mut budget: Budget) raises PolyError -> Int:
    var previous = 0
    var changes = 0
    for p in sequence:
        var sign = _sign_at(p, x, budget)
        if sign:
            changes += Int(previous != 0 and sign != previous)
            previous = sign
    return changes


def _total(sequence: List[DensePoly[QQ]], mut budget: Budget) raises PolyError -> Int:
    var previous_low = 0
    var previous_high = 0
    var low = 0
    var high = 0
    for p in sequence:
        budget.charge(1)
        var positive = p.coefficients[p.degree()].value.sign()
        var negative = -positive if p.degree() % 2 else positive
        low += Int(previous_low != 0 and previous_low != negative)
        high += Int(previous_high != 0 and previous_high != positive)
        previous_low = negative
        previous_high = positive
    return low - high


def _basis(inputs: List[DensePoly[QQ]], mut budget: Budget) raises PolyError -> List[DensePoly[QQ]]:
    var basis = List[DensePoly[QQ]]()
    for f in inputs:
        f.check_invariants(budget)
        if len(inputs) and f.variable != inputs[0].variable:
            raise PolyError(PolyErrorKind.DOMAIN, "final stack requires a common variable")
        if f.degree() <= 0:
            continue
        # This decomposition holds only final factors, not gcd/division proofs.
        # The radical alone would conflate different multiplicities within f.
        var blocks = squarefree(f, budget)
        for block in blocks.factors:
            var todo = block[0].copy()
            var refined = List[DensePoly[QQ]]()
            for b in basis:
                var common = _gcd(b, todo, budget)
                _allocation[DensePoly[QQ]](2, budget)
                if common.degree() <= 0:
                    refined.append(b.copy())
                else:
                    var rest = _exact(b, common, budget)
                    todo = _exact(todo, common, budget)
                    refined.append(common^)
                    if rest.degree() > 0:
                        refined.append(rest^)
            if todo.degree() > 0:
                _allocation[DensePoly[QQ]](1, budget)
                refined.append(todo^)
            basis = refined^
    return basis^


def _bound(f: DensePoly[QQ], mut budget: Budget) raises PolyError -> QQ:
    # Monic Cauchy bound, strictly enlarged to avoid roots at its endpoints.
    var bound = Integer(2)
    for c in f.coefficients:
        budget.charge(1)
        try:
            var value = abs(c.value.numerator()) // c.value.denominator() + Integer(2)
            if value > bound:
                bound = value^
        except e:
            raise PolyError(PolyErrorKind.CAPACITY, String(e))
    return QQ(Rational(bound))


@fieldwise_init
struct _Pending(Copyable, Movable):
    var low: QQ
    var high: QQ
    var left: Int
    var right: Int


def _isolate(basis: List[DensePoly[QQ]], variable: Int64, mut budget: Budget) raises PolyError -> List[_FinalRoot]:
    var product = _constant(variable, QQ.one(NoContext()), budget)
    var sequences = List[List[DensePoly[QQ]]]()
    for b in basis:
        product = product.mul(b, budget)
        _allocation[List[DensePoly[QQ]]](1, budget)
        sequences.append(_sequence(b, budget))
    var chain = _sequence(product, budget)
    var high = _bound(product, budget)
    var low = high.neg(budget)
    _allocation[_Pending](1, budget)
    var pending: List[_Pending] = [_Pending(low.copy(), high.copy(), _variation(chain, low, budget), _variation(chain, high, budget))]
    var roots = List[_FinalRoot]()
    while len(pending):
        budget.charge(1)
        var node = pending.pop()
        var count = node.left - node.right
        if count == 0:
            continue
        if count == 1:
            var owner = -1
            for k in range(len(basis)):
                if _variation(sequences[k], node.low, budget) - _variation(sequences[k], node.high, budget) == 1:
                    owner = k
                    break
            if owner < 0:
                raise PolyError(PolyErrorKind.DOMAIN, "root has no basis member")
            _allocation[_FinalRoot](1, budget)
            roots.append(_FinalRoot(owner, node.low.copy(), node.high.copy()))
            continue
        var middle = _midpoint(node.low, node.high, budget)
        # A rational root is permitted, but this producer keeps open brackets.
        # Avoid making a root an endpoint; finitely many roots ensure progress.
        while _sign_at(product, middle, budget) == 0:
            middle = _midpoint(node.low, middle, budget)
        var variation = _variation(chain, middle, budget)
        _allocation[_Pending](2, budget)
        pending.append(_Pending(middle.copy(), node.high.copy(), variation, node.right))
        pending.append(_Pending(node.low.copy(), middle.copy(), node.left, variation))
    return roots^


def _labels(leading: List[QQ], multiplicities: List[List[Int]], roots: List[_FinalRoot], mut budget: Budget) raises PolyError -> List[List[Int]]:
    var count = _sum(_product(2, len(roots)), 1)
    _allocation[List[Int]](count, budget)
    var labels = List[List[Int]]()
    for _ in range(count):
        _allocation[Int](len(leading), budget)
        labels.append(List[Int](length=len(leading), fill=0))
    _allocation[Int](len(leading), budget)
    var current = List[Int]()
    for c in leading:
        budget.charge(1)
        current.append(c.value.sign())
    _allocation[Int](len(leading), budget)
    labels[count - 1] = current.copy()
    for k in range(len(roots) - 1, -1, -1):
        _allocation[Int](_product(2, len(leading)), budget)
        var section = current.copy()
        for j in range(len(leading)):
            budget.charge(1)
            var exponent = multiplicities[j][roots[k].basis]
            if exponent > 0:
                section[j] = 0
            if exponent % 2:
                current[j] = -current[j]
        labels[2 * k + 1] = section^
        labels[2 * k] = current.copy()
    return labels^


def _construct_stack(inputs: List[DensePoly[QQ]], mut budget: Budget) raises PolyError -> _FinalStack:
    """Untrusted rational-stack producer, retaining final algebra and geometry only."""
    var basis = _basis(inputs, budget)
    var leading = List[QQ]()
    var multiplicities = List[List[Int]]()
    for f in inputs:
        _allocation[QQ](1, budget)
        leading.append(f.coefficients[f.degree()].copy() if not f.is_zero() else QQ.zero(NoContext()))
        _allocation[Int](len(basis), budget)
        var row = List[Int](length=len(basis), fill=0)
        var rest = f.copy()
        if not rest.is_zero():
            for k in range(len(basis)):
                while rest.degree() >= basis[k].degree():
                    var step = divmod(rest, basis[k], budget)
                    if not step.remainder.is_zero():
                        break
                    rest = step.quotient.copy()
                    row[k] += 1
        _allocation[List[Int]](1, budget)
        multiplicities.append(row^)
    var roots = _isolate(basis, inputs[0].variable if len(inputs) else 0, budget)
    var sectors = List[QQ]()
    _allocation[QQ](_sum(len(roots), 1), budget)
    var one = QQ.one(NoContext())
    if len(roots) == 0:
        sectors.append(QQ.zero(NoContext()))
    else:
        sectors.append(roots[0].low.sub(one, budget))
        for i in range(1, len(roots)):
            sectors.append(_midpoint(roots[i - 1].high, roots[i].low, budget))
        sectors.append(roots[len(roots) - 1].high.add(one, budget))
    var signs = _labels(leading, multiplicities, roots, budget)
    return _FinalStack(basis^, leading^, multiplicities^, roots^, sectors^, signs^)


def _check_algebra(inputs: List[DensePoly[QQ]], stack: _FinalStack, mut budget: Budget) raises PolyError -> Bool:
    if len(stack.leading) != len(inputs) or len(stack.multiplicities) != len(inputs):
        return False
    var variable = inputs[0].variable if len(inputs) else 0
    _allocation[Bool](len(stack.basis), budget)
    var used = List[Bool](length=len(stack.basis), fill=False)
    for k in range(len(stack.basis)):
        ref b = stack.basis[k]
        b.check_invariants(budget)
        if b.variable != variable or b.degree() <= 0 or b.coefficients[b.degree()] != QQ.one(NoContext()):
            return False
        if _gcd(b, _derivative_dense(b, budget), budget).degree() != 0:
            return False
        for j in range(k):
            if _gcd(b, stack.basis[j], budget).degree() != 0:
                return False
    for j in range(len(inputs)):
        ref f = inputs[j]
        f.check_invariants(budget)
        if f.variable != variable or len(stack.multiplicities[j]) != len(stack.basis):
            return False
        var expected_lead = f.coefficients[f.degree()].copy() if not f.is_zero() else QQ.zero(NoContext())
        if stack.leading[j] != expected_lead:
            return False
        var remaining = max(0, f.degree())
        for k in range(len(stack.basis)):
            budget.charge(1)
            var exponent = stack.multiplicities[j][k]
            if exponent < 0 or exponent > remaining // stack.basis[k].degree():
                return False
            remaining -= exponent * stack.basis[k].degree()
            if exponent:
                used[k] = True
        if remaining:
            return False
        var product = _constant(variable, stack.leading[j], budget)
        for k in range(len(stack.basis)):
            for _ in range(stack.multiplicities[j][k]):
                product = product.mul(stack.basis[k], budget)
        if product != f:
            return False
    for present in used:
        if not present:
            return False
    return True


def _check_geometry(stack: _FinalStack, mut budget: Budget) raises PolyError -> Bool:
    if len(stack.sectors) != _sum(len(stack.roots), 1) or len(stack.signs) != _sum(_product(2, len(stack.roots)), 1):
        return False
    for row in stack.signs:
        if len(row) != len(stack.leading):
            return False
        for sign in row:
            budget.charge(1)
            if sign < -1 or sign > 1:
                return False
    for k in range(len(stack.roots)):
        ref root = stack.roots[k]
        if root.basis < 0 or root.basis >= len(stack.basis) or _compare(root.low, root.high, budget) > 0:
            return False
        if k and _compare(stack.roots[k - 1].high, root.low, budget) > 0:
            return False
        var exact = _compare(root.low, root.high, budget) == 0
        var before = _compare(stack.sectors[k], root.low, budget)
        var after = _compare(stack.sectors[k + 1], root.high, budget)
        if before > 0 or after < 0 or (exact and (before == 0 or after == 0)):
            return False
    return True


def _check_roots(stack: _FinalStack, mut stats: _SignStats, mut budget: Budget, vectorized: Bool) raises PolyError -> Bool:
    var endpoints = List[QQ]()
    _allocation[QQ](_product(2, len(stack.roots)), budget)
    for root in stack.roots:
        endpoints.append(root.low.copy())
        endpoints.append(root.high.copy())
    for k in range(len(stack.basis)):
        var sequence = _sequence(stack.basis[k], budget)
        _allocation[Int](_product(3, len(endpoints)), budget)
        var variations = List[Int](length=len(endpoints), fill=0)
        var previous = List[Int](length=len(endpoints), fill=0)
        var boundary_signs = List[Int]()
        for i in range(len(sequence)):
            var signs = _signs_at(sequence[i], endpoints, stats, budget, vectorized=vectorized)
            if i == 0:
                boundary_signs = signs.copy()
            for j in range(len(signs)):
                budget.charge(1)
                if signs[j]:
                    variations[j] += Int(previous[j] != 0 and previous[j] != signs[j])
                    previous[j] = signs[j]
        var counted = 0
        for i in range(len(stack.roots)):
            budget.charge(1)
            ref root = stack.roots[i]
            var expected = Int(root.basis == k)
            if _compare(root.low, root.high, budget) == 0:
                if Int(boundary_signs[2 * i] == 0) != expected:
                    return False
            else:
                if boundary_signs[2 * i] == 0 or boundary_signs[2 * i + 1] == 0:
                    return False
                if variations[2 * i] - variations[2 * i + 1] != expected:
                    return False
            counted += expected
        if counted != _total(sequence, budget):
            return False
    return True


def _check_labels(stack: _FinalStack, mut budget: Budget) raises PolyError -> Bool:
    # Check each input's sign from the top down, independently of _labels.
    for j in range(len(stack.leading)):
        var expected = stack.leading[j].value.sign()
        if stack.signs[len(stack.signs) - 1][j] != expected:
            return False
        for k in range(len(stack.roots) - 1, -1, -1):
            budget.charge(1)
            var exponent = stack.multiplicities[j][stack.roots[k].basis]
            if stack.signs[2 * k + 1][j] != (0 if exponent else expected):
                return False
            if exponent & 1:
                expected = -expected
            if stack.signs[2 * k][j] != expected:
                return False
    return True


def _check_final_stack(inputs: List[DensePoly[QQ]], stack: _FinalStack, mut stats: _SignStats, mut budget: Budget, *, vectorized: Bool = True) raises PolyError -> Bool:
    """Verify the entire rational stack from inputs; never trust producer state."""
    try:
        return _check_algebra(inputs, stack, budget) and _check_geometry(stack, budget) and _check_roots(stack, stats, budget, vectorized) and _check_labels(stack, budget)
    except e:
        if e.kind == PolyErrorKind.DOMAIN:
            return False
        raise e^


def _verified_stack(inputs: List[DensePoly[QQ]], mut budget: Budget) raises PolyError -> _FinalStack:
    var result = _construct_stack(inputs, budget)
    var stats = _SignStats()
    if not _check_final_stack(inputs, result, stats, budget):
        raise PolyError(PolyErrorKind.DOMAIN, "final stack verification failed")
    return result^
