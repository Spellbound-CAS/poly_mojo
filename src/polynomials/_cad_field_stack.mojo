"""Private final-data stacks over a shared ordered algebraic field.

Sections remain basis indices and rational isolating brackets. Construction
uses exact Sturm subdivision; verification counts roots from the final data
and at infinity, without replaying subdivisions. No absolute root is built.
"""
from apn_mojo import Integer, Rational, gcd as integer_gcd
from .coefficients import QQ, NoContext
from .dense import DensePoly
from .errors import Budget, PolyError, PolyErrorKind, _allocation, _sum, _product
from ._root_ops import _compare, _midpoint
from ._algebraic_filter import _constant, _gcd_over_field, _dense_interval
from ._cad_sample_field import _SampleField, _check_sample_field, _ksign, _kmul, _kscale, _kmonic, _kderivative, _kdivmod, _kexact, _kevaluate, _kcopy, _ktrim, _kmultiply, _kevaluate_sign
from ._algebraic_filter import _field_mul
from ._images import _denominator
from ._cad_final import _FinalRoot, _Pending
from ._cad_lifting import _simplest
from ._cad_interval_roots import _interval_field_roots, _IntervalRootStats
from ._cad_field_verify import _rational_root_check
from ._cad_root_cache import _RootCountCache, _check_rational_roots
from ._cad_final_signs import _SignStats


@fieldwise_init
struct _FieldStack(Copyable, Movable):
    var basis: List[List[DensePoly[QQ]]]
    var leading: List[Int8]
    var multiplicities: List[List[Int]]
    var roots: List[_FinalRoot]
    var sectors: List[QQ]
    var signs: List[List[Int8]]
    var isolation: _IntervalRootStats


def _kprimitive_positive(f: List[DensePoly[QQ]], mut budget: Budget) raises PolyError -> List[DensePoly[QQ]]:
    """Remove rational content without changing the sign at any real embedding."""
    var denominator = Integer(1)
    for row in f:
        for c in row.coefficients:
            denominator = _denominator(denominator, c, budget)
    var content = Integer(0)
    var integers = List[List[Integer]]()
    for row in f:
        _allocation[Integer](len(row.coefficients), budget)
        var values = List[Integer]()
        for c in row.coefficients:
            budget.charge(1)
            try:
                var value = c.value.numerator() * (denominator // c.value.denominator())
                content = integer_gcd(content, value)
                values.append(value^)
            except e:
                raise PolyError(PolyErrorKind.CAPACITY, String(e))
        _allocation[List[Integer]](1, budget)
        integers.append(values^)
    if not content:
        return _kcopy(f, budget)
    _allocation[DensePoly[QQ]](len(f), budget)
    var result = List[DensePoly[QQ]]()
    for row in integers:
        _allocation[QQ](len(row), budget)
        var values = List[QQ]()
        for value in row:
            budget.charge(1)
            try:
                values.append(QQ(Rational(value // content)))
            except e:
                raise PolyError(PolyErrorKind.CAPACITY, String(e))
        result.append(DensePoly[QQ].from_coefficients(0, NoContext(), values^, budget))
    return result^


def _ksturm_remainder(f: List[DensePoly[QQ]], g: List[DensePoly[QQ]], mut field: _SampleField, mut budget: Budget) raises PolyError -> List[DensePoly[QQ]]:
    """A positive multiple of -rem(f,g), without algebraic field inverses."""
    if len(g) == 0:
        raise PolyError(PolyErrorKind.DOMAIN, "zero Sturm divisor")
    if len(g) == 1:
        return List[DensePoly[QQ]]()
    var remainder = _kcopy(f, budget)
    ref lead = g[len(g) - 1]
    var lead_sign = _ksign(lead, field, budget)
    var correction = -1
    while len(remainder) >= len(g):
        var shift = len(remainder) - len(g)
        var top = remainder[len(remainder) - 1].copy()
        remainder = _kscale(remainder, lead, field, budget)
        for j in range(len(g)):
            remainder[shift + j] = remainder[shift + j].sub(_field_mul(top, g[j], field.modulus, budget), budget)
        _ktrim(remainder)
        # Each elimination multiplies the ordinary remainder by lc(g).
        # Rational primitive normalization contributes only positive factors.
        correction *= lead_sign
        remainder = _kprimitive_positive(remainder, budget)
    if correction < 0:
        for i in range(len(remainder)):
            remainder[i] = remainder[i].neg(budget)
    return remainder^


def _ksequence(f: List[DensePoly[QQ]], mut field: _SampleField, mut budget: Budget) raises PolyError -> List[List[DensePoly[QQ]]]:
    _allocation[List[DensePoly[QQ]]](1, budget)
    var chain: List[List[DensePoly[QQ]]] = [_kprimitive_positive(f, budget)]
    if len(f) > 1:
        _allocation[List[DensePoly[QQ]]](1, budget)
        chain.append(_kprimitive_positive(_kderivative(f, field, budget), budget))
        while True:
            var n = len(chain)
            var remainder = _ksturm_remainder(chain[n - 2], chain[n - 1], field, budget)
            if len(remainder) == 0:
                break
            _allocation[List[DensePoly[QQ]]](1, budget)
            chain.append(remainder^)
    return chain^


def _kvariation(chain: List[List[DensePoly[QQ]]], x: QQ, mut field: _SampleField, mut budget: Budget) raises PolyError -> Int:
    var previous = 0
    var changes = 0
    for f in chain:
        var sign = _kevaluate_sign(f, x, field, budget)
        if sign:
            changes += Int(previous != 0 and previous != sign)
            previous = sign
    return changes


def _ktotal(chain: List[List[DensePoly[QQ]]], mut field: _SampleField, mut budget: Budget) raises PolyError -> Int:
    var previous_low = 0
    var previous_high = 0
    var low = 0
    var high = 0
    for f in chain:
        var positive = _ksign(f[len(f) - 1], field, budget)
        var negative = -positive if (len(f) - 1) % 2 else positive
        low += Int(previous_low != 0 and previous_low != negative)
        high += Int(previous_high != 0 and previous_high != positive)
        previous_low = negative
        previous_high = positive
    return low - high


def _kgcd(a: List[DensePoly[QQ]], b: List[DensePoly[QQ]], field: _SampleField, mut budget: Budget) raises PolyError -> List[DensePoly[QQ]]:
    if len(a) == 0:
        return _kmonic(b, field, budget)
    if len(b) == 0:
        return _kmonic(a, field, budget)
    if len(a) == 1 or len(b) == 1:
        _allocation[DensePoly[QQ]](1, budget)
        return [_constant(QQ.one(NoContext()), budget)]
    if len(a) == 2 and len(b) != 2:
        return _kgcd(b, a, field, budget)
    if len(b) == 2:
        var linear = _kmonic(b, field, budget)
        var root = linear[0].neg(budget)
        var value = DensePoly[QQ](0, NoContext())
        for i in range(len(a) - 1, -1, -1):
            value = _kmultiply(value, root, field.modulus, budget).add(a[i], budget)
        if value.is_zero():
            return linear^
        _allocation[DensePoly[QQ]](1, budget)
        return [_constant(QQ.one(NoContext()), budget)]
    return _gcd_over_field(a, b, field.modulus, budget)


def _kbasis(inputs: List[List[DensePoly[QQ]]], field: _SampleField, mut budget: Budget) raises PolyError -> List[List[DensePoly[QQ]]]:
    var basis = List[List[DensePoly[QQ]]]()
    for f in inputs:
        if len(f) <= 1:
            continue
        var monic = _kmonic(f, field, budget)
        var duplicate = False
        var linear_basis = True
        for b in basis:
            budget.charge(1)
            duplicate = duplicate or b == monic
            linear_basis = linear_basis and len(b) == 2
        if duplicate:
            continue
        # Distinct monic linear polynomials are squarefree and coprime over
        # every field. Their coefficient comparison is already exact.
        if len(monic) == 2 and linear_basis:
            _allocation[List[DensePoly[QQ]]](1, budget)
            basis.append(monic^)
            continue
        var repeated = _kgcd(monic, _kderivative(monic, field, budget), field, budget)
        var remaining = _kexact(monic, repeated, field, budget)
        while len(remaining) > 1:
            var common = _kgcd(remaining, repeated, field, budget)
            var todo = _kexact(remaining, common, field, budget)
            repeated = _kexact(repeated, common, field, budget)
            remaining = common^
            if len(todo) <= 1:
                continue
            var refined = List[List[DensePoly[QQ]]]()
            for k in range(len(basis)):
                var b = List[DensePoly[QQ]]()
                swap(b, basis[k])
                var overlap = _kgcd(b, todo, field, budget)
                _allocation[List[DensePoly[QQ]]](2, budget)
                if len(overlap) <= 1:
                    refined.append(b^)
                else:
                    var rest = _kexact(b, overlap, field, budget)
                    todo = _kexact(todo, overlap, field, budget)
                    refined.append(overlap^)
                    if len(rest) > 1:
                        refined.append(rest^)
            if len(todo) > 1:
                _allocation[List[DensePoly[QQ]]](1, budget)
                refined.append(todo^)
            basis = refined^
    return basis^


def _kbound(f: List[DensePoly[QQ]], field: _SampleField, mut budget: Budget) raises PolyError -> QQ:
    var bound = Integer(2)
    for c in f:
        var low = QQ.zero(NoContext())
        var high = QQ.zero(NoContext())
        if field.embedding:
            var box = _dense_interval(c, field.embedding.value(), budget)
            low = box[0].copy()
            high = box[1].copy()
        elif not c.is_zero():
            low = c.coefficients[0].copy()
            high = low.copy()
        for q in [low.copy(), high.copy()]:
            budget.charge(1)
            try:
                var candidate = abs(q.value.numerator()) // q.value.denominator() + Integer(2)
                if candidate > bound:
                    bound = candidate^
            except e:
                raise PolyError(PolyErrorKind.CAPACITY, String(e))
    return QQ(Rational(bound))


def _kunion_variation(sequences: List[List[List[DensePoly[QQ]]]], x: QQ, mut field: _SampleField, mut budget: Budget) raises PolyError -> Int:
    var result = 0
    for chain in sequences:
        result = _sum(result, _kvariation(chain, x, field, budget))
    return result


def _kisolate(basis: List[List[DensePoly[QQ]]], mut field: _SampleField, mut budget: Budget) raises PolyError -> List[_FinalRoot]:
    # Pairwise coprime basis members have disjoint roots. Sum their Sturm
    # counts without forming a product whose remainder coefficients can swell.
    var sequences = List[List[List[DensePoly[QQ]]]]()
    var high = QQ(Rational(2))
    for b in basis:
        _allocation[List[List[DensePoly[QQ]]]](1, budget)
        sequences.append(_ksequence(b, field, budget))
        var candidate = _kbound(b, field, budget)
        if _compare(candidate, high, budget) > 0:
            high = candidate^
    var low = high.neg(budget)
    _allocation[_Pending](1, budget)
    var pending: List[_Pending] = [_Pending(low.copy(), high.copy(), _kunion_variation(sequences, low, field, budget), _kunion_variation(sequences, high, field, budget))]
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
                if _kvariation(sequences[k], node.low, field, budget) - _kvariation(sequences[k], node.high, field, budget) == 1:
                    owner = k
                    break
            if owner < 0:
                raise PolyError(PolyErrorKind.DOMAIN, "field root has no basis member")
            _allocation[_FinalRoot](1, budget)
            roots.append(_FinalRoot(owner, node.low.copy(), node.high.copy()))
            continue
        var middle = _midpoint(node.low, node.high, budget)
        while True:
            var boundary_root = False
            for b in basis:
                if _kevaluate_sign(b, middle, field, budget) == 0:
                    boundary_root = True
                    break
            if not boundary_root:
                break
            middle = _midpoint(node.low, middle, budget)
        var variation = _kunion_variation(sequences, middle, field, budget)
        _allocation[_Pending](2, budget)
        pending.append(_Pending(middle.copy(), node.high.copy(), variation, node.right))
        pending.append(_Pending(node.low.copy(), middle.copy(), node.left, variation))
    return roots^


def _krefine(mut root: _FinalRoot, basis: List[List[DensePoly[QQ]]], mut field: _SampleField, mut budget: Budget) raises PolyError:
    if root.low == root.high:
        return
    var middle = _midpoint(root.low, root.high, budget)
    var sign = _kevaluate_sign(basis[root.basis], middle, field, budget)
    if not sign:
        root.low = middle.copy()
        root.high = middle^
    elif sign == _kevaluate_sign(basis[root.basis], root.low, field, budget):
        root.low = middle^
    else:
        root.high = middle^


def _kgap(roots: List[_FinalRoot], gap: Int, inner: Bool, mut budget: Budget) raises PolyError -> QQ:
    var low: Optional[QQ] = None
    var high: Optional[QQ] = None
    if gap:
        low = roots[gap - 1].high.copy() if inner else roots[gap - 1].low.copy()
    if gap < len(roots):
        high = roots[gap].low.copy() if inner else roots[gap].high.copy()
    return _simplest(low, high, budget)


def _ksectors(mut roots: List[_FinalRoot], basis: List[List[DensePoly[QQ]]], mut field: _SampleField, mut budget: Budget) raises PolyError -> List[QQ]:
    _allocation[QQ](_sum(len(roots), 1), budget)
    var sectors = List[QQ]()
    for gap in range(len(roots) + 1):
        while True:
            budget.charge(1)
            var outer = _kgap(roots, gap, False, budget)
            var collapsed = False
            # Bisection need not ever hit a rational root, even zero: its
            # position within the initial interval may be a non-dyadic ratio.
            # If the simplest outer candidate is a boundary root, recognize
            # it exactly instead of refining forever towards that boundary.
            for i in range(max(0, gap - 1), min(len(roots), gap + 1)):
                if _compare(roots[i].low, outer, budget) < 0 and _compare(outer, roots[i].high, budget) < 0:
                    if _kevaluate_sign(basis[roots[i].basis], outer, field, budget) == 0:
                        roots[i].low = outer.copy()
                        roots[i].high = outer.copy()
                        collapsed = True
            if collapsed:
                continue
            var separated = gap == 0 or gap == len(roots) or _compare(roots[gap - 1].high, roots[gap].low, budget) < 0
            if separated:
                var candidate = _kgap(roots, gap, True, budget)
                if candidate == outer:
                    sectors.append(candidate^)
                    break
            if gap:
                _krefine(roots[gap - 1], basis, field, budget)
            if gap < len(roots):
                _krefine(roots[gap], basis, field, budget)
    return sectors^


def _construct_field_stack(inputs: List[List[DensePoly[QQ]]], mut field: _SampleField, mut budget: Budget, *, interval: Bool = False, interval_cap: UInt64 = UInt64(2000000)) raises PolyError -> _FieldStack:
    var basis = _kbasis(inputs, field, budget)
    var leading = List[Int8]()
    var multiplicities = List[List[Int]]()
    for f in inputs:
        _allocation[Int8](1, budget)
        var sign = _ksign(f[len(f) - 1], field, budget) if len(f) else 0
        leading.append(Int8(sign))
        _allocation[Int](len(basis), budget)
        var row = List[Int](length=len(basis), fill=0)
        var rest = _kcopy(f, budget)
        for k in range(len(basis)):
            while len(rest) >= len(basis[k]):
                var division = _kdivmod(rest, basis[k], field, budget)
                if len(division[1]):
                    break
                swap(rest, division[0])
                row[k] = _sum(row[k], 1)
        _allocation[List[Int]](1, budget)
        multiplicities.append(row^)
    var isolation = _IntervalRootStats()
    var proposed: Optional[List[_FinalRoot]] = None
    if interval:
        proposed = _interval_field_roots(basis, field, isolation, budget, cap=interval_cap)
    var roots: List[_FinalRoot]
    if proposed:
        roots = proposed.take()
    else:
        isolation.exact_fallbacks = Int(interval)
        roots = _kisolate(basis, field, budget)
    var sectors = _ksectors(roots, basis, field, budget)
    _allocation[List[Int8]](_sum(_product(2, len(roots)), 1), budget)
    var signs = List[List[Int8]]()
    _allocation[Int8](len(leading), budget)
    var current = leading.copy()
    for j in range(len(inputs)):
        if len(inputs[j]) and len(inputs[j]) % 2 == 0:
            current[j] = -current[j]
    _allocation[Int8](len(leading), budget)
    signs.append(current.copy())
    # Build from -infinity; the final verifier independently walks from
    # +infinity after checking exact polynomial product identities.
    for root in roots:
        _allocation[Int8](_product(2, len(leading)), budget)
        var section = current.copy()
        for j in range(len(inputs)):
            budget.charge(1)
            var exponent = multiplicities[j][root.basis]
            if exponent:
                section[j] = 0
            if exponent & 1:
                current[j] = -current[j]
        signs.append(section^)
        signs.append(current.copy())
    return _FieldStack(basis^, leading^, multiplicities^, roots^, sectors^, signs^, isolation^)


def _check_kpoly(f: List[DensePoly[QQ]], field: _SampleField, mut budget: Budget) raises PolyError -> Bool:
    if len(f) and f[len(f) - 1].is_zero():
        return False
    for c in f:
        c.check_invariants(budget)
        if c.variable != 0 or c.degree() >= field.modulus.degree():
            return False
    return True


def _check_field_root_counts(stack: _FieldStack, mut field: _SampleField, mut budget: Budget) raises PolyError -> Bool:
    # Algebra and globally disjoint geometry were checked first. Completeness
    # for every basis member puts every foreign root in its own disjoint
    # bracket, so only the owned endpoints need evaluation for each member.
    for k in range(len(stack.basis)):
        var linear = len(stack.basis[k]) == 2
        var chain = List[List[DensePoly[QQ]]]()
        if not linear:
            chain = _ksequence(stack.basis[k], field, budget)
        var found = 0
        for root in stack.roots:
            if root.basis != k:
                continue
            var low = _kevaluate_sign(stack.basis[k], root.low, field, budget)
            if root.low == root.high:
                if low != 0:
                    return False
            else:
                var high = _kevaluate_sign(stack.basis[k], root.high, field, budget)
                if low == 0 or high == 0:
                    return False
                if linear:
                    if low >= 0 or high <= 0:
                        return False
                elif _kvariation(chain, root.low, field, budget) - _kvariation(chain, root.high, field, budget) != 1:
                    return False
            found += 1
        if found != (1 if linear else _ktotal(chain, field, budget)):
            return False
    return True


def _check_stack_in_field_cached(inputs: List[List[DensePoly[QQ]]], stack: _FieldStack, mut field: _SampleField, mut root_counts: _RootCountCache, mut budget: Budget, *, vectorized: Bool = True) raises PolyError -> Bool:
    """Internal checker; its caller must validate this shared field first."""
    var count = len(stack.roots)
    if len(stack.leading) != len(inputs) or len(stack.multiplicities) != len(inputs) or len(stack.sectors) != _sum(count, 1) or len(stack.signs) != _sum(_product(2, count), 1):
        return False
    _allocation[Bool](len(stack.basis), budget)
    var used = List[Bool](length=len(stack.basis), fill=False)
    var one = _constant(QQ.one(NoContext()), budget)
    for k in range(len(stack.basis)):
        ref b = stack.basis[k]
        if not _check_kpoly(b, field, budget) or len(b) <= 1 or b[len(b) - 1] != one:
            return False
        if len(b) > 2 and len(_kgcd(b, _kderivative(b, field, budget), field, budget)) != 1:
            return False
        for j in range(k):
            if len(b) == 2 and len(stack.basis[j]) == 2:
                if b[0] == stack.basis[j][0]:
                    return False
                continue
            if len(_kgcd(b, stack.basis[j], field, budget)) != 1:
                return False
    for j in range(len(inputs)):
        ref f = inputs[j]
        if not _check_kpoly(f, field, budget) or len(stack.multiplicities[j]) != len(stack.basis):
            return False
        var lead = f[len(f) - 1].copy() if len(f) else DensePoly[QQ](0, NoContext())
        # Recover the leading coefficient from the bound inputs for the
        # exact product identity; the final data only needs its sign.
        if _ksign(lead, field, budget) != Int(stack.leading[j]):
            return False
        var remaining = max(0, len(f) - 1)
        var product: List[DensePoly[QQ]] = [lead^] if len(f) else List[DensePoly[QQ]]()
        for k in range(len(stack.basis)):
            budget.charge(1)
            var exponent = stack.multiplicities[j][k]
            var degree = len(stack.basis[k]) - 1
            if exponent < 0 or exponent > remaining // degree:
                return False
            remaining -= exponent * degree
            if exponent:
                used[k] = True
            for _ in range(exponent):
                product = _kmul(product, stack.basis[k], field, budget)
        if remaining or product != f:
            return False
    for present in used:
        if not present:
            return False
    for row in stack.signs:
        if len(row) != len(inputs):
            return False
    for i in range(count):
        ref root = stack.roots[i]
        if root.basis < 0 or root.basis >= len(stack.basis) or _compare(root.low, root.high, budget) > 0:
            return False
        if i and _compare(stack.roots[i - 1].high, root.low, budget) >= 0:
            return False
        if _compare(stack.sectors[i], root.low, budget) >= 0 or _compare(stack.sectors[i + 1], root.high, budget) <= 0:
            return False
    for gap in range(count + 1):
        if stack.sectors[gap] != _kgap(stack.roots, gap, False, budget):
            return False
    if not field.embedding:
        if not _check_rational_roots(stack.basis, stack.roots, root_counts, budget, vectorized):
            return False
    elif not _check_field_root_counts(stack, field, budget):
        return False
    # Independent top-down sign recurrence, not the producer's label builder.
    for j in range(len(inputs)):
        var expected = Int(stack.leading[j])
        if Int(stack.signs[2 * count][j]) != expected:
            return False
        for k in range(count - 1, -1, -1):
            var exponent = stack.multiplicities[j][stack.roots[k].basis]
            if Int(stack.signs[2 * k + 1][j]) != (0 if exponent else expected):
                return False
            if exponent & 1:
                expected = -expected
            if Int(stack.signs[2 * k][j]) != expected:
                return False
    return True


def _check_stack_in_field(inputs: List[List[DensePoly[QQ]]], stack: _FieldStack, mut field: _SampleField, mut budget: Budget, *, vectorized: Bool = True) raises PolyError -> Bool:
    var root_counts = _RootCountCache()
    return _check_stack_in_field_cached(inputs, stack, field, root_counts, budget, vectorized=vectorized) and root_counts.flush(budget, vectorized)


def _check_field_stack(inputs: List[List[DensePoly[QQ]]], stack: _FieldStack, field: _SampleField, mut budget: Budget, *, vectorized: Bool = True) raises PolyError -> Bool:
    try:
        if not _check_sample_field(field, budget):
            return False
        var checked_field = field.copy()
        return _check_stack_in_field(inputs, stack, checked_field, budget, vectorized=vectorized)
    except e:
        if e.kind == PolyErrorKind.DOMAIN:
            return False
        raise e^
