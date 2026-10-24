"""Exact, pseudo, and ordered multivariate polynomial division.

The loops implement Algorithms V-1/V-4 and the scaling proof of V-2.
"""
from .coefficient_buffer import _view
from .errors import Budget, PolyError, PolyErrorKind, _allocation, _charge_size
from .coefficients import CoefficientRing, CoefficientField, _pow
from ._coefficient_ops import _quotient, _exact
from .dense import DensePoly, _slice, _trim, _word_residues
from .sparse import SparsePoly, Term
from .monomial import Monomial, MonomialOrder
from .certificates import DivisionResult, ListDivision
from ._packed_division import _packed_divide as _heap_divide
from ._sparse_ops import _subtract_multiple, _from_descending
from ._reduction_stream import _ReductionStream
from ._packed_products import _packed_divide
from ._ntt import _ntt_divide_exact
from ._nmod import _Field, _field, _mulm, _from
from ._field_division import _field_divmod, _zz_exact

comptime _WORD_DIVISION = 256
"""Divisor length below which exact GF(p) division runs on word residues instead of transforms."""
from .coefficients import ZZ, ModInt, _submul


def _dense_divide[C: CoefficientRing](f: DensePoly[C], g: DensePoly[C], pseudo: Bool, mut budget: Budget, quotient: Bool = True) raises PolyError -> DivisionResult[C]:
    f._same(g)
    if g.is_zero():
        raise PolyError(PolyErrorKind.DOMAIN, "division by zero polynomial")
    var count = max(0, f.degree() - g.degree() + 1)
    var quotient_count = count if quotient else 0
    _allocation[C](quotient_count, budget)
    var q = List[C](length=quotient_count, fill=C.zero(f.ring))
    var r = _slice(_view(f.coefficients), 0, len(f.coefficients), budget)
    var remaining = count
    var inverse = C.one(f.ring)
    ref lc = g.coefficients[len(g.coefficients) - 1]
    if not pseudo and count:
        inverse = _exact(inverse, lc, budget)
    while len(r) >= len(g.coefficients):
        budget.charge(1)
        var shift = len(r) - len(g.coefficients)
        _charge_size(r[len(r) - 1].size(), budget)
        var factor = r[len(r) - 1].copy()
        if pseudo:
            for i in range(len(q)):
                budget.charge(1)
                q[i] = q[i].mul(lc, budget)
            for i in range(len(r)):
                budget.charge(1)
                r[i] = r[i].mul(lc, budget)
        else:
            factor = factor.mul(inverse, budget)
        if quotient:
            q[shift] = q[shift].add(factor, budget)
        for i in range(len(g.coefficients)):
            budget.charge(1)
            comptime if C == ZZ:
                _submul(rebind[ZZ](r[shift + i]), rebind[ZZ](factor), rebind[ZZ](g.coefficients[i]), budget)
            else:
                r[shift + i] = r[shift + i].sub(factor.mul(g.coefficients[i], budget), budget)
        _trim(r, budget)
        remaining -= 1
    if pseudo and remaining:
        var multiplier = _pow(lc, UInt64(remaining), budget)
        for i in range(len(q)):
            budget.charge(1)
            q[i] = q[i].mul(multiplier, budget)
        for i in range(len(r)):
            budget.charge(1)
            r[i] = r[i].mul(multiplier, budget)
    _trim(q, budget)
    return DivisionResult[C](DensePoly[C]._from_canonical(f.variable, f.ring, q^), DensePoly[C]._from_canonical(f.variable, f.ring, r^))


def _words(p: DensePoly[ModInt], f: _Field) -> List[UInt64]:
    # Residues read run by run, then moved to Montgomery form in place.
    var out = _word_residues(_view(p.coefficients))
    for i in range(len(out)):
        out.unsafe_get(i) = _mulm(out.unsafe_get(i), f.r2, f)
    return out^


def _from_words(words: List[UInt64], template: ModInt, f: _Field, mut budget: Budget) raises PolyError -> List[ModInt]:
    _allocation[ModInt](len(words), budget)
    var out = List[ModInt](capacity=len(words))
    for w in words:
        out.append(template._with_residue(_from(w, f)))
    return out^


def _word_divmod(f: DensePoly[ModInt], g: DensePoly[ModInt], mut budget: Budget) raises PolyError -> Optional[DivisionResult[ModInt]]:
    """divmod over an odd prime field on Montgomery word residues, or None for p = 2."""
    if g.is_zero() or f.is_zero():
        return None
    ref lead = g.coefficients[len(g.coefficients) - 1]
    if lead.modulus % 2 == 0 or lead.modulus < 3 or lead.modulus >= UInt64(1) << 62:
        return None
    var field = _field(lead.modulus)
    _allocation[UInt64](len(f.coefficients) + len(g.coefficients), budget)
    var qr = _field_divmod(_words(f, field), _words(g, field), field, budget)
    var q = DensePoly[ModInt]._from_canonical(f.variable, f.ring, _from_words(qr[0], lead, field, budget))
    var r = DensePoly[ModInt]._from_canonical(f.variable, f.ring, _from_words(qr[1], lead, field, budget))
    return DivisionResult[ModInt](q^, r^)


def divmod[C: CoefficientField](f: DensePoly[C], g: DensePoly[C], mut budget: Budget) raises PolyError -> DivisionResult[C]:
    """Divide univariate polynomials over a field.

    Parameters:
        C: The coefficient field.

    Args:
        f: The dividend.
        g: The nonzero divisor in the same variable and ring.
        budget: The shared work allowance.

    Returns:
        The unique quotient and remainder of degree less than the divisor.

    Raises:
        PolyError: DOMAIN for invalid operands, LIMIT or CAPACITY on work limits.
    """
    comptime if C == ModInt:
        f._same(g)
        var fast = _word_divmod(rebind[DensePoly[ModInt]](f), rebind[DensePoly[ModInt]](g), budget)
        if fast:
            return rebind_var[DivisionResult[C]](fast.take())
    return _dense_divide(f, g, False, budget)


def pseudo_divmod[C: CoefficientRing](f: DensePoly[C], g: DensePoly[C], mut budget: Budget) raises PolyError -> DivisionResult[C]:
    """Divide after multiplying the dividend by the standard leading-coefficient power.

    Parameters:
        C: The coefficient integral domain.

    Args:
        f: The dividend.
        g: The nonzero divisor in the same variable and ring.
        budget: The shared work allowance.

    Returns:
        Quotient and remainder for lc(g)**max(deg(f)-deg(g)+1, 0) times f.

    Raises:
        PolyError: DOMAIN for invalid operands, LIMIT or CAPACITY on work limits.
    """
    return _dense_divide(f, g, True, budget)


def _division_order[C: CoefficientRing](p: SparsePoly[C], others: List[SparsePoly[C]], mut budget: Budget) raises PolyError -> MonomialOrder:
    if p.order.kind != MonomialOrder.KERNEL:
        _charge_size(len(p.order.variables), budget)
        return p.order.copy()
    # KERNEL is a storage order, not a well-order. Use a deterministic lex
    # work order and restore storage order on return.
    var variables = p.variables(budget)
    for other in others:
        budget.charge(1)
        var additions = other.variables(budget)
        for index in additions:
            budget.charge(1)
            var found = False
            for existing in variables:
                budget.charge(1)
                if existing == index:
                    found = True
                    break
            if not found:
                variables.append(index)
    return MonomialOrder(MonomialOrder.LEX, variables^, budget)


def _divide_by_term[C: CoefficientRing](f: SparsePoly[C], g: SparsePoly[C], mut budget: Budget) raises PolyError -> Optional[SparsePoly[C]]:
    """The exact quotient by a one-term divisor, term by term, or None.

    Dividing every monomial by one monomial keeps them distinct and, monomial
    orders being multiplicative, in order; exact quotients of nonzero
    coefficients are nonzero, so the quotient terms are already canonical
    except in the kernel's storage order, which is not multiplicative.
    """
    ref lead = g.terms[0]
    _allocation[Term[C]](len(f.terms), budget)
    var terms = List[Term[C]](capacity=len(f.terms))
    comptime if C == ModInt:
        # One inverse for every coefficient.
        var inverse = rebind[ModInt](lead.coefficient).inverse(budget)
        for t in f.terms:
            budget.charge(1)
            if not lead.monomial.divides(t.monomial, budget):
                return None
            terms.append(Term[C](t.monomial.div(lead.monomial, budget), rebind_var[C](rebind[ModInt](t.coefficient).mul(inverse, budget))))
    else:
        for t in f.terms:
            budget.charge(1)
            if not lead.monomial.divides(t.monomial, budget):
                return None
            var c = _quotient(t.coefficient, lead.coefficient, budget)
            if not c:
                return None
            terms.append(Term[C](t.monomial.div(lead.monomial, budget), c.take()))
    if f.order.kind == MonomialOrder.KERNEL:
        # Storage order compares pair lists, which dividing can reorder.
        return SparsePoly[C].from_terms(terms^, f.order, f.ring, budget)
    var result = SparsePoly[C].zero(f.order, f.ring, budget)
    result.terms = terms^
    return result^


def divide_exact[C: CoefficientRing](f: SparsePoly[C], g: SparsePoly[C], mut budget: Budget) raises PolyError -> Optional[SparsePoly[C]]:
    """Return an exact quotient over the coefficient ring, or no value.

    Parameters:
        C: ZZ or a coefficient field.

    Args:
        f: The dividend.
        g: The nonzero divisor in the same ring and order.
        budget: The shared work allowance.

    Returns:
        The quotient if g divides f, otherwise None.

    Raises:
        PolyError: DOMAIN for invalid operands, UNSUPPORTED for another ring, LIMIT or CAPACITY on work limits.
    """
    f._same(g, budget)
    if g.is_zero():
        raise PolyError(PolyErrorKind.DOMAIN, "division by zero polynomial")
    if len(g.terms) == 1:
        return _divide_by_term(f, g, budget)
    comptime if C == ZZ or C == ModInt:
        if not f.is_zero() and len(g.terms) > 1:
            var fast = _fast_divide_exact(f, g, budget)
            if fast:
                return fast.take()
    var applicable = False
    var packed = _heap_divide(f, g, applicable, budget)
    if applicable:
        return packed^
    var others: List[SparsePoly[C]] = [g._copy(budget)]
    var order = _division_order(f, others, budget)
    var pending = f.reorder(order, budget)
    var divisor = g.reorder(order, budget)
    var terms = List[Term[C]]()
    var leading = divisor.leading(budget)
    while not pending.is_zero():
        budget.charge(1)
        var last = pending.leading(budget)
        if not leading.monomial.divides(last.monomial, budget):
            return None
        var coefficient = _quotient(last.coefficient, leading.coefficient, budget)
        if not coefficient:
            return None
        var term = Term[C](last.monomial.div(leading.monomial, budget), coefficient.take())
        pending = _subtract_multiple(pending, divisor, term, budget)
        _allocation[Term[C]](1, budget)
        terms.append(term^)
    return _from_descending(terms^, divisor, budget).reorder(f.order, budget)


def divide_exact_dense[C: CoefficientRing](f: DensePoly[C], g: DensePoly[C], mut budget: Budget) raises PolyError -> Optional[DensePoly[C]]:
    """Return the exact dense quotient over the coefficient ring, or no value; the dense counterpart of divide_exact.

    Parameters:
        C: ZZ or a coefficient field.

    Args:
        f: The dividend.
        g: The nonzero divisor in the same variable and ring.
        budget: The shared work allowance.

    Returns:
        The quotient if g divides f, otherwise None.

    Raises:
        PolyError: DOMAIN for invalid operands, UNSUPPORTED for another ring, LIMIT or CAPACITY on work limits.
    """
    f._same(g)
    if g.is_zero():
        raise PolyError(PolyErrorKind.DOMAIN, "division by zero polynomial")
    if f.is_zero():
        return DensePoly[C](f.variable, f.ring.copy())
    if len(f.coefficients) < len(g.coefficients):
        return None
    comptime if C == ModInt:
        var fast = _word_divmod(rebind[DensePoly[ModInt]](f), rebind[DensePoly[ModInt]](g), budget)
        if fast:
            if not fast.value().remainder.is_zero():
                return None
            return rebind[DensePoly[C]](fast.value().quotient).copy()
    comptime if C == ZZ:
        if len(g.coefficients) >= 32:
            var q = _ntt_divide_exact(_view(f.coefficients), _view(g.coefficients), budget)
            if q:
                var quotient = q.take()
                if not quotient:
                    return None
                return DensePoly[C]._from_canonical(f.variable, f.ring, quotient.take())
    var order = MonomialOrder(MonomialOrder.LEX, [f.variable], budget)
    var sparse = divide_exact(f.to_sparse(order, budget), g.to_sparse(order, budget), budget)
    if not sparse:
        return None
    return DensePoly[C].from_sparse(sparse.value(), budget, variable=Optional(f.variable))


def _univariate[C: CoefficientRing](p: SparsePoly[C]) -> Optional[Int64]:
    """The single variable of a nonconstant univariate polynomial."""
    var variable = Optional[Int64](None)
    for term in p.terms:
        if len(term.monomial.pairs) > 1:
            return None
        if len(term.monomial.pairs):
            var v = term.monomial.pairs[0][0]
            if variable and variable.value() != v:
                return None
            variable = v
    return variable


def _fast_divide_exact[C: CoefficientRing](f: SparsePoly[C], g: SparsePoly[C], mut budget: Budget) raises PolyError -> Optional[Optional[SparsePoly[C]]]:
    """Exact ZZ or GF(p) division by dense transforms (one variable) or packed chunks (LEX), or None when neither applies."""
    var x = _univariate(f)
    var y = _univariate(g)
    if x and y and x.value() == y.value():
        var a = DensePoly[C].from_sparse(f, budget)
        var b = DensePoly[C].from_sparse(g, budget)
        comptime if C == ModInt:
            # Short divisors divide faster on Montgomery word residues than by transforms.
            if len(b.coefficients) < _WORD_DIVISION:
                var fast = _word_divmod(rebind[DensePoly[ModInt]](a), rebind[DensePoly[ModInt]](b), budget)
                if fast:
                    if not fast.value().remainder.is_zero():
                        return Optional(Optional[SparsePoly[C]](None))
                    return Optional(Optional(rebind[DensePoly[C]](fast.value().quotient).to_sparse(f.order, budget)))
        comptime if C == ZZ:
            if len(a.coefficients) >= len(b.coefficients):
                var small = _zz_exact(_view(a.coefficients), _view(b.coefficients), budget)
                if small:
                    var quotient = small.take()
                    if not quotient:
                        return Optional(Optional[SparsePoly[C]](None))
                    var dense = DensePoly[C]._from_canonical(a.variable, a.ring, quotient.take())
                    return Optional(Optional(dense.to_sparse(f.order, budget)))
        if len(b.coefficients) >= 32 and len(a.coefficients) >= len(b.coefficients):
            var q = _ntt_divide_exact(_view(a.coefficients), _view(b.coefficients), budget)
            if q:
                var quotient = q.take()
                if not quotient:
                    return Optional(Optional[SparsePoly[C]](None))
                var dense = DensePoly[C]._from_canonical(a.variable, a.ring, quotient.take())
                return Optional(Optional(dense.to_sparse(f.order, budget)))
    if f.order.kind == MonomialOrder.LEX and len(f.terms) >= 256:
        return _packed_divide(f, g, budget)
    return None


def divide_by_list[C: CoefficientField](f: SparsePoly[C], divisors: List[SparsePoly[C]], mut budget: Budget) raises PolyError -> ListDivision[C]:
    """Reduce by nonzero divisors in the supplied order under a monomial order.

    Parameters:
        C: The coefficient field.

    Args:
        f: The dividend in LEX, GRLEX, or GREVLEX order.
        divisors: Nonzero divisors in the same ring and order.
        budget: The shared work allowance.

    Returns:
        Quotients and a remainder with no term divisible by a divisor's leading monomial.

    Raises:
        PolyError: DOMAIN for zero divisors, mixed rings/orders, or KERNEL order; LIMIT or CAPACITY on work limits.
    """
    _allocation[Int](len(divisors), budget)
    var indices = List[Int](capacity=len(divisors))
    for i in range(len(divisors)):
        budget.charge(1)
        indices.append(i)
    return _divide_by_indices(f, divisors, indices, True, budget)


struct _ReducerData[C: CoefficientField](Movable):
    var leads: List[Monomial]
    var inverses: List[Self.C]

    def __init__(out self):
        self.leads = List[Monomial]()
        self.inverses = List[Self.C]()

    def append(mut self, g: SparsePoly[Self.C], mut budget: Budget) raises PolyError:
        var leading = g.leading(budget)
        var inverse = leading.coefficient.inverse(budget)
        _allocation[Monomial](1, budget)
        _allocation[Self.C](1, budget)
        var monomial = Monomial()
        swap(monomial, leading.monomial)
        self.leads.append(monomial^)
        self.inverses.append(inverse^)


def _divide_by_indices[C: CoefficientField](f: SparsePoly[C], divisors: List[SparsePoly[C]], indices: List[Int], carry: Bool, mut budget: Budget) raises PolyError -> ListDivision[C]:
    if f.order.kind == MonomialOrder.KERNEL:
        raise PolyError(PolyErrorKind.DOMAIN, "list division requires a monomial order")
    var data = _ReducerData[C]()
    for g in divisors:
        budget.charge(1)
        f._same(g, budget)
        data.append(g, budget)
    return _divide_prepared(f, divisors, indices, carry, data, budget)


def _divide_prepared[C: CoefficientField](f: SparsePoly[C], divisors: List[SparsePoly[C]], indices: List[Int], carry: Bool, data: _ReducerData[C], mut budget: Budget) raises PolyError -> ListDivision[C]:
    # The public wrapper validates ring/order compatibility. The basis producer
    # keeps one metadata entry per immutable arena polynomial and reuses it.
    _allocation[List[Term[C]]](len(indices) if carry else 0, budget)
    var quotient_terms = List[List[Term[C]]]()
    if carry:
        for _ in indices:
            budget.charge(1)
            quotient_terms.append(List[Term[C]]())
    var remainder = List[Term[C]]()
    var pending = f._copy(budget)
    var steps = 0
    while not pending.is_zero():
        if steps == 2 and len(pending.terms) >= 64:
            break
        budget.charge(1)
        var last = pending.leading(budget)
        var reduced = False
        for i in range(len(indices)):
            budget.charge(1)
            if not data.leads[indices[i]].divides(last.monomial, budget):
                continue
            var t = Term[C](last.monomial.div(data.leads[indices[i]], budget), last.coefficient.mul(data.inverses[indices[i]], budget))
            pending = _subtract_multiple(pending, divisors[indices[i]], t, budget)
            if carry:
                _allocation[Term[C]](1, budget)
                quotient_terms[i].append(t^)
            steps = min(2, steps + 1)
            reduced = True
            break
        if not reduced:
            _ = pending.terms.pop()
            _allocation[Term[C]](1, budget)
            remainder.append(last^)
    if not pending.is_zero():
        var stream = _ReductionStream[C](pending^, budget)
        while True:
            var next = stream.next(divisors, budget)
            if not next:
                break
            var last = next.take()
            var reduced = False
            for i in range(len(indices)):
                budget.charge(1)
                if not data.leads[indices[i]].divides(last.monomial, budget):
                    continue
                var t = Term[C](last.monomial.div(data.leads[indices[i]], budget), last.coefficient.mul(data.inverses[indices[i]], budget))
                stream.append(indices[i], t, divisors, budget)
                if carry:
                    _allocation[Term[C]](1, budget)
                    quotient_terms[i].append(t^)
                reduced = True
                break
            if not reduced:
                _allocation[Term[C]](1, budget)
                remainder.append(last^)
    _allocation[SparsePoly[C]](len(quotient_terms), budget)
    var quotients = List[SparsePoly[C]](capacity=len(quotient_terms))
    for i in range(len(quotient_terms)):
        budget.charge(1)
        var terms = List[Term[C]]()
        swap(terms, quotient_terms[i])
        quotients.append(_from_descending(terms^, f, budget))
    return ListDivision[C](quotients^, _from_descending(remainder^, f, budget))
