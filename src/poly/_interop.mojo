"""Checked copies and explicit variable maps at the mutable core boundary."""
from std.builtin.rebind import downcast
from std.builtin.variadics import TypeList
from polynomials import CoefficientRing, DensePoly, SparsePoly, Monomial, MonomialOrder, Term, Budget, PolyError, PolyErrorKind, Modulus
from ._common import _CoreKind, _Ring, _allocation, _charge_size, _product, _Context, _convert, _unwrap, _wrap

comptime _Shape[n: Int] = downcast[Int if n == 1 else Tuple[*TypeList.splat[n, Int]()], Copyable & Deinitable]


def _take_dense[C: CoefficientRing](mut source: DensePoly[C]) -> DensePoly[C]:
    var result = DensePoly[C](source.variable, source.ring.copy())
    swap(result, source)
    return result^


def _take_sparse[C: CoefficientRing](mut source: SparsePoly[C]) -> SparsePoly[C]:
    var result = SparsePoly[C](MonomialOrder(), source.ring.copy())
    swap(result, source)
    return result^


def _context[C: CoefficientRing](ctx: C.Context, mut budget: Budget) raises PolyError -> C.Context:
    comptime k = _CoreKind[C]
    comptime assert k >= 0, "array interoperability supports ZZ, QQ, QQi, and ModInt"
    var context = _Context[k]()
    comptime if k == 3:
        context.observe(rebind[Modulus](ctx).value, budget)
    return rebind_var[C.Context](context.get())


def _dense_copy[C: CoefficientRing](source: DensePoly[C], variable: Int64, mut budget: Budget) raises PolyError -> DensePoly[C]:
    comptime k = _CoreKind[C]
    source.check_invariants(budget)
    var ctx = _context[C](source.ring, budget)
    _allocation[C](len(source.coefficients), budget)
    var values = List[C](capacity=len(source.coefficients))
    for coefficient_index in range(len(source.coefficients)):
        ref c = source.coefficients[coefficient_index]
        _charge_size(c.size(), budget)
        values.append(rebind_var[C](_wrap[k](_convert[k](_unwrap[k](rebind[_Ring[k]](c)), rebind[_Ring[k].Context](ctx), budget))))
    return DensePoly[C].from_coefficients(variable, ctx, values^, budget)


def _variables(width: Int, requested: List[Int64], mut budget: Budget) raises PolyError -> List[Int64]:
    if len(requested) and len(requested) != width:
        raise PolyError(PolyErrorKind.DOMAIN, "variable map must have one core ID per array axis")
    _allocation[Int64](width, budget)
    var result = List[Int64](capacity=width)
    for i in range(width):
        budget.charge(1)
        var value = requested[i] if len(requested) else Int64(i)
        for previous in result:
            budget.charge(1)
            if previous == value:
                raise PolyError(PolyErrorKind.DOMAIN, "variable map contains duplicate IDs")
        result.append(value)
    return result^


def _mapped_order(source: MonomialOrder, mapping: List[Int64], into_core: Bool, mut budget: Budget) raises PolyError -> MonomialOrder:
    # Core fields are mutable, including on a zero polynomial. Validate the
    # stored order itself, not only monomials present in its term list.
    _allocation[Int64](len(source.variables), budget)
    _ = MonomialOrder(source.kind, source.variables.copy(), budget)
    _allocation[Int64](len(mapping), budget)
    var variables = List[Int64](capacity=len(mapping))
    for i in range(len(mapping)):
        budget.charge(1)
        variables.append(mapping[i] if into_core else Int64(i))
    return MonomialOrder(MonomialOrder.LEX, variables^, budget)


def _remap_pairs(monomial: Monomial, mapping: List[Int64], into_core: Bool, mut budget: Budget) raises PolyError -> List[Tuple[Int64, UInt32]]:
    _allocation[Tuple[Int64, UInt32]](len(monomial.pairs), budget)
    var pairs = List[Tuple[Int64, UInt32]](capacity=len(monomial.pairs))
    for entry in monomial.pairs:
        budget.charge(1)
        var target = Int64(0)
        if into_core:
            if entry[0] < 0 or entry[0] >= Int64(len(mapping)):
                raise PolyError(PolyErrorKind.DOMAIN, "array axis is outside its power shape")
            target = mapping[Int(entry[0])]
        else:
            var found = False
            for i in range(len(mapping)):
                budget.charge(1)
                if entry[0] == mapping[i]:
                    target = Int64(i)
                    found = True
                    break
            if not found:
                raise PolyError(PolyErrorKind.DOMAIN, "core polynomial uses a variable absent from the map")
        pairs.append((target, entry[1]))
        var j = len(pairs) - 1
        while j > 0:
            budget.charge(1)
            if pairs[j - 1][0] < pairs[j][0]:
                break
            pairs.swap_elements(j, j - 1)
            j -= 1
    return pairs^


def _sparse_copy[C: CoefficientRing](source: SparsePoly[C], mapping: List[Int64], into_core: Bool, mut budget: Budget) raises PolyError -> SparsePoly[C]:
    comptime k = _CoreKind[C]
    source.check_invariants(budget)
    var ctx = _context[C](source.ring, budget)
    var order = _mapped_order(source.order, mapping, into_core, budget)
    _allocation[Term[C]](_product(len(source.terms), 7), budget)
    var terms = List[Term[C]](capacity=len(source.terms))
    for term in source.terms:
        budget.charge(1)
        var pairs = _remap_pairs(term.monomial, mapping, into_core, budget)
        _charge_size(term.coefficient.size(), budget)
        var coefficient = rebind_var[C](_wrap[k](_convert[k](_unwrap[k](rebind[_Ring[k]](term.coefficient)), rebind[_Ring[k].Context](ctx), budget)))
        terms.append(Term[C](Monomial(pairs^, budget), coefficient^))
    return SparsePoly[C].from_terms(terms^, order, ctx, budget)


def _dense_sparse[C: CoefficientRing](p: DensePoly[C], mut budget: Budget) raises PolyError -> SparsePoly[C]:
    _allocation[Term[C]](_product(len(p.coefficients), 4), budget)
    _allocation[Tuple[Int64, UInt32]](len(p.coefficients), budget)
    _allocation[Int64](1, budget)
    var order = MonomialOrder(MonomialOrder.LEX, [p.variable], budget)
    return p.to_sparse(order, budget)
