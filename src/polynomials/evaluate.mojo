"""Budgeted Horner evaluation, polynomial substitution, and formal derivatives."""
from .errors import Budget, PolyError, PolyErrorKind, _charge_size, _allocation
from .coefficients import CoefficientRing, _pow
from .monomial import Monomial, MonomialOrder
from .sparse import Term, SparsePoly
from .dense import DensePoly
from ._grouping import _coefficient_groups, _CoefficientGroups

def _natural[C: CoefficientRing](n: UInt64, ctx: C.Context, mut budget: Budget) raises PolyError -> C:
    var result = C.zero(ctx)
    var bit = C.one(ctx)
    var value = n
    while value:
        budget.charge(1)
        if value & 1:
            result = result.add(bit, budget)
        value >>= 1
        if value:
            bit = bit.add(bit, budget)
    return result^

@fieldwise_init
struct _HornerTask(Movable):
    var lo: Int
    var hi: Int
    var level: Int
    var exponents: List[UInt32]

def evaluate[C: CoefficientRing](p: SparsePoly[C], point: List[Tuple[Int64, C]], mut budget: Budget) raises PolyError -> C:
    """Evaluate by nested sparse Horner schemes using an explicit work stack.

    Parameters:
        C: The coefficient ring type.

    Args:
        p: The polynomial.
        point: Distinct variable indices and values in the same ring.
        budget: The shared work allowance.

    Returns:
        The exact value in the coefficient ring.

    Raises:
        PolyError: DOMAIN for missing or repeated variables or mixed rings, LIMIT or CAPACITY on work limits.
    """
    var variables = List[Int64]()
    for entry in point:
        budget.charge(1)
        if entry[1].context() != p.ring:
            raise PolyError(PolyErrorKind.DOMAIN, "evaluation point in another ring")
        variables.append(entry[0])
    var order = MonomialOrder(MonomialOrder.LEX, variables^, budget)
    var sorted_poly = p.reorder(order, budget)
    var tasks = List[_HornerTask]()
    var values = List[C]()
    tasks.append(_HornerTask(0, len(sorted_poly.terms), 0, List[UInt32]()))
    while len(tasks):
        budget.charge(1)
        var task = tasks.pop()
        if len(task.exponents):
            var result = values.pop()
            var last = task.exponents[len(task.exponents) - 1]
            for i in range(len(task.exponents) - 1):
                budget.charge(1)
                var e = task.exponents[len(task.exponents) - 2 - i]
                result = result.mul(_pow(point[task.level][1], UInt64(last - e), budget), budget)
                var next_value = values.pop()
                result = result.add(next_value, budget)
                last = e
            if last:
                result = result.mul(_pow(point[task.level][1], UInt64(last), budget), budget)
            values.append(result^)
        elif task.lo == task.hi:
            values.append(C.zero(p.ring))
        elif task.level == len(point):
            _charge_size(sorted_poly.terms[task.lo].coefficient.size(), budget)
            values.append(sorted_poly.terms[task.lo].coefficient.copy())
        else:
            var groups = List[Tuple[Int, Int]]()
            var exponents = List[UInt32]()
            var start = task.lo
            while start < task.hi:
                budget.charge(1)
                var e = sorted_poly.terms[start].monomial.exponent(point[task.level][0], budget)
                var end = start + 1
                while end < task.hi:
                    budget.charge(1)
                    if sorted_poly.terms[end].monomial.exponent(point[task.level][0], budget) != e:
                        break
                    end += 1
                groups.append((start, end))
                exponents.append(e)
                start = end
            tasks.append(_HornerTask(0, 0, task.level, exponents^))
            for i in range(len(groups)):
                budget.charge(1)
                var group = groups[len(groups) - 1 - i]
                tasks.append(_HornerTask(group[0], group[1], task.level + 1, List[UInt32]()))
    return values.pop()

def evaluate_dense[C: CoefficientRing](p: DensePoly[C], x: C, mut budget: Budget) raises PolyError -> C:
    """Evaluate a dense polynomial by Horner's rule.

    Parameters:
        C: The coefficient ring type.

    Args:
        p: The polynomial.
        x: The value in the coefficient ring.
        budget: The shared work allowance.

    Returns:
        The exact value.

    Raises:
        PolyError: DOMAIN for another ring, LIMIT or CAPACITY on work limits.
    """
    if x.context() != p.ring:
        raise PolyError(PolyErrorKind.DOMAIN, "evaluation point in another ring")
    var result = C.zero(p.ring)
    for i in range(len(p.coefficients)):
        budget.charge(1)
        result = result.mul(x, budget).add(p.coefficients[len(p.coefficients) - 1 - i], budget)
    return result^

def substitute[C: CoefficientRing](p: SparsePoly[C], index: Int64, q: SparsePoly[C], mut budget: Budget) raises PolyError -> SparsePoly[C]:
    """Substitute a polynomial for one variable.

    Parameters:
        C: The coefficient ring type.

    Args:
        p: The source polynomial.
        index: The variable to replace.
        q: The replacement in the same ring and order.
        budget: The shared work allowance.

    Returns:
        The canonical substituted polynomial.

    Raises:
        PolyError: DOMAIN for incompatible operands, LIMIT or CAPACITY on work limits.
    """
    p._same(q, budget)
    if len(p.terms) >= 64:
        var groups = _CoefficientGroups(p, index, budget)
        var count = len(groups.starts) - 1
        var result = groups.coefficient(p, count - 1, budget)
        var last = groups.exponent(count - 1)
        for i in range(count - 2, -1, -1):
            var exponent = groups.exponent(i)
            result = result.mul(q.pow(UInt64(last - exponent), budget), budget).add(groups.coefficient(p, i, budget), budget)
            last = exponent
        if last:
            result = result.mul(q.pow(UInt64(last), budget), budget)
        return result^
    var grouped = _coefficient_groups(p, index, False, budget, first=False)
    var result = SparsePoly[C].zero(p.order, p.ring, budget)
    if len(grouped[0]) == 0:
        return result^
    var last = grouped[0][len(grouped[0]) - 1].exponent(index, budget)
    swap(result, grouped[1][len(grouped[1]) - 1])
    # Sparse Horner steps depend on occupied exponents, never the degree box.
    for i in range(len(grouped[0]) - 2, -1, -1):
        budget.charge(1)
        var exponent = grouped[0][i].exponent(index, budget)
        result = result.mul(q.pow(UInt64(last - exponent), budget), budget).add(grouped[1][i], budget)
        last = exponent
    if last:
        result = result.mul(q.pow(UInt64(last), budget), budget)
    return result^


def derivative[C: CoefficientRing](p: SparsePoly[C], index: Int64, mut budget: Budget) raises PolyError -> SparsePoly[C]:
    """Differentiate formally in one variable in any coefficient ring.

    Parameters:
        C: The coefficient ring type.

    Args:
        p: The polynomial.
        index: The variable to differentiate.
        budget: The shared work allowance.

    Returns:
        The canonical formal derivative.

    Raises:
        PolyError: LIMIT or CAPACITY on work limits.
    """
    var terms = List[Term[C]]()
    for t in p.terms:
        budget.charge(1)
        var exponent = t.monomial.exponent(index, budget)
        if exponent == 0:
            continue
        var pairs = List[Tuple[Int64, UInt32]]()
        for pair in t.monomial.pairs:
            budget.charge(1)
            if pair[0] != index:
                pairs.append(pair)
            elif pair[1] > 1:
                pairs.append((index, pair[1] - 1))
        var c = t.coefficient.mul(_natural[C](UInt64(exponent), p.ring, budget), budget)
        if not c.is_zero():
            terms.append(Term[C](Monomial(pairs^, budget), c^))
    return SparsePoly[C].from_terms(terms^, p.order, p.ring, budget)


def _derivative_dense[C: CoefficientRing](f: DensePoly[C], mut budget: Budget) raises PolyError -> DensePoly[C]:
    _allocation[C](max(f.degree(), 0), budget)
    var values = List[C](capacity=max(f.degree(), 0))
    for i in range(1, len(f.coefficients)):
        budget.charge(1)
        values.append(f.coefficients[i].mul(_natural[C](UInt64(i), f.ring, budget), budget))
    return DensePoly[C].from_coefficients(f.variable, f.ring, values^, budget)
