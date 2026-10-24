"""Private exact arithmetic and Lazard specialization in one sample field.

Coefficients are reduced rational vectors, sharing one modulus and one real
embedding per stack. There are no per-coefficient algebraic numbers or operation
histories. A verifier validates each context before using its exact arithmetic.
"""
from apn_mojo import Integer, Rational
from .coefficients import QQ, NoContext, _pow
from .dense import DensePoly
from .sparse import SparsePoly
from .monomial import Monomial
from .evaluate import derivative
from .errors import Budget, PolyError, PolyErrorKind, _allocation, _sum, _product
from .algebraic_types import RealAlgebraic, _check_real
from ._algebraic import _Frame, _refine
from ._algebraic_filter import _constant, _reduce, _field_mul, _field_inverse, _dense_interval, _gcd_over_field
from ._root_ops import _compare
from .factor import factor
from .certificates import check_factorization
from ._unit_bezout import _unit_bezout
from ._dyadic import _of_real, _dense_value, _precision, _point, _zero as _dyadic_zero, _add as _dyadic_add, _mul as _dyadic_mul


@fieldwise_init
struct _SampleField(Copyable, Movable):
    var modulus: DensePoly[QQ]
    var embedding: Optional[RealAlgebraic]
    var variables: List[Int64]
    var coordinates: List[DensePoly[QQ]]


def _sample_field(frame: _Frame, mut budget: Budget) raises PolyError -> _SampleField:
    if not frame.valid:
        raise PolyError(PolyErrorKind.DOMAIN, "sample has no common field")
    var modulus = DensePoly[QQ].from_coefficients(0, NoContext(), [QQ.zero(NoContext()), QQ.one(NoContext())], budget)
    if frame.generator:
        modulus = frame.generator.value().polynomial.copy()
    _allocation[Int64](len(frame.variables), budget)
    _allocation[DensePoly[QQ]](len(frame.representations), budget)
    return _SampleField(modulus^, frame.generator.copy(), frame.variables.copy(), frame.representations.copy())


def _check_sample_field(field: _SampleField, mut budget: Budget, *, modulus_checked: Bool = False) raises PolyError -> Bool:
    """Check the embedding/maps; only a verifier-owned cache may skip factoring."""
    field.modulus.check_invariants(budget)
    if field.modulus.variable != 0 or field.modulus.degree() < 1 or len(field.variables) != len(field.coordinates):
        return False
    if field.modulus.coefficients[field.modulus.degree()] != QQ.one(NoContext()):
        return False
    if field.embedding:
        ref gamma = field.embedding.value()
        if gamma.rational or not _check_real(gamma, budget) or gamma.polynomial != field.modulus:
            return False
        # _check_real proves an isolated root, not irreducibility. Recompute
        # complete rational factorization once for this context, never per value.
        if not modulus_checked:
            var factors = factor(field.modulus, budget)
            if not check_factorization(field.modulus, factors, budget) or len(factors.factors) != 1:
                return False
            if factors.factors[0][1] != 1 or factors.factors[0][0] != field.modulus:
                return False
    else:
        if field.modulus.degree() != 1 or not field.modulus.coefficients[0].is_zero():
            return False
    for i in range(len(field.variables)):
        budget.charge(1)
        if field.variables[i] < 0:
            return False
        for j in range(i):
            budget.charge(1)
            if field.variables[i] == field.variables[j]:
                return False
        ref c = field.coordinates[i]
        c.check_invariants(budget)
        if c.variable != 0 or c.degree() >= field.modulus.degree():
            return False
    return True


def _ksign(value: DensePoly[QQ], mut field: _SampleField, mut budget: Budget) raises PolyError -> Int:
    var reduced = _reduce(value, field.modulus, budget)
    if reduced.is_zero():
        return 0
    if reduced.degree() == 0:
        return reduced.coefficients[0].value.sign()
    if not field.embedding:
        raise PolyError(PolyErrorKind.DOMAIN, "nonrational value without a field embedding")
    while True:
        budget.charge(1)
        ref embedding = field.embedding.value()
        var precision = _precision(embedding.low, embedding.high)
        var enclosure = _dense_value(reduced, _of_real(embedding, precision, budget), precision, budget)
        var sign = enclosure.sign()
        if sign:
            return sign
        # An irreducible modulus and nonzero reduced vector imply nonzero
        # value at this embedding. Outward-rounded dyadic Horner bounds avoid
        # growing the exact interval denominator at every polynomial degree.
        # Refinement increases both embedding accuracy and working precision,
        # so an inconclusive enclosure eventually separates the nonzero value.
        _refine(field.embedding.value(), budget)


def _kmultiply(a: DensePoly[QQ], b: DensePoly[QQ], modulus: DensePoly[QQ], mut budget: Budget) raises PolyError -> DensePoly[QQ]:
    # Multiplying by a rational scalar cannot raise the generator degree.
    # Avoid clearing the entire minimal polynomial when no reduction occurs.
    if a.variable == b.variable and (a.degree() <= 0 or b.degree() <= 0) and max(a.degree(), b.degree()) < modulus.degree():
        return a.mul(b, budget)
    return _field_mul(a, b, modulus, budget)


def _kpow(value: DensePoly[QQ], exponent: UInt32, field: _SampleField, mut budget: Budget) raises PolyError -> DensePoly[QQ]:
    var base = _reduce(value, field.modulus, budget)
    if exponent == 1:
        return base^
    var result = _constant(QQ.one(NoContext()), budget)
    var power = exponent
    while power:
        budget.charge(1)
        if power & 1:
            result = _kmultiply(result, base, field.modulus, budget)
        power >>= 1
        if power:
            base = _kmultiply(base, base, field.modulus, budget)
    return result^


@fieldwise_init
struct _KTerm(Copyable, Movable):
    var monomial: Monomial
    var value: DensePoly[QQ]


def _partial_field(f: SparsePoly[QQ], field: _SampleField, prefix: Int, mut budget: Budget) raises PolyError -> List[_KTerm]:
    """Substitute exactly a coordinate prefix; retain every other variable."""
    if prefix < 0 or prefix > len(field.variables):
        raise PolyError(PolyErrorKind.DOMAIN, "invalid field specialization prefix")
    var terms = List[_KTerm]()
    for term in f.terms:
        budget.charge(1)
        var scalar = term.coefficient.copy()
        var algebraic: Optional[DensePoly[QQ]] = None
        var monomial = Monomial()
        for pair in term.monomial.pairs:
            var index = -1
            for i in range(prefix):
                budget.charge(1)
                if field.variables[i] == pair[0]:
                    index = i
                    break
            if index >= 0:
                ref coordinate = field.coordinates[index]
                if coordinate.is_zero():
                    scalar = QQ.zero(NoContext())
                    break
                if coordinate.degree() == 0:
                    var power = coordinate.coefficients[0].copy() if pair[1] == 1 else _pow(coordinate.coefficients[0], UInt64(pair[1]), budget)
                    scalar = scalar.mul(power, budget)
                else:
                    var power = _kpow(coordinate, pair[1], field, budget)
                    if algebraic:
                        algebraic = _kmultiply(algebraic.value(), power, field.modulus, budget)
                    else:
                        algebraic = power^
            else:
                _allocation[Tuple[Int64, UInt32]](1, budget)
                monomial.pairs.append(pair)
        if scalar.is_zero():
            continue
        var value = _constant(scalar, budget)
        if algebraic:
            value = _kmultiply(value, algebraic.value(), field.modulus, budget)
        if value.is_zero():
            continue
        var found = -1
        for i in range(len(terms)):
            budget.charge(UInt64(_sum(len(terms[i].monomial.pairs), len(monomial.pairs))))
            if terms[i].monomial == monomial:
                found = i
                break
        if found >= 0:
            terms[found].value = terms[found].value.add(value, budget)
        else:
            _allocation[_KTerm](1, budget)
            terms.append(_KTerm(monomial^, value^))
    # Compact surviving terms in their existing storage. They are uniquely
    # owned here, and cancellation must retain the order of the survivors.
    var kept = 0
    for i in range(len(terms)):
        budget.charge(1)
        if not terms[i].value.is_zero():
            if kept != i:
                # List indices share an origin in Mojo; move through a local
                # so the two mutable borrows cannot alias.
                var surviving = _KTerm(Monomial(), DensePoly[QQ](0, NoContext()))
                swap(surviving, terms[i])
                swap(surviving, terms[kept])
            kept += 1
    while len(terms) > kept:
        _ = terms.pop()
    return terms^


def _specialize_field(f: SparsePoly[QQ], variable: Int64, field: _SampleField, mut budget: Budget) raises PolyError -> List[DensePoly[QQ]]:
    var terms = _partial_field(f, field, len(field.variables), budget)
    var size = 0
    for term in terms:
        for pair in term.monomial.pairs:
            budget.charge(1)
            if pair[0] != variable:
                raise PolyError(PolyErrorKind.DOMAIN, "missing specialization coordinate")
        size = max(size, _sum(Int(term.monomial.exponent(variable, budget)), 1))
    _allocation[DensePoly[QQ]](size, budget)
    var result = List[DensePoly[QQ]](length=size, fill=DensePoly[QQ](0, NoContext()))
    for i in range(len(terms)):
        ref term = terms[i]
        var index = Int(term.monomial.exponent(variable, budget))
        if result[index].is_zero():
            swap(result[index], term.value)
        else:
            result[index] = result[index].add(term.value, budget)
    return result^


def _lazard_field(f: SparsePoly[QQ], variable: Int64, field: _SampleField, mut budget: Budget) raises PolyError -> Tuple[List[Int], List[DensePoly[QQ]]]:
    if f.is_zero():
        raise PolyError(PolyErrorKind.DOMAIN, "Lazard evaluation requires a nonzero polynomial")
    for assigned in field.variables:
        budget.charge(1)
        if assigned == variable:
            raise PolyError(PolyErrorKind.DOMAIN, "lifting variable is already assigned")
    # A nonzero full specialization implies every preceding specialization
    # is nonzero. Most stacks therefore need one substitution, not a separate
    # substitution of every coordinate prefix followed by the full one.
    var direct = _specialize_field(f, variable, field, budget)
    if len(direct):
        _allocation[Int](len(field.variables), budget)
        return (List[Int](length=len(field.variables), fill=0), direct^)
    var current = f._copy(budget)
    _allocation[Int](len(field.variables), budget)
    var orders = List[Int]()
    for i in range(len(field.variables)):
        var order = 0
        while len(_partial_field(current, field, i + 1, budget)) == 0:
            current = derivative(current, field.variables[i], budget)
            order = _sum(order, 1)
            if current.is_zero():
                raise PolyError(PolyErrorKind.DOMAIN, "Lazard derivative chain unexpectedly vanished")
        orders.append(order)
    return (orders^, _specialize_field(current, variable, field, budget))


def _ktrim(mut f: List[DensePoly[QQ]]):
    while len(f) and f[len(f) - 1].is_zero():
        _ = f.pop()


def _kcopy(f: List[DensePoly[QQ]], mut budget: Budget) raises PolyError -> List[DensePoly[QQ]]:
    _allocation[DensePoly[QQ]](len(f), budget)
    return f.copy()


def _kmul(f: List[DensePoly[QQ]], g: List[DensePoly[QQ]], field: _SampleField, mut budget: Budget) raises PolyError -> List[DensePoly[QQ]]:
    if len(f) == 0 or len(g) == 0:
        return List[DensePoly[QQ]]()
    var size = _sum(len(f), len(g)) - 1
    _allocation[DensePoly[QQ]](size, budget)
    var result = List[DensePoly[QQ]](length=size, fill=DensePoly[QQ](0, NoContext()))
    for i in range(len(f)):
        for j in range(len(g)):
            result[i + j] = result[i + j].add(_kmultiply(f[i], g[j], field.modulus, budget), budget)
    _ktrim(result)
    return result^


def _kscale(f: List[DensePoly[QQ]], scalar: DensePoly[QQ], field: _SampleField, mut budget: Budget) raises PolyError -> List[DensePoly[QQ]]:
    if scalar.is_zero():
        return List[DensePoly[QQ]]()
    if scalar.degree() == 0 and scalar.coefficients[0] == QQ.one(NoContext()):
        return _kcopy(f, budget)
    _allocation[DensePoly[QQ]](len(f), budget)
    var result = List[DensePoly[QQ]]()
    for c in f:
        result.append(_kmultiply(c, scalar, field.modulus, budget))
    _ktrim(result)
    return result^


def _kinverse(value: DensePoly[QQ], field: _SampleField, mut budget: Budget) raises PolyError -> DensePoly[QQ]:
    if value.degree() == 0:
        return _constant(value.coefficients[0].inverse(budget), budget)
    if field.modulus.degree() >= 4:
        # Reconstruct just the inverse via modular arithmetic. The existing
        # kernel verifies u*value + v*modulus == 1 exactly before returning;
        # only u survives as the field value, with no Euclidean history.
        var modular = _unit_bezout(value, field.modulus, budget)
        if modular:
            return _reduce(modular.value()[0], field.modulus, budget)
    return _field_inverse(value, field.modulus, budget)


def _kmonic(f: List[DensePoly[QQ]], field: _SampleField, mut budget: Budget) raises PolyError -> List[DensePoly[QQ]]:
    if len(f) == 0:
        return List[DensePoly[QQ]]()
    if f[len(f) - 1].degree() == 0 and f[len(f) - 1].coefficients[0] == QQ.one(NoContext()):
        return _kcopy(f, budget)
    return _kscale(f, _kinverse(f[len(f) - 1], field, budget), field, budget)


def _kderivative(f: List[DensePoly[QQ]], field: _SampleField, mut budget: Budget) raises PolyError -> List[DensePoly[QQ]]:
    _allocation[DensePoly[QQ]](max(0, len(f) - 1), budget)
    var result = List[DensePoly[QQ]]()
    for i in range(1, len(f)):
        result.append(_kmultiply(f[i], _constant(QQ(Rational(Integer(i))), budget), field.modulus, budget))
    _ktrim(result)
    return result^


def _kdivmod(f: List[DensePoly[QQ]], g: List[DensePoly[QQ]], field: _SampleField, mut budget: Budget) raises PolyError -> Tuple[List[DensePoly[QQ]], List[DensePoly[QQ]]]:
    if len(g) == 0:
        raise PolyError(PolyErrorKind.DOMAIN, "zero field polynomial divisor")
    if len(f) < len(g):
        return (List[DensePoly[QQ]](), _kcopy(f, budget))
    if len(g) == 1:
        return (_kscale(f, _kinverse(g[0], field, budget), field, budget), List[DensePoly[QQ]]())
    var remainder = _kcopy(f, budget)
    var size = max(0, len(f) - len(g) + 1)
    _allocation[DensePoly[QQ]](size, budget)
    var quotient = List[DensePoly[QQ]](length=size, fill=DensePoly[QQ](0, NoContext()))
    var inverse = _kinverse(g[len(g) - 1], field, budget)
    while len(remainder) >= len(g):
        var shift = len(remainder) - len(g)
        var top = _kmultiply(remainder[len(remainder) - 1], inverse, field.modulus, budget)
        for j in range(len(g)):
            remainder[shift + j] = remainder[shift + j].sub(_kmultiply(top, g[j], field.modulus, budget), budget)
        quotient[shift] = top^
        _ktrim(remainder)
    _ktrim(quotient)
    return (quotient^, remainder^)


def _kexact(f: List[DensePoly[QQ]], g: List[DensePoly[QQ]], field: _SampleField, mut budget: Budget) raises PolyError -> List[DensePoly[QQ]]:
    var result = _kdivmod(f, g, field, budget)
    if len(result[1]):
        raise PolyError(PolyErrorKind.DOMAIN, "inexact field polynomial division")
    var quotient = List[DensePoly[QQ]]()
    swap(quotient, result[0])
    return quotient^


def _kevaluate(f: List[DensePoly[QQ]], x: QQ, field: _SampleField, mut budget: Budget) raises PolyError -> DensePoly[QQ]:
    var value = DensePoly[QQ](0, NoContext())
    var argument = _constant(x, budget)
    for i in range(len(f) - 1, -1, -1):
        value = _kmultiply(value, argument, field.modulus, budget).add(f[i], budget)
    return value^


def _kevaluate_sign(f: List[DensePoly[QQ]], x: QQ, mut field: _SampleField, mut budget: Budget) raises PolyError -> Int:
    """Filter a rational endpoint directly, without first forming f(x) in K."""
    if field.embedding:
        ref gamma = field.embedding.value()
        var precision = max(_precision(gamma.low, gamma.high), _precision(x, x))
        var generator = _of_real(gamma, precision, budget)
        var argument = _point(x, precision, budget)
        var value = _dyadic_zero()
        for i in range(len(f) - 1, -1, -1):
            value = _dyadic_add(_dyadic_mul(value, argument, precision, budget), _dense_value(f[i], generator, precision, budget), precision, budget)
        var sign = value.sign()
        if sign:
            return sign
    # This also recognizes exact zeros; an interval containing zero never
    # supplies a zero decision. Only the reduced exact field value can do so.
    return _ksign(_kevaluate(f, x, field, budget), field, budget)
