"""Complete algebraic factorization by bounded Kronecker recombination.

The producer splits a worklist by exact division. The checker consumes retained
univariate factorizations and independently exhausts divisor multiplicities;
it never calls the multivariate producer or re-solves an image.
"""
from polynomials import Budget, PolyError, PolyErrorKind, CoefficientField, CoefficientRing, DensePoly, SparsePoly, Term, Monomial, Factorization, MultivariateFactorization, divide_exact, squarefree_multivariate, content, MonomialOrder
from ._common import _allocation, _charge_size, _sum, _product, _grow
from ._algebraic_value import Real
from . import _algebraic_factor as _uni
from . import _complex_roots as _roots


def _copy[C: CoefficientRing](value: C, mut budget: Budget) raises PolyError -> C:
    _charge_size(value.size(), budget)
    return value.copy()


def _normalize[C: CoefficientField](p: SparsePoly[C], mut budget: Budget) raises PolyError -> SparsePoly[C]:
    return p.scale(p.terms[len(p.terms) - 1].coefficient.inverse(budget), budget)


def _compare[C: CoefficientRing](a: SparsePoly[C], b: SparsePoly[C], mut budget: Budget) raises PolyError -> Int:
    var da = a.total_degree(budget)
    var db = b.total_degree(budget)
    if da != db:
        return -1 if da < db else 1
    for offset in range(min(len(a.terms), len(b.terms))):
        budget.charge(1)
        ref x = a.terms[len(a.terms) - 1 - offset]
        ref y = b.terms[len(b.terms) - 1 - offset]
        var relation = a.order.compare(x.monomial, y.monomial, budget)
        if relation:
            return relation
        relation = _roots._compare(_uni._complex(x.coefficient), _uni._complex(y.coefficient), budget)
        if relation:
            return relation
    return -1 if len(a.terms) < len(b.terms) else 1 if len(a.terms) > len(b.terms) else 0


def _simple[C: CoefficientRing](p: SparsePoly[C], mut budget: Budget) raises PolyError -> Bool:
    # A*x_j+B(other variables), with A a nonzero scalar, has no nonunit
    # divisor in the parameter ring and is linear over its fraction field.
    for term in p.terms:
        budget.charge(1)
        if len(term.monomial.pairs) != 1 or term.monomial.pairs[0][1] != 1:
            continue
        var axis = term.monomial.pairs[0][0]
        var count = 0
        for other in p.terms:
            budget.charge(1)
            if other.monomial.exponent(axis, budget):
                count += 1
        if count == 1:
            return True
    return False


def _box[C: CoefficientRing](p: SparsePoly[C], mut budget: Budget) raises PolyError -> Tuple[List[Int64], List[Int], Int]:
    var variables = p.variables(budget)
    _allocation[Int](len(variables), budget)
    var radices = List[Int](capacity=len(variables))
    var size = 1
    for axis in variables:
        budget.charge(1)
        var radix = _sum(p.degree_in(axis, budget), 1)
        radices.append(radix)
        size = _product(size, radix)
    return variables^, radices^, size


def _encode[C: CoefficientRing](p: SparsePoly[C], variables: List[Int64], radices: List[Int], size: Int, mut budget: Budget) raises PolyError -> DensePoly[C]:
    # Account for the entire degree box before any dense allocation.
    _allocation[C](size, budget)
    var coefficients = List[C](length=size, fill=C.zero(p.ring))
    for term in p.terms:
        budget.charge(1)
        var index = 0
        for j in range(len(variables) - 1, -1, -1):
            budget.charge(1)
            index = _sum(_product(index, radices[j]), Int(UInt64(term.monomial.exponent(variables[j], budget))))
        coefficients[index] = _copy(term.coefficient, budget)
    return DensePoly[C].from_coefficients(0, p.ring, coefficients^, budget)


def _decode[C: CoefficientField](p: DensePoly[C], variables: List[Int64], radices: List[Int], template: SparsePoly[C], mut budget: Budget) raises PolyError -> SparsePoly[C]:
    var terms = List[Term[C]]()
    for i in range(len(p.coefficients)):
        budget.charge(1)
        if p.coefficients[i].is_zero():
            continue
        _allocation[Tuple[Int64, UInt32]](len(variables), budget)
        var pairs = List[Tuple[Int64, UInt32]](capacity=len(variables))
        var power = i
        for j in range(len(variables)):
            budget.charge(1)
            var exponent = power % radices[j]
            power //= radices[j]
            if exponent:
                pairs.append((variables[j], UInt32(exponent)))
        if power:
            raise PolyError(PolyErrorKind.DOMAIN, "algebraic factor image exceeds its degree box")
        _grow(terms, budget)
        terms.append(Term[C](Monomial(pairs^, budget), _copy(p.coefficients[i], budget)))
    return _normalize(SparsePoly[C].from_terms(terms^, template.order, template.ring, budget), budget)


def _candidate[C: CoefficientRing](proof: Factorization[C], counts: List[Int], degree: Int, mut budget: Budget) raises PolyError -> Optional[DensePoly[C]]:
    var selected_degree = 0
    for i in range(len(counts)):
        budget.charge(1)
        selected_degree = _sum(selected_degree, _product(counts[i], proof.factors[i][0].degree()))
    # Of two proper complementary factors, one has image degree at most d/2.
    if not selected_degree or selected_degree > degree // 2:
        return None
    _allocation[C](1, budget)
    var value = DensePoly[C].from_coefficients(0, proof.unit.context(), [C.one(proof.unit.context())], budget)
    for i in range(len(counts)):
        budget.charge(1)
        if counts[i]:
            value = value.mul(proof.factors[i][0].pow(UInt64(counts[i]), budget), budget)
    return value^


def _find_divisor[C: CoefficientField](p: SparsePoly[C], variables: List[Int64], radices: List[Int], image: DensePoly[C], proof: Factorization[C], mut budget: Budget) raises PolyError -> Optional[SparsePoly[C]]:
    _allocation[Int](len(proof.factors), budget)
    var counts = List[Int](length=len(proof.factors), fill=0)
    var degree = p.total_degree(budget)
    while True:
        budget.charge(1)
        # Ascending multiplicity vectors: identical image factors are never
        # expanded into exponentially many duplicate occurrence subsets.
        var j = len(counts) - 1
        while j >= 0 and counts[j] == proof.factors[j][1]:
            budget.charge(1)
            counts[j] = 0
            j -= 1
        if j < 0:
            return None
        counts[j] += 1
        var selected = _candidate(proof, counts, image.degree(), budget)
        if not selected:
            continue
        var divisor = _decode(selected.value(), variables, radices, p, budget)
        var d = divisor.total_degree(budget)
        if d > 0 and d < degree and divide_exact(p, divisor, budget):
            return divisor^


def _no_divisor[C: CoefficientField](p: SparsePoly[C], variables: List[Int64], radices: List[Int], image: DensePoly[C], proof: Factorization[C], mut budget: Budget) raises PolyError -> Bool:
    # Independent descending traversal. This checker consumes a verified
    # complete image factorization, never a producer success flag or search log.
    _allocation[Int](len(proof.factors), budget)
    var counts = List[Int](capacity=len(proof.factors))
    for entry in proof.factors:
        budget.charge(1)
        counts.append(entry[1])
    var degree = p.total_degree(budget)
    while True:
        budget.charge(1)
        var position = 0
        while position < len(counts) and counts[position] == 0:
            budget.charge(1)
            counts[position] = proof.factors[position][1]
            position += 1
        if position == len(counts):
            return True
        counts[position] -= 1
        var selected = _candidate(proof, counts, image.degree(), budget)
        if selected:
            var divisor = _decode(selected.value(), variables, radices, p, budget)
            var d = divisor.total_degree(budget)
            if d > 0 and d < degree and divide_exact(p, divisor, budget):
                return False


def _irreducible[C: CoefficientField](p: SparsePoly[C], image_proof: Optional[Factorization[C]], mut budget: Budget) raises PolyError -> Bool:
    if _simple(p, budget):
        return not image_proof
    var axis = _uni._axis(p, budget)
    if axis[0]:
        if image_proof:
            return False
        comptime maximum = 2 if C == Real else 1
        if p.degree_in(axis[1], budget) > maximum:
            return False
        var dense = DensePoly[C].from_sparse(p, budget, variable=axis[1])
        _allocation[Tuple[DensePoly[C], Int]](1, budget)
        _allocation[C](len(dense.coefficients), budget)
        return _uni._check(dense, Factorization[C](C.one(p.ring), [(dense.copy(), 1)]), budget)
    if not image_proof:
        return False
    var box = _box(p, budget)
    var image = _encode(p, box[0], box[1], box[2], budget)
    if not _uni._check(image, image_proof.value(), budget):
        return False
    return _no_divisor(p, box[0], box[1], image, image_proof.value(), budget)


def _valid[C: CoefficientRing](p: SparsePoly[C], mut budget: Budget) raises PolyError -> Bool:
    p.check_invariants(budget)
    for term in p.terms:
        budget.charge(1)
        if not _roots._valid(_uni._complex(term.coefficient), budget):
            return False
    return True


def _check[C: CoefficientField](p: SparsePoly[C], proof: MultivariateFactorization[C], images: List[Optional[Factorization[C]]], mut budget: Budget) raises PolyError -> Bool:
    if not _valid(p, budget) or not _roots._valid(_uni._complex(proof.unit), budget):
        return False
    var degree = p.total_degree(budget)
    if degree <= 0:
        return not len(proof.factors) and not len(images) and (proof.unit.is_zero() if p.is_zero() else proof.unit.sub(p.terms[0].coefficient, budget).is_zero())
    if proof.unit.is_zero() or len(proof.factors) > degree or len(images) != len(proof.factors):
        return False
    var total = 0
    for i in range(len(proof.factors)):
        budget.charge(1)
        ref q = proof.factors[i][0]
        var multiplicity = proof.factors[i][1]
        if not _valid(q, budget):
            return False
        _charge_size(_sum(len(p.order.variables), len(q.order.variables)), budget)
        if q.order != p.order:
            return False
        var d = q.total_degree(budget)
        if d < 1 or multiplicity <= 0 or multiplicity > (degree - total) // d:
            return False
        total += d * multiplicity
        # Degree bounds reject malicious huge proofs before encoding or powers.
        for axis in q.variables(budget):
            budget.charge(1)
            if q.degree_in(axis, budget) > p.degree_in(axis, budget) // multiplicity:
                return False
        if not q.terms[len(q.terms) - 1].coefficient.sub(C.one(p.ring), budget).is_zero():
            return False
        if i and _compare(proof.factors[i - 1][0], q, budget) >= 0:
            return False
    if total != degree:
        return False
    var product = SparsePoly[C].constant(_copy(proof.unit, budget), p.order, budget)
    for entry in proof.factors:
        budget.charge(1)
        product = product.mul(entry[0].pow(UInt64(entry[1]), budget), budget)
    if not product.sub(p, budget).is_zero():
        return False
    for i in range(len(proof.factors)):
        budget.charge(1)
        if not _irreducible(proof.factors[i][0], images[i], budget):
            return False
    return True


def _factor[C: CoefficientField](p: SparsePoly[C], mut budget: Budget) raises PolyError -> Tuple[MultivariateFactorization[C], List[Optional[Factorization[C]]]]:
    if not _valid(p, budget):
        raise PolyError(PolyErrorKind.DOMAIN, "invalid algebraic coefficient witness")
    var factors = List[Tuple[SparsePoly[C], Int]]()
    var images = List[Optional[Factorization[C]]]()
    var remaining = p.reorder(p.order, budget)
    # Remove coordinate factors before squarefree decomposition and dense images.
    for axis in p.variables(budget):
        budget.charge(1)
        var exponent = p.degree_in(axis, budget)
        for term in p.terms:
            budget.charge(1)
            exponent = min(exponent, Int(UInt64(term.monomial.exponent(axis, budget))))
        if exponent:
            var coordinate = SparsePoly[C].variable(axis, p.order, p.ring, budget)
            var quotient = divide_exact(remaining, coordinate.pow(UInt64(exponent), budget), budget)
            if not quotient:
                raise PolyError(PolyErrorKind.DOMAIN, "coordinate factor reconstruction failed")
            remaining = quotient.take()
            _grow(factors, budget)
            _grow(images, budget)
            factors.append((coordinate^, exponent))
            images.append(None)
    var blocks = squarefree_multivariate(remaining, budget)
    _allocation[Tuple[SparsePoly[C], Int]](len(blocks.factors), budget)
    var pending = List[Tuple[SparsePoly[C], Int]](capacity=len(blocks.factors))
    for entry in blocks.factors:
        budget.charge(1)
        pending.append((entry[0].reorder(p.order, budget), entry[1]))
    while len(pending):
        budget.charge(1)
        var task = pending.pop()
        var value = _normalize(task[0], budget)
        var multiplicity = task[1]
        if _simple(value, budget):
            _grow(factors, budget)
            _grow(images, budget)
            factors.append((value^, multiplicity))
            images.append(None)
            continue
        var axis = _uni._axis(value, budget)
        if axis[0]:
            var split = _uni._factor(DensePoly[C].from_sparse(value, budget, variable=axis[1]), budget)
            for entry in split.factors:
                budget.charge(1)
                _grow(factors, budget)
                _grow(images, budget)
                factors.append((_normalize(entry[0].to_sparse(p.order, budget), budget), _product(multiplicity, entry[1])))
                images.append(None)
            continue
        # Factors independent of the first variable are coefficient content.
        # Extract them before constructing an unnecessarily large root problem.
        var variables = value.variables(budget)
        _allocation[Int64](len(variables), budget)
        var order = MonomialOrder(MonomialOrder.LEX, variables.copy(), budget)
        var common = content(value.reorder(order, budget), variables[0], budget).reorder(p.order, budget)
        if common.total_degree(budget) > 0:
            var quotient = divide_exact(value, common, budget)
            if not quotient or quotient.value().total_degree(budget) <= 0:
                raise PolyError(PolyErrorKind.DOMAIN, "algebraic coefficient content is not a proper divisor")
            _grow(pending, budget)
            pending.append((common^, multiplicity))
            _grow(pending, budget)
            pending.append((quotient.take(), multiplicity))
            continue
        var box = _box(value, budget)
        var image = _encode(value, box[0], box[1], box[2], budget)
        var split = _uni._factor(image, budget)
        var divisor = _find_divisor(value, box[0], box[1], image, split, budget)
        if divisor:
            var q = divisor.take()
            var complement = divide_exact(value, q, budget)
            if not complement:
                raise PolyError(PolyErrorKind.DOMAIN, "algebraic factor reconstruction failed")
            _grow(pending, budget)
            pending.append((q^, multiplicity))
            _grow(pending, budget)
            pending.append((complement.take(), multiplicity))
        else:
            _grow(factors, budget)
            _grow(images, budget)
            factors.append((value^, multiplicity))
            images.append(split^)
    for i in range(1, len(factors)):
        var j = i
        while j and _compare(factors[j][0], factors[j - 1][0], budget) < 0:
            budget.charge(1)
            factors.swap_elements(j, j - 1)
            images.swap_elements(j, j - 1)
            j -= 1
    # Recover the scalar by exact division. KERNEL is a storage order whose
    # greatest terms need not multiply to the greatest term of the product.
    remaining = p.reorder(p.order, budget)
    for entry in factors:
        budget.charge(1)
        var quotient = divide_exact(remaining, entry[0].pow(UInt64(entry[1]), budget), budget)
        if not quotient:
            raise PolyError(PolyErrorKind.DOMAIN, "algebraic factors do not reconstruct the source")
        remaining = quotient.take()
    if remaining.total_degree(budget) > 0:
        raise PolyError(PolyErrorKind.DOMAIN, "algebraic factorization left a nonconstant cofactor")
    var unit = C.zero(p.ring) if remaining.is_zero() else _copy(remaining.terms[0].coefficient, budget)
    var proof = MultivariateFactorization[C](unit^, factors^)
    if not _check(p, proof, images, budget):
        raise PolyError(PolyErrorKind.DOMAIN, "algebraic multivariate factorization failed independent verification")
    return proof^, images^
