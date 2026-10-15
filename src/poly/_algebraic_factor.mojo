"""AA and QQbar factorization dispatch and univariate irreducibility checks."""
from polynomials import Budget, PolyError, PolyErrorKind, NoContext, CoefficientRing, DensePoly, SparsePoly, Factorization, MultivariateFactorization
from ._algebraic_value import Real, Complex
from . import _complex_roots as _roots
from ._common import _allocation, _charge_size


def _complex[C: CoefficientRing](value: C) -> Complex:
    comptime if C == Real:
        return Complex(rebind[Real](value).copy(), Real.zero(NoContext()))
    else:
        comptime assert C == Complex, "algebraic factorization requires AA or QQbar"
        return rebind[Complex](value).copy()


def _compare[C: CoefficientRing](a: DensePoly[C], b: DensePoly[C], mut budget: Budget) raises PolyError -> Int:
    budget.charge(1)
    if a.degree() != b.degree():
        return -1 if a.degree() < b.degree() else 1
    for i in range(a.degree(), -1, -1):
        budget.charge(1)
        var relation = _roots._compare(_complex(a.coefficients[i]), _complex(b.coefficients[i]), budget)
        if relation:
            return relation
    return 0


def _sort[C: CoefficientRing](mut factors: List[Tuple[DensePoly[C], Int]], mut budget: Budget) raises PolyError:
    for i in range(1, len(factors)):
        budget.charge(1)
        var j = i
        while j > 0 and _compare(factors[j][0], factors[j - 1][0], budget) < 0:
            budget.charge(1)
            factors.swap_elements(j, j - 1)
            j -= 1


def _check[C: CoefficientRing](p: DensePoly[C], proof: Factorization[C], mut budget: Budget) raises PolyError -> Bool:
    # No solving or factoring of the source polynomial: degree and discriminant
    # establish irreducibility, then the product establishes completeness.
    # Scalar normalization and arithmetic still verify their own evidence.
    p.check_invariants(budget)
    if not _roots._valid(_complex(proof.unit), budget):
        return False
    for coefficient_index in range(len(p.coefficients)):
        ref c = p.coefficients[coefficient_index]
        if not _roots._valid(_complex(c), budget):
            return False
    if p.degree() <= 0:
        if len(proof.factors):
            return False
        return proof.unit.is_zero() if p.is_zero() else proof.unit.sub(p.coefficients[0], budget).is_zero()
    if not proof.unit.sub(p.coefficients[p.degree()], budget).is_zero() or len(proof.factors) > p.degree():
        return False
    var total = 0
    for i in range(len(proof.factors)):
        budget.charge(1)
        ref q = proof.factors[i][0]
        var m = proof.factors[i][1]
        q.check_invariants(budget)
        var d = q.degree()
        comptime maximum = 2 if C == Real else 1
        if q.variable != p.variable or d < 1 or d > maximum or m <= 0 or m > (p.degree() - total) // d:
            return False
        total += d * m
        for coefficient_index in range(len(q.coefficients)):
            ref c = q.coefficients[coefficient_index]
            if not _roots._valid(_complex(c), budget):
                return False
        if not q.coefficients[d].sub(C.one(p.ring), budget).is_zero():
            return False
        comptime if C == Real:
            if d == 2:
                var twice = q.coefficients[0].add(q.coefficients[0], budget)
                var discriminant = q.coefficients[1].mul(q.coefficients[1], budget).sub(twice.add(twice, budget), budget)
                if _roots._compare_real(rebind[Real](discriminant), Real.zero(NoContext()), budget) >= 0:
                    return False
        if i and _compare(proof.factors[i - 1][0], q, budget) >= 0:
            return False
    if total != p.degree():
        return False
    _allocation[C](1, budget)
    _charge_size(proof.unit.size(), budget)
    var coefficients: List[C] = [proof.unit.copy()]
    var product = DensePoly[C].from_coefficients(p.variable, p.ring, coefficients^, budget)
    for entry in proof.factors:
        budget.charge(1)
        product = product.mul(entry[0].pow(UInt64(entry[1]), budget), budget)
    if len(product.coefficients) != len(p.coefficients):
        return False
    for i in range(len(p.coefficients)):
        budget.charge(1)
        if not product.coefficients[i].sub(p.coefficients[i], budget).is_zero():
            return False
    return True


def _factor[C: CoefficientRing](p: DensePoly[C], mut budget: Budget) raises PolyError -> Factorization[C]:
    comptime assert C == Real or C == Complex, "algebraic factorization requires AA or QQbar"
    p.check_invariants(budget)
    var unit = C.zero(p.ring)
    if not p.is_zero():
        _charge_size(p.coefficients[p.degree()].size(), budget)
        unit = p.coefficients[p.degree()].copy()
    var factors = List[Tuple[DensePoly[C], Int]]()
    if p.degree() > 0:
        _allocation[Complex](len(p.coefficients), budget)
        var coefficients = List[Complex](capacity=len(p.coefficients))
        for coefficient_index in range(len(p.coefficients)):
            ref c = p.coefficients[coefficient_index]
            _charge_size(c.size(), budget)
            coefficients.append(_complex(c))
        var roots = _roots._solve(DensePoly[Complex].from_coefficients(0, NoContext(), coefficients^, budget), budget)
        _allocation[Tuple[DensePoly[C], Int]](len(roots), budget)
        factors.reserve(len(roots))
        for i in range(len(roots)):
            budget.charge(1)
            var root = roots[i]
            var entries: List[C]
            comptime if C == Complex:
                _allocation[C](2, budget)
                entries = [rebind_var[C](root.neg(budget)), C.one(p.ring)]
            else:
                var re = root.real()
                var im = root.imag()
                var sign = _roots._compare_real(im, Real.zero(NoContext()), budget)
                if sign > 0:
                    continue
                if sign == 0:
                    _allocation[C](2, budget)
                    entries = [rebind_var[C](re.neg(budget)), C.one(p.ring)]
                else:
                    # Check the partner explicitly: roots with the same real
                    # part can contain several different conjugate pairs.
                    var partner = Complex(re.copy(), im.neg(budget))
                    var found = False
                    for j in range(i + 1, len(roots)):
                        budget.charge(1)
                        if _roots._compare(roots[j], partner, budget) == 0:
                            found = roots.multiplicity(j) == roots.multiplicity(i)
                            break
                    if not found:
                        raise PolyError(PolyErrorKind.DOMAIN, "real polynomial is missing a conjugate root with matching multiplicity")
                    _allocation[C](3, budget)
                    entries = [rebind_var[C](re.mul(re, budget).add(im.mul(im, budget), budget)), rebind_var[C](re.add(re, budget).neg(budget)), C.one(p.ring)]
            factors.append((DensePoly[C].from_coefficients(p.variable, p.ring, entries^, budget), roots.multiplicity(i)))
    _sort(factors, budget)
    var proof = Factorization[C](unit^, factors^)
    if not _check(p, proof, budget):
        raise PolyError(PolyErrorKind.DOMAIN, "algebraic factorization failed irreducibility or complete reconstruction")
    return proof^


def _axis[C: CoefficientRing](p: SparsePoly[C], mut budget: Budget) raises PolyError -> Tuple[Bool, Int64]:
    var axis = Int64(0)
    var found = False
    for term in p.terms:
        budget.charge(1)
        if len(term.monomial.pairs) > 1:
            return False, axis
        if len(term.monomial.pairs):
            var current = term.monomial.pairs[0][0]
            if found and current != axis:
                return False, axis
            axis = current
            found = True
    return True, axis


def _factor_sparse[C: CoefficientRing](p: SparsePoly[C], mut budget: Budget) raises PolyError -> MultivariateFactorization[C]:
    var active = _axis(p, budget)
    if not active[0]:
        var complete = _factor_sparse_evidence(p, budget)
        var result = MultivariateFactorization[C](C.zero(p.ring), List[Tuple[SparsePoly[C], Int]]())
        swap(result, complete[0])
        return result^
    var univariate = _factor(DensePoly[C].from_sparse(p, budget, variable=active[1]), budget)
    _allocation[Tuple[SparsePoly[C], Int]](len(univariate.factors), budget)
    var factors = List[Tuple[SparsePoly[C], Int]](capacity=len(univariate.factors))
    for entry in univariate.factors:
        budget.charge(1)
        factors.append((entry[0].to_sparse(p.order, budget), entry[1]))
    _charge_size(univariate.unit.size(), budget)
    return MultivariateFactorization[C](univariate.unit.copy(), factors^)


def _check_sparse[C: CoefficientRing](p: SparsePoly[C], proof: MultivariateFactorization[C], mut budget: Budget) raises PolyError -> Bool:
    p.check_invariants(budget)
    var active = _axis(p, budget)
    if not active[0]:
        return False
    var source = DensePoly[C].from_sparse(p, budget, variable=active[1])
    if len(proof.factors) > max(0, source.degree()):
        return False
    _allocation[Tuple[DensePoly[C], Int]](len(proof.factors), budget)
    var factors = List[Tuple[DensePoly[C], Int]](capacity=len(proof.factors))
    for entry in proof.factors:
        budget.charge(1)
        entry[0].check_invariants(budget)
        var variable = _axis(entry[0], budget)
        if not variable[0] or (len(entry[0].terms) and variable[1] != active[1]):
            return False
        # Reject impossible factor degrees before any dense proof allocation.
        comptime maximum = 2 if C == Real else 1
        for term in entry[0].terms:
            budget.charge(1)
            if len(term.monomial.pairs) and term.monomial.pairs[0][1] > UInt32(maximum):
                return False
        factors.append((DensePoly[C].from_sparse(entry[0], budget, variable=active[1]), entry[1]))
    _charge_size(proof.unit.size(), budget)
    return _check(source, Factorization[C](proof.unit.copy(), factors^), budget)


def _factor_sparse_evidence[C: CoefficientRing](p: SparsePoly[C], mut budget: Budget) raises PolyError -> Tuple[MultivariateFactorization[C], List[Optional[Factorization[C]]]]:
    if _axis(p, budget)[0]:
        return _factor_sparse(p, budget), List[Optional[Factorization[C]]]()
    from . import _algebraic_multivariate
    comptime if C == Real:
        return rebind_var[Tuple[MultivariateFactorization[C], List[Optional[Factorization[C]]]]](_algebraic_multivariate._factor(rebind[SparsePoly[Real]](p), budget))
    else:
        return rebind_var[Tuple[MultivariateFactorization[C], List[Optional[Factorization[C]]]]](_algebraic_multivariate._factor(rebind[SparsePoly[Complex]](p), budget))


def _check_sparse_evidence[C: CoefficientRing](p: SparsePoly[C], proof: MultivariateFactorization[C], images: List[Optional[Factorization[C]]], mut budget: Budget) raises PolyError -> Bool:
    if _axis(p, budget)[0]:
        return len(images) == 0 and _check_sparse(p, proof, budget)
    from . import _algebraic_multivariate
    comptime if C == Real:
        return _algebraic_multivariate._check(rebind[SparsePoly[Real]](p), rebind[MultivariateFactorization[Real]](proof), rebind[List[Optional[Factorization[Real]]]](images), budget)
    else:
        return _algebraic_multivariate._check(rebind[SparsePoly[Complex]](p), rebind[MultivariateFactorization[Complex]](proof), rebind[List[Optional[Factorization[Complex]]]](images), budget)
