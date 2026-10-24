"""Private Brown–McCallum eligibility tests for Lazard projection.

An empty coefficient variety permits omission of the trailing coefficient.
The unit-ideal claim is checked through exact generator identities, not just
Buchberger's criterion (which alone accepts the spurious basis {1}). Analysis
is capped and charged to the caller; no witness is stored in a CAD result.

Reference: Davenport et al., arXiv:2302.05813, section 2.7, Theorem 2.1,
Hypothesis 1 and Definition 11. Coordinate eliminants for finite nullification
cut every exceptional coordinate into sections, so their common zeros are
zero-dimensional base cells without a new stored point-witness type.
"""
from .coefficients import QQ
from .sparse import SparsePoly
from .monomial import MonomialOrder
from .errors import Budget, PolyError, PolyErrorKind, _allocation
from ._views import _lc, _is_unit
from ._grouping import _coefficient_groups
from .groebner import GroebnerBasis, groebner, check_groebner_basis
from .certificates import MembershipCertificate, check_membership
from .gcd import gcd_value

comptime _BROWN_ANALYSIS_LIMIT = UInt64(1000000)


def _check_coefficient_basis(generators: List[SparsePoly[QQ]], basis: GroebnerBasis[QQ], mut budget: Budget) raises PolyError -> Bool:
    """Establish both ideal containments, not merely input reduction to zero."""
    if not basis.cofactors or len(basis.cofactors.value()) != len(basis.basis):
        return False
    for i in range(len(basis.basis)):
        budget.charge(1)
        if basis.basis[i].is_zero():
            return False
        var row = List[SparsePoly[QQ]]()
        _allocation[SparsePoly[QQ]](len(basis.cofactors.value()[i]), budget)
        for c in basis.cofactors.value()[i]:
            row.append(c.reorder(basis.order, budget))
        if not check_membership(basis.basis[i], generators, MembershipCertificate[QQ](row^), budget):
            return False
    return check_groebner_basis(generators, basis, budget)


def _coefficient_basis(coefficients: List[SparsePoly[QQ]], variables: List[Int64], mut budget: Budget, cap: UInt64, *, checking: Bool = True) raises PolyError -> Optional[GroebnerBasis[QQ]]:
    var available = budget.remaining()
    var allowance = min(cap, available)
    var local = Budget(allowance)
    try:
        _allocation[Int64](len(variables), local)
        var order = MonomialOrder(MonomialOrder.LEX, variables.copy(), local)
        var generators = List[SparsePoly[QQ]]()
        _allocation[SparsePoly[QQ]](len(coefficients), local)
        for c in coefficients:
            generators.append(c.reorder(order, local))
        var basis = groebner(generators, order, local, track_cofactors=checking, certify=checking)
        if checking and not _check_coefficient_basis(generators, basis, local):
            raise PolyError(PolyErrorKind.DOMAIN, "coefficient ideal failed independent verification")
        budget.charge(allowance - local.remaining())
        return basis^
    except e:
        budget.charge(allowance - local.remaining())
        if e.kind == PolyErrorKind.LIMIT and cap < available:
            return None
        raise e^


def _omit_trailing(f: SparsePoly[QQ], variable: Int64, lower_variables: List[Int64], mut budget: Budget, *, cap: UInt64 = _BROWN_ANALYSIS_LIMIT) raises PolyError -> Bool:
    f.check_invariants(budget)
    if _is_unit(_lc(f, variable, budget), budget):
        return True
    # Group only nonzero coefficients. A sparse polynomial's degree may be
    # enormous; eligibility must not allocate a dense degree-sized array.
    var grouped = _coefficient_groups(f, variable, False, budget)
    var coefficients = List[SparsePoly[QQ]]()
    swap(coefficients, grouped[1])
    # Any nonzero constant coefficient proves that simultaneous vanishing is
    # impossible. This includes the constant-leading-coefficient fast path.
    for c in coefficients:
        if _is_unit(c, budget):
            return True
    var basis = _coefficient_basis(coefficients, lower_variables, budget, cap)
    if not basis:
        return False
    for b in basis.value().basis:
        if _is_unit(b, budget):
            return True
    return False


def _zero_dimensional(basis: GroebnerBasis[QQ], variables: List[Int64], mut budget: Budget) raises PolyError -> Bool:
    # Every ambient coordinate needs a pure power in the leading ideal.
    # A variable absent from all coefficients is still an unconstrained axis.
    for variable in variables:
        var found = False
        for f in basis.basis:
            budget.charge(1)
            var term = f.leading(budget)
            ref leading = term.monomial
            if len(leading.pairs) == 1 and leading.pairs[0][0] == variable:
                found = True
                break
        if not found:
            return False
    return True


def _coordinate_eliminant(basis: GroebnerBasis[QQ], variable: Int64, mut budget: Budget) raises PolyError -> Optional[SparsePoly[QQ]]:
    for f in basis.basis:
        var used = f.variables(budget)
        if len(used) == 1 and used[0] == variable:
            # The checked generator identity makes this an ideal member;
            # monicity prevents a leading coefficient from vanishing later.
            var normalized = f.scale(f.leading(budget).coefficient.inverse(budget), budget)
            return normalized.reorder(MonomialOrder(), budget)
    return None


def _finite_eliminants(f: SparsePoly[QQ], variable: Int64, lower_variables: List[Int64], mut budget: Budget, *, cap: UInt64 = _BROWN_ANALYSIS_LIMIT, checking: Bool = True) raises PolyError -> Optional[List[SparsePoly[QQ]]]:
    """None keeps tc; empty means no nullification; otherwise one eliminant/axis."""
    f.check_invariants(budget)
    if _is_unit(_lc(f, variable, budget), budget):
        return List[SparsePoly[QQ]]()
    var grouped = _coefficient_groups(f, variable, False, budget)
    var coefficients = List[SparsePoly[QQ]]()
    swap(coefficients, grouped[1])
    for c in coefficients:
        if _is_unit(c, budget):
            return List[SparsePoly[QQ]]()
    if len(coefficients) == 0:
        return None
    var available = budget.remaining()
    var allowance = min(cap, available)
    var local = Budget(allowance)
    try:
        if len(lower_variables) == 1:
            # In Q[x], the coefficient ideal is generated by their exact gcd.
            # Recompute it independently in the checker. Multivariate
            # Buchberger/cofactor construction can swell enormously here even
            # when the gcd is 1; no such evidence is needed in a PID.
            var generator = coefficients[0]._copy(local)
            for i in range(1, len(coefficients)):
                generator = gcd_value(generator, coefficients[i], local)
                if _is_unit(generator, local):
                    break
            var result = List[SparsePoly[QQ]]()
            if not _is_unit(generator, local):
                _allocation[SparsePoly[QQ]](1, local)
                result.append(generator.scale(generator.leading(local).coefficient.inverse(local), local).reorder(MonomialOrder(), local))
            budget.charge(allowance - local.remaining())
            return result^
        var checked = _coefficient_basis(coefficients, lower_variables, local, UInt64.MAX, checking=checking)
        if not checked:
            raise PolyError(PolyErrorKind.DOMAIN, "missing coefficient basis")
        var basis = checked.take()
        for b in basis.basis:
            if _is_unit(b, local):
                budget.charge(allowance - local.remaining())
                return List[SparsePoly[QQ]]()
        if not _zero_dimensional(basis, lower_variables, local):
            budget.charge(allowance - local.remaining())
            return None
        _allocation[SparsePoly[QQ]](len(lower_variables), local)
        var result = List[SparsePoly[QQ]]()
        for v in lower_variables:
            var eliminant = _coordinate_eliminant(basis, v, local)
            if not eliminant:
                # Lex elimination retains the last variable. Rotate the order
                # only when the first basis has no univariate member for v.
                _allocation[Int64](len(lower_variables), local)
                var rotated = List[Int64]()
                for w in lower_variables:
                    if w != v:
                        rotated.append(w)
                rotated.append(v)
                var other = _coefficient_basis(coefficients, rotated, local, UInt64.MAX, checking=checking)
                if not other:
                    raise PolyError(PolyErrorKind.DOMAIN, "missing elimination basis")
                eliminant = _coordinate_eliminant(other.value(), v, local)
            if not eliminant:
                raise PolyError(PolyErrorKind.DOMAIN, "zero-dimensional ideal lacks a coordinate eliminant")
            result.append(eliminant.take())
        budget.charge(allowance - local.remaining())
        return result^
    except e:
        budget.charge(allowance - local.remaining())
        if e.kind == PolyErrorKind.LIMIT and cap < available:
            return None
        raise e^
