"""Private projection values with no determinant or coprimality histories.

The producer retains normalized factors only. Final verification recomputes
the projection obligations, checks factor coverage by exact division, and
establishes squarefreeness/primitive coprimality by exact gcds. Brown producer
heuristics are untrusted; omitted contributions require checked ideal analysis.
"""
from .coefficients import QQ, NoContext
from .sparse import SparsePoly
from .monomial import MonomialOrder
from .errors import Budget, PolyError, PolyErrorKind, _allocation, _sum
from .evaluate import derivative
from .division import divide_exact
from .gcd import gcd_value
from .resultant import resultant_value
from .multivariate_factor import factor_multivariate
from ._views import _coefficient, _lc, _is_unit
from ._cad_brown import _finite_eliminants, _BROWN_ANALYSIS_LIMIT
from .cad_projection import _level
from .certificates import _equal


@fieldwise_init
struct _FinalProjection(Copyable, Movable):
    var variables: List[Int64]
    var brown: Bool
    var levels: List[List[SparsePoly[QQ]]]


def _projection_variables(variables: List[Int64], mut budget: Budget) raises PolyError -> Bool:
    for i in range(len(variables)):
        budget.charge(1)
        if variables[i] < 0:
            return False
        for j in range(i):
            budget.charge(1)
            if variables[i] == variables[j]:
                return False
    return True


def _collect_values(f: SparsePoly[QQ], maximum: Int, mut result: _FinalProjection, mut budget: Budget) raises PolyError:
    # This factorization consists only of final factors/multiplicities; it
    # carries no search, determinant, or gcd/cofactor transcripts.
    var factors = factor_multivariate(f.reorder(MonomialOrder(), budget), budget)
    for entry in factors.factors:
        var level = _level(entry[0], result.variables, budget)
        if level < 0 or level >= maximum:
            raise PolyError(PolyErrorKind.DOMAIN, "invalid projection factor level")
        var found = False
        for known in result.levels[level]:
            if _equal(entry[0], known, budget):
                found = True
                break
        if not found:
            _allocation[SparsePoly[QQ]](1, budget)
            result.levels[level].append(entry[0]._copy(budget))


def _prefer_values(tc: SparsePoly[QQ], eliminants: List[SparsePoly[QQ]], projection: _FinalProjection, lower: List[Int64], mut budget: Budget) raises PolyError -> Bool:
    var available = budget.remaining()
    var allowance = min(_BROWN_ANALYSIS_LIMIT, available)
    var local = Budget(allowance)
    try:
        var counts = List[Int]()
        for choice in range(2):
            var inputs = List[SparsePoly[QQ]]()
            for level in range(len(lower)):
                for f in projection.levels[level]:
                    _allocation[SparsePoly[QQ]](1, local)
                    inputs.append(f._copy(local))
            if choice == 0:
                _allocation[SparsePoly[QQ]](1, local)
                inputs.append(tc._copy(local))
            else:
                for f in eliminants:
                    _allocation[SparsePoly[QQ]](1, local)
                    inputs.append(f._copy(local))
            var probe = _project_values(inputs, lower, local, brown=False)
            var count = 0
            for level in probe.levels:
                local.charge(1)
                count = _sum(count, len(level))
            _allocation[Int](1, local)
            counts.append(count)
        budget.charge(allowance - local.remaining())
        return counts[1] < counts[0]
    except e:
        budget.charge(allowance - local.remaining())
        if e.kind == PolyErrorKind.LIMIT and _BROWN_ANALYSIS_LIMIT < available:
            return False
        raise e^


def _project_values(inputs: List[SparsePoly[QQ]], variables: List[Int64], mut budget: Budget, *, brown: Bool = True) raises PolyError -> _FinalProjection:
    if not _projection_variables(variables, budget):
        raise PolyError(PolyErrorKind.DOMAIN, "invalid projection variable order")
    _allocation[Int64](len(variables), budget)
    _allocation[List[SparsePoly[QQ]]](len(variables), budget)
    var result = _FinalProjection(variables.copy(), brown, List[List[SparsePoly[QQ]]](length=len(variables), fill=List[SparsePoly[QQ]]()))
    for f in inputs:
        _ = _level(f, variables, budget)
        _collect_values(f, len(variables), result, budget)
    for level in range(len(variables) - 1, 0, -1):
        var variable = variables[level]
        var size = len(result.levels[level])
        for i in range(size):
            var f = result.levels[level][i]._copy(budget)
            _collect_values(_lc(f, variable, budget), level, result, budget)
            _collect_values(resultant_value(f, derivative(f, variable, budget), variable, budget), level, result, budget)
            for j in range(i):
                var g = result.levels[level][j]._copy(budget)
                _collect_values(resultant_value(f, g, variable, budget), level, result, budget)
        _allocation[Int64](level, budget)
        var lower = List[Int64]()
        for k in range(level):
            lower.append(variables[k])
        for i in range(size):
            var f = result.levels[level][i]._copy(budget)
            var tc = _coefficient(f, variable, 0, budget)
            if brown:
                var replacement = _finite_eliminants(f, variable, lower, budget, checking=False)
                if replacement:
                    if len(replacement.value()) == 0:
                        continue
                    if _prefer_values(tc, replacement.value(), result, lower, budget):
                        for eliminant in replacement.value():
                            _collect_values(eliminant, level, result, budget)
                        continue
            _collect_values(tc, level, result, budget)
    return result^


def _covered(f: SparsePoly[QQ], maximum: Int, projection: _FinalProjection, mut budget: Budget) raises PolyError -> Bool:
    var rest = f.reorder(MonomialOrder(), budget)
    if rest.is_zero():
        return True
    for level in range(maximum):
        for factor in projection.levels[level]:
            while not _is_unit(rest, budget):
                var divided = divide_exact(rest, factor, budget)
                if not divided:
                    break
                rest = divided.take()
    return _is_unit(rest, budget)


def _check_projection_values(inputs: List[SparsePoly[QQ]], projection: _FinalProjection, mut budget: Budget) raises PolyError -> Bool:
    try:
        if not _projection_variables(projection.variables, budget) or len(projection.levels) != len(projection.variables):
            return False
        for level in range(len(projection.levels)):
            var variable = projection.variables[level]
            for i in range(len(projection.levels[level])):
                ref f = projection.levels[level][i]
                f.check_invariants(budget)
                if f.order != MonomialOrder() or _level(f, projection.variables, budget) != level:
                    return False
                if f.terms[len(f.terms) - 1].coefficient != QQ.one(NoContext()):
                    return False
                if not _is_unit(gcd_value(f, derivative(f, variable, budget), budget), budget):
                    return False
                for j in range(i):
                    if not _is_unit(gcd_value(f, projection.levels[level][j], budget), budget):
                        return False
        for f in inputs:
            _ = _level(f, projection.variables, budget)
            if not _covered(f, len(projection.levels), projection, budget):
                return False
        for level in range(len(projection.levels) - 1, 0, -1):
            var variable = projection.variables[level]
            _allocation[Int64](level, budget)
            var lower = List[Int64]()
            for k in range(level):
                lower.append(projection.variables[k])
            for i in range(len(projection.levels[level])):
                ref f = projection.levels[level][i]
                if not _covered(_lc(f, variable, budget), level, projection, budget):
                    return False
                if not _covered(resultant_value(f, derivative(f, variable, budget), variable, budget), level, projection, budget):
                    return False
                for j in range(i):
                    if not _covered(resultant_value(f, projection.levels[level][j], variable, budget), level, projection, budget):
                        return False
                if _covered(_coefficient(f, variable, 0, budget), level, projection, budget):
                    continue
                if not projection.brown:
                    return False
                # Verify only mathematical coverage. The producer's choice
                # heuristic and local cost estimate are not trusted premises.
                var replacements = _finite_eliminants(f, variable, lower, budget)
                if not replacements:
                    return False
                for eliminant in replacements.value():
                    if not _covered(eliminant, level, projection, budget):
                        return False
        return True
    except e:
        if e.kind == PolyErrorKind.DOMAIN:
            return False
        raise e^
