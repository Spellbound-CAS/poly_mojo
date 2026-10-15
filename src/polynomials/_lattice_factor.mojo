"""Van Hoeij trace-knapsack recombination with exact lattice reduction.

Fresh implementation from Mark van Hoeij (2002), sections 2.1 and 2.2,
DOI 10.1006/jnth.2001.2763. Scaled traces handle nonmonic inputs.
"""
from apn_mojo import Integer
from .coefficients import ZZ, QQ, ModInt, NoContext, _work
from .dense import DensePoly
from .errors import Budget, PolyError, PolyErrorKind, _allocation, _product, _sum
from ._coefficient_ops import _exact, _copy as _zcopy
from ._dense_ops import _copy, _constant
from ._hensel import _HenselState, _hensel_start, _hensel_extend, _truncate, _primitive, _trial_divide
from .primes import _symmetric
from ._lattice import _lll, _q
from ._root_ops import _compare, _copy as _copy_q


def _greater(a: ZZ, b: ZZ, mut budget: Budget) raises PolyError -> Bool:
    _work(a.size(), b.size(), False, budget)
    return a.value > b.value


def _traces(f: DensePoly[ZZ], count: Int, modulus: ZZ, mut budget: Budget) raises PolyError -> List[ZZ]:
    # Newton sums of the roots of the monic modular factor.
    _allocation[ZZ](_sum(count, 1), budget)
    var sums: List[ZZ] = [ZZ.zero(NoContext())]
    for k in range(1, _sum(count, 1)):
        budget.charge(1)
        var total = ZZ.zero(NoContext())
        for j in range(1, min(k, _sum(f.degree(), 1))):
            budget.charge(1)
            total = total.add(f.coefficients[f.degree() - j].mul(sums[k - j], budget), budget)
            total = ZZ(_symmetric(total.value, modulus.value, budget))
        if k <= f.degree():
            total = total.add(ZZ(Integer(k)).mul(f.coefficients[f.degree() - k], budget), budget)
        sums.append(ZZ(_symmetric(total.neg(budget).value, modulus.value, budget)))
    return sums^


def _projection_rref(basis: List[List[ZZ]], rows: Int, columns: Int, mut budget: Budget) raises PolyError -> List[List[QQ]]:
    _allocation[List[QQ]](rows, budget)
    _allocation[QQ](_product(rows, columns), budget)
    var result = List[List[QQ]]()
    for i in range(rows):
        budget.charge(1)
        var row = List[QQ]()
        for j in range(columns):
            budget.charge(1)
            row.append(_q(basis[i][j], budget))
        result.append(row^)
    var rank = 0
    for column in range(columns):
        budget.charge(1)
        var pivot = rank
        while pivot < rows and result[pivot][column].is_zero():
            budget.charge(1)
            pivot += 1
        if pivot == rows:
            continue
        if pivot != rank:
            var temporary = List[QQ]()
            swap(temporary, result[pivot])
            swap(temporary, result[rank])
            swap(temporary, result[pivot])
        var inverse = result[rank][column].inverse(budget)
        for j in range(column, columns):
            budget.charge(1)
            result[rank][j] = result[rank][j].mul(inverse, budget)
        for i in range(rows):
            budget.charge(1)
            if i == rank or result[i][column].is_zero():
                continue
            var multiplier = _copy_q(result[i][column], budget)
            for j in range(column, columns):
                budget.charge(1)
                result[i][j] = result[i][j].sub(multiplier.mul(result[rank][j], budget), budget)
        rank += 1
        if rank == rows:
            break
    while len(result) > rank:
        budget.charge(1)
        _ = result.pop()
    return result^


def _partition(rows: List[List[QQ]], count: Int, mut budget: Budget) raises PolyError -> Bool:
    var one = QQ.one(NoContext())
    for j in range(count):
        budget.charge(1)
        var occurrences = 0
        for row in rows:
            budget.charge(1)
            if row[j].is_zero():
                continue
            if _compare(row[j], one, budget) != 0:
                return False
            occurrences += 1
        if occurrences != 1:
            return False
    return True


def _van_hoeij(f: DensePoly[ZZ], images: List[DensePoly[ModInt]], reconstruction_bound: ZZ, mut budget: Budget) raises PolyError -> List[DensePoly[ZZ]]:
    var state = _hensel_start(images, budget)
    return _van_hoeij_lifted(f, images, reconstruction_bound, state, budget)


def _van_hoeij_lifted(f: DensePoly[ZZ], images: List[DensePoly[ModInt]], reconstruction_bound: ZZ, mut state: _HenselState, mut budget: Budget) raises PolyError -> List[DensePoly[ZZ]]:
    var count = len(images)
    if count < 2 or f.degree() < 2:
        raise PolyError(PolyErrorKind.DOMAIN, "lattice recombination requires at least two modular factors")
    var p = ZZ(Integer(Int(images[0].ring.value)))
    var root_bound = ZZ.zero(NoContext())
    for coefficient_index in range(len(f.coefficients)):
        ref coefficient = f.coefficients[coefficient_index]
        budget.charge(1)
        var magnitude = coefficient.neg(budget) if coefficient.value.sign() < 0 else _zcopy(coefficient, budget)
        if _greater(magnitude, root_bound, budget):
            root_bound = magnitude^
    ref lc = f.coefficients[f.degree()]
    var leading = lc.neg(budget) if lc.value.sign() < 0 else _zcopy(lc, budget)
    root_bound = root_bound.add(leading, budget)
    var traces_count = min(f.degree(), 2)
    var padding = _zcopy(p, budget)
    while True:
        budget.charge(1)
        # For a rational factor, lc(f)^k * Tr_k is an integer with
        # absolute value at most degree(f)*(|lc(f)|+max |f_i|)^k.
        _allocation[ZZ](traces_count, budget)
        var cuts = List[ZZ]()
        var power = ZZ.one(NoContext())
        var largest_cut = ZZ.one(NoContext())
        for _ in range(traces_count):
            budget.charge(1)
            power = power.mul(root_bound, budget)
            var bound = power.mul(ZZ(Integer(_product(2, f.degree()))), budget)
            var cut = ZZ.one(NoContext())
            while not _greater(cut, bound, budget):
                budget.charge(1)
                cut = cut.mul(p, budget)
            largest_cut = _zcopy(cut, budget)
            cuts.append(cut^)
        var target = largest_cut.mul(padding, budget)
        if _greater(reconstruction_bound, target, budget):
            target = _zcopy(reconstruction_bound, budget)
        _hensel_extend(f, images, target, state, budget)
        ref lifted = state.lifted
        ref modulus = state.modulus
        var dimension = _sum(count, traces_count)
        _allocation[List[ZZ]](dimension, budget)
        _allocation[ZZ](_product(dimension, dimension), budget)
        var lattice = List[List[ZZ]](length=dimension, fill=List[ZZ](length=dimension, fill=ZZ.zero(NoContext())))
        for i in range(count):
            budget.charge(1)
            lattice[i][i] = ZZ.one(NoContext())
            var sums = _traces(lifted[i], traces_count, modulus, budget)
            var scale = ZZ.one(NoContext())
            for j in range(traces_count):
                budget.charge(1)
                scale = ZZ(_symmetric(scale.mul(leading, budget).value, modulus.value, budget))
                var trace = ZZ(_symmetric(sums[j + 1].mul(scale, budget).value, modulus.value, budget))
                var low = ZZ(_symmetric(trace.value, cuts[j].value, budget))
                lattice[i][count + j] = _exact(trace.sub(low, budget), cuts[j], budget)
        for j in range(traces_count):
            budget.charge(1)
            lattice[count + j][count + j] = _exact(modulus, cuts[j], budget)
        var reduced = _lll[False](lattice, budget)
        # Each true factor's indicator has norm squared <= n+s*n^2/4.
        # The looser integer bound avoids rounding. Dropping only a suffix
        # of larger Gram-Schmidt norms retains every such short vector.
        var bound_squared = _q(ZZ(Integer(_sum(count, _product(traces_count, _product(count, count))))), budget)
        var retained = dimension
        while retained and _compare(reduced.norms[retained - 1], bound_squared, budget) > 0:
            budget.charge(1)
            retained -= 1
        var partition = _projection_rref(reduced.basis, retained, count, budget)
        if _partition(partition, count, budget):
            var remaining = _copy(f, budget)
            var factors = List[DensePoly[ZZ]]()
            var valid = True
            for row in partition:
                budget.charge(1)
                var candidate = _constant(f, leading, budget)
                for j in range(count):
                    budget.charge(1)
                    if not row[j].is_zero():
                        candidate = _truncate(candidate.mul(lifted[j], budget), modulus, budget)
                candidate = _primitive(candidate, budget)
                var quotient = _trial_divide(remaining, candidate, budget)
                if candidate.degree() < 1 or not quotient:
                    valid = False
                    break
                remaining = quotient.take()
                _allocation[DensePoly[ZZ]](1, budget)
                factors.append(candidate^)
            if valid and remaining.degree() == 0:
                return factors^
        if traces_count < f.degree():
            traces_count = min(f.degree(), _product(2, traces_count))
        else:
            padding = padding.mul(p, budget)
