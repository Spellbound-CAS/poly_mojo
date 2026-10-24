"""Bivariate rational resultants with their cofactors by evaluation and interpolation modulo primes.

For f, g in Q[y][x] of positive degrees n = deg_x f and m = deg_x g with a
nonzero resultant r, there are unique u, v with u f + v g = r, deg_x u < m
and deg_x v < n: the subresultant sequence of `resultant` returns exactly
them. r itself comes exactly from the integer evaluation-interpolation
resultant after clearing denominators. Modulo a prime and at a point y0 where
both leading coefficients and r(y0) are nonzero, u(y0) and v(y0) are r(y0)
times the modular unit Bezout cofactors of f(x, y0) and g(x, y0). Each
coefficient has degree at most n deg_y g + m deg_y f in y, so that many
points plus one interpolate it; primes are combined by the Chinese remainder
theorem and coefficients recovered as fractions. A candidate is returned only
once u f + v g = r holds exactly, which by uniqueness makes it the pair
`resultant` computes; otherwise the caller falls back to it.
"""
from apn_mojo import Integer, Rational, isqrt, lcm as integer_lcm
from .errors import Budget, PolyError, PolyErrorKind, _allocation, _sum, _product
from .coefficients import ZZ, QQ, NoContext
from .sparse import SparsePoly, Term
from .monomial import Monomial, MonomialOrder
from .certificates import ResultantResult
from .primes import _residue
from ._rational_groebner import _next_prime, _reconstruct
from ._multivariate_resultant import _integer_resultant
from ._unit_bezout import _inverse, _modular_bezout
from ._modular_simd import _Prime32, _inverse32, _mod_add, _mod_sub
from ._word_simd import _WORD_LANES

comptime _COFACTOR_PRIMES = 1024
"""Primes tried before the caller falls back to the subresultant sequence."""
comptime _COFACTOR_UNLUCKY = 4
"""Consecutive primes without enough evaluation points before giving up."""


def _common_denominator(p: SparsePoly[QQ], mut budget: Budget) raises PolyError -> Integer:
    var common = Integer(1)
    for term in p.terms:
        budget.charge(1)
        try:
            common = integer_lcm(common, term.coefficient.value.denominator())
        except e:
            raise PolyError(PolyErrorKind.CAPACITY, String(e))
    return common^


def _integral(p: SparsePoly[QQ], common: Integer, mut budget: Budget) raises PolyError -> SparsePoly[ZZ]:
    """common times p, with integer coefficients."""
    _allocation[Term[ZZ]](len(p.terms), budget)
    var terms = List[Term[ZZ]](capacity=len(p.terms))
    for term in p.terms:
        budget.charge(1)
        try:
            terms.append(Term[ZZ](term.monomial.copy(), ZZ((term.coefficient.value * Rational(common)).to_integer_exact())))
        except e:
            raise PolyError(PolyErrorKind.CAPACITY, String(e))
    var result = SparsePoly[ZZ].zero(p.order, NoContext(), budget)
    result.terms = terms^
    return result^


def _rational(p: SparsePoly[ZZ], scale: Integer, order: MonomialOrder, mut budget: Budget) raises PolyError -> SparsePoly[QQ]:
    """p divided by scale, over QQ."""
    _allocation[Term[QQ]](len(p.terms), budget)
    var terms = List[Term[QQ]](capacity=len(p.terms))
    for term in p.terms:
        budget.charge(1)
        try:
            terms.append(Term[QQ](term.monomial.copy(), QQ(Rational(term.coefficient.value.copy(), scale.copy()))))
        except e:
            raise PolyError(PolyErrorKind.CAPACITY, String(e))
    return SparsePoly[QQ].from_terms(terms^, order, NoContext(), budget)


def _grid(p: SparsePoly[QQ], variable: Int64, other: Int64, rows: Int, columns: Int) -> List[List[Rational]]:
    """Coefficients of x^i y^j at [i][j]."""
    var grid = List[List[Rational]](capacity=rows)
    for _ in range(rows):
        grid.append(List[Rational](length=columns, fill=Rational(0)))
    for term in p.terms:
        var i = 0
        var j = 0
        for pair in term.monomial.pairs:
            if pair[0] == variable:
                i = Int(UInt64(pair[1]))
            elif pair[0] == other:
                j = Int(UInt64(pair[1]))
        grid[i][j] = term.coefficient.value.copy()
    return grid^


def _grid_image(grid: List[List[Rational]], p: UInt64, mut budget: Budget) raises PolyError -> Optional[List[List[UInt64]]]:
    """The grid modulo p, or None when p divides a denominator."""
    var result = List[List[UInt64]](capacity=len(grid))
    for row in grid:
        var image = List[UInt64](capacity=len(row))
        for value in row:
            budget.charge(1)
            if value.numerator() == Integer(0):
                image.append(0)
                continue
            var d = _residue(value.denominator(), p, budget)
            if d == 0:
                return None
            image.append(UInt64((UInt128(_residue(value.numerator(), p, budget)) * UInt128(_inverse(d, p))) % UInt128(p)))
        result.append(image^)
    return result^


def _horner(values: List[UInt64], y: UInt64, p: UInt64) -> UInt64:
    if p < (UInt64(1) << 32):
        var field = _Prime32.of(p)
        var result = UInt64(0)
        for k in range(len(values) - 1, -1, -1):
            result = _mod_add[1](field.mul[1](result, y), values[k], p)
        return result
    var result = UInt64(0)
    for k in range(len(values) - 1, -1, -1):
        result = UInt64((UInt128(result) * UInt128(y) + UInt128(values[k])) % UInt128(p))
    return result


def _specialize(grid: List[List[UInt64]], y: UInt64, p: UInt64) -> List[UInt64]:
    var result = List[UInt64](capacity=len(grid))
    for row in grid:
        result.append(_horner(row, y, p))
    while len(result) and result[len(result) - 1] == 0:
        _ = result.pop()
    return result^


def _interpolate(points: List[UInt64], values: List[UInt64], inverses: List[List[UInt64]], p: UInt64) -> List[UInt64]:
    """Newton interpolation; inverses[level-1][k-level] = 1/(points[k]-points[k-level])."""
    if p < (UInt64(1) << 32):
        return _interpolate32(points, values, inverses, _Prime32.of(p))
    var n = len(points)
    var c = values.copy()
    for level in range(1, n):
        for k in range(n - 1, level - 1, -1):
            var difference = (c[k] + p - c[k - 1]) % p
            c[k] = UInt64((UInt128(difference) * UInt128(inverses[level - 1][k - level])) % UInt128(p))
    # Horner in the Newton basis: result = c[n-1]; result = result (y - points[k]) + c[k].
    var result = List[UInt64](length=n, fill=0)
    result[0] = c[n - 1]
    var length = 1
    for k in range(n - 2, -1, -1):
        # result *= (y - points[k])
        var shifted = List[UInt64](length=n, fill=0)
        for t in range(length):
            shifted[t + 1] = (shifted[t + 1] + result[t]) % p
            var product = UInt64((UInt128(result[t]) * UInt128(points[k])) % UInt128(p))
            shifted[t] = (shifted[t] + p - product) % p
        length += 1
        shifted[0] = (shifted[0] + c[k]) % p
        result = shifted^
    return result^


def _interpolate32[vectorized: Bool = True](points: List[UInt64], values: List[UInt64], inverses: List[List[UInt64]], field: _Prime32) -> List[UInt64]:
    """Contiguous divided differences and descending in-place Newton expansion."""
    var n = len(points)
    var c = values.copy()
    var pc = c.unsafe_ptr()
    var p = field.p
    comptime W = _WORD_LANES
    for level in range(1, n):
        var end = n
        var pi = inverses[level - 1].unsafe_ptr()
        comptime if vectorized:
            while end - level >= W:
                var start = end - W
                # Load both overlapping input blocks before writing. Moving
                # from high to low indices preserves every old neighbour.
                var difference = _mod_sub[W](pc.unsafe_offset(start).unsafe_load[width=W](), pc.unsafe_offset(start - 1).unsafe_load[width=W](), p)
                pc.unsafe_offset(start).unsafe_store[width=W](field.mul[W](difference, pi.unsafe_offset(start - level).unsafe_load[width=W]()))
                end -= W
        for k in range(end - 1, level - 1, -1):
            c[k] = field.mul[1](_mod_sub[1](c[k], c[k - 1], p), inverses[level - 1][k - level])
    var result = List[UInt64](length=n, fill=0)
    result[0] = c[n - 1]
    var data = result.unsafe_ptr()
    var length = 1
    for k in range(n - 2, -1, -1):
        result[length] = result[length - 1]
        var end = length
        comptime if vectorized:
            while end - 1 >= W:
                var start = end - W
                var old = data.unsafe_offset(start).unsafe_load[width=W]()
                var previous = data.unsafe_offset(start - 1).unsafe_load[width=W]()
                data.unsafe_offset(start).unsafe_store[width=W](_mod_sub[W](previous, field.mul[W](old, SIMD[DType.uint64, W](points[k])), p))
                end -= W
        for t in range(end - 1, 0, -1):
            result[t] = _mod_sub[1](result[t - 1], field.mul[1](result[t], points[k]), p)
        result[0] = _mod_sub[1](c[k], field.mul[1](result[0], points[k]), p)
        length += 1
    return result^


def _bivariate_resultant(f: SparsePoly[QQ], g: SparsePoly[QQ], variable: Int64, variables: List[Int64], mut budget: Budget) raises PolyError -> Optional[ResultantResult[QQ]]:
    """`resultant(f, g, variable)` for inputs in `variable` and one other, or None.

    None means only that this method did not find it, including a zero
    resultant; the caller decides by the subresultant sequence.
    """
    if len(variables) != 2:
        return None
    var other = variables[1] if variables[0] == variable else variables[0]
    var n = f.degree_in(variable, budget)
    var m = g.degree_in(variable, budget)
    if n <= 0 or m <= 0:
        return None
    var a = max(f.degree_in(other, budget), 0)
    var b = max(g.degree_in(other, budget), 0)
    var c = _common_denominator(f, budget)
    var d = _common_denominator(g, budget)
    var integer = _integer_resultant(_integral(f, c, budget), _integral(g, d, budget), variable, variables, budget)
    if not integer or integer.value().is_zero():
        return None
    var scale: Integer
    try:
        scale = (c ** m) * (d ** n)
    except e:
        raise PolyError(PolyErrorKind.CAPACITY, String(e))
    var r = _rational(integer.value(), scale, f.order, budget)
    var bound_y = _sum(_product(n, b), _product(m, a))
    var columns = bound_y + 1
    _allocation[Rational](_sum(_product(n + 1, a + 1), _product(m + 1, b + 1)), budget)
    var f_grid = _grid(f, variable, other, n + 1, a + 1)
    var g_grid = _grid(g, variable, other, m + 1, b + 1)
    var r_column = List[Rational](length=r.degree_in(other, budget) + 1, fill=Rational(0))
    for term in r.terms:
        var j = 0
        for pair in term.monomial.pairs:
            if pair[0] == other:
                j = Int(UInt64(pair[1]))
        r_column[j] = term.coefficient.value.copy()
    var r_grid = List[List[Rational]]()
    r_grid.append(r_column^)
    # The cofactors' degree in y is usually that of the resultant, far below
    # the Sylvester bound: try it first, refuted by two extra points per prime.
    var degree = r.degree_in(other, budget)
    if degree < bound_y:
        var guess = _cofactors(f, g, r, f_grid, g_grid, r_grid, variable, other, n, m, degree + 1, True, budget)
        if guess:
            return ResultantResult[QQ](r^, guess.take())
    var full = _cofactors(f, g, r, f_grid, g_grid, r_grid, variable, other, n, m, columns, False, budget)
    if full:
        return ResultantResult[QQ](r^, full.take())
    return None


def _cofactors(f: SparsePoly[QQ], g: SparsePoly[QQ], r: SparsePoly[QQ], f_grid: List[List[Rational]], g_grid: List[List[Rational]], r_grid: List[List[Rational]], variable: Int64, other: Int64, n: Int, m: Int, columns: Int, checked: Bool, mut budget: Budget) raises PolyError -> Optional[Tuple[SparsePoly[QQ], SparsePoly[QQ]]]:
    """u and v with `columns` coefficients in y each, or None; `checked` samples two more points to refute the degree."""
    # Residues of u's m and v's n coefficient polynomials, lifted by Garner's step.
    var series = m + n
    var count = _product(series, columns)
    _allocation[Integer](count, budget)
    var values = List[Integer](length=count, fill=Integer(0))
    var modulus = Integer(1)
    var primes = 0
    var unlucky = 0
    var attempt = 1
    var p = (UInt64(1) << 31) + 1
    for _ in range(_COFACTOR_PRIMES):
        p = _next_prime(p, budget)
        var fp = _grid_image(f_grid, p, budget)
        var gp = _grid_image(g_grid, p, budget)
        var rp = _grid_image(r_grid, p, budget)
        if not fp or not gp or not rp:
            continue
        var needed = columns + 2 if checked else columns
        var points = List[UInt64](capacity=needed)
        var samples = List[List[UInt64]](capacity=series)
        for _ in range(series):
            samples.append(List[UInt64](capacity=needed))
        var y = UInt64(0)
        var limit = _sum(_product(3, needed), 16)
        while len(points) < needed and Int(y) < limit:
            budget.charge(UInt64(series))
            var fx = _specialize(fp.value(), y, p)
            var gx = _specialize(gp.value(), y, p)
            var ry = _horner(rp.value()[0], y, p)
            if len(fx) != n + 1 or len(gx) != m + 1 or ry == 0:
                y += 1
                continue
            var pair = _modular_bezout(fx, gx, p)
            if not pair:
                y += 1
                continue
            ref up = pair.value()[0]
            ref vp = pair.value()[1]
            for k in range(m):
                var value = up[k] if k < len(up) else UInt64(0)
                samples[k].append(UInt64((UInt128(value) * UInt128(ry)) % UInt128(p)))
            for k in range(n):
                var value = vp[k] if k < len(vp) else UInt64(0)
                samples[m + k].append(UInt64((UInt128(value) * UInt128(ry)) % UInt128(p)))
            points.append(y)
            y += 1
        if len(points) < needed:
            unlucky += 1
            if unlucky >= _COFACTOR_UNLUCKY:
                return None
            continue
        unlucky = 0
        _allocation[UInt64](_product(columns, columns), budget)
        var inverses = List[List[UInt64]](capacity=columns)
        # Store each divided-difference level contiguously, shared by all
        # coefficient series. Independent denominators invert in SIMD batches.
        for level in range(1, columns):
            var row = List[UInt64](capacity=columns - level)
            for k in range(level, columns):
                row.append((points[k] + p - points[k - level]) % p)
            if p < (UInt64(1) << 32):
                _inverse32(row, _Prime32.of(p))
            else:
                for k in range(len(row)):
                    row[k] = _inverse(row[k], p)
            inverses.append(row^)
        budget.charge(UInt64(_product(count, 2)))
        var inverse = _inverse(_residue(modulus, p, budget), p)
        var extra = List[UInt64]()
        if checked:
            extra = [points[columns], points[columns + 1]]
            _ = points.pop()
            _ = points.pop()
        for s in range(series):
            var check = List[UInt64]()
            if checked:
                check = [samples[s][columns], samples[s][columns + 1]]
                _ = samples[s].pop()
                _ = samples[s].pop()
            var coefficients = _interpolate(points, samples[s], inverses, p)
            for e in range(len(extra)):
                if _horner(coefficients, extra[e], p) != check[e]:
                    return None
            for j in range(columns):
                var k = s * columns + j
                var x = _residue(values[k], p, budget)
                var t = UInt64((UInt128((coefficients[j] + p - x) % p) * UInt128(inverse)) % UInt128(p))
                if t:
                    try:
                        values[k] = values[k] + modulus * Integer(Int(t))
                    except e:
                        raise PolyError(PolyErrorKind.CAPACITY, String(e))
        try:
            modulus = modulus * Integer(Int(p))
        except e:
            raise PolyError(PolyErrorKind.CAPACITY, String(e))
        primes += 1
        if primes < attempt:
            continue
        attempt = 2 * primes
        var bound: Integer
        try:
            bound = isqrt(modulus // Integer(2))
        except e:
            raise PolyError(PolyErrorKind.CAPACITY, String(e))
        var u_terms = List[Term[QQ]]()
        var v_terms = List[Term[QQ]]()
        var complete = True
        for k in range(count):
            if values[k] == Integer(0):
                continue
            var fraction = _reconstruct(values[k], modulus, bound, budget)
            if not fraction:
                complete = False
                break
            var s = k // columns
            var j = k % columns
            var pairs = List[Tuple[Int64, UInt32]]()
            var power = s if s < m else s - m
            if power:
                pairs.append((variable, UInt32(power)))
            if j:
                pairs.append((other, UInt32(j)))
            if variable > other and len(pairs) == 2:
                pairs = [pairs[1], pairs[0]]
            var term = Term[QQ](Monomial(pairs^, budget), QQ(fraction.take()))
            if s < m:
                u_terms.append(term^)
            else:
                v_terms.append(term^)
        if not complete:
            continue
        var u = SparsePoly[QQ].from_terms(u_terms^, f.order, NoContext(), budget)
        var v = SparsePoly[QQ].from_terms(v_terms^, f.order, NoContext(), budget)
        if u.mul(f, budget).add(v.mul(g, budget), budget) == r:
            return (u^, v^)
    return None
