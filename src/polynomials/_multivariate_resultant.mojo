"""Multivariate resultants by dense evaluation and interpolation modulo primes.

Res_x(f, g) has degree at most n·deg_y(f) + m·deg_y(g) in each remaining
variable y, where m = deg_x f and n = deg_x g, and at most tdeg(f)·tdeg(g)
in all of them together (the Bezout bound). Modulo a prime the
remaining variables are specialized on the grid 0 … D_y one variable at a
time (Horner on dense arrays), univariate resultants come from the word-size
recurrence, and Newton interpolation on the points 0 … D_y rebuilds each
variable. Where a leading coefficient vanishes at a point, the Sylvester
resultant of the formal degrees follows from the actual-degree one:
Res_{m,n} = (-1)^((m-m')n) lc(g)^(m-m') Res_{m',n} and
Res_{m,n} = lc(f)^(n-n') Res_{m,n'}; it is zero when both degrees drop.
Over ZZ the images are combined by CRT until the product of primes exceeds
twice (sum_i |f_i|_1^2)^(n/2) (sum_j |g_j|_1^2)^(m/2), with f_i, g_j the
coefficients in x: a bound on every coefficient of the resultant by
Hadamard's inequality on the unit torus of the remaining variables.
"""
from apn_mojo import Integer
from std.collections import Dict
from .errors import Budget, PolyError, _allocation, _charge_size, _sum, _product
from .coefficients import CoefficientRing, ZZ, ModInt, NoContext
from .sparse import SparsePoly, Term
from .monomial import Monomial
from .primes import prime, _residue
from ._nmod_euclid import _Field, _resultant_words, _resultant_in, _trim_words, _inverses as _inverse_words
from ._univariate_modular import _Reconstruction

# Specialization grids above this many points use the subresultant PRS.
comptime _GRID_LIMIT = 1 << 22


struct _Layout(Movable):
    """Dense array axes: the remaining variables in order, then the eliminated one innermost."""
    var others: List[Int64]
    var f_extents: List[Int]
    var g_extents: List[Int]
    var bounds: List[Int]
    var m: Int
    var n: Int

    def __init__(out self):
        self.others = List[Int64]()
        self.f_extents = List[Int]()
        self.g_extents = List[Int]()
        self.bounds = List[Int]()
        self.m = 0
        self.n = 0


def _size(extents: List[Int], start: Int) -> Int:
    var size = 1
    for i in range(start, len(extents)):
        size *= extents[i]
    return size


def _layout[C: CoefficientRing](f: SparsePoly[C], g: SparsePoly[C], variable: Int64, variables: List[Int64], mut budget: Budget) raises PolyError -> Optional[_Layout]:
    var layout = _Layout()
    layout.m = f.degree_in(variable, budget)
    layout.n = g.degree_in(variable, budget)
    # Homogenizing gives the Bezout bound tdeg(f) tdeg(g) on the resultant's
    # total degree, so on its degree in each remaining variable as well.
    var total = _product(_total_degree(f, budget), _total_degree(g, budget))
    var points = 1
    for v in variables:
        budget.charge(1)
        if v == variable:
            continue
        var df = max(f.degree_in(v, budget), 0)
        var dg = max(g.degree_in(v, budget), 0)
        var bound = min(_sum(_product(layout.n, df), _product(layout.m, dg)), total)
        points = _product(points, _sum(bound, 1))
        if points > _GRID_LIMIT:
            return None
        layout.others.append(v)
        layout.f_extents.append(df + 1)
        layout.g_extents.append(dg + 1)
        layout.bounds.append(bound)
    layout.f_extents.append(layout.m + 1)
    layout.g_extents.append(layout.n + 1)
    if _product(_size(layout.f_extents, 0), 1) > _GRID_LIMIT or _size(layout.g_extents, 0) > _GRID_LIMIT:
        return None
    return layout^


def _total_degree[C: CoefficientRing](p: SparsePoly[C], mut budget: Budget) raises PolyError -> Int:
    var top = 0
    for t in p.terms:
        budget.charge(1)
        var d = 0
        for pair in t.monomial.pairs:
            d += Int(pair[1])
        top = max(top, d)
    return top


def _offset(m: Monomial, variable: Int64, layout: _Layout, extents: List[Int]) -> Int:
    var x = 0
    for pair in m.pairs:
        if pair[0] == variable:
            x = Int(pair[1])
    var offset = x
    var stride = extents[len(extents) - 1]
    for axis in range(len(layout.others) - 1, -1, -1):
        var e = 0
        for pair in m.pairs:
            if pair[0] == layout.others[axis]:
                e = Int(pair[1])
        offset += e * stride
        stride *= extents[axis]
    return offset


def _dense_image(p: SparsePoly[ZZ], variable: Int64, layout: _Layout, extents: List[Int], field: _Field, mut budget: Budget) raises PolyError -> List[UInt64]:
    var size = _size(extents, 0)
    _allocation[UInt64](size, budget)
    var out = List[UInt64](length=size, fill=0)
    for t in p.terms:
        budget.charge(1)
        _charge_size(len(t.monomial.pairs), budget)
        out[_offset(t.monomial, variable, layout, extents)] = field.enter(_residue(t.coefficient.value, field.p, budget))
    return out^


def _dense_image_mod(p: SparsePoly[ModInt], variable: Int64, layout: _Layout, extents: List[Int], field: _Field, mut budget: Budget) raises PolyError -> List[UInt64]:
    var size = _size(extents, 0)
    _allocation[UInt64](size, budget)
    var out = List[UInt64](length=size, fill=0)
    for t in p.terms:
        budget.charge(1)
        _charge_size(len(t.monomial.pairs), budget)
        out[_offset(t.monomial, variable, layout, extents)] = field.enter(t.coefficient.value)
    return out^


def _horner(a: List[UInt64], extent: Int, point: UInt64, field: _Field, mut budget: Budget) raises PolyError -> List[UInt64]:
    """Specialize the outermost axis of a dense array at a Montgomery-form point."""
    var block = len(a) // extent
    _allocation[UInt64](block, budget)
    _charge_size(len(a), budget)
    var out = List[UInt64](capacity=block)
    var top = (extent - 1) * block
    for i in range(block):
        out.append(a[top + i])
    for k in range(extent - 2, -1, -1):
        var base = k * block
        for i in range(block):
            out[i] = field.add(field.mul(out[i], point), a[base + i])
    return out^


def _horner_into(mut out: List[UInt64], a: List[UInt64], extent: Int, point: UInt64, field: _Field, mut budget: Budget) raises PolyError:
    """_horner into a buffer of the block's length."""
    var block = len(out)
    _charge_size(len(a), budget)
    var top = (extent - 1) * block
    for i in range(block):
        out.unsafe_get(i) = a.unsafe_get(top + i)
    for k in range(extent - 2, -1, -1):
        var base = k * block
        for i in range(block):
            out.unsafe_get(i) = field.add(field.mul(out.unsafe_get(i), point), a.unsafe_get(base + i))


def _leaf_in(f: List[UInt64], g: List[UInt64], mut fw: List[UInt64], mut gw: List[UInt64], m: Int, n: Int, field: _Field, mut budget: Budget) raises PolyError -> Tuple[UInt64, UInt64]:
    """_leaf of copies of f and g made in the work buffers fw and gw, as a numerator over a nonzero denominator."""
    fw.clear()
    gw.clear()
    for v in f:
        fw.append(v)
    for v in g:
        gw.append(v)
    _trim_words(fw)
    _trim_words(gw)
    if not len(fw) or not len(gw):
        return (UInt64(0), field.one())
    var df = len(fw) - 1
    var dg = len(gw) - 1
    if df < m and dg < n:
        return (UInt64(0), field.one())
    var lc_f = fw[df]
    var lc_g = gw[dg]
    var fraction = _resultant_in(fw, gw, field, budget)
    var r = fraction[0]
    if df < m:
        r = field.mul(r, field.power(lc_g, m - df))
        if ((m - df) * n) & 1:
            r = field.sub(0, r)
    if dg < n:
        r = field.mul(r, field.power(lc_f, n - dg))
    return (r, fraction[1])


def _leaf(var f: List[UInt64], var g: List[UInt64], m: Int, n: Int, field: _Field, mut budget: Budget) raises PolyError -> UInt64:
    _trim_words(f)
    _trim_words(g)
    if not len(f) or not len(g):
        return 0
    var df = len(f) - 1
    var dg = len(g) - 1
    if df < m and dg < n:
        return 0
    var lc_f = f[df]
    var lc_g = g[dg]
    var r = _resultant_words(f^, g^, field, budget)
    if df < m:
        r = field.mul(r, field.power(lc_g, m - df))
        if ((m - df) * n) & 1:
            r = field.sub(0, r)
    if dg < n:
        r = field.mul(r, field.power(lc_f, n - dg))
    return r


def _image(f: List[UInt64], g: List[UInt64], level: Int, layout: _Layout, inverses: List[UInt64], field: _Field, mut budget: Budget) raises PolyError -> List[UInt64]:
    """The resultant modulo p as a dense array over the remaining variables from `level` on."""
    if level == len(layout.others):
        budget.charge(1)
        return [_leaf(f.copy(), g.copy(), layout.m, layout.n, field, budget)]
    var bound = layout.bounds[level]
    var block = 1
    for axis in range(level + 1, len(layout.others)):
        block *= layout.bounds[axis] + 1
    var count = (bound + 1) * block
    _allocation[UInt64](count, budget)
    var values = List[UInt64](length=count, fill=0)
    if level + 1 == len(layout.others):
        # The last variable: specialize into reused buffers and take each
        # univariate resultant in place, with no allocation per point.
        var f_block = len(f) // layout.f_extents[level]
        var g_block = len(g) // layout.g_extents[level]
        _allocation[UInt64](_sum(_product(2, f_block), _product(2, g_block)), budget)
        var fa = List[UInt64](length=f_block, fill=0)
        var ga = List[UInt64](length=g_block, fill=0)
        var fw = List[UInt64](capacity=f_block)
        var gw = List[UInt64](capacity=g_block)
        var denominators = List[UInt64](capacity=bound + 1)
        for a in range(bound + 1):
            budget.charge(1)
            var point = field.enter(UInt64(a))
            _horner_into(fa, f, layout.f_extents[level], point, field, budget)
            _horner_into(ga, g, layout.g_extents[level], point, field, budget)
            budget.charge(1)
            var fraction = _leaf_in(fa, ga, fw, gw, layout.m, layout.n, field, budget)
            values[a] = fraction[0]
            denominators.append(fraction[1])
        # One inversion for the whole row.
        var inverted = _inverse_words(denominators, field)
        for a in range(bound + 1):
            values[a] = field.mul(values[a], inverted[a])
    else:
        for a in range(bound + 1):
            budget.charge(1)
            var point = field.enter(UInt64(a))
            var fa = _horner(f, layout.f_extents[level], point, field, budget)
            var ga = _horner(g, layout.g_extents[level], point, field, budget)
            var r = _image(fa, ga, level + 1, layout, inverses, field, budget)
            for t in range(block):
                values[a * block + t] = r[t]
    # Newton divided differences on the points 0 … bound: x_i - x_(i-j) = j.
    _charge_size(_product(_product(bound + 1, bound + 1), block), budget)
    for j in range(1, bound + 1):
        for i in range(bound, j - 1, -1):
            for t in range(block):
                values[i * block + t] = field.mul(field.sub(values[i * block + t], values[(i - 1) * block + t]), inverses[j])
    # Expand the Newton form c_0 + c_1 y + c_2 y(y - 1) + … into monomial coefficients.
    var out = List[UInt64](length=count, fill=0)
    for t in range(block):
        out[t] = values[bound * block + t]
    var length = 1
    for j in range(bound - 1, -1, -1):
        var shift = field.enter(UInt64(j))
        # out = out·(y - j) + c_j over the first `length + 1` coefficients.
        for t in range(block):
            var carry = UInt64(0)
            for i in range(length + 1):
                var current = out[i * block + t] if i < length else UInt64(0)
                out[i * block + t] = field.sub(carry, field.mul(current, shift))
                carry = current
            out[t] = field.add(out[t], values[j * block + t])
        length += 1
    return out^


def _inverses(bound: Int, field: _Field) raises PolyError -> List[UInt64]:
    var out = List[UInt64](length=bound + 1, fill=0)
    for j in range(1, bound + 1):
        out[j] = field.inverse(field.enter(UInt64(j)))
    return out^


def _terms[C: CoefficientRing](values: List[C], layout: _Layout, exemplar: SparsePoly[C], mut budget: Budget) raises PolyError -> SparsePoly[C]:
    var terms = List[Term[C]]()
    for index in range(len(values)):
        budget.charge(1)
        if values[index].is_zero():
            continue
        var pairs = List[Tuple[Int64, UInt32]]()
        var rest = index
        var exponents = List[Int](length=len(layout.others), fill=0)
        for axis in range(len(layout.others) - 1, -1, -1):
            exponents[axis] = rest % (layout.bounds[axis] + 1)
            rest //= layout.bounds[axis] + 1
        for axis in range(len(layout.others)):
            if exponents[axis]:
                pairs.append((layout.others[axis], UInt32(exponents[axis])))
        _allocation[Term[C]](1, budget)
        terms.append(Term[C](Monomial(pairs^, budget), values[index].copy()))
    return SparsePoly[C].from_terms(terms^, exemplar.order, exemplar.ring, budget)


def _max_bound(layout: _Layout) -> Int:
    var top = 1
    for b in layout.bounds:
        top = max(top, b)
    return top


def _row_squares(p: SparsePoly[ZZ], variable: Int64, mut budget: Budget) raises PolyError -> ZZ:
    """The sum over k of |p_k|_1^2, p_k the coefficient of variable^k."""
    var norms = Dict[UInt32, ZZ]()
    for t in p.terms:
        budget.charge(1)
        var k = UInt32(0)
        for pair in t.monomial.pairs:
            if pair[0] == variable:
                k = pair[1]
        var c = t.coefficient.copy()
        if c.value.sign() < 0:
            c = c.neg(budget)
        var old = norms.get(k, ZZ(Integer(0)))
        norms[k] = old.add(c, budget)
    var total = ZZ(Integer(0))
    for entry in norms.items():
        budget.charge(1)
        total = total.add(entry.value.mul(entry.value, budget), budget)
    return total^


def _power(base: ZZ, k: Int, mut budget: Budget) raises PolyError -> ZZ:
    var result = ZZ(Integer(1))
    for _ in range(k):
        budget.charge(1)
        result = result.mul(base, budget)
    return result^


def _integer_resultant(f: SparsePoly[ZZ], g: SparsePoly[ZZ], variable: Int64, variables: List[Int64], mut budget: Budget) raises PolyError -> Optional[SparsePoly[ZZ]]:
    """Res_x(f, g) for positive degrees in x and at least one other variable, or None for a large grid."""
    var maybe = _layout(f, g, variable, variables, budget)
    if not maybe:
        return None
    var layout = maybe.take()
    # Each coefficient of the resultant is at most its maximum modulus on the
    # unit torus |y_i| = 1, where every Sylvester entry f_i(y) has modulus at
    # most |f_i|_1; Hadamard bounds the determinant there by
    # (sum_i |f_i|_1^2)^(n/2) (sum_j |g_j|_1^2)^(m/2).
    # The primes suffice once their product M exceeds twice that bound, that is
    # once M^2 > 4 S_f^n S_g^m, compared exactly.
    var target = ZZ(Integer(4)).mul(_power(_row_squares(f, variable, budget), layout.n, budget), budget).mul(_power(_row_squares(g, variable, budget), layout.m, budget), budget)
    var modulus = ZZ(Integer(1))
    var count = 1
    for b in layout.bounds:
        count *= b + 1
    var values = _Reconstruction(count, budget)
    var index = 0
    while not (modulus.mul(modulus, budget).value > target.value):
        budget.charge(1)
        var p = prime(index, budget)
        index = _sum(index, 1)
        var field = _Field.of(p)
        var image = _image(_dense_image(f, variable, layout, layout.f_extents, field, budget), _dense_image(g, variable, layout, layout.g_extents, field, budget), 0, layout, _inverses(_max_bound(layout), field), field, budget)
        for i in range(len(image)):
            image[i] = field.leave(image[i])
        _ = values.merge(image, p, budget)
        modulus = modulus.mul(ZZ(Integer(p)), budget)
    _allocation[ZZ](count, budget)
    var coefficients = List[ZZ](capacity=count)
    for i in range(count):
        coefficients.append(ZZ(values.values[i].copy()))
    return _terms(coefficients, layout, f, budget)


def _field_resultant(f: SparsePoly[ModInt], g: SparsePoly[ModInt], variable: Int64, variables: List[Int64], mut budget: Budget) raises PolyError -> Optional[SparsePoly[ModInt]]:
    """Res_x(f, g) over an odd prime field with more points than any degree bound, or None."""
    var maybe = _layout(f, g, variable, variables, budget)
    if not maybe:
        return None
    var layout = maybe.take()
    var p = f.ring.value
    if p <= 2 or UInt64(_max_bound(layout)) >= p:
        return None
    var field = _Field.of(p)
    var image = _image(_dense_image_mod(f, variable, layout, layout.f_extents, field, budget), _dense_image_mod(g, variable, layout, layout.g_extents, field, budget), 0, layout, _inverses(_max_bound(layout), field), field, budget)
    _allocation[ModInt](len(image), budget)
    var coefficients = List[ModInt](capacity=len(image))
    for v in image:
        coefficients.append(ModInt(field.leave(v), f.ring))
    return _terms(coefficients, layout, f, budget)
