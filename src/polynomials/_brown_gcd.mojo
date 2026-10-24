"""Brown's dense modular gcd for multivariate integer polynomials.

Modulo each prime the inputs are dense word arrays over their exponent
boxes, the last variable outermost. One level of the recursion removes the
contents in that variable y, takes gamma = gcd of the inputs' leading
coefficients (LEX over the other variables) in Z_p[y], specializes y at
points where neither leading coefficient vanishes, recurses, and keeps the
images whose leading monomial is smallest (larger ones come from unlucky
points). Images scaled by gamma(point) are interpolated in y from
deg(gamma) + d + 1 points, where d bounds deg_y of the gcd by one univariate
image in y with the other variables specialized. The primitive part of the
interpolant times the content gcd is the image gcd, made monic in LEX.

Over ZZ the images, scaled by the gcd of the inputs' LEX leading integer
coefficients, are combined by symmetric CRT until a prime changes nothing;
the primitive candidate is accepted only when divide_exact shows that it
divides both inputs. Anything unexpected returns None, and the caller uses
another algorithm.
"""
from apn_mojo import Integer
from .errors import Budget, PolyError, _allocation, _charge_size, _sum, _product
from .coefficients import ZZ, NoContext
from .sparse import SparsePoly, Term
from .monomial import Monomial
from ._gcd_verify import _divides_both
from .primes import prime, _residue
from ._coefficient_ops import _gcd_coefficient, _exact
from ._nmod_euclid import _Field, _gcd_words, _trim_words, _inverses
from ._univariate_modular import _Reconstruction, _Words

# Exponent boxes above this many cells, or more primes than this, give up.
comptime _CELLS = 1 << 22
comptime _PRIMES = 64


struct _Grid(Copyable, Movable):
    """A dense array; extents are outermost first, the first variable innermost."""
    var data: List[UInt64]
    var extents: List[Int]

    def __init__(out self, var data: List[UInt64], var extents: List[Int]):
        self.data = data^
        self.extents = extents^


def _size(extents: List[Int], start: Int) -> Int:
    var size = 1
    for i in range(start, len(extents)):
        size *= extents[i]
    return size


def _decode(index: Int, extents: List[Int]) -> List[Int]:
    """Exponents of a cell, innermost (first) variable first."""
    var out = List[Int](capacity=len(extents))
    var rest = index
    for axis in range(len(extents) - 1, -1, -1):
        out.append(rest % extents[axis])
        rest //= extents[axis]
    return out^


def _encode(exponents: List[Int], extents: List[Int]) -> Int:
    """The cell of exponents (first variable first), or -1 outside the box."""
    var index = 0
    var n = len(extents)
    for axis in range(n):
        var e = exponents[n - 1 - axis]
        if e >= extents[axis]:
            return -1
        index = index * extents[axis] + e
    return index


def _lead_from(g: _Grid, first: Int) -> Int:
    """The cell of axes first … n-1 holding the LEX-leading nonzero entry, or -1.

    Cells are visited from the largest key down: the innermost axis is the
    most significant digit and axis `first` the least, so the first nonzero
    cell met is the leading one. With first = 1 a cell is the coefficient
    vector in the outermost variable.
    """
    var n = len(g.extents)
    var block = _size(g.extents, first)
    var outer = len(g.data) // block if block else 0
    var strides = List[Int](length=n, fill=1)
    for axis in range(n - 2, first - 1, -1):
        strides[axis] = strides[axis + 1] * g.extents[axis + 1]
    var digits = List[Int](length=n, fill=0)
    var index = 0
    for axis in range(first, n):
        digits[axis] = g.extents[axis] - 1
        index += digits[axis] * strides[axis]
    while True:
        for j in range(outer):
            if g.data.unsafe_get(j * block + index) != 0:
                return index
        var axis = first
        while axis < n and digits[axis] == 0:
            digits[axis] = g.extents[axis] - 1
            index += digits[axis] * strides[axis]
            axis += 1
        if axis == n:
            return -1
        digits[axis] -= 1
        index -= strides[axis]


def _lead_index(g: _Grid) -> Int:
    """The LEX-leading nonzero cell, or -1."""
    return _lead_from(g, 0)


def _less(a: List[Int], b: List[Int]) -> Bool:
    """LEX comparison with the first variable most significant."""
    for i in range(len(a)):
        if a[i] != b[i]:
            return a[i] < b[i]
    return False


def _vector(g: _Grid, t: Int) -> List[UInt64]:
    """The coefficients in the outermost variable at inner cell t."""
    var block = _size(g.extents, 1)
    var out = List[UInt64](capacity=g.extents[0])
    for j in range(g.extents[0]):
        out.append(g.data[j * block + t])
    _trim_words(out)
    return out^


def _evaluate(a: List[UInt64], x: UInt64, field: _Field) -> UInt64:
    var value = UInt64(0)
    for i in range(len(a) - 1, -1, -1):
        value = field.add(field.mul(value, x), a[i])
    return value


def _product_words(a: List[UInt64], b: List[UInt64], field: _Field) -> List[UInt64]:
    if not len(a) or not len(b):
        return List[UInt64]()
    var out = List[UInt64](length=len(a) + len(b) - 1, fill=0)
    for i in range(len(a)):
        for j in range(len(b)):
            out[i + j] = field.add(out[i + j], field.mul(a[i], b[j]))
    return out^


def _quotient_words(a: List[UInt64], b: List[UInt64], field: _Field) raises PolyError -> List[UInt64]:
    """a / b for b dividing a."""
    var r = a.copy()
    _trim_words(r)
    var db = len(b) - 1
    if len(r) <= db:
        return List[UInt64]()
    var inverse = field.inverse(b[db])
    var q = List[UInt64](length=len(r) - db, fill=0)
    for k in range(len(r) - 1, db - 1, -1):
        var c = field.mul(r[k], inverse)
        q[k - db] = c
        if c == 0:
            continue
        for i in range(db + 1):
            r[k - db + i] = field.sub(r[k - db + i], field.mul(c, b[i]))
    return q^


def _content(g: _Grid, field: _Field, mut budget: Budget) raises PolyError -> List[UInt64]:
    var block = _size(g.extents, 1)
    var c = List[UInt64]()
    for t in range(block):
        budget.charge(1)
        var v = _vector(g, t)
        if not len(v):
            continue
        c = _gcd_words(c^, v^, field, budget)
        if len(c) == 1:
            break
    return c^


def _divide_vectors(g: _Grid, c: List[UInt64], field: _Field, mut budget: Budget) raises PolyError -> _Grid:
    if len(c) <= 1:
        return g.copy()
    var block = _size(g.extents, 1)
    var outer = g.extents[0] - (len(c) - 1)
    _allocation[UInt64](outer * block, budget)
    var data = List[UInt64](length=outer * block, fill=0)
    for t in range(block):
        budget.charge(1)
        var q = _quotient_words(_vector(g, t), c, field)
        for j in range(len(q)):
            data[j * block + t] = q[j]
    var extents = g.extents.copy()
    extents[0] = outer
    return _Grid(data^, extents^)


def _inner_lead(g: _Grid) -> Int:
    """The inner cell of the LEX-largest nonzero coefficient vector."""
    return _lead_from(g, 1)


def _lead(g: _Grid) -> List[Int]:
    """Exponents of the LEX-leading cell, first variable first."""
    var index = _lead_index(g)
    if index < 0:
        return List[Int]()
    return _decode(index, g.extents)


def _horner_pair(g: _Grid, x: UInt64, field: _Field, mut budget: Budget) raises PolyError -> List[_Grid]:
    """Specialize the outermost variable at x and at -x in one pass.

    With g = E(y^2) + y O(y^2) over the outermost variable y, the images are
    E(x^2) + x O(x^2) and E(x^2) - x O(x^2).
    """
    var block = _size(g.extents, 1)
    var n = g.extents[0]
    _charge_size(len(g.data), budget)
    _allocation[UInt64](_product(4, block), budget)
    var even = List[UInt64](length=block, fill=0)
    var odd = List[UInt64](length=block, fill=0)
    var square = field.mul(x, x)
    var source = g.data.unsafe_ptr()
    var pe = even.unsafe_ptr()
    var po = odd.unsafe_ptr()
    for k in range(n - 1, -1, -1):
        var row = source.unsafe_offset(k * block)
        var target = pe if k % 2 == 0 else po
        for i in range(block):
            var cell = target.unsafe_offset(i)
            cell.unsafe_store(field.add(field.mul(cell.unsafe_load(), square), row.unsafe_offset(i).unsafe_load()))
    var at_plus = List[UInt64](capacity=block)
    var at_minus = List[UInt64](capacity=block)
    for i in range(block):
        var t = field.mul(odd[i], x)
        at_plus.append(field.add(even[i], t))
        at_minus.append(field.sub(even[i], t))
    var extents = List[Int]()
    for i in range(1, len(g.extents)):
        extents.append(g.extents[i])
    var out = List[_Grid](capacity=2)
    out.append(_Grid(at_plus^, extents.copy()))
    out.append(_Grid(at_minus^, extents^))
    return out^


def _outer_degree(g: _Grid) -> Int:
    var block = _size(g.extents, 1)
    for j in range(g.extents[0] - 1, -1, -1):
        for t in range(block):
            if g.data[j * block + t] != 0:
                return j
    return -1


def _outer_image(g: _Grid, attempt: Int, field: _Field, mut budget: Budget) raises PolyError -> List[UInt64]:
    """The univariate polynomial in the outermost variable with the others specialized.

    The inner axes collapse innermost first by Horner over contiguous runs, so
    the cost is one pass over the cells and no cell is decoded.
    """
    _charge_size(len(g.data), budget)
    var data = g.data.copy()
    for axis in range(len(g.extents) - 1, 0, -1):
        var extent = g.extents[axis]
        var point = field.enter(UInt64(3 + axis + attempt * len(g.extents)))
        var groups = len(data) // extent
        var out = List[UInt64](capacity=groups)
        for i in range(groups):
            var base = i * extent
            var value = UInt64(0)
            for k in range(extent - 1, -1, -1):
                value = field.add(field.mul(value, point), data[base + k])
            out.append(value)
        data = out^
    return data^


def _brown(a: _Grid, b: _Grid, field: _Field, mut budget: Budget) raises PolyError -> Optional[_Grid]:
    """The LEX-monic gcd of two nonzero grids over the same variables, or None."""
    budget.charge(1)
    if len(a.extents) == 1:
        var g = _gcd_words(a.data.copy(), b.data.copy(), field, budget)
        var n = len(g)
        return _Grid(g^, [n])
    var content_a = _content(a, field, budget)
    var content_b = _content(b, field, budget)
    var common = _gcd_words(content_a.copy(), content_b.copy(), field, budget)
    var pa = _divide_vectors(a, content_a, field, budget)
    var pb = _divide_vectors(b, content_b, field, budget)
    var lc_a = _vector(pa, _inner_lead(pa))
    var lc_b = _vector(pb, _inner_lead(pb))
    var gamma = _gcd_words(lc_a.copy(), lc_b.copy(), field, budget)
    # Bound deg_y gcd(pa, pb) by one image in y, or by the smaller degree.
    var da = _outer_degree(pa)
    var db = _outer_degree(pb)
    var bound = min(da, db)
    for attempt in range(3):
        var ia = _outer_image(pa, attempt, field, budget)
        var ib = _outer_image(pb, attempt, field, budget)
        if ia[da] != 0 and ib[db] != 0:
            bound = len(_gcd_words(ia^, ib^, field, budget)) - 1
            break
    var need = bound + len(gamma)
    var n = len(a.extents)
    var inner = List[Int]()
    for i in range(1, n):
        inner.append(min(pa.extents[i], pb.extents[i]))
    var block = _size(inner, 0)
    _allocation[UInt64](_product(need, block), budget)
    var values = List[UInt64](length=need * block, fill=0)
    var points = List[UInt64](length=need, fill=0)
    var best = List[Int]()
    var count = 0
    var x = 0
    while count < need:
        budget.charge(1)
        x += 1
        if UInt64(x) >= field.p or x > 4 * need + 32:
            return None
        # The points x and -x share one pass over each grid.
        var plus = field.enter(UInt64(x))
        var minus = field.sub(0, plus)
        var good_plus = _evaluate(lc_a, plus, field) != 0 and _evaluate(lc_b, plus, field) != 0
        var good_minus = _evaluate(lc_a, minus, field) != 0 and _evaluate(lc_b, minus, field) != 0
        if not good_plus and not good_minus:
            continue
        var images_a = _horner_pair(pa, plus, field, budget)
        var images_b = _horner_pair(pb, plus, field, budget)
        for side in range(2):
            if count >= need:
                break
            if not (good_plus if side == 0 else good_minus):
                continue
            var point = plus if side == 0 else minus
            var image = _brown(images_a[side], images_b[side], field, budget)
            if not image:
                return None
            ref g = image.value()
            var lead = _lead(g)
            if count == 0 or _less(lead, best):
                best = lead^
                count = 0
            elif _less(best, lead):
                continue
            var scale = _evaluate(gamma, point, field)
            var base = count * block
            for t in range(block):
                values[base + t] = 0
            for index in range(len(g.data)):
                if g.data[index] == 0:
                    continue
                var cell = _recode(index, g.extents, inner)
                if cell < 0:
                    return None
                values[base + cell] = field.mul(g.data[index], scale)
            points[count] = point
            count += 1
    # Newton divided differences and expansion, all inner cells at once.
    _charge_size(_product(_product(need, need), block), budget)
    # Rows are contiguous blocks; every index below stays within need * block.
    var differences = List[UInt64](capacity=need * (need - 1) // 2)
    for j in range(1, need):
        for i in range(need - 1, j - 1, -1):
            differences.append(field.sub(points[i], points[i - j]))
    var inverses = _inverses(differences, field)
    var next = 0
    var pv = values.unsafe_ptr()
    for j in range(1, need):
        for i in range(need - 1, j - 1, -1):
            var inverse = inverses[next]
            next += 1
            var row = pv.unsafe_offset(i * block)
            var previous = pv.unsafe_offset((i - 1) * block)
            for t in range(block):
                var cell = row.unsafe_offset(t)
                cell.unsafe_store(field.mul(field.sub(cell.unsafe_load(), previous.unsafe_offset(t).unsafe_load()), inverse))
    var expanded = List[UInt64](length=need * block, fill=0)
    var pe = expanded.unsafe_ptr()
    for t in range(block):
        pe.unsafe_offset(t).unsafe_store(pv.unsafe_offset((need - 1) * block + t).unsafe_load())
    var length = 1
    for j in range(need - 2, -1, -1):
        # Multiply by (y - points[j]) in place, top row first: row i becomes
        # row i - 1 minus points[j] times row i, with rows -1 and length zero.
        var x = points[j]
        var top = pe.unsafe_offset(length * block)
        var below = pe.unsafe_offset((length - 1) * block)
        for t in range(block):
            top.unsafe_offset(t).unsafe_store(below.unsafe_offset(t).unsafe_load())
        for i in range(length - 1, 0, -1):
            var row = pe.unsafe_offset(i * block)
            var lower = pe.unsafe_offset((i - 1) * block)
            for t in range(block):
                var cell = row.unsafe_offset(t)
                cell.unsafe_store(field.sub(lower.unsafe_offset(t).unsafe_load(), field.mul(cell.unsafe_load(), x)))
        var source = pv.unsafe_offset(j * block)
        for t in range(block):
            var cell = pe.unsafe_offset(t)
            cell.unsafe_store(field.add(field.sub(0, field.mul(cell.unsafe_load(), x)), source.unsafe_offset(t).unsafe_load()))
        length += 1
    var extents = List[Int]()
    extents.append(need)
    for e in inner:
        extents.append(e)
    var interpolant = _Grid(expanded^, extents^)
    var primitive = _divide_vectors(interpolant, _content(interpolant, field, budget), field, budget)
    # Multiply back the content gcd in y.
    var pblock = _size(primitive.extents, 1)
    var outer = primitive.extents[0] + len(common) - 1
    var data = List[UInt64](length=outer * pblock, fill=0)
    for t in range(pblock):
        var v = _product_words(_vector(primitive, t), common, field)
        for j in range(len(v)):
            data[j * pblock + t] = v[j]
    var result_extents = primitive.extents.copy()
    result_extents[0] = outer
    var result = _Grid(data^, result_extents^)
    var top = _lead_index(result)
    var inverse = field.inverse(result.data[top])
    for i in range(len(result.data)):
        result.data[i] = field.mul(result.data[i], inverse)
    return result^


def _box(tops: List[Int], mut budget: Budget) raises PolyError -> List[Int]:
    """Extents outermost first from per-variable degrees: the last variable first, the first variable last."""
    _allocation[Int](len(tops), budget)
    var out = List[Int](capacity=len(tops))
    for i in range(len(tops) - 1, -1, -1):
        out.append(tops[i] + 1)
    return out^


def _cells(p: SparsePoly[ZZ], variables: List[Int64], extents: List[Int], mut budget: Budget) raises PolyError -> List[Int]:
    """Each term's cell; the extents (from _box) hold every exponent."""
    var n = len(extents)
    var weights = List[Int](length=n, fill=1)
    for axis in range(n - 2, -1, -1):
        weights[axis] = weights[axis + 1] * extents[axis + 1]
    _allocation[Int](len(p.terms), budget)
    var out = List[Int](capacity=len(p.terms))
    for t in p.terms:
        budget.charge(1)
        var index = 0
        var position = 0
        for pair in t.monomial.pairs:
            while variables[position] != pair[0]:
                position += 1
            index += Int(pair[1]) * weights[n - 1 - position]
        out.append(index)
    return out^


def _lead_term(p: SparsePoly[ZZ], variables: List[Int64], extents: List[Int]) -> Int:
    """The position in p.terms of the LEX-leading term; the extents (from _box) hold every exponent."""
    var n = len(extents)
    # Mixed-radix LEX weights by variable position, the first variable most significant.
    var weights = List[Int](length=n, fill=1)
    for position in range(n - 2, -1, -1):
        weights[position] = weights[position + 1] * extents[n - 2 - position]
    var best = 0
    var best_key = -1
    for i in range(len(p.terms)):
        var key = 0
        var position = 0
        for pair in p.terms[i].monomial.pairs:
            while variables[position] != pair[0]:
                position += 1
            key += Int(pair[1]) * weights[position]
        if key > best_key:
            best = i
            best_key = key
    return best


def _image(words: _Words, cells: List[Int], extents: List[Int], field: _Field, mut budget: Budget) raises PolyError -> _Grid:
    """The grid of p modulo field.p, from p's coefficient words (flattened once for every prime)."""
    var size = _size(extents, 0)
    _allocation[UInt64](size, budget)
    var data = List[UInt64](length=size, fill=0)
    var residues = words.residues(field, budget)
    for i in range(len(residues)):
        budget.charge(1)
        data[cells[i]] = residues[i]
    return _Grid(data^, extents.copy())


def _recode(index: Int, source: List[Int], target: List[Int]) -> Int:
    """The cell of `target` with the exponents of cell `index` of `source` (axes outermost first), or -1 outside it."""
    var rest = index
    var out = 0
    var weight = 1
    for axis in range(len(source) - 1, -1, -1):
        var e = rest % source[axis]
        rest //= source[axis]
        if e >= target[axis]:
            return -1
        out += e * weight
        weight *= target[axis]
    return out


def _integer_content(p: SparsePoly[ZZ], mut budget: Budget) raises PolyError -> ZZ:
    var c = ZZ(Integer(0))
    for t in p.terms:
        c = _gcd_coefficient(c, t.coefficient, budget)
        if c.value == Integer(1):
            break
    return c^


def _divided(p: SparsePoly[ZZ], c: ZZ, mut budget: Budget) raises PolyError -> SparsePoly[ZZ]:
    var terms = List[Term[ZZ]](capacity=len(p.terms))
    for t in p.terms:
        terms.append(Term[ZZ](t.monomial.copy(), _exact(t.coefficient, c, budget)))
    var out = SparsePoly[ZZ].zero(p.order, NoContext(), budget)
    out.terms = terms^
    return out^


def _primitive(var p: SparsePoly[ZZ], mut budget: Budget) raises PolyError -> SparsePoly[ZZ]:
    var c = _integer_content(p, budget)
    if c.value == Integer(1):
        return p^
    return _divided(p, c, budget)


def _brown_gcd(f: SparsePoly[ZZ], h: SparsePoly[ZZ], variables: List[Int64], tops_f: List[Int], tops_h: List[Int], mut budget: Budget) raises PolyError -> Optional[SparsePoly[ZZ]]:
    """The gcd up to sign of two nonzero integer polynomials in two or more variables, or None."""
    var ef = _box(tops_f, budget)
    var eh = _box(tops_h, budget)
    var cells_f = 1
    var cells_h = 1
    for i in range(len(ef)):
        cells_f = _product(cells_f, ef[i])
        cells_h = _product(cells_h, eh[i])
        if cells_f > _CELLS or cells_h > _CELLS:
            return None
    var cf = _integer_content(f, budget)
    var ch = _integer_content(h, budget)
    var common = _gcd_coefficient(cf, ch, budget)
    # Primitive inputs, the common case, are used in place.
    var result: Optional[SparsePoly[ZZ]]
    if cf.value == Integer(1) and ch.value == Integer(1):
        result = _brown_primitive(f, h, ef, eh, variables, budget)
    elif cf.value == Integer(1):
        result = _brown_primitive(f, _divided(h, ch, budget), ef, eh, variables, budget)
    elif ch.value == Integer(1):
        result = _brown_primitive(_divided(f, cf, budget), h, ef, eh, variables, budget)
    else:
        result = _brown_primitive(_divided(f, cf, budget), _divided(h, ch, budget), ef, eh, variables, budget)
    if not result:
        return None
    return result.take().scale(common, budget)


def _brown_primitive(pf: SparsePoly[ZZ], ph: SparsePoly[ZZ], ef: List[Int], eh: List[Int], variables: List[Int64], mut budget: Budget) raises PolyError -> Optional[SparsePoly[ZZ]]:
    """The primitive gcd up to sign of two primitive polynomials with boxes ef and eh, or None."""
    var where_f = _cells(pf, variables, ef, budget)
    var where_h = _cells(ph, variables, eh, budget)
    ref lf = pf.terms[_lead_term(pf, variables, ef)].coefficient
    ref lh = ph.terms[_lead_term(ph, variables, eh)].coefficient
    var gamma = _gcd_coefficient(lf, lh, budget)
    var words_f = _Words(pf, budget)
    var words_h = _Words(ph, budget)
    var shape = List[Int]()
    for i in range(len(ef)):
        shape.append(min(ef[i], eh[i]))
    var count = _size(shape, 0)
    var values = _Reconstruction(0, budget)
    var best = List[Int]()
    var primes = 0
    var failures = 0
    for index in range(_PRIMES):
        budget.charge(1)
        var p = prime(index, budget)
        if _residue(lf.value, p, budget) == 0 or _residue(lh.value, p, budget) == 0:
            continue
        var field = _Field.of(p)
        var image = _brown(_image(words_f, where_f, ef, field, budget), _image(words_h, where_h, eh, field, budget), field, budget)
        if not image:
            return None
        ref g = image.value()
        var lead = _lead(g)
        if primes == 0 or _less(lead, best):
            best = lead^
            values = _Reconstruction(count, budget)
            primes = 0
        elif _less(best, lead):
            continue
        var scale = field.enter(_residue(gamma.value, p, budget))
        var residues = List[UInt64](length=count, fill=0)
        for i in range(len(g.data)):
            if g.data[i] == 0:
                continue
            var cell = _recode(i, g.extents, shape)
            if cell < 0:
                return None
            residues[cell] = field.leave(field.mul(g.data[i], scale))
        var changed = values.merge(residues, p, budget)
        primes += 1
        # Small coefficients are exact after one prime: try the first image, then
        # only images that a further prime left unchanged.
        if primes != 1 and changed:
            continue
        var terms = List[Term[ZZ]]()
        for cell in range(count):
            if values.values[cell] == Integer(0):
                continue
            var pairs = List[Tuple[Int64, UInt32]]()
            var rest = cell
            for axis in range(len(shape) - 1, -1, -1):
                var e = rest % shape[axis]
                rest //= shape[axis]
                if e:
                    pairs.append((variables[len(shape) - 1 - axis], UInt32(e)))
            terms.append(Term[ZZ](Monomial(pairs^, budget), ZZ(values.values[cell].copy())))
        var candidate = _primitive(SparsePoly[ZZ].from_terms(terms^, pf.order, NoContext(), budget), budget)
        if _divides_both(pf, ph, candidate, budget):
            return candidate^
        if primes == 1:
            continue
        failures += 1
        if failures >= 3:
            return None
    return None
