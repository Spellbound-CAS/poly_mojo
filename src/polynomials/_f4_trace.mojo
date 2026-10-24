"""F4 over word primes recording a trace of every new polynomial.

The run follows `_f4._f4_run`: the same pair selection, matrices and
insertion order, with the reductions done row by row so that every pivot
multiple taken out is recorded. Trace polynomial i, for i below the number
of generators, is generator i; each later one is a step, a list of
(source, monomial, coefficient) terms that sum to it:

- a generator top-reduced by the elements so far, made monic;
- a pair or reducer row reduced by the pivots, made monic;
- a new pivot reduced by the new pivots right of its leading column;
- a basis element's tail reduced by the final matrix's rows.

A matrix row is a monomial times an element, whose trace number is known;
a pivot found earlier in the round has its own number. Only top reduction
runs on the generators: insertion needs leading monomials that no other
divides, and the reduced basis is the same whichever path reaches it.
"""
from .coefficients import ModInt, Modulus
from .errors import Budget, PolyError, _allocation
from ._packed_groebner import _Engine, _Packed, _Term, _compare, _div, _divides, _equal, _packed
from ._f4 import _Matrix, _Rows, _Store, _raw_hash, _reciprocal, _residue, _scatter, _times, _insertion, _inserted
from ._f4 import _LANES, _Lanes, _LaneRows, _Primes, _axpy_lanes, _mod_lanes
from ._f4 import _BLOCK, _Block, _Field, _Lane, _axpy_block, _residues
from .errors import _product


struct _Traced[W: Int](Movable):
    """The steps so far, flattened, and each step's start."""
    var sources: List[Int]
    var monomials: List[_Packed[Self.W]]
    var coefficients: List[UInt64]
    var starts: List[Int]
    var generators: Int

    def __init__(out self, generators: Int):
        self.sources = List[Int]()
        self.monomials = List[_Packed[Self.W]]()
        self.coefficients = List[UInt64]()
        self.starts = [0]
        self.generators = generators

    def add(mut self, source: Int, monomial: _Packed[Self.W], coefficient: UInt64):
        self.sources.append(source)
        self.monomials.append(monomial)
        self.coefficients.append(coefficient)

    def close(mut self) -> Int:
        """End the current step; its trace number."""
        self.starts.append(len(self.sources))
        return self.generators + len(self.starts) - 2

    def scale_open(mut self, factor: UInt64, prime: UInt64):
        """Multiply the terms of the step still open by factor."""
        for k in range(self.starts[len(self.starts) - 1], len(self.sources)):
            self.coefficients[k] = UInt64((UInt128(self.coefficients[k]) * UInt128(factor)) % UInt128(prime))


def _one[W: Int]() -> _Packed[W]:
    return _packed(SIMD[DType.uint32, W](0))


def _sub_multiple[W: Int](p: List[_Term[ModInt, W]], q: List[_Term[ModInt, W]], m: _Packed[W], c: UInt64, kind: UInt8, ring: Modulus) -> List[_Term[ModInt, W]]:
    """p - c m q for increasing term lists, leading term last."""
    var prime = ring.value
    var result = List[_Term[ModInt, W]](capacity=len(p) + len(q))
    var i = 0
    var j = 0
    while i < len(p) or j < len(q):
        var take_p: Bool
        var mq = _times(q[j].monomial, m) if j < len(q) else _one[W]()
        if i >= len(p):
            take_p = False
        elif j >= len(q):
            take_p = True
        else:
            var cmp = _compare(kind, p[i].monomial, mq)
            if cmp == 0:
                var v = (p[i].coefficient.value + prime - UInt64((UInt128(c) * UInt128(q[j].coefficient.value)) % UInt128(prime))) % prime
                if v != 0:
                    result.append(_Term[ModInt, W](p[i].monomial, ModInt(v, ring)))
                i += 1
                j += 1
                continue
            take_p = cmp < 0
        if take_p:
            result.append(_Term[ModInt, W](p[i].monomial, p[i].coefficient.copy()))
            i += 1
        else:
            var v = UInt64((UInt128(prime - c) * UInt128(q[j].coefficient.value)) % UInt128(prime))
            if v != 0:
                result.append(_Term[ModInt, W](mq, ModInt(v, ring)))
            j += 1
    return result^


def _reduce_recording[W: Int](mut dense: List[UInt64], start: Int, pivot: List[Int], matrix: _Matrix[W], values: List[List[UInt32]], mut found: _Rows, prime: UInt64, mut factors: List[Tuple[Int, UInt64]], mut budget: Budget) raises PolyError:
    """As `_f4._reduce_row_in` with entries below p^2, recording each pivot q taken v times as (q, v).

    Residues come by Barrett's method and runs of zero words are skipped
    eight at a time, as there.
    """
    var square = prime * prime
    var reciprocal = _reciprocal(prime)
    var ncols = len(dense)
    var rows = len(matrix.sources)
    var data = dense.unsafe_ptr()
    budget.charge(UInt64(ncols - start))
    var c = start
    while c < ncols:
        if c & 7 == 0 and c + 8 <= ncols:
            var block = data.unsafe_offset(c).unsafe_load[width=8]()
            if not block.ne(0).reduce_or():
                c += 8
                continue
        var v = data.unsafe_offset(c).unsafe_load()
        if v == 0:
            c += 1
            continue
        data.unsafe_offset(c).unsafe_store(UInt64(0))
        v = _residue(v, prime, reciprocal)
        if v == 0:
            c += 1
            continue
        var q = pivot.unsafe_get(c)
        if q < 0:
            found.columns.append(Int32(c))
            found.values.append(UInt32(v))
            c += 1
            continue
        factors.append((q, v))
        var g = prime - v
        if q >= rows:
            var r = q - rows
            var s = found.starts.unsafe_get(r)
            var n = found.starts.unsafe_get(r + 1) - s
            budget.charge(UInt64(n))
            var columns = found.columns.unsafe_ptr().unsafe_offset(s)
            var entries = found.values.unsafe_ptr().unsafe_offset(s)
            for k in range(1, n):
                var at = data.unsafe_offset(Int(columns.unsafe_offset(k).unsafe_load()))
                var x = at.unsafe_load() + g * UInt64(entries.unsafe_offset(k).unsafe_load())
                at.unsafe_store(x - square if x >= square else x)
        else:
            var s = matrix.starts.unsafe_get(q)
            var n = matrix.starts.unsafe_get(q + 1) - s
            var columns = matrix.entries.unsafe_ptr().unsafe_offset(s)
            var entries = values[matrix.sources.unsafe_get(q)].unsafe_ptr()
            budget.charge(UInt64(n))
            for k in range(1, n):
                var at = data.unsafe_offset(Int(columns.unsafe_offset(k).unsafe_load()))
                var x = at.unsafe_load() + g * UInt64(entries.unsafe_offset(k).unsafe_load())
                at.unsafe_store(x - square if x >= square else x)
        c += 1


comptime _FOUND = 1 << 30
"""Pivot codes from this on name a round's found rows; below it, its kept matrix rows."""

comptime _REPLAY_ENTRIES = 1 << 27
"""Matrix entries a recording keeps at most; past it the run is not replayed in lanes."""

comptime _BLOCK_WEIGHT = 4
"""Multiples per pivot, over a round's first `_BLOCK` rows, from which `_traced_round` reduces the rest in blocks: rewriting a block's multiples costs about what that many row-by-row multiples do."""


struct _ReplayRound(Movable):
    """One matrix of a traced run, kept for `_lanes_trace`.

    Rows are the matrix's pivots, then, in a round, the rows that gave new
    pivots, in turn; the final matrix keeps all its rows, the first `basis`
    being the basis rows. An operation is a kept row's reduction, then each
    new pivot's reduction by those right of it in `created` order, or in the
    final matrix a basis row's: its factors are the pivots it took, as codes
    (a kept row's number, or `_FOUND` plus a found row's), in column order,
    and its terms start at trace slot `slots[op]`, the row itself first.
    """
    var final: Bool
    var ncols: Int
    var sources: List[Int]
    var starts: List[Int]
    var entries: List[Int32]
    var pivots: Int
    var basis: Int
    var leads: List[Int32]
    var created: List[Int]
    var factors: List[Int32]
    var columns: List[Int32]
    var fstarts: List[Int]
    var slots: List[Int]
    var results: List[Int]
    var supports: List[List[Int32]]

    def __init__(out self, final: Bool, ncols: Int):
        self.final = final
        self.ncols = ncols
        self.sources = List[Int]()
        self.starts = [0]
        self.entries = List[Int32]()
        self.pivots = 0
        self.basis = 0
        self.leads = List[Int32]()
        self.created = List[Int]()
        self.factors = List[Int32]()
        self.columns = List[Int32]()
        self.fstarts = [0]
        self.slots = List[Int]()
        self.results = List[Int]()
        self.supports = List[List[Int32]]()

    def keep[W: Int](mut self, matrix: _Matrix[W], r: Int):
        for k in range(matrix.starts[r], matrix.starts[r + 1]):
            self.entries.append(matrix.entries[k])
        self.sources.append(matrix.sources[r])
        self.starts.append(len(self.entries))

    def operation(mut self, slot: Int, factors: List[Tuple[Int, UInt64]], codes: List[Int], pivot_columns: List[Int32], found: _Rows, rows: Int):
        """Record one reduction: its first slot and its factors' codes and columns."""
        self.slots.append(slot)
        for f in factors:
            var q = f[0]
            if q >= rows:
                self.factors.append(Int32(_FOUND + q - rows))
                self.columns.append(found.columns[found.starts[q - rows]])
            else:
                self.factors.append(Int32(codes[q]))
                self.columns.append(pivot_columns[q])
        self.fstarts.append(len(self.factors))


struct _Replay[W: Int](Movable):
    """What `_lanes_trace` needs of a traced run: its generators' elements and its matrices."""
    var inputs: List[List[_Packed[Self.W]]]
    """Each element the generators leave, its monomials in increasing order."""
    var generator_steps: Int
    """Trace steps the generators take, before the first round's."""
    var rounds: List[_ReplayRound]
    var entries: Int
    var complete: Bool
    """False once the recording passed `_REPLAY_ENTRIES` and stopped."""

    def __init__(out self):
        self.inputs = List[List[_Packed[Self.W]]]()
        self.generator_steps = 0
        self.rounds = List[_ReplayRound]()
        self.entries = 0
        self.complete = True


struct _BlockRecord(Movable):
    """What `_reduce_block_recording` records of each lane: the pivots it took, as (column, pivot, multiple), and the earlier lanes it took, as (lane, column, multiple of that lane's state there)."""
    var taken: List[List[Tuple[Int, Int, UInt64]]]
    var lanes: List[List[Tuple[Int, Int, UInt64]]]
    var events: Int
    """Pivot multiples taken, each by every lane nonzero there at once."""

    def __init__(out self):
        self.events = 0
        self.taken = List[List[Tuple[Int, Int, UInt64]]]()
        self.lanes = List[List[Tuple[Int, Int, UInt64]]]()
        for _ in range(_BLOCK):
            self.taken.append(List[Tuple[Int, Int, UInt64]]())
            self.lanes.append(List[Tuple[Int, Int, UInt64]]())

    def clear(mut self):
        self.events = 0
        for l in range(_BLOCK):
            self.taken[l].clear()
            self.lanes[l].clear()


def _reduce_block_recording[small: Bool, W: Int](mut dense: List[UInt64], rows: List[Int], first: Int, count: Int, pivot: List[Int], matrix: _Matrix[W], values: List[List[UInt32]], found: _Rows, mut lanes: List[_Lane], field: _Field, ring: Modulus, mut record: _BlockRecord, mut budget: Budget) raises PolyError:
    """`_f4._reduce_block`, recording in `record` every multiple each lane takes: of a pivot, or of an earlier lane's state where that lane leads."""
    var ncols = len(pivot)
    var rows_in_matrix = len(matrix.sources)
    var data = dense.unsafe_ptr()
    var start = ncols
    record.clear()
    for l in range(count):
        var r = rows[first + l]
        var s = matrix.starts[r]
        ref v = values[matrix.sources[r]]
        for k in range(matrix.starts[r + 1] - s):
            data.unsafe_offset(_BLOCK * Int(matrix.entries.unsafe_get(s + k)) + l).unsafe_store(UInt64(v.unsafe_get(k)))
        start = min(start, Int(matrix.entries[s]))
        lanes[l].columns.clear()
        lanes[l].values.clear()
    var led = SIMD[DType.bool, _BLOCK](fill=False)
    budget.charge(UInt64(ncols - start))
    for c in range(start, ncols):
        var at = data.unsafe_offset(_BLOCK * c)
        var v = at.unsafe_load[width=_BLOCK]()
        if not v.ne(0).reduce_or():
            continue
        at.unsafe_store[width=_BLOCK](_Block(0))
        v = _residues(v, field)
        var nonzero = v.ne(0)
        if not nonzero.reduce_or():
            continue
        var q = pivot.unsafe_get(c)
        if q >= 0:
            record.events += 1
            for l in range(count):
                if nonzero[l]:
                    record.taken[l].append((c, q, v[l]))
            var factor = nonzero.select(field.primes - v, _Block(0)).cast[DType.uint32]()
            var weight = UInt64(nonzero.cast[DType.uint8]().reduce_add())
            if q >= rows_in_matrix:
                var r = q - rows_in_matrix
                var s = found.starts.unsafe_get(r)
                var n = found.starts.unsafe_get(r + 1) - s
                budget.charge(UInt64(n) * weight)
                _axpy_block[small](data, found.columns.unsafe_ptr().unsafe_offset(s), found.values.unsafe_ptr().unsafe_offset(s), n, factor, field.squares)
            else:
                var s = matrix.starts.unsafe_get(q)
                var n = matrix.starts.unsafe_get(q + 1) - s
                budget.charge(UInt64(n) * weight)
                _axpy_block[small](data, matrix.entries.unsafe_ptr().unsafe_offset(s), values[matrix.sources.unsafe_get(q)].unsafe_ptr(), n, factor, field.squares)
            continue
        var leader = -1
        for l in range(count):
            if nonzero[l] and not led[l]:
                leader = l
                break
        for l in range(count):
            if nonzero[l] and (leader < 0 or l <= leader):
                lanes[l].columns.append(Int32(c))
                lanes[l].values.append(UInt32(v[l]))
        if leader < 0:
            continue
        led[leader] = True
        # The later lanes subtract their multiple of the leader from here on.
        var inverse = ModInt(v[leader], ring).inverse(budget).value
        var factors = SIMD[DType.uint32, _BLOCK](0)
        var any = False
        for l in range(leader + 1, count):
            if nonzero[l]:
                var multiple = _residue(v[l] * inverse, field.prime, field.reciprocal)
                record.lanes[l].append((leader, c, multiple))
                factors[l] = UInt32(field.prime - multiple)
                any = True
        if not any:
            continue
        budget.charge(UInt64(ncols - c))
        var factor = factors.cast[DType.uint64]()
        for k in range(c + 1, ncols):
            var to = data.unsafe_offset(_BLOCK * k)
            var x = to.unsafe_offset(leader).unsafe_load()
            if x == 0:
                continue
            var y = to.unsafe_load[width=_BLOCK]() + factor * _Block(UInt64(UInt32(_residue(x, field.prime, field.reciprocal))))
            comptime if small:
                to.unsafe_store[width=_BLOCK](y)
            else:
                to.unsafe_store[width=_BLOCK](y.ge(field.squares).select(y - field.squares, y))


def _merged(a: List[Tuple[Int, Int, UInt64]], s: UInt64, b: List[Tuple[Int, Int, UInt64]], prime: UInt64, mut budget: Budget) raises PolyError -> List[Tuple[Int, Int, UInt64]]:
    """a + s b for sparse rows of (column, pivot, coefficient) by increasing column."""
    budget.charge(UInt64(len(a) + len(b)))
    var result = List[Tuple[Int, Int, UInt64]](capacity=len(a) + len(b))
    var i = 0
    var j = 0
    while i < len(a) or j < len(b):
        if j == len(b) or (i < len(a) and a[i][0] < b[j][0]):
            result.append(a[i])
            i += 1
            continue
        var add = UInt64((UInt128(s) * UInt128(b[j][2])) % UInt128(prime))
        if i < len(a) and a[i][0] == b[j][0]:
            result.append((a[i][0], a[i][1], (a[i][2] + add) % prime))
            i += 1
        else:
            result.append((b[j][0], b[j][1], add))
        j += 1
    return result^


def _block_factors(record: _BlockRecord, lanes: List[_Lane], count: Int, found_base: Int, rows: Int, prime: UInt64, mut budget: Budget) raises PolyError -> List[List[Tuple[Int, UInt64]]]:
    """Each nonzero lane's pivots and multiples as reducing its row alone, after the earlier rows, records them: (pivot, multiple) by increasing column.

    A lane i that leads in column c_i left its state D_i there in the lanes
    after it; it ends as U_i = D_i - T_i, with T_i its later multiples of
    pivots and of the states of the lanes that lead after c_i. Row j alone
    takes U_i, as the found row it becomes, where the block takes D_i: so
    R_j - U_j = sum of j's pivot multiples + sum over the lanes i it took of
    its multiple of (U_i + T_i). Pivots and found rows lead in distinct
    columns, so these coefficients are the ones the row alone records, a
    found row's as the multiple of its monic form (times its leading value).
    Sparse rows hold (column, pivot, coefficient); U_i stands in its leading
    column with pivot -1 - i until the multiple of its monic form is taken.
    """
    var leads = List[Int](length=count, fill=-1)
    var lead_values = List[UInt64](length=count, fill=0)
    var numbers = List[Int](length=count, fill=-1)
    var next = found_base
    for l in range(count):
        if len(lanes[l].columns):
            leads[l] = Int(lanes[l].columns[0])
            lead_values[l] = UInt64(lanes[l].values[0])
            numbers[l] = next
            next += 1
    var tails = List[List[Tuple[Int, Int, UInt64]]]()
    var result = List[List[Tuple[Int, UInt64]]]()
    for l in range(count):
        ref taken = record.taken[l]
        var all = List[Tuple[Int, Int, UInt64]](capacity=len(taken))
        var tail = List[Tuple[Int, Int, UInt64]]()
        for t in taken:
            all.append(t)
            if leads[l] >= 0 and t[0] > leads[l]:
                tail.append(t)
        for t in record.lanes[l]:
            var i = t[0]
            var c = t[1]
            # s (U_i + T_i), U_i in its leading column.
            var part = List[Tuple[Int, Int, UInt64]](capacity=len(tails[i]) + 1)
            part.append((c, -1 - i, UInt64(1)))
            for e in tails[i]:
                part.append(e)
            all = _merged(all, t[2], part, prime, budget)
            if leads[l] >= 0 and c > leads[l]:
                tail = _merged(tail, t[2], part, prime, budget)
        tails.append(tail^)
        var factors = List[Tuple[Int, UInt64]]()
        if leads[l] >= 0:
            factors.reserve(len(all))
            for e in all:
                if e[2] == 0:
                    continue
                if e[1] < 0:
                    var lane = -1 - e[1]
                    factors.append((rows + numbers[lane], UInt64((UInt128(e[2]) * UInt128(lead_values[lane])) % UInt128(prime))))
                else:
                    factors.append((e[1], e[2]))
        result.append(factors^)
    return result^

struct _Round[W: Int](Movable):
    """Trace numbers of the engine's elements and of the current round's found rows."""
    var ids: List[Int]
    var found: List[Int]

    def __init__(out self):
        self.ids = List[Int]()
        self.found = List[Int]()


def _record_factors[W: Int](mut traced: _Traced[W], factors: List[Tuple[Int, UInt64]], matrix: _Matrix[W], order: List[Int], engine: _Engine[ModInt, W, False], numbers: _Round[W], prime: UInt64, mut budget: Budget) raises PolyError:
    """Each recorded (q, v) as minus v times pivot q: a matrix row's element and multiplier, or a found row."""
    var rows = len(matrix.sources)
    for f in factors:
        budget.charge(1)
        var q = f[0]
        var minus = prime - f[1]
        if q >= rows:
            traced.add(numbers.found[q - rows], _one[W](), minus)
        else:
            var source = matrix.sources[q]
            var lead = matrix.columns.monomials[order[Int(matrix.entries[matrix.starts[q]])]]
            traced.add(numbers.ids[source], _div(lead, engine.leads[source]), minus)


def _traced_round[W: Int, record: Bool](mut engine: _Engine[ModInt, W, False], store: _Store[W], prime: UInt64, ring: Modulus, mut traced: _Traced[W], mut numbers: _Round[W], mut replay: _Replay[W], mut budget: Budget) raises PolyError -> List[List[_Term[ModInt, W]]]:
    """`_f4._round`, row by row, recording every new row as a step; the new elements, monic, in insertion order, their numbers appended to `numbers.ids`. With `record`, the round's matrix goes to `replay`."""
    var degree = engine.pairs[0].lcm.degree
    var matrix = _Matrix[W]()
    var seen = Dict[Int, Bool]()
    var stride = len(engine.polys)
    while len(engine.pairs) and engine.pairs[0].lcm.degree == degree:
        var pair = engine.select(budget)
        matrix.columns.reserve(1, budget)
        var lead = Int(matrix.columns.find(pair.lcm, _raw_hash(pair.lcm)))
        for source in [pair.first, pair.second]:
            budget.charge(1)
            var key = lead * stride + source
            if key in seen:
                continue
            seen[key] = True
            matrix.add_row(engine, store, source, _div(pair.lcm, engine.leads[source]), budget)
    var pair_rows = len(matrix.sources)
    matrix.preprocess(engine, store, budget)
    var nrows = len(matrix.sources)
    var order = matrix.number(engine.kind, budget)
    var ncols = len(order)
    _allocation[Int](ncols, budget)
    var pivot = List[Int](length=ncols, fill=-1)
    for r in range(pair_rows, nrows):
        pivot[Int(matrix.entries[matrix.starts[r]])] = r
    var reduce = List[Int]()
    for r in range(pair_rows):
        var c = Int(matrix.entries[matrix.starts[r]])
        if pivot[c] < 0:
            pivot[c] = r
        else:
            reduce.append(r)
    # With `record`: the pivots' numbers among the kept rows, each pivot's
    # column, and the round's operations.
    var codes = List[Int]()
    var pivot_columns = List[Int32]()
    var kept = List[Int]()
    var round = _ReplayRound(False, ncols)
    comptime if record:
        if replay.complete:
            codes = List[Int](length=nrows, fill=-1)
            pivot_columns = List[Int32](length=nrows, fill=-1)
            for r in range(nrows):
                var c = Int(matrix.entries[matrix.starts[r]])
                pivot_columns[r] = Int32(c)
                if pivot[c] == r:
                    codes[r] = len(round.sources)
                    round.keep(matrix, r)
            round.pivots = len(round.sources)
    var found = _Rows()
    numbers.found = List[Int]()
    var reciprocal = _reciprocal(prime)
    _allocation[UInt64](ncols, budget)
    var dense = List[UInt64](length=ncols, fill=0)
    var factors = List[Tuple[Int, UInt64]]()
    # Rows `_BLOCK` at a time, as `_f4._eliminate` reduces them, each lane's
    # multiples then rewritten as reducing its row alone records them
    # (`_block_factors`): each nonzero remainder becomes a monic pivot. Both
    # ways record the same factors, so a round may switch: its first `_BLOCK`
    # rows go row by row, and the rest in blocks only when those rows took
    # each pivot they used `_BLOCK_WEIGHT` times or more on average; with
    # fewer lanes per multiple the rewriting costs more than blocks save.
    var field = _Field(prime)
    var block = List[UInt64]()
    var lanes = List[_Lane](capacity=_BLOCK)
    var block_record = _BlockRecord()
    var blocked = False
    var taken = 0
    var distinct = Dict[Int, Bool]()
    var first = 0
    while first < len(reduce):
        if first == _BLOCK:
            blocked = taken >= _BLOCK_WEIGHT * len(distinct)
            if blocked:
                break
        var r = reduce[first]
        first += 1
        _scatter(dense, matrix, store.values, r, 0)
        var before = len(found.columns)
        factors.clear()
        _reduce_recording(dense, Int(matrix.entries[matrix.starts[r]]), pivot, matrix, store.values, found, prime, factors, budget)
        if first <= _BLOCK:
            taken += len(factors)
            for f in factors:
                distinct[f[0]] = True
        if len(found.columns) == before:
            continue
        var inverse = ModInt(UInt64(found.values[before]), ring).inverse(budget).value
        for k in range(before, len(found.columns)):
            found.values[k] = UInt32(_residue(UInt64(found.values[k]) * inverse, prime, reciprocal))
        var source = matrix.sources[r]
        var lead = matrix.columns.monomials[order[Int(matrix.entries[matrix.starts[r]])]]
        comptime if record:
            if replay.complete:
                kept.append(r)
                round.leads.append(found.columns[before])
                round.operation(len(traced.sources), factors, codes, pivot_columns, found, nrows)
        traced.add(numbers.ids[source], _div(lead, engine.leads[source]), 1)
        _record_factors(traced, factors, matrix, order, engine, numbers, prime, budget)
        traced.scale_open(inverse, prime)
        numbers.found.append(traced.close())
        pivot[Int(found.columns[before])] = nrows + len(found.starts) - 1
        found.starts.append(len(found.columns))
    if blocked:
        _allocation[UInt64](_product(ncols, _BLOCK), budget)
        block = List[UInt64](length=_BLOCK * ncols, fill=0)
        for _ in range(_BLOCK):
            lanes.append(_Lane())
    while first < len(reduce):
        var count = min(_BLOCK, len(reduce) - first)
        if prime < 65536:
            _reduce_block_recording[True](block, reduce, first, count, pivot, matrix, store.values, found, lanes, field, ring, block_record, budget)
        else:
            _reduce_block_recording[False](block, reduce, first, count, pivot, matrix, store.values, found, lanes, field, ring, block_record, budget)
        var lane_factors = _block_factors(block_record, lanes, count, len(found.starts) - 1, nrows, prime, budget)
        for l in range(count):
            ref lane = lanes[l]
            if len(lane.columns) == 0:
                continue
            var r = reduce[first + l]
            _allocation[UInt32](2 * len(lane.columns), budget)
            var before = len(found.columns)
            var inverse = ModInt(UInt64(lane.values[0]), ring).inverse(budget).value
            for k in range(len(lane.columns)):
                found.columns.append(lane.columns[k])
                found.values.append(UInt32(_residue(UInt64(lane.values[k]) * inverse, prime, reciprocal)))
            factors.clear()
            for f in lane_factors[l]:
                factors.append(f)
            var source = matrix.sources[r]
            var lead = matrix.columns.monomials[order[Int(matrix.entries[matrix.starts[r]])]]
            comptime if record:
                if replay.complete:
                    kept.append(r)
                    round.leads.append(found.columns[before])
                    round.operation(len(traced.sources), factors, codes, pivot_columns, found, nrows)
            traced.add(numbers.ids[source], _div(lead, engine.leads[source]), 1)
            _record_factors(traced, factors, matrix, order, engine, numbers, prime, budget)
            traced.scale_open(inverse, prime)
            numbers.found.append(traced.close())
            pivot[Int(found.columns[before])] = nrows + len(found.starts) - 1
            found.starts.append(len(found.columns))
        first += count
    comptime if record:
        if replay.complete:
            for r in kept:
                round.keep(matrix, r)
    # The new pivots by decreasing leading column, each reduced by those right of it.
    var created = List[Int]()
    for i in range(len(found.starts) - 1):
        var j = len(created)
        created.append(i)
        while j > 0 and found.columns[found.starts[created[j - 1]]] < found.columns[found.starts[i]]:
            created[j] = created[j - 1]
            j -= 1
        created[j] = i
    var results = List[Int](capacity=len(created))
    for i in created:
        var lead = Int(found.columns[found.starts[i]])
        for k in range(found.starts[i] + 1, found.starts[i + 1]):
            dense[Int(found.columns[k])] = UInt64(found.values[k])
        found.columns.append(Int32(lead))
        found.values.append(UInt32(1))
        factors.clear()
        _reduce_recording(dense, lead + 1, pivot, matrix, store.values, found, prime, factors, budget)
        comptime if record:
            if replay.complete:
                round.created.append(i)
                round.operation(len(traced.sources), factors, codes, pivot_columns, found, nrows)
        traced.add(numbers.found[i], _one[W](), 1)
        _record_factors(traced, factors, matrix, order, engine, numbers, prime, budget)
        numbers.found.append(traced.close())
        pivot[lead] = nrows + len(found.starts) - 1
        results.append(len(found.starts) - 1)
        found.starts.append(len(found.columns))
    var inserted = _insertion(found, results, matrix.columns.monomials, order, budget)
    results = _inserted(results, inserted)
    comptime if record:
        if replay.complete:
            for i in results:
                round.results.append(i)
                var support = List[Int32](capacity=found.starts[i + 1] - found.starts[i])
                for k in range(found.starts[i], found.starts[i + 1]):
                    support.append(found.columns[k])
                round.supports.append(support^)
            replay.entries += len(round.entries)
            if replay.entries > _REPLAY_ENTRIES:
                replay.complete = False
                replay.rounds = List[_ReplayRound]()
            else:
                replay.rounds.append(round^)
    var fresh = List[List[_Term[ModInt, W]]]()
    for i in results:
        var n = found.starts[i + 1] - found.starts[i]
        _allocation[_Term[ModInt, W]](n, budget)
        var terms = List[_Term[ModInt, W]](capacity=n)
        for k in range(found.starts[i + 1] - 1, found.starts[i] - 1, -1):
            terms.append(_Term[ModInt, W](matrix.columns.monomials[order[Int(found.columns[k])]], ModInt(UInt64(found.values[k]), ring)))
        fresh.append(terms^)
        numbers.ids.append(numbers.found[i])
    return fresh^


struct _TracedBasis[W: Int](Movable):
    var basis: List[List[_Term[ModInt, Self.W]]]
    var traced: _Traced[Self.W]
    var sources: List[Int]

    def __init__(out self, var basis: List[List[_Term[ModInt, Self.W]]], var traced: _Traced[Self.W], var sources: List[Int]):
        self.basis = basis^
        self.traced = traced^
        self.sources = sources^


def _f4_traced[W: Int](generators: List[List[_Term[ModInt, W]]], kind: UInt8, ring: Modulus, mut budget: Budget) raises PolyError -> _TracedBasis[W]:
    """The reduced monic basis, in decreasing leading-monomial order, with the trace of every polynomial the run makes."""
    var replay = _Replay[W]()
    return _f4_traced_run[W, False](generators, kind, ring, replay, budget)


def _f4_recorded[W: Int](generators: List[List[_Term[ModInt, W]]], kind: UInt8, ring: Modulus, mut replay: _Replay[W], mut budget: Budget) raises PolyError -> _TracedBasis[W]:
    """As `_f4_traced`, keeping in `replay` what `_lanes_trace` needs to trace other primes the same way."""
    replay = _Replay[W]()
    return _f4_traced_run[W, True](generators, kind, ring, replay, budget)


def _traced_generators[W: Int](generators: List[List[_Term[ModInt, W]]], kind: UInt8, ring: Modulus, mut engine: _Engine[ModInt, W, False], mut traced: _Traced[W], mut numbers: _Round[W], mut budget: Budget) raises PolyError:
    """The generators, each top-reduced by the elements so far, made monic and inserted, recording a step each."""
    var prime = ring.value
    for i in range(len(generators)):
        budget.charge(1)
        if len(generators[i]) == 0:
            continue
        var p = generators[i].copy()
        traced.add(i, _one[W](), 1)
        while len(p):
            var lead = p[len(p) - 1].monomial
            var divisor = -1
            for e in engine.active:
                budget.charge(1)
                if _divides(engine.leads[e], lead):
                    divisor = e
                    break
            if divisor < 0:
                break
            var c = p[len(p) - 1].coefficient.value
            var m = _div(lead, engine.leads[divisor])
            budget.charge(UInt64(len(engine.polys[divisor])))
            p = _sub_multiple(p, engine.polys[divisor], m, c, kind, ring)
            traced.add(numbers.ids[divisor], m, prime - c)
        if len(p) == 0:
            # A generator in the ideal of the earlier ones: its step is never used.
            _ = traced.close()
            continue
        traced.scale_open(ModInt(p[len(p) - 1].coefficient.value, ring).inverse(budget).value, prime)
        numbers.ids.append(traced.close())
        engine.insert(p^, budget)


def _f4_traced_run[W: Int, record: Bool](generators: List[List[_Term[ModInt, W]]], kind: UInt8, ring: Modulus, mut replay: _Replay[W], mut budget: Budget) raises PolyError -> _TracedBasis[W]:
    var prime = ring.value
    var traced = _Traced[W](len(generators))
    var numbers = _Round[W]()
    var engine = _Engine[ModInt, W, False](kind, ring)
    _traced_generators(generators, kind, ring, engine, traced, numbers, budget)
    comptime if record:
        for p in engine.polys:
            var monomials = List[_Packed[W]](capacity=len(p))
            for t in p:
                monomials.append(t.monomial)
            replay.inputs.append(monomials^)
        replay.generator_steps = len(traced.starts) - 1
    var store = _Store[W]()
    while len(engine.pairs):
        budget.charge(1)
        engine.sort_active(budget)
        store.extend(engine, budget)
        store.release(engine)
        var fresh = _traced_round[W, record](engine, store, prime, ring, traced, numbers, replay, budget)
        for k in range(len(fresh)):
            var p = List[_Term[ModInt, W]]()
            swap(p, fresh[k])
            engine.insert(p^, budget, heap=False)
        engine.heapify_pairs(budget)
    store.extend(engine, budget)
    store.release(engine)
    # The reduced basis: each active element's tail reduced by the final matrix.
    engine.sort_active(budget)
    var matrix = _Matrix[W]()
    for position in range(len(engine.active) - 1, -1, -1):
        matrix.add_row(engine, store, engine.active[position], _one[W](), budget)
    var basis_rows = len(matrix.sources)
    matrix.preprocess(engine, store, budget)
    var nrows = len(matrix.sources)
    var order = matrix.number(engine.kind, budget)
    var ncols = len(order)
    _allocation[Int](ncols, budget)
    var pivot = List[Int](length=ncols, fill=-1)
    for r in range(nrows):
        pivot[Int(matrix.entries[matrix.starts[r]])] = r
    var round = _ReplayRound(True, ncols)
    var codes = List[Int]()
    var pivot_columns = List[Int32]()
    comptime if record:
        if replay.complete:
            for r in range(nrows):
                codes.append(r)
                pivot_columns.append(matrix.entries[matrix.starts[r]])
                round.keep(matrix, r)
            round.pivots = nrows
            round.basis = basis_rows
    var found = _Rows()
    numbers.found = List[Int]()
    _allocation[UInt64](ncols, budget)
    var dense = List[UInt64](length=ncols, fill=0)
    var factors = List[Tuple[Int, UInt64]]()
    var sources = List[Int](capacity=basis_rows)
    for r in range(basis_rows):
        var first = Int(matrix.entries[matrix.starts[r]])
        _scatter(dense, matrix, store.values, r, 1)
        found.columns.append(Int32(first))
        found.values.append(store.values[matrix.sources[r]][0])
        factors.clear()
        _reduce_recording(dense, first + 1, pivot, matrix, store.values, found, prime, factors, budget)
        found.starts.append(len(found.columns))
        comptime if record:
            if replay.complete:
                round.operation(len(traced.sources), factors, codes, pivot_columns, found, nrows)
        traced.add(numbers.ids[matrix.sources[r]], _one[W](), 1)
        _record_factors(traced, factors, matrix, order, engine, numbers, prime, budget)
        sources.append(traced.close())
    comptime if record:
        if replay.complete:
            replay.entries += len(round.entries)
            if replay.entries > _REPLAY_ENTRIES:
                replay.complete = False
                replay.rounds = List[_ReplayRound]()
            else:
                replay.rounds.append(round^)
    var basis = List[List[_Term[ModInt, W]]](capacity=basis_rows)
    for r in range(basis_rows):
        var n = found.starts[r + 1] - found.starts[r]
        _allocation[_Term[ModInt, W]](n, budget)
        var terms = List[_Term[ModInt, W]](capacity=n)
        for k in range(found.starts[r + 1] - 1, found.starts[r] - 1, -1):
            terms.append(_Term[ModInt, W](matrix.columns.monomials[order[Int(found.columns[k])]], ModInt(UInt64(found.values[k]), ring)))
        basis.append(terms^)
    return _TracedBasis[W](basis^, traced^, sources^)


struct _LaneTrace(Movable):
    """The trace coefficients of `_LANES` primes: slot k's residue modulo lane l's prime is values[_LANES * k + l]; lanes 0 to live - 1 have them."""
    var values: List[UInt32]
    var live: Int

    def __init__(out self, slots: Int):
        self.values = List[UInt32](length=_LANES * slots, fill=0)
        self.live = 0


@always_inline
def _store_slot(mut trace: _LaneTrace, slot: Int, v: _Lanes):
    trace.values.unsafe_ptr().unsafe_offset(_LANES * slot).unsafe_store[width=_LANES](v.cast[DType.uint32]())


def _scatter_row(mut dense: List[UInt64], round: _ReplayRound, values: List[List[UInt32]], row: Int, skip: Int):
    """Write a kept row's residues in every lane, after its first `skip` entries, into the dense accumulator."""
    var s = round.starts[row]
    var source = values[round.sources[row]].unsafe_ptr()
    var data = dense.unsafe_ptr()
    for k in range(skip, round.starts[row + 1] - s):
        var at = data.unsafe_offset(_LANES * Int(round.entries.unsafe_get(s + k)))
        at.unsafe_store[width=_LANES](source.unsafe_offset(_LANES * k).unsafe_load[width=_LANES]().cast[DType.uint64]())


def _reduce_traced_lanes(mut dense: List[UInt64], start: Int, pivot: List[Int], round: _ReplayRound, values: List[List[UInt32]], mut found: _LaneRows, mut lanes: _Primes, op: Int, mut trace: _LaneTrace, mut budget: Budget) raises PolyError:
    """As `_reduce_recording` in every lane, taking the pivots operation `op` of the recorded run took.

    Each recorded factor's slot gets minus the entry it cancels, zero in a
    lane where the entry is zero. A lane with a nonzero entry where the
    recorded run took no pivot departs; entries in columns with no pivot go
    to `found`.
    """
    var ncols = round.ncols
    var data = dense.unsafe_ptr()
    var live = lanes.mask()
    var f = round.fstarts[op]
    var end = round.fstarts[op + 1]
    var slot = round.slots[op] + 1 - f
    budget.charge(UInt64(_LANES * (ncols - start)))
    for c in range(start, ncols):
        var at = data.unsafe_offset(_LANES * c)
        var v = at.unsafe_load[width=_LANES]()
        if not v.ne(0).reduce_or():
            continue
        at.unsafe_store[width=_LANES](_Lanes(0))
        v = _mod_lanes(v, lanes) & live
        if not v.ne(0).reduce_or():
            continue
        var q = pivot.unsafe_get(c)
        if q < 0:
            found.columns.append(Int32(c))
            var n = len(found.values)
            found.values.resize(n + _LANES, 0)
            found.values.unsafe_ptr().unsafe_offset(n).unsafe_store[width=_LANES](v.cast[DType.uint32]())
            continue
        while f < end and Int(round.columns.unsafe_get(f)) < c:
            _store_slot(trace, slot + f, _Lanes(0))
            f += 1
        if f == end or Int(round.columns.unsafe_get(f)) != c:
            # The recorded run had a zero here: the lanes that have not depart.
            for l in range(lanes.live):
                if v[l] != 0:
                    lanes.depart(l)
                    break
            live = lanes.mask()
            continue
        var minus = v.eq(0).select(_Lanes(0), lanes.primes - v)
        _store_slot(trace, slot + f, minus)
        f += 1
        var factor = minus.cast[DType.uint32]()
        if q >= _FOUND:
            var r = q - _FOUND
            var s = found.starts.unsafe_get(r)
            var n = found.starts.unsafe_get(r + 1) - s
            budget.charge(UInt64(_LANES * n))
            _axpy_lanes(data, found.columns.unsafe_ptr().unsafe_offset(s), found.values.unsafe_ptr().unsafe_offset(_LANES * s), n, factor, lanes.squares)
        else:
            var s = round.starts.unsafe_get(q)
            var n = round.starts.unsafe_get(q + 1) - s
            budget.charge(UInt64(_LANES * n))
            _axpy_lanes(data, round.entries.unsafe_ptr().unsafe_offset(s), values[round.sources.unsafe_get(q)].unsafe_ptr(), n, factor, lanes.squares)
    while f < end:
        _store_slot(trace, slot + f, _Lanes(0))
        f += 1


def _lane_results(found: _LaneRows, round: _ReplayRound, mut values: List[List[UInt32]], mut lanes: _Primes, mut budget: Budget) raises PolyError:
    """Append the round's new elements' residues on their recorded supports; a lane departs when one leaves its support or loses its leading entry."""
    for t in range(len(round.results)):
        var i = round.results[t]
        ref support = round.supports[t]
        _allocation[UInt32](_product(len(support), _LANES), budget)
        var residues = List[UInt32](length=_LANES * len(support), fill=0)
        var j = 0
        for k in range(found.starts[i], found.starts[i + 1]):
            budget.charge(UInt64(_LANES))
            while j < len(support) and support[j] < found.columns[k]:
                j += 1
            if j == len(support) or support[j] != found.columns[k]:
                for l in range(lanes.live):
                    if found.values[_LANES * k + l] != 0:
                        lanes.depart(l)
                        break
                continue
            for l in range(_LANES):
                residues[_LANES * j + l] = found.values[_LANES * k + l]
        for l in range(lanes.live):
            if residues[l] == 0:
                lanes.depart(l)
                break
        values.append(residues^)


def _lanes_trace[W: Int](replay: _Replay[W], run: _TracedBasis[W], generators: List[List[List[_Term[ModInt, W]]]], kind: UInt8, rings: List[Modulus], mut budget: Budget) raises PolyError -> _LaneTrace:
    """The trace coefficients of `_LANES` primes below 2^31 at once, one per SIMD lane, over the structure of the recorded run.

    Each lane's generators are top-reduced as the recorded run's were, and
    must take the same steps; the rounds then follow the recorded matrices:
    each operation takes the recorded pivots, in the recorded order, so
    every lane's trace has the recorded sources and monomials, with a zero
    coefficient where the lane's entry vanishes. A lane departs, with the
    lanes after it, where its run leaves the recorded one otherwise; rows
    that reduced to zero in the recorded run are not reduced again.
    """
    ref traced = run.traced
    var slots = len(traced.sources)
    var trace = _LaneTrace(slots)
    var lanes = _Primes(rings)
    var values = List[List[UInt32]]()
    var generator_slots = traced.starts[replay.generator_steps]
    for l in range(_LANES):
        var engine = _Engine[ModInt, W, False](kind, rings[l])
        var own = _Traced[W](len(generators[l]))
        var numbers = _Round[W]()
        _traced_generators(generators[l], kind, rings[l], engine, own, numbers, budget)
        var follows = len(own.starts) == replay.generator_steps + 1 and len(own.sources) == generator_slots and len(engine.polys) == len(replay.inputs)
        if follows:
            for k in range(len(own.starts)):
                if own.starts[k] != traced.starts[k]:
                    follows = False
                    break
        if follows:
            for k in range(generator_slots):
                budget.charge(1)
                if own.sources[k] != traced.sources[k] or not _equal(own.monomials[k], traced.monomials[k]):
                    follows = False
                    break
        if follows and l == 0:
            for i in range(len(replay.inputs)):
                _allocation[UInt32](_product(len(replay.inputs[i]), _LANES), budget)
                values.append(List[UInt32](length=_LANES * len(replay.inputs[i]), fill=0))
        if follows:
            for i in range(len(engine.polys)):
                ref p = engine.polys[i]
                ref monomials = replay.inputs[i]
                var n = len(monomials)
                var j = 0
                for t in p:
                    budget.charge(1)
                    while j < n and not _equal(monomials[j], t.monomial):
                        j += 1
                    if j == n:
                        follows = False
                        break
                    values[i][_LANES * (n - 1 - j) + l] = UInt32(t.coefficient.value)
                    j += 1
                if not follows:
                    break
        if not follows:
            lanes.depart(l)
            break
        for k in range(generator_slots):
            trace.values[_LANES * k + l] = UInt32(own.coefficients[k])
    if lanes.live == 0:
        return trace^
    var width = 0
    for s in range(len(replay.rounds)):
        width = max(width, replay.rounds[s].ncols)
    _allocation[UInt64](_product(width, _LANES), budget)
    var dense = List[UInt64](length=_LANES * width, fill=0)
    for s in range(len(replay.rounds)):
        budget.charge(1)
        ref round = replay.rounds[s]
        var ncols = round.ncols
        var nrows = len(round.sources)
        _allocation[Int](ncols, budget)
        var pivot = List[Int](length=ncols, fill=-1)
        for r in range(round.pivots):
            pivot[Int(round.entries[round.starts[r]])] = r
        var found = _LaneRows()
        if round.final:
            for r in range(round.basis):
                var first = Int(round.entries[round.starts[r]])
                _scatter_row(dense, round, values, r, 1)
                _store_slot(trace, round.slots[r], _Lanes(1))
                _reduce_traced_lanes(dense, first + 1, pivot, round, values, found, lanes, r, trace, budget)
                found = _LaneRows()
                if lanes.live == 0:
                    break
            trace.live = lanes.live
            return trace^
        var kept = nrows - round.pivots
        for i in range(kept):
            var r = round.pivots + i
            _scatter_row(dense, round, values, r, 0)
            var before = len(found.columns)
            _reduce_traced_lanes(dense, Int(round.entries[round.starts[r]]), pivot, round, values, found, lanes, i, trace, budget)
            _allocation[UInt32](_product(len(found.columns) - before, 2 * _LANES), budget)
            # Each lane that follows the run leads in the recorded column.
            var lead = round.leads[i]
            var at = before
            while at < len(found.columns) and found.columns[at] < lead:
                at += 1
            var present = at < len(found.columns) and found.columns[at] == lead
            for l in range(lanes.live):
                var leads_there = present and found.values[_LANES * at + l] != 0
                for k in range(before, at):
                    if found.values[_LANES * k + l] != 0:
                        leads_there = False
                if not leads_there:
                    lanes.depart(l)
                    break
            if lanes.live == 0:
                return trace^
            if at > before:
                var n = len(found.columns) - at
                for k in range(n):
                    found.columns[before + k] = found.columns[at + k]
                    for l in range(_LANES):
                        found.values[_LANES * (before + k) + l] = found.values[_LANES * (at + k) + l]
                found.columns.resize(before + n, 0)
                found.values.resize(_LANES * (before + n), 0)
            # Monic in every lane that follows, and the step scaled to match.
            var inverse = _Lanes(0)
            for l in range(lanes.live):
                inverse[l] = ModInt(UInt64(found.values[_LANES * before + l]), lanes.rings[l]).inverse(budget).value
            var p = found.values.unsafe_ptr()
            for k in range(before, len(found.columns)):
                var x = p.unsafe_offset(_LANES * k).unsafe_load[width=_LANES]().cast[DType.uint64]()
                p.unsafe_offset(_LANES * k).unsafe_store[width=_LANES](_mod_lanes(x * inverse, lanes).cast[DType.uint32]())
            _store_slot(trace, round.slots[i], _Lanes(1))
            var first = round.slots[i]
            var last = first + 1 + round.fstarts[i + 1] - round.fstarts[i]
            var t = trace.values.unsafe_ptr()
            budget.charge(UInt64(last - first))
            for k in range(first, last):
                var x = t.unsafe_offset(_LANES * k).unsafe_load[width=_LANES]().cast[DType.uint64]()
                t.unsafe_offset(_LANES * k).unsafe_store[width=_LANES](_mod_lanes(x * inverse, lanes).cast[DType.uint32]())
            pivot[Int(lead)] = _FOUND + len(found.starts) - 1
            found.starts.append(len(found.columns))
        # The new pivots by decreasing leading column, each reduced by those right of it.
        var data = dense.unsafe_ptr()
        for j in range(len(round.created)):
            var i = round.created[j]
            var op = kept + j
            var lead = Int(found.columns[found.starts[i]])
            for k in range(found.starts[i] + 1, found.starts[i + 1]):
                data.unsafe_offset(_LANES * Int(found.columns[k])).unsafe_store[width=_LANES](found.values.unsafe_ptr().unsafe_offset(_LANES * k).unsafe_load[width=_LANES]().cast[DType.uint64]())
            var before = len(found.columns)
            found.columns.append(Int32(lead))
            for _ in range(_LANES):
                found.values.append(UInt32(1))
            _store_slot(trace, round.slots[op], _Lanes(1))
            _reduce_traced_lanes(dense, lead + 1, pivot, round, values, found, lanes, op, trace, budget)
            _allocation[UInt32](_product(len(found.columns) - before, 2 * _LANES), budget)
            pivot[lead] = _FOUND + len(found.starts) - 1
            found.starts.append(len(found.columns))
        if lanes.live == 0:
            return trace^
        _lane_results(found, round, values, lanes, budget)
        if lanes.live == 0:
            return trace^
    trace.live = lanes.live
    return trace^
