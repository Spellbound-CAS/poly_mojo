"""Fixed-coordinate scopes survive checked string and fragmented archives."""
from poly import cad, dense, algebraic, io
from poly._archive_meta import _Archive
from poly._algebraic_io import _Writer, _Sink
from poly._result_io import _tag
from polynomials import Budget, PolyError, PolyErrorKind, PolynomialTextSource
from std.testing import assert_true, assert_equal


@fieldwise_init
struct Source(PolynomialTextSource):
    var text: String
    var position: Int
    def read(mut self, mut budget: Budget) raises PolyError -> Optional[String]:
        budget.charge(1)
        if self.position == self.text.byte_length():
            return None
        var chunk = String(self.text[byte=self.position:self.position+1])
        self.position += 1
        return chunk^


def unchecked[T: _Archive](value: T) raises PolyError -> String:
    var budget = Budget()
    var sink = _Sink(String())
    var writer = _Writer()
    writer.line("POLYNOMIALS-VALUE 1", sink, budget)
    writer.line(_tag[T](budget), sink, budget)
    value._archive_write(writer, sink, budget, 0)
    writer.flush(sink, budget)
    var result = String()
    swap(result, sink.text)
    return result^


def roundtrip[T: _Archive](value: T) raises:
    var text = io.dumps(value)
    var restored = io.loads[T](text)
    var budget = Budget()
    assert_true(restored._archive_check(budget))
    assert_equal(io.dumps(restored), text)
    var fragmented = io.load[T](Source(text.copy(), 0))
    assert_true(fragmented._archive_check(budget))
    assert_equal(io.dumps(fragmented), text)


def main() raises:
    var p = [(0, -1), (1, 1)]
    var result = cad.decompose([p.copy()], fixed=[(0, 2)])
    roundtrip(result)
    roundtrip(cad.solve([cad.positive(p)], fixed=[(0, 2)]))
    roundtrip(cad.solve([cad.negative(p)], fixed=[(0, 2)]))
    roundtrip(cad.cover([cad.positive(p)], fixed=[(0, 2)]))
    roundtrip(cad.cover([cad.negative(p)], fixed=[(0, 2)]))
    var a = dense.root_between([-2, 0, 1], 1, 2)
    var aa = [(0, algebraic.negative(a)), (1, algebraic.asreal(1))]
    roundtrip(cad.decompose([aa.copy()], fixed=[(0, a.copy())]))
    var restored = io.loads[type_of(result)](io.dumps(result))
    assert_equal(cad.fixed_coordinates(restored)[0][0], 0)
    assert_true(cad.fixed_coordinates(restored)[0][1] == algebraic.asreal(2))
    assert_equal(cad.cells(restored)[0].signs, [1])
    var text = io.dumps(result)
    var budget = Budget()
    _ = io.load[type_of(result)](Source(text.copy(), 0), budget)
    var spent = UInt64.MAX - budget.remaining()
    var exact_budget = Budget(spent)
    _ = io.load[type_of(result)](Source(text.copy(), 0), exact_budget)
    assert_equal(exact_budget.remaining(), UInt64(0))
    var rejected = False
    try:
        _ = io.load[type_of(result)](Source(text.copy(), 0), max_steps=spent-1)
    except e:
        rejected = e.kind == PolyErrorKind.LIMIT
    assert_true(rejected)
    for malformed in [String(text, "0\n"), String(text[byte=0:text.byte_length()-1])]:
        rejected = False
        try:
            _ = io.loads[type_of(result)](malformed)
        except e:
            rejected = e.kind == PolyErrorKind.ENCODING
        assert_true(rejected)
    # Bypass only the writer check: loading still binds the selected conjugate.
    var wrong = cad.decompose([p.copy()], fixed=[(0, a.copy())])
    wrong._data[].fixed[0][1] = algebraic.negative(a)
    rejected = False
    try:
        _ = io.loads[type_of(wrong)](unchecked(wrong))
    except e:
        rejected = e.kind == PolyErrorKind.ENCODING
    assert_true(rejected)
    rejected = False
    try:
        _ = io.dumps(wrong)
    except e:
        rejected = e.kind == PolyErrorKind.ENCODING
    assert_true(rejected)
    # The empty scope retains the original archive layout.
    var ordinary = cad.decompose([p.copy()])
    var explicit_empty = cad.decompose([p.copy()], fixed=())
    assert_equal(io.dumps(ordinary), io.dumps(explicit_empty))
    roundtrip(ordinary)
    print("fixed CAD archive checks passed")
