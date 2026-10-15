"""Archive transports reject corrupt evidence and preserve budget failures."""
from poly import dense, sparse, roots, io, exact, batch
from poly._archive_meta import _Archive
from poly._algebraic_io import _Writer, _Sink
from poly._result_io import _tag
from polynomials import Budget, PolyError, PolyErrorKind, PolynomialTextSource, ProjectionResultantCertificate
from std.testing import assert_true, assert_equal

@fieldwise_init
struct Source(PolynomialTextSource):
    var text: String
    var position: Int
    def read(mut self, mut budget: Budget) raises PolyError -> Optional[String]:
        budget.charge(1)
        if self.position == self.text.byte_length(): return None
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

def main() raises:
    var value = dense.count_roots([-2, 0, 1], -2, 2)
    var text = io.dumps(value)
    var restored = io.load[type_of(value)](Source(text.copy(), 0))
    assert_true(roots.check(restored))
    assert_equal(restored.count(), 2)
    var work = Budget()
    _ = io.load[type_of(value)](Source(text.copy(), 0), work)
    var spent = UInt64.MAX - work.remaining()
    var exact_budget = Budget(spent)
    _ = io.load[type_of(value)](Source(text.copy(), 0), exact_budget)
    assert_equal(exact_budget.remaining(), UInt64(0))
    var failed = False
    try:
        _ = io.load[type_of(value)](Source(text.copy(), 0), max_steps=spent-1)
    except e:
        failed = e.kind == PolyErrorKind.LIMIT
    assert_true(failed)
    failed = False
    try:
        _ = io.load[type_of(value)](Source(text.copy(), 0), max_field_bytes=8)
    except e:
        failed = e.kind == PolyErrorKind.CAPACITY
    assert_true(failed)
    for malformed in [String(text, "0\n"), String(text[byte=0:text.byte_length()-1])]:
        failed = False
        try:
            _ = io.loads[type_of(value)](malformed)
        except e:
            failed = e.kind == PolyErrorKind.ENCODING
        assert_true(failed)
    var values = batch.asarray([value.copy(), value.copy()], shape=List[Int]([1, 2]))
    var saved = io.loads[type_of(values)](io.dumps(values))
    assert_equal(saved.shape(), List[Int]([1, 2]))
    assert_true(roots.check(saved[1]))
    # Bypass only the writer gate to test the independent reader gate.
    value._data[].proof.count += 1
    assert_true(not roots.check(value))
    failed = False
    try:
        _ = io.loads[type_of(value)](unchecked(value))
    except e:
        failed = e.kind == PolyErrorKind.ENCODING
    assert_true(failed)
    failed = False
    try:
        _ = io.dumps(value)
    except e:
        failed = e.kind == PolyErrorKind.ENCODING
    assert_true(failed)
    var division = dense.divide_checked([1, 2, 1], [1, 1])
    division._data[].proof.remainder.coefficients.append(division._data[].proof.quotient.coefficients[0].copy())
    assert_true(not exact.check(division))
    failed = False
    try:
        _ = io.loads[type_of(division)](unchecked(division))
    except e:
        failed = e.kind == PolyErrorKind.ENCODING
    assert_true(failed)
    var shaped = sparse.resultant_certified([((1, 0), 1)], [((0, 1), 1)], axis=0)
    var shaped_text = io.dumps(shaped)
    failed = False
    try:
        _ = io.loads[exact.Result[1, Tuple[Int, Int, Int], True, ProjectionResultantCertificate]](shaped_text)
    except e:
        failed = e.kind == PolyErrorKind.ENCODING
    assert_true(failed)
    print("completion archive stream checks passed")
