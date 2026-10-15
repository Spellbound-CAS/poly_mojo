"""Version-one line records shared by exact and numerical interchange."""
from .errors import Budget, PolyError, PolyErrorKind, _allocation, _charge_size, _sum
from .monomial import Monomial, MonomialOrder
from ._text_ops import _unsigned


struct _TextWriter(Movable):
    var buffer: String

    def __init__(out self):
        self.buffer = String()

    def add(mut self, var line: String, mut budget: Budget) raises PolyError:
        var added = _sum(line.byte_length(), 1)
        _ = _sum(self.buffer.byte_length(), _sum(added, 1))
        _charge_size(added, budget)
        self.buffer.write(line, "\n")

    def monomial(mut self, value: Monomial, mut budget: Budget) raises PolyError:
        self.add(String(len(value.pairs)), budget)
        for pair in value.pairs:
            self.add(String(pair[0]), budget)
            self.add(String(pair[1]), budget)

    def finish(mut self, mut budget: Budget) raises PolyError -> String:
        budget.charge(1)
        var result = String()
        swap(result, self.buffer)
        return result^


@fieldwise_init
struct _TextReader(Movable):
    var text: String
    var position: Int

    def line(mut self, mut budget: Budget) raises PolyError -> String:
        var start = self.position
        var bytes = self.text.as_bytes()
        while self.position < len(bytes) and bytes[self.position] != 10:
            budget.charge(1)
            self.position += 1
        if self.position == len(bytes):
            raise PolyError(PolyErrorKind.ENCODING, "unterminated interchange field")
        var end = self.position
        self.position += 1
        _charge_size(end - start, budget)
        return String(self.text[byte=start:end])

    def expect(mut self, expected: String, mut budget: Budget) raises PolyError:
        var actual = self.line(budget)
        _charge_size(_sum(actual.byte_length(), expected.byte_length()), budget)
        if actual != expected:
            raise PolyError(PolyErrorKind.ENCODING, "unexpected interchange header")

    def number(mut self, maximum: UInt64, mut budget: Budget) raises PolyError -> UInt64:
        return _number(self.line(budget), maximum, budget)

    def identifier(mut self, mut budget: Budget) raises PolyError -> Int64:
        return _identifier(self.line(budget), budget)

    def count(mut self, mut budget: Budget) raises PolyError -> Int:
        # Every remaining item needs at least one byte; reject oversized counts
        # before allocating from untrusted metadata.
        return Int(self.number(UInt64(self.text.byte_length() - self.position), budget))


@fieldwise_init
struct _TextTerm(Copyable, Movable):
    var monomial: Monomial
    var coefficient: String


@fieldwise_init
struct _RecordReader(Movable):
    var reader: _TextReader
    var dense: Bool
    var domain: String
    var variable: Int64
    var order: MonomialOrder
    var count: Int

    @always_inline
    def monomial(mut self, mut budget: Budget) raises PolyError -> Monomial:
        var pairs = self.reader.count(budget)
        _allocation[Tuple[Int64, UInt32]](pairs, budget)
        var result = Monomial()
        for j in range(pairs):
            var identifier = self.reader.identifier(budget)
            var exponent = UInt32(self.reader.number(UInt64(UInt32.MAX), budget))
            if exponent == 0 or (j > 0 and identifier <= result.pairs[j - 1][0]):
                raise PolyError(PolyErrorKind.ENCODING, "noncanonical serialized monomial")
            result.pairs.append((identifier, exponent))
        try:
            self.order.validate(result, budget)
        except e:
            if e.kind == PolyErrorKind.DOMAIN:
                raise PolyError(PolyErrorKind.ENCODING, "serialized variable absent from monomial order")
            raise e^
        return result^

    def increasing(self, previous: Monomial, current: Monomial, mut budget: Budget) raises PolyError:
        if self.order.compare(previous, current, budget) >= 0:
            raise PolyError(PolyErrorKind.ENCODING, "nonincreasing serialized terms")

    def finish(self) raises PolyError:
        if self.reader.position != self.reader.text.byte_length():
            raise PolyError(PolyErrorKind.ENCODING, "trailing polynomial interchange data")


def _read_header(text: String, mut budget: Budget) raises PolyError -> _RecordReader:
    _charge_size(text.byte_length(), budget)
    var reader = _TextReader(text.copy(), 0)
    reader.expect("POLYNOMIALS 1", budget)
    var representation = reader.line(budget)
    if representation != "DENSE" and representation != "SPARSE":
        raise PolyError(PolyErrorKind.ENCODING, "unknown polynomial representation")
    var domain = reader.line(budget)
    var variable = Int64(0)
    var order = MonomialOrder()
    if representation == "DENSE":
        variable = reader.identifier(budget)
    else:
        var kind = UInt8(reader.number(3, budget))
        var count = reader.count(budget)
        _allocation[Int64](count, budget)
        var variables = List[Int64](capacity=count)
        for _ in range(count):
            variables.append(reader.identifier(budget))
        try:
            order = MonomialOrder(kind, variables^, budget)
        except e:
            if e.kind == PolyErrorKind.DOMAIN:
                raise PolyError(PolyErrorKind.ENCODING, "invalid serialized monomial order")
            raise e^
    var count = reader.count(budget)
    return _RecordReader(reader^, representation == "DENSE", domain^, variable, order^, count)


def _write_header(dense: Bool, var domain: String, variable: Int64, order: MonomialOrder, count: Int, mut budget: Budget) raises PolyError -> _TextWriter:
    if not dense:
        _charge_size(len(order.variables), budget)
        _ = MonomialOrder(order.kind, order.variables.copy(), budget)
    var writer = _TextWriter()
    writer.add("POLYNOMIALS 1", budget)
    writer.add("DENSE" if dense else "SPARSE", budget)
    writer.add(domain^, budget)
    if dense:
        writer.add(String(variable), budget)
    else:
        writer.add(String(order.kind), budget)
        writer.add(String(len(order.variables)), budget)
        for identifier in order.variables:
            writer.add(String(identifier), budget)
    writer.add(String(count), budget)
    return writer^


def _number(text: String, maximum: UInt64, mut budget: Budget) raises PolyError -> UInt64:
    var value = _unsigned(text, maximum, budget)
    if text != String(value):
        raise PolyError(PolyErrorKind.ENCODING, "noncanonical integer metadata")
    return value


def _identifier(text: String, mut budget: Budget) raises PolyError -> Int64:
    var value: Int64
    if text.startswith("-"):
        var magnitude = _unsigned(String(text[byte=1:]), UInt64(1) << 63, budget)
        value = Int64.MIN if magnitude == UInt64(1) << 63 else -Int64(magnitude)
    else:
        value = Int64(_unsigned(text, UInt64(Int64.MAX), budget))
    if text != String(value):
        raise PolyError(PolyErrorKind.ENCODING, "noncanonical variable identifier")
    return value
