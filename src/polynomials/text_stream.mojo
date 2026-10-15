"""Bounded text transport for canonical version-one polynomial interchange."""
from .errors import Budget, PolyError, PolyErrorKind, _allocation, _charge_size
from .monomial import Monomial, MonomialOrder
from ._interchange import _TextWriter, _number, _identifier


trait PolynomialTextSource(Movable, Deinitable):
    """Supply arbitrary UTF-8 chunks of one canonical polynomial record."""

    def read(mut self, mut budget: Budget) raises PolyError -> Optional[String]:
        """Read the next owned text chunk.

        Args:
            budget: The shared transport allowance.

        Returns:
            A chunk, possibly empty, or None at the end of input.

        Raises:
            PolyError: ENCODING for transport failure or LIMIT on exhaustion.
        """
        ...


trait PolynomialTextSink:
    """Consume canonical text chunks without retaining the complete record."""

    def write(mut self, var chunk: String, mut budget: Budget) raises PolyError:
        """Consume the next owned text chunk in order.

        Args:
            chunk: The next portion of the version-one record.
            budget: The shared transport allowance.

        Raises:
            PolyError: ENCODING for transport failure or LIMIT on exhaustion.
        """
        ...


def _flush[S: PolynomialTextSink](mut writer: _TextWriter, mut sink: S, target: Int, mut budget: Budget, final: Bool = False) raises PolyError:
    if writer.buffer.byte_length() and (final or writer.buffer.byte_length() >= target):
        sink.write(writer.finish(budget), budget)


@fieldwise_init
struct _ChunkReader[S: PolynomialTextSource](Movable):
    var source: Self.S
    var chunk: String
    var position: Int
    var maximum: Int

    def line(mut self, mut budget: Budget) raises PolyError -> String:
        var result = String()
        while True:
            budget.charge(1)
            if self.position == self.chunk.byte_length():
                var next = self.source.read(budget)
                if not next:
                    raise PolyError(PolyErrorKind.ENCODING, "unterminated interchange field")
                self.chunk = next.take()
                self.position = 0
            var start = self.position
            var bytes = self.chunk.as_bytes()
            while self.position < len(bytes) and bytes[self.position] != 10:
                budget.charge(1)
                self.position += 1
            var length = self.position - start
            if length > self.maximum - result.byte_length():
                raise PolyError(PolyErrorKind.CAPACITY, "interchange field exceeds byte limit")
            _charge_size(length, budget)
            result.write(self.chunk[byte=start:self.position])
            if self.position < self.chunk.byte_length():
                self.position += 1
                return result^

    def number(mut self, maximum: UInt64, mut budget: Budget) raises PolyError -> UInt64:
        return _number(self.line(budget), maximum, budget)

    def identifier(mut self, mut budget: Budget) raises PolyError -> Int64:
        return _identifier(self.line(budget), budget)

    def count(mut self, mut budget: Budget) raises PolyError -> Int:
        return Int(self.number(UInt64(Int.MAX), budget))

    def finish(mut self, mut budget: Budget) raises PolyError:
        while self.position == self.chunk.byte_length():
            budget.charge(1)
            var next = self.source.read(budget)
            if not next:
                return
            self.chunk = next.take()
            self.position = 0
        raise PolyError(PolyErrorKind.ENCODING, "trailing polynomial interchange data")


@fieldwise_init
struct _ChunkRecord[S: PolynomialTextSource](Movable):
    var reader: _ChunkReader[Self.S]
    var order: MonomialOrder
    var variable: Int64
    var count: Int

    def monomial(mut self, mut budget: Budget) raises PolyError -> Monomial:
        var count = self.reader.count(budget)
        var pairs = List[Tuple[Int64, UInt32]]()
        for i in range(count):
            var identifier = self.reader.identifier(budget)
            var exponent = UInt32(self.reader.number(UInt64(UInt32.MAX), budget))
            if exponent == 0 or (i and identifier <= pairs[i - 1][0]):
                raise PolyError(PolyErrorKind.ENCODING, "noncanonical serialized monomial")
            _allocation[Tuple[Int64, UInt32]](1, budget)
            pairs.append((identifier, exponent))
        var result = Monomial()
        result.pairs = pairs^
        try:
            self.order.validate(result, budget)
        except e:
            if e.kind == PolyErrorKind.DOMAIN:
                raise PolyError(PolyErrorKind.ENCODING, "serialized variable absent from monomial order")
            raise e^
        return result^


def _chunk_header[S: PolynomialTextSource](var source: S, dense: Bool, domain: String, maximum: Int, mut budget: Budget) raises PolyError -> _ChunkRecord[S]:
    if maximum < 1:
        raise PolyError(PolyErrorKind.DOMAIN, "interchange field limit must be positive")
    var reader = _ChunkReader[S](source^, String(), 0, maximum)
    if reader.line(budget) != "POLYNOMIALS 1" or reader.line(budget) != ("DENSE" if dense else "SPARSE") or reader.line(budget) != domain:
        raise PolyError(PolyErrorKind.ENCODING, "incompatible polynomial interchange header")
    var variable = Int64(0)
    var order = MonomialOrder()
    if dense:
        variable = reader.identifier(budget)
    else:
        var kind = UInt8(reader.number(3, budget))
        var count = reader.count(budget)
        var variables = List[Int64]()
        for _ in range(count):
            var identifier = reader.identifier(budget)
            _allocation[Int64](1, budget)
            variables.append(identifier)
        try:
            order = MonomialOrder(kind, variables^, budget)
        except e:
            if e.kind == PolyErrorKind.DOMAIN:
                raise PolyError(PolyErrorKind.ENCODING, "invalid serialized monomial order")
            raise e^
    var count = reader.count(budget)
    return _ChunkRecord[S](reader^, order^, variable, count)
