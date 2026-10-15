"""Versioned algebraic array records with incremental transport and checked witnesses."""
from apn_mojo import Rational, ConversionLimits
from polynomials import Budget, PolyError, PolyErrorKind, CoefficientRing, NoContext, QQ, DensePoly, SparsePoly, Term, Monomial, MonomialOrder, RealAlgebraic, SturmSequence, check_real_algebraic, PolynomialTextSource, PolynomialTextSink
from ._common import _allocation, _charge_size, _sum, _product, _grow
from ._algebraic_value import Real, Complex, _plain, _copy_witness


def _bad(message: String) -> PolyError:
    return PolyError(PolyErrorKind.ENCODING, message)


def _work(n: Int, mut budget: Budget) raises PolyError -> Int:
    var extent = _sum(n, 1)
    var work = _product(extent, extent)
    _charge_size(work, budget)
    return _product(work, 64)


def _qtext(value: QQ, mut budget: Budget) raises PolyError -> String:
    _ = _work(value.size(), budget)
    try:
        return value.value.to_string()
    except e:
        raise _bad(String(e))


def _unsigned(text: String, maximum: UInt64, mut budget: Budget) raises PolyError -> UInt64:
    _charge_size(text.byte_length(), budget)
    if not text.byte_length() or (text.byte_length() > 1 and text.as_bytes()[0] == 48):
        raise _bad("noncanonical unsigned metadata")
    var result = UInt64(0)
    for byte in text.as_bytes():
        if byte < 48 or byte > 57:
            raise _bad("invalid unsigned metadata")
        var digit = UInt64(byte - 48)
        if digit > maximum or result > (maximum - digit) // 10:
            raise PolyError(PolyErrorKind.CAPACITY, "record integer exceeds capacity")
        result = result * 10 + digit
    return result


@fieldwise_init
struct _Reader[S: PolynomialTextSource](Movable):
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
                    raise _bad("unterminated record field")
                self.chunk = next.take()
                self.position = 0
            var start = self.position
            var bytes = self.chunk.as_bytes()
            while self.position < len(bytes) and bytes[self.position] != 10:
                budget.charge(1)
                self.position += 1
            var length = self.position - start
            if length > self.maximum - result.byte_length():
                raise PolyError(PolyErrorKind.CAPACITY, "record field exceeds byte limit")
            _charge_size(length, budget)
            _ = _sum(_sum(result.byte_length(), length), 1)
            result.write(self.chunk[byte=start:self.position])
            if self.position < self.chunk.byte_length():
                self.position += 1
                return result^

    def number(mut self, maximum: UInt64, mut budget: Budget) raises PolyError -> UInt64:
        return _unsigned(self.line(budget), maximum, budget)

    def count(mut self, mut budget: Budget) raises PolyError -> Int:
        return Int(self.number(UInt64(Int.MAX), budget))

    def identifier(mut self, mut budget: Budget) raises PolyError -> Int64:
        var text = self.line(budget)
        var value: Int64
        if text.startswith("-"):
            var magnitude = _unsigned(String(text[byte=1:]), UInt64(1) << 63, budget)
            if not magnitude:
                raise _bad("negative zero variable identifier")
            value = Int64.MIN if magnitude == UInt64(1) << 63 else -Int64(magnitude)
        else:
            value = Int64(_unsigned(text, UInt64(Int64.MAX), budget))
        return value

    def rational(mut self, mut budget: Budget) raises PolyError -> QQ:
        var text = self.line(budget)
        var capacity = _work(text.byte_length(), budget)
        for byte in text.as_bytes():
            if (byte < 48 or byte > 57) and byte != 45 and byte != 47:
                raise _bad("expected rational witness field")
        var value: Rational
        try:
            value = Rational.parse(text, limits=ConversionLimits(max_input_bytes=text.byte_length(), max_digits=text.byte_length(), max_allocated_bytes=capacity))
        except e:
            raise _bad(String(e))
        var q = QQ(value^)
        if text != _qtext(q, budget):
            raise _bad("noncanonical rational witness field")
        return q^

    def polynomial(mut self, mut budget: Budget) raises PolyError -> DensePoly[QQ]:
        var count = self.count(budget)
        var coefficients = List[QQ]()
        for _ in range(count):
            budget.charge(1)
            var c = self.rational(budget)
            _grow(coefficients, budget)
            coefficients.append(c^)
        if count and coefficients[count - 1].is_zero():
            raise _bad("trailing zero in witness polynomial")
        return DensePoly[QQ].from_coefficients(0, NoContext(), coefficients^, budget)

    def real(mut self, mut budget: Budget) raises PolyError -> Real:
        var rational = Bool(self.number(1, budget))
        var low = self.rational(budget)
        var high = self.rational(budget)
        var polynomial = self.polynomial(budget)
        var polynomials = List[DensePoly[QQ]]()
        var count = self.count(budget)
        for _ in range(count):
            budget.charge(1)
            var p = self.polynomial(budget)
            _grow(polynomials, budget)
            polynomials.append(p^)
        var quotients = List[DensePoly[QQ]]()
        count = self.count(budget)
        for _ in range(count):
            budget.charge(1)
            var p = self.polynomial(budget)
            _grow(quotients, budget)
            quotients.append(p^)
        var witness = RealAlgebraic(polynomial^, low^, high^, SturmSequence(polynomials^, quotients^), rational)
        if not check_real_algebraic(witness, budget):
            raise _bad("invalid algebraic root witness")
        # Preserve the supplied witness; compute the normal equality cache only
        # after the independent Sturm/evaluation check has accepted its value.
        return _plain(witness^, budget)

    def scalar[C: CoefficientRing](mut self, mut budget: Budget) raises PolyError -> C:
        var real = self.real(budget)
        comptime if C == Real:
            return rebind_var[C](real^)
        else:
            return rebind_var[C](Complex(real^, self.real(budget)))

    def finish(mut self, mut budget: Budget) raises PolyError:
        while self.position == self.chunk.byte_length():
            budget.charge(1)
            var next = self.source.read(budget)
            if not next:
                return
            self.chunk = next.take()
            self.position = 0
        raise _bad("trailing record data")


struct _Writer(Movable):
    var buffer: String

    def __init__(out self):
        self.buffer = String()

    def flush[S: PolynomialTextSink](mut self, mut sink: S, mut budget: Budget) raises PolyError:
        if self.buffer.byte_length():
            var chunk = String()
            swap(chunk, self.buffer)
            sink.write(chunk^, budget)

    def line[S: PolynomialTextSink](mut self, var text: String, mut sink: S, mut budget: Budget) raises PolyError:
        var added = _sum(text.byte_length(), 1)
        _charge_size(added, budget)
        _ = _sum(_sum(self.buffer.byte_length(), added), 1)
        self.buffer.write(text, "\n")
        if self.buffer.byte_length() >= 4096:
            self.flush(sink, budget)

    def polynomial[S: PolynomialTextSink](mut self, p: DensePoly[QQ], mut sink: S, mut budget: Budget) raises PolyError:
        p.check_invariants(budget)
        if p.variable != 0:
            raise PolyError(PolyErrorKind.DOMAIN, "algebraic witness variable must be zero")
        self.line(String(len(p.coefficients)), sink, budget)
        for coefficient_index in range(len(p.coefficients)):
            ref c = p.coefficients[coefficient_index]
            self.line(_qtext(c, budget), sink, budget)

    def real[S: PolynomialTextSink](mut self, value: Real, mut sink: S, mut budget: Budget) raises PolyError:
        ref witness = value.witness()
        if not check_real_algebraic(witness, budget):
            raise PolyError(PolyErrorKind.DOMAIN, "invalid algebraic root witness")
        var normalized = _plain(_copy_witness(witness, budget), budget)
        if value != normalized or value.is_zero() != normalized.is_zero():
            raise PolyError(PolyErrorKind.DOMAIN, "invalid algebraic equality cache")
        self.line(String(Int(witness.rational)), sink, budget)
        self.line(_qtext(witness.low, budget), sink, budget)
        self.line(_qtext(witness.high, budget), sink, budget)
        self.polynomial(witness.polynomial, sink, budget)
        self.line(String(len(witness.sequence.polynomials)), sink, budget)
        for p in witness.sequence.polynomials:
            self.polynomial(p, sink, budget)
        self.line(String(len(witness.sequence.quotients)), sink, budget)
        for p in witness.sequence.quotients:
            self.polynomial(p, sink, budget)

    def scalar[C: CoefficientRing, S: PolynomialTextSink](mut self, value: C, mut sink: S, mut budget: Budget) raises PolyError:
        comptime if C == Real:
            self.real(rebind[Real](value), sink, budget)
        else:
            self.real(rebind[Complex](value).real(), sink, budget)
            self.real(rebind[Complex](value).imag(), sink, budget)


def _header[C: CoefficientRing, S: PolynomialTextSink](dense: Bool, mut writer: _Writer, mut sink: S, mut budget: Budget) raises PolyError:
    writer.line("POLYNOMIALS-ALGEBRAIC 1", sink, budget)
    writer.line("DENSE" if dense else "SPARSE", sink, budget)
    writer.line("AA" if C == Real else "QQbar", sink, budget)


def _reader[C: CoefficientRing, S: PolynomialTextSource](var source: S, dense: Bool, maximum: Int, mut budget: Budget) raises PolyError -> _Reader[S]:
    if maximum < 1:
        raise PolyError(PolyErrorKind.DOMAIN, "field limit must be positive")
    var reader = _Reader[S](source^, String(), 0, maximum)
    if reader.line(budget) != "POLYNOMIALS-ALGEBRAIC 1" or reader.line(budget) != ("DENSE" if dense else "SPARSE") or reader.line(budget) != ("AA" if C == Real else "QQbar"):
        raise _bad("incompatible algebraic array record")
    return reader^


def _dump_dense[C: CoefficientRing, S: PolynomialTextSink](p: DensePoly[C], mut sink: S, mut budget: Budget) raises PolyError:
    p.check_invariants(budget)
    var writer = _Writer()
    _header[C](True, writer, sink, budget)
    writer.line(String(p.variable), sink, budget)
    writer.line(String(len(p.coefficients)), sink, budget)
    for coefficient_index in range(len(p.coefficients)):
        ref c = p.coefficients[coefficient_index]
        writer.scalar(c, sink, budget)
    writer.flush(sink, budget)


def _dump_sparse[C: CoefficientRing, S: PolynomialTextSink](p: SparsePoly[C], mut sink: S, mut budget: Budget) raises PolyError:
    p.check_invariants(budget)
    _allocation[Int64](len(p.order.variables), budget)
    _ = MonomialOrder(p.order.kind, p.order.variables.copy(), budget)
    var writer = _Writer()
    _header[C](False, writer, sink, budget)
    writer.line(String(p.order.kind), sink, budget)
    writer.line(String(len(p.order.variables)), sink, budget)
    for axis in p.order.variables:
        writer.line(String(axis), sink, budget)
    writer.line(String(len(p.terms)), sink, budget)
    for term in p.terms:
        writer.line(String(len(term.monomial.pairs)), sink, budget)
        for pair in term.monomial.pairs:
            writer.line(String(pair[0]), sink, budget)
            writer.line(String(pair[1]), sink, budget)
        writer.scalar(term.coefficient, sink, budget)
    writer.flush(sink, budget)


def _load_dense[C: CoefficientRing, S: PolynomialTextSource](var source: S, mut budget: Budget, maximum: Int) raises PolyError -> DensePoly[C]:
    try:
        var reader = _reader[C](source^, True, maximum, budget)
        var variable = reader.identifier(budget)
        var count = reader.count(budget)
        var coefficients = List[C]()
        for _ in range(count):
            budget.charge(1)
            var value = reader.scalar[C](budget)
            _grow(coefficients, budget)
            coefficients.append(value^)
        reader.finish(budget)
        if count and coefficients[count - 1].is_zero():
            raise _bad("trailing zero algebraic dense coefficient")
        return DensePoly[C].from_coefficients(variable, rebind_var[C.Context](NoContext()), coefficients^, budget)
    except e:
        if e.kind == PolyErrorKind.DOMAIN and maximum > 0:
            raise _bad(e.message)
        raise e^


def _load_sparse[C: CoefficientRing, S: PolynomialTextSource](var source: S, mut budget: Budget, maximum: Int) raises PolyError -> SparsePoly[C]:
    try:
        var reader = _reader[C](source^, False, maximum, budget)
        var kind = UInt8(reader.number(3, budget))
        var variables = List[Int64]()
        var count = reader.count(budget)
        for _ in range(count):
            budget.charge(1)
            var axis = reader.identifier(budget)
            _grow(variables, budget)
            variables.append(axis)
        var order = MonomialOrder(kind, variables^, budget)
        count = reader.count(budget)
        var result = SparsePoly[C](order^, rebind_var[C.Context](NoContext()))
        for i in range(count):
            budget.charge(1)
            var pairs = List[Tuple[Int64, UInt32]]()
            var length = reader.count(budget)
            for j in range(length):
                budget.charge(1)
                var axis = reader.identifier(budget)
                var power = UInt32(reader.number(UInt64(UInt32.MAX), budget))
                if not power or (j and axis <= pairs[j - 1][0]):
                    raise _bad("noncanonical algebraic monomial")
                _grow(pairs, budget)
                pairs.append((axis, power))
            var monomial = Monomial(pairs^, budget)
            result.order.validate(monomial, budget)
            if i and result.order.compare(result.terms[i - 1].monomial, monomial, budget) >= 0:
                raise _bad("nonincreasing algebraic sparse terms")
            var coefficient = reader.scalar[C](budget)
            if coefficient.is_zero():
                raise _bad("zero algebraic sparse coefficient")
            _grow(result.terms, budget)
            result.terms.append(Term[C](monomial^, coefficient^))
        reader.finish(budget)
        return result^
    except e:
        if e.kind == PolyErrorKind.DOMAIN and maximum > 0:
            raise _bad(e.message)
        raise e^


@fieldwise_init
struct _Source(PolynomialTextSource):
    var text: Optional[String]

    def read(mut self, mut budget: Budget) raises PolyError -> Optional[String]:
        budget.charge(1)
        var result: Optional[String] = None
        swap(result, self.text)
        return result^


@fieldwise_init
struct _Sink(PolynomialTextSink):
    var text: String

    def write(mut self, var chunk: String, mut budget: Budget) raises PolyError:
        _charge_size(chunk.byte_length(), budget)
        _ = _sum(_sum(self.text.byte_length(), chunk.byte_length()), 1)
        self.text.write(chunk)


def _dumps_dense[C: CoefficientRing](p: DensePoly[C], mut budget: Budget) raises PolyError -> String:
    var sink = _Sink(String())
    _dump_dense(p, sink, budget)
    var result = String()
    swap(result, sink.text)
    return result^


def _dumps_sparse[C: CoefficientRing](p: SparsePoly[C], mut budget: Budget) raises PolyError -> String:
    var sink = _Sink(String())
    _dump_sparse(p, sink, budget)
    var result = String()
    swap(result, sink.text)
    return result^


def _loads_dense[C: CoefficientRing](text: String, mut budget: Budget) raises PolyError -> DensePoly[C]:
    _charge_size(text.byte_length(), budget)
    return _load_dense[C](_Source(Optional(text.copy())), budget, max(1, text.byte_length()))


def _loads_sparse[C: CoefficientRing](text: String, mut budget: Budget) raises PolyError -> SparsePoly[C]:
    _charge_size(text.byte_length(), budget)
    return _load_sparse[C](_Source(Optional(text.copy())), budget, max(1, text.byte_length()))
