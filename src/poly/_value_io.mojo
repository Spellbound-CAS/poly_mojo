"""Typed, incremental value records built from checked scalar fields."""
from std.builtin.rebind import downcast
from std.memory import ArcPointer
from apn_mojo import Integer, Rational, ExactComplex, Float, Complex as NumericComplex, Ball, ComplexBall, ConversionLimits
from polynomials import Budget, PolyError, PolyErrorKind, PolynomialTextSource, PolynomialTextSink, ModInt, Modulus, QQ, DensePoly, SparsePoly, Term, Monomial, NumericDensePoly, NumericSparsePoly, NumericTerm, TruncatedSeries, BasisPolynomial, NumericSeries, NumericBasisPolynomial
from . import dense, sparse, series, basis, rational, batch
from . import _numeric as _num
from ._algebraic_value import Real, Complex
from ._algebraic_io import _Reader, _Writer, _bad, _work, _qtext
from polynomials import CoefficientBuffer
from ._common import _Kind, _Scalar, _Ring, _Context, _wrap, _unwrap, _size, _allocation, _charge_size, _grow, _sum
from ._archive_meta import _Archive
from . import _result_io
from ._io_meta import _Record, _Dense, _Sparse, _Series, _Basis, _Fraction, _Batch

comptime _Numeric[T: AnyType] = T == Float or T == NumericComplex or T == Ball or T == ComplexBall
comptime _Supported[T: AnyType] = T == Int or T == Integer or T == Rational or T == ExactComplex or T == ModInt or T == Real or T == Complex or _Numeric[T]


def _schema[T: Copyable & Deinitable](mut budget: Budget) raises PolyError -> String:
    budget.charge(1)
    comptime if conforms_to(T, _Archive):
        return _result_io._tag[T](budget)
    elif conforms_to(T, _Record):
        comptime V = downcast[T, _Record]
        comptime C = V._Item
        comptime k = _Kind[C]
        comptime D = _num._StoredDomain[C]
        comptime if V._Format == _Batch:
            comptime assert T == batch.Array[C], "unsupported batch value type"
        elif V._Format == _Fraction:
            comptime if k >= 0:
                comptime assert T == rational.Function[V._Power, k], "unsupported rational value type"
            else:
                comptime assert T == rational.NumericFunction[V._Power, D], "unsupported numerical fraction type"
        else:
            comptime assert k >= 0 or _Numeric[C], "unsupported coefficient value type"
            comptime if V._Format == _Dense:
                comptime if k >= 0:
                    comptime assert T == dense.Array[k], "unsupported dense value type"
                else:
                    comptime assert T == dense.NumericArray[D], "unsupported dense value type"
            elif V._Format == _Sparse:
                comptime if k >= 0:
                    comptime assert T == sparse.Array[V._Power, k], "unsupported sparse value type"
                else:
                    comptime assert T == sparse.NumericArray[V._Power, D], "unsupported sparse value type"
            elif V._Format == _Series:
                comptime if k >= 0:
                    comptime assert T == series.Array[k], "unsupported series value type"
                else:
                    comptime assert T == series.NumericArray[D], "unsupported series value type"
            elif V._Format == _Basis:
                comptime if k >= 0:
                    comptime assert T == basis.Array[k], "unsupported basis value type"
                else:
                    comptime assert T == basis.NumericArray[D], "unsupported basis value type"
            else:
                comptime assert False, "unsupported value representation"
        var item = _schema[V._Item](budget)
        _charge_size(_sum(item.byte_length(), 64), budget)
        comptime if V._Format == _Batch:
            return String("BATCH(", item, ")")
        elif V._Format == _Dense:
            return String("DENSE(", item, ")")
        elif V._Format == _Series:
            return String("SERIES(", item, ")")
        elif V._Format == _Basis:
            return String("BASIS(", item, ")")
        else:
            return String("SPARSE(" if V._Format == _Sparse else "RATIONAL(", sparse._width[V._Power](), ",", item, ")")
    else:
        comptime assert _Supported[T], "poly.io supports prepared polynomial values, batches, and supported coefficient scalars"
        comptime if T == Int: return "Int"
        elif T == Integer: return "ZZ"
        elif T == Rational: return "QQ"
        elif T == ExactComplex: return "QQi"
        elif T == ModInt: return "GF"
        elif T == Real: return "AA"
        elif T == Complex: return "QQbar"
        elif T == Float: return "Float"
        elif T == NumericComplex: return "Complex"
        elif T == Ball: return "Ball"
        else: return "ComplexBall"


def _numeric_text[T: Copyable & Deinitable](value: T, mut budget: Budget) raises PolyError -> String:
    var box = _num._bounds(value, budget)
    var size = _sum(_sum(box.real_low.size(), box.real_high.size()), _sum(box.imag_low.size(), box.imag_high.size()))
    _ = _work(size, budget)
    try:
        comptime if T == Float:
            return rebind[Float](value).to_json()
        elif T == NumericComplex:
            return rebind[NumericComplex](value).to_json()
        elif T == Ball:
            return rebind[Ball](value).to_json()
        else:
            return rebind[ComplexBall](value).to_json()
    except e:
        raise _bad(String(e))


def _write_scalar[T: Copyable & Deinitable, S: PolynomialTextSink](value: T, mut writer: _Writer, mut sink: S, mut budget: Budget) raises PolyError:
    _charge_size(_size(value), budget)
    comptime if T == Real or T == Complex:
        writer.scalar(rebind[Real if T == Real else Complex](value), sink, budget)
    elif T == ModInt:
        ref v = rebind[ModInt](value)
        if v.value >= v.modulus:
            raise PolyError(PolyErrorKind.DOMAIN, "noncanonical modular coefficient")
        writer.line(String(v.value), sink, budget)
    elif _Numeric[T]:
        writer.line(_numeric_text(value, budget), sink, budget)
    elif T == Int:
        writer.line(String(rebind[Int](value)), sink, budget)
    elif T == Integer:
        writer.line(_qtext(QQ(Rational(rebind[Integer](value))), budget), sink, budget)
    elif T == Rational:
        writer.line(_qtext(QQ(rebind[Rational](value)), budget), sink, budget)
    else:
        comptime assert T == ExactComplex, "poly.io supports prepared polynomial values, batches, and supported coefficient scalars"
        writer.line(_qtext(QQ(rebind[ExactComplex](value).real()), budget), sink, budget)
        writer.line(_qtext(QQ(rebind[ExactComplex](value).imag()), budget), sink, budget)


def _read_scalar[T: Copyable & Deinitable, S: PolynomialTextSource](mut reader: _Reader[S], modulus: UInt64, mut budget: Budget) raises PolyError -> T:
    comptime if T == Real or T == Complex:
        return rebind_var[T](reader.scalar[Real if T == Real else Complex](budget))
    elif T == ModInt:
        var ctx = Modulus(modulus, budget)
        var value = reader.number(UInt64.MAX, budget)
        if value >= modulus:
            raise _bad("noncanonical modular coefficient")
        return rebind_var[T](ModInt(value, ctx))
    elif T == Int:
        return rebind_var[T](Int(reader.identifier(budget)))
    elif T == Integer:
        var q = reader.rational(budget)
        if q.value.denominator() != 1:
            raise _bad("expected an integer coefficient")
        return rebind_var[T](q.value.numerator().copy())
    elif T == Rational:
        var q = reader.rational(budget)
        return rebind_var[T](q.value.copy())
    elif T == ExactComplex:
        var real = reader.rational(budget)
        var imag = reader.rational(budget)
        return rebind_var[T](ExactComplex(real.value, imag.value))
    else:
        comptime assert _Numeric[T], "poly.io supports prepared polynomial values, batches, and supported coefficient scalars"
        var text = reader.line(budget)
        var capacity = _work(text.byte_length(), budget)
        var value: T
        try:
            var limits = ConversionLimits(max_input_bytes=text.byte_length(), max_digits=text.byte_length(), max_allocated_bytes=capacity)
            comptime if T == Float:
                value = rebind_var[T](Float.from_json(text, limits=limits))
            elif T == NumericComplex:
                value = rebind_var[T](NumericComplex.from_json(text, limits=limits))
            elif T == Ball:
                value = rebind_var[T](Ball.from_json(text, limits=limits))
            else:
                value = rebind_var[T](ComplexBall.from_json(text, limits=limits))
        except e:
            raise _bad(String(e))
        if _numeric_text(value, budget) != text:
            raise _bad("noncanonical numerical coefficient")
        return value^


def _zero[T: Copyable & Deinitable](value: T, mut budget: Budget) raises PolyError -> Bool:
    comptime if _Kind[T] >= 0:
        _charge_size(_size(value), budget)
        return _wrap[_Kind[T]](rebind[_Scalar[_Kind[T]]](value)).is_zero()
    else:
        var box = _num._bounds(value, budget)
        return box.real_low.is_zero() and box.real_high.is_zero() and box.imag_low.is_zero() and box.imag_high.is_zero()


def _read_values[T: Copyable & Deinitable, S: PolynomialTextSource](mut reader: _Reader[S], modulus: UInt64, mut budget: Budget) raises PolyError -> List[T]:
    var count = reader.count(budget)
    var result = List[T]()
    for _ in range(count):
        budget.charge(1)
        var value = _read_scalar[T](reader, modulus, budget)
        _grow(result, budget)
        result.append(value^)
    return result^


def _write_values[T: Copyable & Deinitable, S: PolynomialTextSink](values: List[T], mut writer: _Writer, mut sink: S, mut budget: Budget) raises PolyError:
    writer.line(String(len(values)), sink, budget)
    for value in values:
        _write_scalar(value, writer, sink, budget)


def _write_core_values[k: Int, S: PolynomialTextSink](values: List[_Ring[k]], mut writer: _Writer, mut sink: S, mut budget: Budget) raises PolyError:
    writer.line(String(len(values)), sink, budget)
    for value in values:
        _charge_size(value.size(), budget)
        _write_scalar(_unwrap[k](value), writer, sink, budget)


def _write_core_values[k: Int, S: PolynomialTextSink](values: CoefficientBuffer[_Ring[k]], mut writer: _Writer, mut sink: S, mut budget: Budget) raises PolyError:
    writer.line(String(len(values)), sink, budget)
    for index in range(len(values)):
        ref value = values[index]
        _charge_size(value.size(), budget)
        _write_scalar(_unwrap[k](value), writer, sink, budget)


def _core_values[k: Int](var values: List[_Scalar[k]], mut budget: Budget) raises PolyError -> List[_Ring[k]]:
    var result = List[_Ring[k]]()
    for value in values:
        _charge_size(_size(value), budget)
        _grow(result, budget)
        result.append(_wrap[k](value))
    return result^


def _write_sparse[P: Copyable & Deinitable, k: Int, S: PolynomialTextSink](p: SparsePoly[_Ring[k]], mut writer: _Writer, mut sink: S, mut budget: Budget) raises PolyError:
    p.check_invariants(budget)
    var expected = sparse._order[P](budget)
    if p.order != expected:
        raise PolyError(PolyErrorKind.DOMAIN, "array requires canonical axis order")
    writer.line(String(len(p.terms)), sink, budget)
    for i in range(len(p.terms)):
        budget.charge(1)
        ref term = p.terms[i]
        term.monomial.check_invariants(budget)
        p.order.validate(term.monomial, budget)
        if i and p.order.compare(p.terms[i - 1].monomial, term.monomial, budget) >= 0:
            raise PolyError(PolyErrorKind.DOMAIN, "nonincreasing sparse terms")
        if _zero(_unwrap[k](term.coefficient), budget):
            raise PolyError(PolyErrorKind.DOMAIN, "zero sparse coefficient")
        writer.line(String(len(term.monomial.pairs)), sink, budget)
        for pair in term.monomial.pairs:
            writer.line(String(pair[0]), sink, budget)
            writer.line(String(pair[1]), sink, budget)
        _write_scalar(_unwrap[k](term.coefficient), writer, sink, budget)


def _write_dense[k: Int, S: PolynomialTextSink](p: DensePoly[_Ring[k]], mut writer: _Writer, mut sink: S, mut budget: Budget) raises PolyError:
    p.check_invariants(budget)
    if p.variable != 0 or (len(p.coefficients) and _zero(_unwrap[k](p.coefficients[len(p.coefficients) - 1]), budget)):
        raise PolyError(PolyErrorKind.DOMAIN, "noncanonical dense array")
    _write_core_values[k](p.coefficients, writer, sink, budget)


def _write_sparse[P: Copyable & Deinitable, D: _num.NumericDomain, S: PolynomialTextSink](p: NumericSparsePoly[D], mut writer: _Writer, mut sink: S, mut budget: Budget) raises PolyError:
    var expected = sparse._order[P](budget)
    if p.order != expected:
        raise PolyError(PolyErrorKind.DOMAIN, "array requires canonical axis order")
    writer.line(String(len(p.terms)), sink, budget)
    for i in range(len(p.terms)):
        budget.charge(1)
        ref term = p.terms[i]
        term.monomial.check_invariants(budget)
        p.order.validate(term.monomial, budget)
        if i and p.order.compare(p.terms[i - 1].monomial, term.monomial, budget) >= 0:
            raise PolyError(PolyErrorKind.DOMAIN, "nonincreasing sparse terms")
        if _zero(term.coefficient, budget):
            raise PolyError(PolyErrorKind.DOMAIN, "zero sparse coefficient")
        writer.line(String(len(term.monomial.pairs)), sink, budget)
        for pair in term.monomial.pairs:
            writer.line(String(pair[0]), sink, budget)
            writer.line(String(pair[1]), sink, budget)
        _write_scalar(term.coefficient, writer, sink, budget)


def _write_dense[D: _num.NumericDomain, S: PolynomialTextSink](p: NumericDensePoly[D], mut writer: _Writer, mut sink: S, mut budget: Budget) raises PolyError:
    if p.variable != 0 or (len(p.coefficients) and _zero(p.coefficients[len(p.coefficients) - 1], budget)):
        raise PolyError(PolyErrorKind.DOMAIN, "noncanonical dense array")
    _write_values(p.coefficients, writer, sink, budget)


def _write[T: Copyable & Deinitable, S: PolynomialTextSink](value: T, mut writer: _Writer, mut sink: S, mut budget: Budget) raises PolyError:
    writer.line(_schema[T](budget), sink, budget)
    comptime if conforms_to(T, _Archive):
        _result_io._write(rebind[downcast[T, _Archive]](value), writer, sink, budget)
    elif conforms_to(T, _Record):
        comptime V = downcast[T, _Record]
        comptime C = V._Item
        comptime F = V._Format
        comptime k = _Kind[C]
        comptime D = _num._StoredDomain[C]
        comptime if F == _Batch:
            ref v = rebind[batch.Array[C]](value)
            writer.line(String(v.ndim()), sink, budget)
            for axis in v.shape():
                writer.line(String(axis), sink, budget)
            writer.line(String(len(v)), sink, budget)
            for item in v:
                _write(item, writer, sink, budget)
        else:
            comptime if k == 3:
                from ._common import _WithModulus
                var modulus = rebind[downcast[T, _WithModulus]](value).modulus()
                _ = Modulus(modulus, budget)
                writer.line(String(modulus), sink, budget)
            comptime if F == _Fraction:
                comptime if k >= 0:
                    ref v = rebind[rational.Function[V._Power, k]](value)
                    if not rational.check(v, budget):
                        raise PolyError(PolyErrorKind.DOMAIN, "invalid rational-function evidence")
                    _write(rational.numerator(v, budget), writer, sink, budget)
                    _write(rational.denominator(v, budget), writer, sink, budget)
                else:
                    ref v = rebind[rational.NumericFunction[V._Power, D]](value)
                    _ = rational._nf._make(v._numerator.copy(), v._denominator.copy(), budget)
                    _write(v._numerator, writer, sink, budget)
                    _write(v._denominator, writer, sink, budget)
            elif F == _Sparse:
                comptime if k >= 0:
                    _write_sparse[V._Power, k](rebind[sparse.Array[V._Power, k]](value)._core[], writer, sink, budget)
                else:
                    _write_sparse[V._Power, D](rebind[sparse.NumericArray[V._Power, D]](value)._core[], writer, sink, budget)
            elif F == _Basis:
                comptime if k >= 0:
                    ref v = rebind[basis.Array[k]](value)
                    ref p = v._core[]
                    if p.kind > 2 or p.variable != 0:
                        raise PolyError(PolyErrorKind.DOMAIN, "invalid basis metadata")
                    basis._shape(p.kind, len(p.coefficients), len(p.nodes))
                    writer.line(v.basis(), sink, budget)
                    for c in p.coefficients:
                        budget.charge(1)
                        if c.context() != p.ring:
                            raise PolyError(PolyErrorKind.DOMAIN, "mixed basis coefficient contexts")
                    for c in p.nodes:
                        budget.charge(1)
                        if c.context() != p.ring:
                            raise PolyError(PolyErrorKind.DOMAIN, "mixed basis node contexts")
                    _write_core_values[k](rebind[List[_Ring[k]]](p.coefficients), writer, sink, budget)
                    _write_core_values[k](rebind[List[_Ring[k]]](p.nodes), writer, sink, budget)
                else:
                    ref v = rebind[basis.NumericArray[D]](value)
                    ref p = v._core[]
                    if p.kind > 2 or p.variable != 0:
                        raise PolyError(PolyErrorKind.DOMAIN, "invalid basis metadata")
                    basis._shape(p.kind, len(p.coefficients), len(p.nodes))
                    writer.line(v.basis(), sink, budget)
                    _write_values(p.coefficients, writer, sink, budget)
                    _write_values(p.nodes, writer, sink, budget)
            elif F == _Series:
                comptime if k >= 0:
                    ref v = rebind[series.Array[k]](value)
                    ref p = v._core[].polynomial
                    if v.order() < len(p.coefficients) or v.order() < 0:
                        raise PolyError(PolyErrorKind.DOMAIN, "invalid series order")
                    writer.line(String(v.order()), sink, budget)
                    _write_dense[k](rebind[DensePoly[_Ring[k]]](p), writer, sink, budget)
                else:
                    ref v = rebind[series.NumericArray[D]](value)
                    ref p = v._core[].polynomial
                    if v.order() < len(p.coefficients) or v.order() < 0:
                        raise PolyError(PolyErrorKind.DOMAIN, "invalid series order")
                    writer.line(String(v.order()), sink, budget)
                    if len(v._zero[]) != 1 or not _zero(v._zero[][0], budget):
                        raise PolyError(PolyErrorKind.DOMAIN, "invalid numerical series padding")
                    _write_scalar(v._zero[][0], writer, sink, budget)
                    _write_dense(p, writer, sink, budget)
            elif F == _Dense:
                comptime if k >= 0:
                    ref v = rebind[dense.Array[k]](value)
                    _write_dense[k](v._core[], writer, sink, budget)
                else:
                    ref v = rebind[dense.NumericArray[D]](value)
                    _write_dense(v._core[], writer, sink, budget)
    else:
        comptime if T == ModInt:
            _ = Modulus(rebind[ModInt](value).modulus, budget)
            writer.line(String(rebind[ModInt](value).modulus), sink, budget)
        _write_scalar(value, writer, sink, budget)


def _read[T: Copyable & Deinitable, S: PolynomialTextSource](mut reader: _Reader[S], mut budget: Budget) raises PolyError -> T:
    if reader.line(budget) != _schema[T](budget):
        raise _bad("value record has an incompatible type or shape")
    comptime if conforms_to(T, _Archive):
        return rebind_var[T](_result_io._read[downcast[T, _Archive]](reader, budget))
    elif conforms_to(T, _Record):
        comptime V = downcast[T, _Record]
        comptime C = V._Item
        comptime F = V._Format
        comptime k = _Kind[C]
        comptime D = _num._StoredDomain[C]
        comptime if F == _Batch:
            var shape = List[Int]()
            var rank = reader.count(budget)
            for _ in range(rank):
                var extent = reader.count(budget)
                _grow(shape, budget)
                shape.append(extent)
            var layout = batch._layout(shape, budget)
            if reader.count(budget) != layout.count:
                raise _bad("batch shape does not match its item count")
            var items = List[C]()
            for _ in range(layout.count):
                var item = _read[C](reader, budget)
                _grow(items, budget)
                items.append(item^)
            return rebind_var[T](batch._take(items^, layout^, budget))
        else:
            var modulus = UInt64(0)
            comptime if k == 3:
                modulus = reader.number(UInt64.MAX, budget)
                _ = Modulus(modulus, budget)
            comptime if F == _Fraction:
                comptime if k >= 0:
                    comptime A = sparse.Array[V._Power, k]
                    var p = _read[A](reader, budget)
                    var q = _read[A](reader, budget)
                    comptime if k == 3:
                        if p.modulus() != modulus or q.modulus() != modulus:
                            raise _bad("fraction modulus disagrees with its coefficient arrays")
                    var result = rational._normalize[V._Power, k](p._core[].reorder(p._core[].order, budget), q._core[].reorder(q._core[].order, budget), budget)
                    var a = rational.numerator(result, budget)
                    var b = rational.denominator(result, budget)
                    if not a._core[].sub(p._core[], budget).is_zero() or not b._core[].sub(q._core[], budget).is_zero():
                        raise _bad("rational function is not reduced with a monic denominator")
                    # Normalization proves the saved pair is canonical. Retain the
                    # supplied witnesses rather than its equivalent arithmetic outputs.
                    result._data[].value.numerator = rebind_var[SparsePoly[rational._Field[k]]](p._core[].reorder(p._core[].order, budget))
                    result._data[].value.denominator = rebind_var[SparsePoly[rational._Field[k]]](q._core[].reorder(q._core[].order, budget))
                    if not rational.check(result, budget):
                        raise _bad("invalid restored rational-function evidence")
                    return rebind_var[T](result^)
                else:
                    comptime A = sparse.NumericArray[V._Power, D]
                    var p = _read[A](reader, budget)
                    var q = _read[A](reader, budget)
                    return rebind_var[T](rational._nf._make(p^, q^, budget))
            elif F == _Sparse:
                comptime if k >= 0:
                    var order = sparse._order[V._Power](budget)
                    var count = reader.count(budget)
                    var ctx = _Context[k]()
                    ctx.observe(modulus, budget)
                    var p = SparsePoly[_Ring[k]](order.copy(), ctx.get())
                    for i in range(count):
                        var length = reader.count(budget)
                        if length > sparse._width[V._Power]():
                            raise _bad("too many sparse axes")
                        var pairs = List[Tuple[Int64, UInt32]]()
                        for j in range(length):
                            var axis = reader.identifier(budget)
                            var power = UInt32(reader.number(UInt64(UInt32.MAX), budget))
                            if axis < 0 or axis >= Int64(sparse._width[V._Power]()) or not power or (j and axis <= pairs[j - 1][0]):
                                raise _bad("noncanonical sparse power")
                            _grow(pairs, budget)
                            pairs.append((axis, power))
                        var monomial = Monomial(pairs^, budget)
                        if i and order.compare(p.terms[i - 1].monomial, monomial, budget) >= 0:
                            raise _bad("nonincreasing sparse terms")
                        var value = _read_scalar[C](reader, modulus, budget)
                        if _zero(value, budget):
                            raise _bad("zero sparse coefficient")
                        _grow(p.terms, budget)
                        p.terms.append(Term[_Ring[k]](monomial^, _wrap[k](rebind[_Scalar[k]](value))))
                    return rebind_var[T](sparse.Array[V._Power, k](p^, budget))
                else:
                    var order = sparse._order[V._Power](budget)
                    var count = reader.count(budget)
                    var p = NumericSparsePoly[D](order.copy())
                    for i in range(count):
                        var length = reader.count(budget)
                        if length > sparse._width[V._Power]():
                            raise _bad("too many sparse axes")
                        var pairs = List[Tuple[Int64, UInt32]]()
                        for j in range(length):
                            var axis = reader.identifier(budget)
                            var power = UInt32(reader.number(UInt64(UInt32.MAX), budget))
                            if axis < 0 or axis >= Int64(sparse._width[V._Power]()) or not power or (j and axis <= pairs[j - 1][0]):
                                raise _bad("noncanonical sparse power")
                            _grow(pairs, budget)
                            pairs.append((axis, power))
                        var monomial = Monomial(pairs^, budget)
                        if i and order.compare(p.terms[i - 1].monomial, monomial, budget) >= 0:
                            raise _bad("nonincreasing sparse terms")
                        var value = _read_scalar[C](reader, modulus, budget)
                        if _zero(value, budget):
                            raise _bad("zero sparse coefficient")
                        _grow(p.terms, budget)
                        p.terms.append(NumericTerm[D](monomial^, rebind_var[D.Scalar](value^)))
                    return rebind_var[T](sparse.NumericArray[V._Power, D](p^, budget))
            elif F == _Basis:
                var kind = basis._kind(reader.line(budget))
                var values = _read_values[C](reader, modulus, budget)
                var nodes = _read_values[C](reader, modulus, budget)
                basis._shape(kind, len(values), len(nodes))
                comptime if k >= 0:
                    var ctx = _Context[k]()
                    ctx.observe(modulus, budget)
                    var core = BasisPolynomial[basis._Field[k]](0, kind, rebind_var[basis._Field[k].Context](ctx.get()), rebind_var[List[basis._Field[k]]](_core_values[k](rebind_var[List[_Scalar[k]]](values^), budget)), rebind_var[List[basis._Field[k]]](_core_values[k](rebind_var[List[_Scalar[k]]](nodes^), budget)))
                    return rebind_var[T](basis._wrap_basis[k](core^, budget))
                else:
                    _allocation[NumericBasisPolynomial[D]](1, budget)
                    var core = NumericBasisPolynomial[D](0, kind, rebind_var[List[D.Scalar]](values^), rebind_var[List[D.Scalar]](nodes^))
                    return rebind_var[T](basis.NumericArray[D](ArcPointer(core^)))
            else:
                var order: Int
                var zeros = List[C]()
                comptime if F == _Series:
                    order = reader.count(budget)
                    comptime if k < 0:
                        var zero = _read_scalar[C](reader, 0, budget)
                        if not _zero(zero, budget):
                            raise _bad("numerical series padding must be exact zero")
                        _grow(zeros, budget)
                        zeros.append(zero^)
                else:
                    order = 0
                var values = _read_values[C](reader, modulus, budget)
                if len(values) and _zero(values[len(values) - 1], budget):
                    raise _bad("trailing zero polynomial coefficient")
                comptime if F == _Series:
                    if len(values) > order:
                        raise _bad("series coefficients exceed its known order")
                comptime if k >= 0:
                    var ctx = _Context[k]()
                    ctx.observe(modulus, budget)
                    var p = DensePoly[_Ring[k]].from_coefficients(0, ctx.get(), _core_values[k](rebind_var[List[_Scalar[k]]](values^), budget), budget)
                    comptime if F == _Series:
                        return rebind_var[T](series._wrap_series[k](TruncatedSeries[series._Field[k]](rebind_var[DensePoly[series._Field[k]]](p^), order), budget))
                    else:
                        return rebind_var[T](dense.Array[k](p^, budget))
                else:
                    var p = NumericDensePoly[D].from_coefficients(0, rebind_var[List[D.Scalar]](values^), budget)
                    comptime if F == _Series:
                        _allocation[NumericSeries[D]](1, budget)
                        _allocation[List[D.Scalar]](1, budget)
                        _allocation[D.Scalar](1, budget)
                        return rebind_var[T](series.NumericArray[D](ArcPointer(NumericSeries[D](p^, order)), ArcPointer(rebind_var[List[D.Scalar]](zeros^))))
                    else:
                        return rebind_var[T](dense.NumericArray[D](p^, budget))
    else:
        var modulus = UInt64(0)
        comptime if T == ModInt:
            modulus = reader.number(UInt64.MAX, budget)
            _ = Modulus(modulus, budget)
        return _read_scalar[T](reader, modulus, budget)
