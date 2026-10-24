"""Private version-one compact CAD payload, with bounded exact decoding.

Integer fields use an eight-byte value offset by one (the only negative value
in this schema is -1). Counts bound allocations by remaining input bytes.
Rationals use the canonical kernel rational codec; polynomials retain their
monomial order. Search counters are diagnostic and are not serialized.
The public archive boundary must independently verify the decoded claim.
"""
from .coefficients import QQ, NoContext
from .dense import DensePoly
from .sparse import SparsePoly, Term
from .monomial import Monomial, MonomialOrder
from .algebraic_types import RealAlgebraic
from .errors import Budget, PolyError, PolyErrorKind, _allocation, _charge_size
from .codec import _write_word, _read_word, _rational_bytes, _read_rational
from ._cad_final import _FinalRoot
from ._cad_sample_field import _SampleField
from ._cad_field_stack import _FieldStack
from ._cad_interval_roots import _IntervalRootStats
from ._cad_projection_final import _FinalProjection
from ._cad_tree_final import _TreeStack, _FinalCAD


def _bad_final(message: String) -> PolyError:
    return PolyError(PolyErrorKind.ENCODING, message)


def _read_count(data: Span[UInt8, _], mut position: Int, mut budget: Budget) raises PolyError -> Int:
    var count = _read_word(data, position, 8, budget)
    if count > UInt64(len(data) - position):
        raise _bad_final("compact CAD count exceeds remaining payload")
    return Int(count)


def _put_final[T: Copyable & Deinitable](value: T, mut data: List[UInt8], mut budget: Budget) raises PolyError:
    budget.charge(1)
    comptime if T == Bool:
        _write_word(data, UInt64(rebind[Bool](value)), 1, budget)
    elif T == Int or T == Int64 or T == Int8 or T == UInt8 or T == UInt32:
        var number: Int
        comptime if T == Int:
            number = rebind[Int](value)
        elif T == Int64:
            number = Int(rebind[Int64](value))
        elif T == Int8:
            number = Int(rebind[Int8](value))
        elif T == UInt8:
            number = Int(rebind[UInt8](value))
        else:
            number = Int(rebind[UInt32](value))
        if number < -1:
            raise _bad_final("compact CAD integer is below -1")
        _write_word(data, UInt64(number) + 1 if number >= 0 else 0, 8, budget)
    elif T == QQ:
        var bytes = _rational_bytes(rebind[QQ](value).value, budget)
        _write_word(data, UInt64(len(bytes)), 8, budget)
        _allocation[UInt8](len(bytes), budget)
        for byte in bytes:
            data.append(byte)
    elif T == DensePoly[QQ]:
        ref p = rebind[DensePoly[QQ]](value)
        p.check_invariants(budget)
        _put_final(p.variable, data, budget)
        _write_word(data, UInt64(len(p.coefficients)), 8, budget)
        for c in p.coefficients:
            _put_final(c, data, budget)
    elif T == SparsePoly[QQ]:
        ref p = rebind[SparsePoly[QQ]](value)
        p.check_invariants(budget)
        _put_final(p.order.kind, data, budget)
        _put_final(p.order.variables, data, budget)
        _put_final(p.terms, data, budget)
    elif T == Tuple[Int64, UInt32]:
        ref pair = rebind[Tuple[Int64, UInt32]](value)
        _put_final(pair[0], data, budget)
        _put_final(pair[1], data, budget)
    elif T == Optional[RealAlgebraic]:
        ref optional = rebind[Optional[RealAlgebraic]](value)
        _put_final(Bool(optional), data, budget)
        if optional:
            _put_final(optional.value(), data, budget)
    elif T == List[QQ]:
        ref values = rebind[List[QQ]](value)
        _write_word(data, UInt64(len(values)), 8, budget)
        for item in values:
            _put_final(item, data, budget)
    elif T == List[Int64]:
        ref values = rebind[List[Int64]](value)
        _write_word(data, UInt64(len(values)), 8, budget)
        for item in values:
            _put_final(item, data, budget)
    elif T == List[Term[QQ]]:
        ref values = rebind[List[Term[QQ]]](value)
        _write_word(data, UInt64(len(values)), 8, budget)
        for item in values:
            _put_final(item, data, budget)
    elif T == List[Tuple[Int64, UInt32]]:
        ref values = rebind[List[Tuple[Int64, UInt32]]](value)
        _write_word(data, UInt64(len(values)), 8, budget)
        for item in values:
            _put_final(item, data, budget)
    elif T == List[DensePoly[QQ]]:
        ref values = rebind[List[DensePoly[QQ]]](value)
        _write_word(data, UInt64(len(values)), 8, budget)
        for item in values:
            _put_final(item, data, budget)
    elif T == List[List[DensePoly[QQ]]]:
        ref values = rebind[List[List[DensePoly[QQ]]]](value)
        _write_word(data, UInt64(len(values)), 8, budget)
        for item in values:
            _put_final(item, data, budget)
    elif T == List[Int8]:
        ref values = rebind[List[Int8]](value)
        _write_word(data, UInt64(len(values)), 8, budget)
        for item in values:
            _put_final(item, data, budget)
    elif T == List[List[Int]]:
        ref values = rebind[List[List[Int]]](value)
        _write_word(data, UInt64(len(values)), 8, budget)
        for item in values:
            _put_final(item, data, budget)
    elif T == List[Int]:
        ref values = rebind[List[Int]](value)
        _write_word(data, UInt64(len(values)), 8, budget)
        for item in values:
            _put_final(item, data, budget)
    elif T == List[_FinalRoot]:
        ref values = rebind[List[_FinalRoot]](value)
        _write_word(data, UInt64(len(values)), 8, budget)
        for item in values:
            _put_final(item, data, budget)
    elif T == List[List[Int8]]:
        ref values = rebind[List[List[Int8]]](value)
        _write_word(data, UInt64(len(values)), 8, budget)
        for item in values:
            _put_final(item, data, budget)
    elif T == List[List[SparsePoly[QQ]]]:
        ref values = rebind[List[List[SparsePoly[QQ]]]](value)
        _write_word(data, UInt64(len(values)), 8, budget)
        for item in values:
            _put_final(item, data, budget)
    elif T == List[SparsePoly[QQ]]:
        ref values = rebind[List[SparsePoly[QQ]]](value)
        _write_word(data, UInt64(len(values)), 8, budget)
        for item in values:
            _put_final(item, data, budget)
    elif T == List[Bool]:
        ref values = rebind[List[Bool]](value)
        _write_word(data, UInt64(len(values)), 8, budget)
        for item in values:
            _put_final(item, data, budget)
    elif T == List[_TreeStack]:
        ref values = rebind[List[_TreeStack]](value)
        _write_word(data, UInt64(len(values)), 8, budget)
        for item in values:
            _put_final(item, data, budget)
    elif T == Monomial:
        ref record = rebind[Monomial](value)
        _put_final(record.pairs, data, budget)
    elif T == Term[QQ]:
        ref record = rebind[Term[QQ]](value)
        _put_final(record.monomial, data, budget)
        _put_final(record.coefficient, data, budget)
    elif T == RealAlgebraic:
        ref record = rebind[RealAlgebraic](value)
        _put_final(record.polynomial, data, budget)
        _put_final(record.low, data, budget)
        _put_final(record.high, data, budget)
        _put_final(record.rational, data, budget)
    elif T == _FinalRoot:
        ref record = rebind[_FinalRoot](value)
        _put_final(record.basis, data, budget)
        _put_final(record.low, data, budget)
        _put_final(record.high, data, budget)
    elif T == _SampleField:
        ref record = rebind[_SampleField](value)
        _put_final(record.modulus, data, budget)
        _put_final(record.embedding, data, budget)
        _put_final(record.variables, data, budget)
        _put_final(record.coordinates, data, budget)
    elif T == _FieldStack:
        ref record = rebind[_FieldStack](value)
        _put_final(record.basis, data, budget)
        _put_final(record.leading, data, budget)
        _put_final(record.multiplicities, data, budget)
        _put_final(record.roots, data, budget)
        _put_final(record.sectors, data, budget)
        _put_final(record.signs, data, budget)
    elif T == _FinalProjection:
        ref record = rebind[_FinalProjection](value)
        _put_final(record.variables, data, budget)
        _put_final(record.brown, data, budget)
        _put_final(record.levels, data, budget)
    elif T == _TreeStack:
        ref record = rebind[_TreeStack](value)
        _put_final(record.parent, data, budget)
        _put_final(record.slot, data, budget)
        _put_final(record.level, data, budget)
        _put_final(record.field, data, budget)
        _put_final(record.parent_generator, data, budget)
        _put_final(record.orders, data, budget)
        _put_final(record.stack, data, budget)
        _put_final(record.children, data, budget)
        _put_final(record.truth, data, budget)
    elif T == _FinalCAD:
        ref record = rebind[_FinalCAD](value)
        _put_final(record.inputs, data, budget)
        _put_final(record.masks, data, budget)
        _put_final(record.projection, data, budget)
        _put_final(record.stacks, data, budget)
        _put_final(record.satisfiable, data, budget)
    else:
        comptime assert False, "unsupported compact CAD field"


def _get_final[T: Copyable](data: Span[UInt8, _], mut position: Int, mut budget: Budget) raises PolyError -> T:
    budget.charge(1)
    comptime if T == Bool:
        var value = _read_word(data, position, 1, budget)
        if value > 1:
            raise _bad_final("invalid compact CAD Boolean")
        return rebind_var[T](Bool(value))
    elif T == Int or T == Int64 or T == Int8 or T == UInt8 or T == UInt32:
        var encoded = _read_word(data, position, 8, budget)
        if encoded > UInt64(Int.MAX) + 1:
            raise _bad_final("compact CAD integer exceeds native range")
        var number = Int(encoded - 1) if encoded else -1
        comptime if T == Int8:
            if number > 127:
                raise _bad_final("compact CAD signed byte exceeds range")
        elif T == UInt8:
            if number < 0 or number > 255:
                raise _bad_final("compact CAD byte exceeds range")
        elif T == UInt32:
            if number < 0 or UInt64(number) > UInt64(UInt32.MAX):
                raise _bad_final("compact CAD exponent exceeds range")
        comptime if T == Int:
            return rebind_var[T](number)
        elif T == Int64:
            return rebind_var[T](Int64(number))
        elif T == Int8:
            return rebind_var[T](Int8(number))
        elif T == UInt8:
            return rebind_var[T](UInt8(number))
        else:
            return rebind_var[T](UInt32(number))
    elif T == QQ:
        var count = _read_count(data, position, budget)
        var value = _read_rational(data[position:position + count], budget)
        position += count
        return rebind_var[T](value^)
    elif T == DensePoly[QQ]:
        var variable = _get_final[Int64](data, position, budget)
        var values = _get_final[List[QQ]](data, position, budget)
        if len(values) and values[len(values) - 1].is_zero():
            raise _bad_final("noncanonical compact CAD dense polynomial")
        return rebind_var[T](DensePoly[QQ].from_coefficients(variable, NoContext(), values^, budget))
    elif T == SparsePoly[QQ]:
        var kind = _get_final[UInt8](data, position, budget)
        var variables = _get_final[List[Int64]](data, position, budget)
        var order = MonomialOrder(kind, variables^, budget)
        var values = _get_final[List[Term[QQ]]](data, position, budget)
        var p = SparsePoly[QQ](order^, NoContext())
        p.terms = values^
        p.check_invariants(budget)
        return rebind_var[T](p^)
    elif T == Tuple[Int64, UInt32]:
        var axis = _get_final[Int64](data, position, budget)
        var power = _get_final[UInt32](data, position, budget)
        return rebind_var[T]((axis, power))
    elif T == Optional[RealAlgebraic]:
        var present = _get_final[Bool](data, position, budget)
        var value = Optional[RealAlgebraic]()
        if present:
            value = _get_final[RealAlgebraic](data, position, budget)
        return rebind_var[T](value^)
    elif T == List[QQ]:
        var count = _read_count(data, position, budget)
        _allocation[QQ](count, budget)
        var values = List[QQ](capacity=count)
        for _ in range(count):
            values.append(_get_final[QQ](data, position, budget))
        return rebind_var[T](values^)
    elif T == List[Int64]:
        var count = _read_count(data, position, budget)
        _allocation[Int64](count, budget)
        var values = List[Int64](capacity=count)
        for _ in range(count):
            values.append(_get_final[Int64](data, position, budget))
        return rebind_var[T](values^)
    elif T == List[Term[QQ]]:
        var count = _read_count(data, position, budget)
        _allocation[Term[QQ]](count, budget)
        var values = List[Term[QQ]](capacity=count)
        for _ in range(count):
            values.append(_get_final[Term[QQ]](data, position, budget))
        return rebind_var[T](values^)
    elif T == List[Tuple[Int64, UInt32]]:
        var count = _read_count(data, position, budget)
        _allocation[Tuple[Int64, UInt32]](count, budget)
        var values = List[Tuple[Int64, UInt32]](capacity=count)
        for _ in range(count):
            values.append(_get_final[Tuple[Int64, UInt32]](data, position, budget))
        return rebind_var[T](values^)
    elif T == List[DensePoly[QQ]]:
        var count = _read_count(data, position, budget)
        _allocation[DensePoly[QQ]](count, budget)
        var values = List[DensePoly[QQ]](capacity=count)
        for _ in range(count):
            values.append(_get_final[DensePoly[QQ]](data, position, budget))
        return rebind_var[T](values^)
    elif T == List[List[DensePoly[QQ]]]:
        var count = _read_count(data, position, budget)
        _allocation[List[DensePoly[QQ]]](count, budget)
        var values = List[List[DensePoly[QQ]]](capacity=count)
        for _ in range(count):
            values.append(_get_final[List[DensePoly[QQ]]](data, position, budget))
        return rebind_var[T](values^)
    elif T == List[Int8]:
        var count = _read_count(data, position, budget)
        _allocation[Int8](count, budget)
        var values = List[Int8](capacity=count)
        for _ in range(count):
            values.append(_get_final[Int8](data, position, budget))
        return rebind_var[T](values^)
    elif T == List[List[Int]]:
        var count = _read_count(data, position, budget)
        _allocation[List[Int]](count, budget)
        var values = List[List[Int]](capacity=count)
        for _ in range(count):
            values.append(_get_final[List[Int]](data, position, budget))
        return rebind_var[T](values^)
    elif T == List[Int]:
        var count = _read_count(data, position, budget)
        _allocation[Int](count, budget)
        var values = List[Int](capacity=count)
        for _ in range(count):
            values.append(_get_final[Int](data, position, budget))
        return rebind_var[T](values^)
    elif T == List[_FinalRoot]:
        var count = _read_count(data, position, budget)
        _allocation[_FinalRoot](count, budget)
        var values = List[_FinalRoot](capacity=count)
        for _ in range(count):
            values.append(_get_final[_FinalRoot](data, position, budget))
        return rebind_var[T](values^)
    elif T == List[List[Int8]]:
        var count = _read_count(data, position, budget)
        _allocation[List[Int8]](count, budget)
        var values = List[List[Int8]](capacity=count)
        for _ in range(count):
            values.append(_get_final[List[Int8]](data, position, budget))
        return rebind_var[T](values^)
    elif T == List[List[SparsePoly[QQ]]]:
        var count = _read_count(data, position, budget)
        _allocation[List[SparsePoly[QQ]]](count, budget)
        var values = List[List[SparsePoly[QQ]]](capacity=count)
        for _ in range(count):
            values.append(_get_final[List[SparsePoly[QQ]]](data, position, budget))
        return rebind_var[T](values^)
    elif T == List[SparsePoly[QQ]]:
        var count = _read_count(data, position, budget)
        _allocation[SparsePoly[QQ]](count, budget)
        var values = List[SparsePoly[QQ]](capacity=count)
        for _ in range(count):
            values.append(_get_final[SparsePoly[QQ]](data, position, budget))
        return rebind_var[T](values^)
    elif T == List[Bool]:
        var count = _read_count(data, position, budget)
        _allocation[Bool](count, budget)
        var values = List[Bool](capacity=count)
        for _ in range(count):
            values.append(_get_final[Bool](data, position, budget))
        return rebind_var[T](values^)
    elif T == List[_TreeStack]:
        var count = _read_count(data, position, budget)
        _allocation[_TreeStack](count, budget)
        var values = List[_TreeStack](capacity=count)
        for _ in range(count):
            values.append(_get_final[_TreeStack](data, position, budget))
        return rebind_var[T](values^)
    elif T == Monomial:
        var field0 = _get_final[List[Tuple[Int64, UInt32]]](data, position, budget)
        return rebind_var[T](Monomial(field0^, budget))
    elif T == Term[QQ]:
        var field0 = _get_final[Monomial](data, position, budget)
        var field1 = _get_final[QQ](data, position, budget)
        return rebind_var[T](Term[QQ](field0^, field1^))
    elif T == RealAlgebraic:
        var field0 = _get_final[DensePoly[QQ]](data, position, budget)
        var field1 = _get_final[QQ](data, position, budget)
        var field2 = _get_final[QQ](data, position, budget)
        var field3 = _get_final[Bool](data, position, budget)
        return rebind_var[T](RealAlgebraic(field0^, field1^, field2^, field3))
    elif T == _FinalRoot:
        var field0 = _get_final[Int](data, position, budget)
        var field1 = _get_final[QQ](data, position, budget)
        var field2 = _get_final[QQ](data, position, budget)
        return rebind_var[T](_FinalRoot(field0, field1^, field2^))
    elif T == _SampleField:
        var field0 = _get_final[DensePoly[QQ]](data, position, budget)
        var field1 = _get_final[Optional[RealAlgebraic]](data, position, budget)
        var field2 = _get_final[List[Int64]](data, position, budget)
        var field3 = _get_final[List[DensePoly[QQ]]](data, position, budget)
        return rebind_var[T](_SampleField(field0^, field1^, field2^, field3^))
    elif T == _FieldStack:
        var field0 = _get_final[List[List[DensePoly[QQ]]]](data, position, budget)
        var field1 = _get_final[List[Int8]](data, position, budget)
        var field2 = _get_final[List[List[Int]]](data, position, budget)
        var field3 = _get_final[List[_FinalRoot]](data, position, budget)
        var field4 = _get_final[List[QQ]](data, position, budget)
        var field5 = _get_final[List[List[Int8]]](data, position, budget)
        return rebind_var[T](_FieldStack(field0^, field1^, field2^, field3^, field4^, field5^, _IntervalRootStats()))
    elif T == _FinalProjection:
        var field0 = _get_final[List[Int64]](data, position, budget)
        var field1 = _get_final[Bool](data, position, budget)
        var field2 = _get_final[List[List[SparsePoly[QQ]]]](data, position, budget)
        return rebind_var[T](_FinalProjection(field0^, field1, field2^))
    elif T == _TreeStack:
        var field0 = _get_final[Int](data, position, budget)
        var field1 = _get_final[Int](data, position, budget)
        var field2 = _get_final[Int](data, position, budget)
        var field3 = _get_final[_SampleField](data, position, budget)
        var field4 = _get_final[DensePoly[QQ]](data, position, budget)
        var field5 = _get_final[List[List[Int]]](data, position, budget)
        var field6 = _get_final[_FieldStack](data, position, budget)
        var field7 = _get_final[List[Int]](data, position, budget)
        var field8 = _get_final[List[Bool]](data, position, budget)
        return rebind_var[T](_TreeStack(field0, field1, field2, field3^, field4^, field5^, field6^, field7^, field8^))
    elif T == _FinalCAD:
        var field0 = _get_final[List[SparsePoly[QQ]]](data, position, budget)
        var field1 = _get_final[List[Int]](data, position, budget)
        var field2 = _get_final[_FinalProjection](data, position, budget)
        var field3 = _get_final[List[_TreeStack]](data, position, budget)
        var field4 = _get_final[Bool](data, position, budget)
        return rebind_var[T](_FinalCAD(field0^, field1^, field2^, field3^, field4, 0, 0, 0))
    else:
        comptime assert False, "unsupported compact CAD field"


def _encode_final(result: _FinalCAD, mut budget: Budget) raises PolyError -> List[UInt8]:
    var data = List[UInt8]()
    # ASCII QCF1, independent of the enclosing facade archive version.
    _write_word(data, 0x31464351, 4, budget)
    _put_final(result, data, budget)
    return data^


def _decode_final(data: Span[UInt8, _], mut budget: Budget) raises PolyError -> _FinalCAD:
    var position = 0
    if _read_word(data, position, 4, budget) != 0x31464351:
        raise _bad_final("unsupported compact CAD payload")
    var result = _get_final[_FinalCAD](data, position, budget)
    if position != len(data):
        raise _bad_final("trailing compact CAD bytes")
    return result^
