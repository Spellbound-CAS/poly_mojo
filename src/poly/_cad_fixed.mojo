"""Exact user-coordinate fibers over the public CADFiber boundary."""
from std.iter import Iterable, iter, next, StopIteration
from std.memory import ArcPointer
from std.builtin.rebind import downcast
from polynomials import Budget, PolyError, PolyErrorKind, QQ, SparsePoly, MonomialOrder, RealAlgebraic, CADFiber, CADSink, CADStreamChecker, CADProjectionCertificate, CADStreamEvent, CADStreamReceipt
from . import sparse, algebraic, _cad_algebraic as _aa
from .algebraic import Real, _copy_witness
from ._common import _allocation, _grow, _sum


trait _Pair:
    def _append_fixed(self, mut output: List[Tuple[Int, Real]], width: Int, mut budget: Budget) raises PolyError: ...


__extension Tuple(_Pair):
    def _append_fixed(self, mut output: List[Tuple[Int, Real]], width: Int, mut budget: Budget) raises PolyError:
        comptime assert Self.Ts.length == 2, "fixed CAD entries must be (axis, value) pairs"
        comptime assert Self.Ts[0] == Int, "fixed CAD axes must be Int"
        var axis = rebind[Int](self[0])
        budget.charge(1)
        if axis < 0 or axis >= width:
            raise PolyError(PolyErrorKind.DOMAIN, "fixed CAD axis is outside the polynomial shape")
        for entry in output:
            budget.charge(1)
            if entry[0] == axis:
                raise PolyError(PolyErrorKind.DOMAIN, "duplicate fixed CAD axis")
        comptime C = downcast[Self.Ts[1], Copyable & Deinitable]
        var value = _aa._validate(algebraic.asreal(rebind[C](self[1]), budget), budget)
        _grow(output, budget)
        output.append((axis, value^))


def _append[T: AnyType](pair: T, mut output: List[Tuple[Int, Real]], width: Int, mut budget: Budget) raises PolyError:
    comptime assert conforms_to(T, _Pair), "fixed CAD entries must be (axis, value) pairs"
    rebind[downcast[T, _Pair]](pair)._append_fixed(output, width, budget)


trait _TupleBindings:
    def _collect_fixed(self, mut output: List[Tuple[Int, Real]], width: Int, mut budget: Budget) raises PolyError: ...


__extension Tuple(_TupleBindings):
    def _collect_fixed(self, mut output: List[Tuple[Int, Real]], width: Int, mut budget: Budget) raises PolyError:
        comptime for i in range(Self.Ts.length):
            _append(self[i], output, width, budget)


def _collect[F: AnyType](fixed: F, width: Int, mut budget: Budget) raises PolyError -> List[Tuple[Int, Real]]:
    var output = List[Tuple[Int, Real]]()
    budget.charge(1)
    comptime if conforms_to(F, _TupleBindings):
        rebind[downcast[F, _TupleBindings]](fixed)._collect_fixed(output, width, budget)
    elif conforms_to(F, Iterable):
        comptime I = downcast[F, Iterable]
        var cursor = iter(rebind[I](fixed))
        while True:
            budget.charge(1)
            var item: downcast[type_of(cursor).Element, Movable & Deinitable]
            try:
                item = rebind_var[downcast[type_of(cursor).Element, Movable & Deinitable]](next(cursor))
            except StopIteration:
                break
            _append(item, output, width, budget)
    else:
        comptime assert False, "fixed must be an iterable or tuple of axis/value pairs"
    return output^


def _order[P: Copyable & Deinitable](mut fixed: List[Tuple[Int, Real]], axes: List[Int], mut budget: Budget) raises PolyError -> List[Int64]:
    var requested = sparse._axes[P](axes, budget)
    _allocation[Tuple[Int, Real]](len(fixed), budget)
    var ordered = List[Tuple[Int, Real]](capacity=len(fixed))
    _allocation[Int64](len(requested), budget)
    var variables = List[Int64](capacity=len(requested))
    # The explicit order is authoritative; otherwise use increasing fixed axes.
    for i in range(len(requested)):
        for entry in fixed:
            budget.charge(1)
            if requested[i] == Int64(entry[0]):
                if len(axes) and i >= len(fixed):
                    raise PolyError(PolyErrorKind.DOMAIN, "fixed CAD axes must lead the explicit traversal order")
                ordered.append((entry[0], entry[1].copy()))
                variables.append(requested[i])
    for axis in requested:
        var bound = False
        for entry in fixed:
            budget.charge(1)
            bound = bound or axis == Int64(entry[0])
        if not bound:
            variables.append(axis)
    fixed = ordered^
    return variables^


def _defining(axis: Int, value: Real, mut budget: Budget) raises PolyError -> SparsePoly[QQ]:
    # The caller has independently validated this equality cache. Use its
    # minimal polynomial so equivalent witnesses prepare the same CAD input.
    var p = value._canonical[].polynomial.to_sparse(MonomialOrder(), budget)
    for i in range(len(p.terms)):
        budget.charge(1)
        if len(p.terms[i].monomial.pairs):
            p.terms[i].monomial.pairs[0][0] = Int64(axis)
    return p.reorder(MonomialOrder(), budget)


def _augment(mut inputs: List[SparsePoly[QQ]], fixed: List[Tuple[Int, Real]], mut budget: Budget) raises PolyError:
    for entry in fixed:
        var p = _defining(entry[0], entry[1], budget)
        _grow(inputs, budget)
        inputs.append(p^)


def _masks(mut masks: List[UInt8], fixed: List[Tuple[Int, Real]], mut budget: Budget) raises PolyError:
    for _ in range(len(fixed)):
        _grow(masks, budget)
        masks.append(2)


def _fiber(embedding: Optional[ArcPointer[_aa._Embedding]], fixed: List[Tuple[Int, Real]], mut budget: Budget) raises PolyError -> CADFiber:
    return CADFiber(_prefix(embedding, fixed, budget), budget)


def _prefix(embedding: Optional[ArcPointer[_aa._Embedding]], fixed: List[Tuple[Int, Real]], mut budget: Budget) raises PolyError -> List[RealAlgebraic]:
    var count = len(embedding.value()[].roots) if embedding else 0
    _allocation[RealAlgebraic](_sum(count, len(fixed)), budget)
    var values = List[RealAlgebraic](capacity=_sum(count, len(fixed)))
    if embedding:
        for root in embedding.value()[].roots:
            values.append(_copy_witness(root.witness(), budget))
    for entry in fixed:
        values.append(_copy_witness(entry[1].witness(), budget))
    return values^


def _check_axes(variables: List[Int64], embedding: Optional[ArcPointer[_aa._Embedding]], fixed: List[Tuple[Int, Real]], width: Int, mut budget: Budget) raises PolyError -> Bool:
    if len(fixed) > width or not _aa._check_axes(variables, embedding, width, budget):
        return False
    var offset = len(embedding.value()[].roots) if embedding else 0
    for i in range(len(fixed)):
        budget.charge(1)
        if variables[offset + i] != Int64(fixed[i][0]):
            return False
    return True


def _check(inputs: List[SparsePoly[QQ]], embedding: Optional[ArcPointer[_aa._Embedding]], fixed: List[Tuple[Int, Real]], variables: List[Int64], width: Int, mut budget: Budget) raises PolyError -> Bool:
    if len(fixed) > len(inputs) or not _check_axes(variables, embedding, fixed, width, budget):
        return False
    try:
        if not _aa._check(inputs, embedding, width, budget, tail=len(fixed)):
            return False
        for i in range(len(fixed)):
            var value = _aa._validate(fixed[i][1], budget)
            var p = _defining(fixed[i][0], value, budget)
            if not inputs[len(inputs) - len(fixed) + i].sub(p, budget).is_zero():
                return False
        return True
    except e:
        if e.kind == PolyErrorKind.DOMAIN:
            return False
        raise e^


def _check_masks(inputs: List[SparsePoly[QQ]], masks: List[UInt8], embedding: Optional[ArcPointer[_aa._Embedding]], fixed: List[Tuple[Int, Real]], mut budget: Budget) raises PolyError -> Bool:
    if len(masks) != len(inputs) or len(fixed) > len(masks) or not _aa._check_masks(masks, embedding, budget, tail=len(fixed)):
        return False
    for i in range(len(masks) - len(fixed), len(masks)):
        budget.charge(1)
        if masks[i] != 2:
            return False
    return True


@fieldwise_init
struct _Checker(CADSink, Movable):
    var core: CADStreamChecker
    var variables: List[Int64]
    var failed: Bool

    def begin(mut self, projection: CADProjectionCertificate, mut budget: Budget) raises PolyError:
        if self.failed:
            raise PolyError(PolyErrorKind.DOMAIN, "failed fixed-coordinate CAD checker")
        self.failed = True
        if len(projection.variables) != len(self.variables):
            raise PolyError(PolyErrorKind.DOMAIN, "fixed-coordinate CAD stream has the wrong axes")
        for i in range(len(self.variables)):
            budget.charge(1)
            if projection.variables[i] != self.variables[i]:
                raise PolyError(PolyErrorKind.DOMAIN, "fixed-coordinate CAD stream has the wrong axis order")
        self.core.begin(projection, budget)
        self.failed = False

    def accept(mut self, var event: CADStreamEvent, mut budget: Budget) raises PolyError -> Bool:
        if self.failed:
            return False
        self.failed = True
        var valid = self.core.accept(event^, budget)
        self.failed = not valid
        return valid

    def write(mut self, var event: CADStreamEvent, mut budget: Budget) raises PolyError:
        if not self.accept(event^, budget):
            raise PolyError(PolyErrorKind.DOMAIN, "invalid fixed-coordinate CAD stream event")

    def finish(self, receipt: CADStreamReceipt, mut budget: Budget) raises PolyError -> Bool:
        return not self.failed and self.core.finish(receipt, budget)
