"""Checked real-algebraic coefficient embeddings into rational CAD problems."""
from std.memory import ArcPointer
from polynomials import QQ, NoContext, SparsePoly, Term, Monomial, MonomialOrder, Budget, PolyError, PolyErrorKind, RealAlgebraic, check_real_algebraic, CADCell
from ._algebraic_value import Real, _plain, _copy_witness
from ._common import _allocation, _grow, _sum, _product, _charge_size


@fieldwise_init
struct _Embedding(Copyable, Movable):
    var originals: List[SparsePoly[Real]]
    var roots: List[Real]
    var width: Int


@fieldwise_init
struct _Prepared(Movable):
    var inputs: List[SparsePoly[QQ]]
    var embedding: Optional[ArcPointer[_Embedding]]


def _validate(value: Real, mut budget: Budget) raises PolyError -> Real:
    if not check_real_algebraic(value.witness(), budget):
        raise PolyError(PolyErrorKind.DOMAIN, "invalid CAD exact-real witness")
    var normalized = _plain(_copy_witness(value.witness(), budget), budget)
    _charge_size(_sum(value.size(), normalized.size()), budget)
    if normalized != value or normalized.is_zero() != value.is_zero():
        raise PolyError(PolyErrorKind.DOMAIN, "invalid CAD exact-real equality cache")
    return normalized^


def _encode(var originals: List[SparsePoly[Real]], width: Int, mut budget: Budget) raises PolyError -> _Prepared:
    if width < 0:
        raise PolyError(PolyErrorKind.DOMAIN, "negative CAD axis count")
    var roots = List[Real]()
    var inputs = List[SparsePoly[QQ]]()
    for p in originals:
        p.check_invariants(budget)
        var terms = List[Term[QQ]]()
        for term in p.terms:
            var value = _validate(term.coefficient, budget)
            _allocation[Tuple[Int64, UInt32]](_sum(len(term.monomial.pairs), 1), budget)
            var pairs = List[Tuple[Int64, UInt32]](capacity=_sum(len(term.monomial.pairs), 1))
            for pair in term.monomial.pairs:
                budget.charge(1)
                if pair[0] < 0 or pair[0] >= Int64(width):
                    raise PolyError(PolyErrorKind.DOMAIN, "CAD coefficient polynomial has an unmapped axis")
                pairs.append(pair)
            var coefficient = QQ.one(NoContext())
            ref minimal = value._canonical[].polynomial
            if minimal.degree() == 1:
                coefficient = minimal.coefficients[0].neg(budget)
            else:
                var index = -1
                for i in range(len(roots)):
                    _charge_size(_sum(value.size(), roots[i].size()), budget)
                    if roots[i] == value:
                        index = i
                        break
                    var opposite = roots[i].neg(budget)
                    _charge_size(_sum(value.size(), opposite.size()), budget)
                    if opposite == value:
                        index = i
                        coefficient = QQ.one(NoContext()).neg(budget)
                        break
                if index < 0:
                    index = len(roots)
                    _grow(roots, budget)
                    roots.append(value.copy())
                pairs.append((Int64(_sum(width, index)), UInt32(1)))
            _grow(terms, budget)
            terms.append(Term[QQ](Monomial(pairs^, budget), coefficient^))
        _grow(inputs, budget)
        inputs.append(SparsePoly[QQ].from_terms(terms^, MonomialOrder(), NoContext(), budget))
    for i in range(len(roots)):
        ref root = roots[i].witness()
        var axis = Int64(_sum(width, i))
        var p = root.polynomial.to_sparse(MonomialOrder(), budget)
        # Witness polynomials use variable zero; map it to the auxiliary axis.
        for j in range(len(p.terms)):
            if len(p.terms[j].monomial.pairs):
                p.terms[j].monomial.pairs[0][0] = axis
        p = p.reorder(MonomialOrder(), budget)
        _grow(inputs, budget)
        inputs.append(p^)
        var x = SparsePoly[QQ].variable(axis, MonomialOrder(), NoContext(), budget)
        _grow(inputs, budget)
        inputs.append(x.sub(SparsePoly[QQ].constant(root.low, MonomialOrder(), budget), budget))
        _grow(inputs, budget)
        inputs.append(SparsePoly[QQ].constant(root.high, MonomialOrder(), budget).sub(x, budget))
    _allocation[_Embedding](1, budget)
    return _Prepared(inputs^, ArcPointer(_Embedding(originals^, roots^, width)))


def _axes(variables: List[Int64], embedding: Optional[ArcPointer[_Embedding]], mut budget: Budget) raises PolyError -> List[Int64]:
    var count = len(embedding.value()[].roots) if embedding else 0
    _allocation[Int64](_sum(len(variables), count), budget)
    var result = List[Int64](capacity=_sum(len(variables), count))
    if embedding:
        ref data = embedding.value()[]
        for i in range(count):
            result.append(Int64(_sum(data.width, i)))
    for axis in variables:
        budget.charge(1)
        result.append(axis)
    return result^


def _masks(mut masks: List[UInt8], embedding: Optional[ArcPointer[_Embedding]], mut budget: Budget) raises PolyError:
    if embedding:
        for _ in range(len(embedding.value()[].roots)):
            _grow(masks, budget)
            masks.append(2)
            _grow(masks, budget)
            masks.append(4)
            _grow(masks, budget)
            masks.append(4)


def _query(var query: List[RealAlgebraic], embedding: Optional[ArcPointer[_Embedding]], mut budget: Budget) raises PolyError -> List[RealAlgebraic]:
    if not embedding:
        return query^
    ref data = embedding.value()[]
    _allocation[RealAlgebraic](_sum(len(query), len(data.roots)), budget)
    var result = List[RealAlgebraic](capacity=_sum(len(query), len(data.roots)))
    for root in data.roots:
        result.append(_copy_witness(root.witness(), budget))
    for point in query:
        result.append(_copy_witness(point, budget))
    return result^


def _check(inputs: List[SparsePoly[QQ]], embedding: Optional[ArcPointer[_Embedding]], width: Int, mut budget: Budget, *, tail: Int = 0) raises PolyError -> Bool:
    if not embedding:
        return True
    ref data = embedding.value()[]
    if data.width != width:
        return False
    _allocation[SparsePoly[Real]](len(data.originals), budget)
    var originals = List[SparsePoly[Real]](capacity=len(data.originals))
    for p in data.originals:
        originals.append(p.reorder(p.order, budget))
    var rebuilt = _encode(originals^, width, budget)
    ref roots = rebuilt.embedding.value()[].roots
    if len(inputs) - tail != len(rebuilt.inputs) or len(roots) != len(data.roots):
        return False
    for i in range(len(roots)):
        var validated = _validate(data.roots[i], budget)
        _charge_size(_sum(roots[i].size(), validated.size()), budget)
        if roots[i] != validated:
            return False
    for i in range(len(rebuilt.inputs)):
        if not inputs[i].sub(rebuilt.inputs[i], budget).is_zero():
            return False
    return True


def _check_axes(variables: List[Int64], embedding: Optional[ArcPointer[_Embedding]], width: Int, mut budget: Budget) raises PolyError -> Bool:
    var count = len(embedding.value()[].roots) if embedding else 0
    if len(variables) != _sum(width, count):
        return False
    for i in range(count):
        budget.charge(1)
        if variables[i] != Int64(_sum(width, i)):
            return False
    for i in range(count, len(variables)):
        budget.charge(1)
        if variables[i] < 0 or variables[i] >= Int64(width):
            return False
        for j in range(count, i):
            budget.charge(1)
            if variables[i] == variables[j]:
                return False
    return True


def _check_masks(masks: List[UInt8], embedding: Optional[ArcPointer[_Embedding]], mut budget: Budget, *, tail: Int = 0) raises PolyError -> Bool:
    if not embedding:
        return True
    ref data = embedding.value()[]
    if len(masks) - tail != _sum(len(data.originals), _product(3, len(data.roots))):
        return False
    for i in range(len(data.roots)):
        budget.charge(3)
        var offset = len(data.originals) + 3 * i
        if masks[offset] != 2 or masks[offset+1] != 4 or masks[offset+2] != 4:
            return False
    return True


def _on_fiber(cell: CADCell, embedding: Optional[ArcPointer[_Embedding]], mut budget: Budget, *, tail: Int = 0) raises PolyError -> Bool:
    if not embedding:
        return True
    ref data = embedding.value()[]
    if len(cell.signs) - tail != _sum(len(data.originals), _product(3, len(data.roots))):
        raise PolyError(PolyErrorKind.DOMAIN, "CAD cell has the wrong sign count")
    for i in range(len(data.roots)):
        budget.charge(3)
        var offset = len(data.originals) + 3 * i
        if cell.signs[offset].relation != 0 or cell.signs[offset+1].relation != 1 or cell.signs[offset+2].relation != 1:
            return False
    return True


def _specialize(p: SparsePoly[QQ], embedding: Optional[ArcPointer[_Embedding]], mut budget: Budget) raises PolyError -> SparsePoly[Real]:
    if not embedding:
        raise PolyError(PolyErrorKind.DOMAIN, "missing algebraic CAD embedding")
    ref data = embedding.value()[]
    var terms = List[Term[Real]]()
    for term in p.terms:
        var value = _plain_real(term.coefficient, budget)
        var powers = List[Tuple[Int64, UInt32]]()
        for pair in term.monomial.pairs:
            budget.charge(1)
            if pair[0] < 0:
                raise PolyError(PolyErrorKind.DOMAIN, "negative CAD variable")
            if pair[0] < Int64(data.width):
                _grow(powers, budget)
                powers.append(pair)
            else:
                var i = Int(pair[0]) - data.width
                if i >= len(data.roots):
                    raise PolyError(PolyErrorKind.DOMAIN, "unknown CAD auxiliary variable")
                var power = Real.one(NoContext())
                var base = data.roots[i].copy()
                _charge_size(base.size(), budget)
                var exponent = pair[1]
                while exponent:
                    if exponent & 1:
                        power = power.mul(base, budget)
                    exponent >>= 1
                    if exponent:
                        base = base.mul(base, budget)
                value = value.mul(power, budget)
        _grow(terms, budget)
        terms.append(Term[Real](Monomial(powers^, budget), value^))
    return SparsePoly[Real].from_terms(terms^, MonomialOrder(), NoContext(), budget)


def _plain_real(q: QQ, mut budget: Budget) raises PolyError -> Real:
    from polynomials import real_rational
    return _plain(real_rational(q, budget), budget)


def _assignments(point: List[Tuple[Int64, RealAlgebraic]], embedding: Optional[ArcPointer[_Embedding]], mut budget: Budget) raises PolyError -> List[Tuple[Int64, RealAlgebraic]]:
    var count = len(embedding.value()[].roots) if embedding else 0
    _allocation[Tuple[Int64, RealAlgebraic]](_sum(count, len(point)), budget)
    var result = List[Tuple[Int64, RealAlgebraic]](capacity=_sum(count, len(point)))
    if embedding:
        ref data = embedding.value()[]
        for i in range(count):
            result.append((Int64(_sum(data.width, i)), _copy_witness(data.roots[i].witness(), budget)))
    for entry in point:
        result.append((entry[0], _copy_witness(entry[1], budget)))
    return result^


def _same_root(value: RealAlgebraic, root: Real, mut budget: Budget) raises PolyError -> Bool:
    if not check_real_algebraic(value, budget):
        return False
    var normalized = _plain(_copy_witness(value, budget), budget)
    _charge_size(_sum(normalized.size(), root.size()), budget)
    return normalized == root


def _check_query(point: List[RealAlgebraic], embedding: Optional[ArcPointer[_Embedding]], mut budget: Budget) raises PolyError -> Bool:
    if not embedding:
        return True
    ref data = embedding.value()[]
    if len(point) != _sum(data.width, len(data.roots)):
        return False
    for i in range(len(data.roots)):
        if not _same_root(point[i], data.roots[i], budget):
            return False
    return True


def _check_lazard(polynomial: SparsePoly[QQ], point: List[Tuple[Int64, RealAlgebraic]], embedding: Optional[ArcPointer[_Embedding]], width: Int, mut budget: Budget) raises PolyError -> Bool:
    if not embedding:
        return True
    ref data = embedding.value()[]
    if data.width != width or len(data.originals) != 1 or len(point) < len(data.roots):
        return False
    _allocation[SparsePoly[Real]](1, budget)
    var original = data.originals[0].reorder(data.originals[0].order, budget)
    var originals = List[SparsePoly[Real]]()
    originals.append(original^)
    var prepared = _encode(originals^, width, budget)
    if not _check(prepared.inputs, embedding, width, budget):
        return False
    if not polynomial.sub(prepared.inputs[0], budget).is_zero():
        return False
    for i in range(len(data.roots)):
        if point[i][0] != Int64(_sum(width, i)) or not _same_root(point[i][1], data.roots[i], budget):
            return False
    for i in range(len(data.roots), len(point)):
        budget.charge(1)
        if point[i][0] < 0 or point[i][0] >= Int64(width):
            return False
    return True


from polynomials import CADFiber, CADSink, CADStreamChecker, CADProjectionCertificate, CADStreamEvent, CADStreamReceipt


@fieldwise_init
struct _Checker(CADSink, Movable):
    var core: CADStreamChecker
    var embedding: Optional[ArcPointer[_Embedding]]
    var width: Int
    var failed: Bool

    def begin(mut self, projection: CADProjectionCertificate, mut budget: Budget) raises PolyError:
        if self.failed:
            raise PolyError(PolyErrorKind.DOMAIN, "failed algebraic CAD stream checker")
        self.failed = True
        if not _check_axes(projection.variables, self.embedding, self.width, budget):
            raise PolyError(PolyErrorKind.DOMAIN, "algebraic CAD stream must fix coefficient axes first")
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
            raise PolyError(PolyErrorKind.DOMAIN, "invalid algebraic CAD stream event")

    def finish(self, receipt: CADStreamReceipt, mut budget: Budget) raises PolyError -> Bool:
        return not self.failed and self.core.finish(receipt, budget)


def _fiber(embedding: Optional[ArcPointer[_Embedding]], mut budget: Budget) raises PolyError -> CADFiber:
    var fixed = List[RealAlgebraic]()
    if embedding:
        ref roots = embedding.value()[].roots
        _allocation[RealAlgebraic](len(roots), budget)
        fixed = List[RealAlgebraic](capacity=len(roots))
        for root in roots:
            fixed.append(_copy_witness(root.witness(), budget))
    return CADFiber(fixed, budget)
