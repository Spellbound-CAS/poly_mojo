"""Exact complex roots by real-coordinate elimination and checked reconstruction."""
from ._archive_meta import _Archive
from ._algebraic_io import _Reader, _Writer
from polynomials import PolynomialTextSource, PolynomialTextSink, BallDomain, QQ
from ._common import _grow
from std.memory import ArcPointer
from std.iter import Iterable, Iterator, StopIteration
from apn_mojo import ExactComplex
from polynomials import Budget, PolyError, PolyErrorKind, NoContext, QQi, CoefficientField, DensePoly, SparsePoly, Term, Monomial, MonomialOrder
from polynomials import resultant, divmod, squarefree, isolate_real_roots, real_root, real_rational, algebraic_roots, RealAlgebraic, algebraic_compare, check_algebraic_comparison, check_real_algebraic
from ._algebraic_value import Real, Complex, _plain, _copy_witness
from ._common import _allocation, _sum, _charge_size
from ._interop import _dense_copy, _take_dense


@fieldwise_init
struct ComplexRootsCertificate(Copyable, Movable):
    """Exact root witnesses and positive multiplicities for full reconstruction."""
    var roots: List[Complex]
    """Distinct roots in lexicographic order of real and imaginary components."""
    var multiplicities: List[Int]
    """One positive multiplicity per root; their sum is the source degree."""


@fieldwise_init
struct _Data(Copyable, Movable):
    var polynomial: DensePoly[Complex]
    var proof: ComplexRootsCertificate


@fieldwise_init
struct ComplexRoots( _Archive, Copyable, Movable, Iterable, Sized):
    """All distinct exact complex roots, with multiplicities and completeness evidence.

    Indexing and iteration return algebraic.Complex values. Copies share the
    immutable source and certificate. Constants yield an empty result.
    """
    var _data: ArcPointer[_Data]
    comptime IteratorType[iterable_mut: Bool, //, iterable_origin: Origin[mut=iterable_mut]]: Iterator = _Iterator[ImmOrigin(iterable_origin)]
    """A borrowing iterator over exact root values.

    Parameters:
        iterable_mut: Whether the source origin permits mutation.
        iterable_origin: The source lifetime.
    """

    def __len__(self) -> Int:
        """Return the number of distinct roots.

        Returns:
            The root count, excluding multiplicity.
        """
        return len(self._data[].proof.roots)

    def __getitem__(self, index: Int) -> Complex:
        """Return an exact root at a valid zero-based index.

        Args:
            index: The position in lexicographic component order.

        Returns:
            An immutable value sharing its component witnesses.
        """
        return self._data[].proof.roots[index].copy()

    def __iter__(ref self) -> Self.IteratorType[origin_of(self)]:
        """Iterate over distinct exact roots in deterministic order.

        Returns:
            A borrowing iterator over the result.
        """
        return _Iterator[ImmOrigin(origin_of(self))](Pointer(to=self).as_imm(), 0)

    def multiplicity(self, index: Int) -> Int:
        """Return the multiplicity of the root at a valid index.

        Args:
            index: The root's zero-based position.

        Returns:
            A positive algebraic multiplicity.
        """
        return self._data[].proof.multiplicities[index]

    def evidence(self) -> ref[ImmOrigin(origin_of(self._data[].proof))] ComplexRootsCertificate:
        """Borrow the immutable roots and multiplicities certificate.

        Returns:
            Evidence verified by roots.check or roots.check_complex.
        """
        return self._data[].proof



    @staticmethod
    def _archive_shape() -> Int:
        return -1

    def _archive_check(self, mut budget: Budget) raises PolyError -> Bool:
        return _check(self._data[].polynomial, self._data[].proof, budget)

    def _archive_write[S: PolynomialTextSink](self, mut writer: _Writer, mut sink: S, mut budget: Budget, depth: Int) raises PolyError:
        from . import _result_io
        comptime _ArchiveRing = Complex
        comptime _ArchiveDomain = BallDomain
        comptime _ArchiveProof = Int
        _result_io._depth(depth, budget)
        _result_io._write_node[_Data, _ArchiveRing, _ArchiveDomain, _ArchiveProof](self._data[], writer, sink, budget, depth + 1)

    @staticmethod
    def _archive_read[S: PolynomialTextSource](mut reader: _Reader[S], mut budget: Budget, depth: Int) raises PolyError -> Self:
        from . import _result_io
        comptime _ArchiveRing = Complex
        comptime _ArchiveDomain = BallDomain
        comptime _ArchiveProof = Int
        _result_io._depth(depth, budget)
        var field701 = _result_io._read_node[_Data, _ArchiveRing, _ArchiveDomain, _ArchiveProof](reader, budget, depth + 1)
        _allocation[_Data](1, budget)
        var field700 = ArcPointer(field701^)
        return Self(field700^)

@fieldwise_init
struct _Iterator[origin: ImmOrigin](Iterator):
    comptime Element = Complex
    var source: Pointer[ComplexRoots, Self.origin]
    var position: Int

    def __next__(mut self) raises StopIteration -> Complex:
        if self.position == len(self.source[]):
            raise StopIteration()
        var value = self.source[][self.position]
        self.position += 1
        return value^


def _compare_real(a: Real, b: Real, mut budget: Budget) raises PolyError -> Int:
    var proof = algebraic_compare(a.witness(), b.witness(), budget)
    if not check_algebraic_comparison(a.witness(), b.witness(), proof, budget):
        raise PolyError(PolyErrorKind.DOMAIN, "invalid algebraic ordering evidence")
    return proof.relation


def _compare(a: Complex, b: Complex, mut budget: Budget) raises PolyError -> Int:
    var relation = _compare_real(a.real(), b.real(), budget)
    return relation if relation else _compare_real(a.imag(), b.imag(), budget)


def _real[C: CoefficientField](c: C, mut budget: Budget) raises PolyError -> Real:
    comptime if C == QQ:
        return _plain(real_rational(rebind[QQ](c), budget), budget)
    else:
        comptime assert C == Real, "complex root coordinates require QQ or AA"
        _charge_size(c.size(), budget)
        return rebind[Real](c).copy()


def _sort(mut proof: ComplexRootsCertificate, mut budget: Budget) raises PolyError:
    for i in range(1, len(proof.roots)):
        budget.charge(1)
        var j = i
        while j > 0 and _compare(proof.roots[j], proof.roots[j - 1], budget) < 0:
            budget.charge(1)
            var saved = proof.roots[j].copy()
            proof.roots[j] = proof.roots[j - 1].copy()
            proof.roots[j - 1] = saved^
            var saved_multiplicity = proof.multiplicities[j]
            proof.multiplicities[j] = proof.multiplicities[j - 1]
            proof.multiplicities[j - 1] = saved_multiplicity
            j -= 1


def _real_roots[C: CoefficientField](p: DensePoly[C], mut budget: Budget) raises PolyError -> List[Real]:
    var values = List[Real]()
    if p.degree() == 0:
        return values^
    comptime if C == QQ:
        ref q = rebind[DensePoly[QQ]](p)
        var isolated = isolate_real_roots(q, budget)
        _allocation[Real](_sum(len(isolated.exact_roots), len(isolated.intervals)), budget)
        values.reserve(_sum(len(isolated.exact_roots), len(isolated.intervals)))
        for point in isolated.exact_roots:
            values.append(_plain(real_rational(point, budget), budget))
        for interval in isolated.intervals:
            values.append(_plain(real_root(q, interval.low, interval.high, budget), budget))
    else:
        comptime assert C == Real, "complex root coordinates require QQ or AA"
        _allocation[RealAlgebraic](len(p.coefficients), budget)
        var coefficients = List[RealAlgebraic](capacity=len(p.coefficients))
        for coefficient_index in range(len(p.coefficients)):
            ref c = p.coefficients[coefficient_index]
            coefficients.append(_copy_witness(rebind[Real](c).witness(), budget))
        var isolated = algebraic_roots(coefficients, budget)
        _allocation[Real](len(isolated.roots), budget)
        values.reserve(len(isolated.roots))
        for root in isolated.roots:
            values.append(_plain(_copy_witness(root, budget), budget))
    return values^


def _constant[C: CoefficientField](c: C, mut budget: Budget) raises PolyError -> SparsePoly[C]:
    _allocation[Term[C]](1, budget)
    _charge_size(c.size(), budget)
    var terms: List[Term[C]] = [Term[C](Monomial(), c.copy())]
    return SparsePoly[C].from_terms(terms^, MonomialOrder(), c.context(), budget)


def _coordinate[C: CoefficientField](axis: Int64, ctx: C.Context, mut budget: Budget) raises PolyError -> SparsePoly[C]:
    _allocation[Tuple[Int64, UInt32]](1, budget)
    _allocation[Term[C]](1, budget)
    var pairs: List[Tuple[Int64, UInt32]] = [(axis, UInt32(1))]
    var terms: List[Term[C]] = [Term[C](Monomial(pairs^, budget), C.one(ctx))]
    return SparsePoly[C].from_terms(terms^, MonomialOrder(), ctx, budget)


def _parts[C: CoefficientField](reals: List[C], imags: List[C], ctx: C.Context, mut budget: Budget) raises PolyError -> Tuple[SparsePoly[C], SparsePoly[C]]:
    var x = _coordinate[C](0, ctx, budget)
    var y = _coordinate[C](1, ctx, budget)
    var re = SparsePoly[C].zero(MonomialOrder(), ctx, budget)
    var im = SparsePoly[C].zero(MonomialOrder(), ctx, budget)
    for offset in range(len(reals)):
        budget.charge(1)
        var index = len(reals) - 1 - offset
        var next_re = re.mul(x, budget).sub(im.mul(y, budget), budget).add(_constant(reals[index], budget), budget)
        var next_im = re.mul(y, budget).add(im.mul(x, budget), budget).add(_constant(imags[index], budget), budget)
        re = next_re^
        im = next_im^
    return (re^, im^)


def _fiber[C: CoefficientField](p: SparsePoly[C], point: Real, mut budget: Budget) raises PolyError -> DensePoly[Real]:
    var degree = p.degree_in(1, budget)
    var extent = 0 if degree < 0 else _sum(degree, 1)
    _allocation[Real](extent, budget)
    var coefficients = List[Real](length=extent, fill=Real.zero(NoContext()))
    var x_degree = p.degree_in(0, budget)
    _allocation[Real](_sum(max(x_degree, 0), 1), budget)
    var powers = List[Real](capacity=_sum(max(x_degree, 0), 1))
    powers.append(Real.one(NoContext()))
    for _ in range(max(x_degree, 0)):
        powers.append(powers[len(powers) - 1].mul(point, budget))
    for term in p.terms:
        budget.charge(1)
        var xp = Int(term.monomial.exponent(0, budget))
        var yp = Int(term.monomial.exponent(1, budget))
        coefficients[yp] = coefficients[yp].add(_real(term.coefficient, budget).mul(powers[xp], budget), budget)
    return DensePoly[Real].from_coefficients(1, NoContext(), coefficients^, budget)


def _gcd(a: DensePoly[Real], b: DensePoly[Real], mut budget: Budget) raises PolyError -> DensePoly[Real]:
    var f = _dense_copy(a, 1, budget)
    var g = _dense_copy(b, 1, budget)
    while not g.is_zero():
        budget.charge(1)
        var step = divmod(f, g, budget)
        f = g^
        g = _take_dense(step.remainder)
    if not f.is_zero():
        return f.scale(f.coefficients[f.degree()].inverse(budget), budget)
    raise PolyError(PolyErrorKind.DOMAIN, "complex root specialization vanished identically")


def _append(mut proof: ComplexRootsCertificate, var value: Complex, multiplicity: Int, mut budget: Budget) raises PolyError:
    _grow(proof.roots, budget)
    _grow(proof.multiplicities, budget)
    proof.roots.append(value^)
    proof.multiplicities.append(multiplicity)


def _block[C: CoefficientField, S: ExactComplexRootSink](reals: List[C], imags: List[C], ctx: C.Context, multiplicity: Int, mut sink: S, mut checker: ExactComplexRootChecker, mut budget: Budget) raises PolyError:
    if len(reals) == 2:
        var c = Complex(_real(reals[0], budget), _real(imags[0], budget))
        var lead = Complex(_real(reals[1], budget), _real(imags[1], budget))
        _emit(sink, checker, c.neg(budget).mul(lead.inverse(budget), budget), multiplicity, budget)
        return
    var parts = _parts(reals, imags, ctx, budget)
    var elimination = resultant(parts[0], parts[1], 1, budget)
    var projection = DensePoly[C].from_sparse(elimination.resultant, budget, variable=Int64(0))
    if projection.is_zero():
        raise PolyError(PolyErrorKind.DOMAIN, "complex root projection vanished")
    var real_parts = _real_roots(projection, budget)
    for real in real_parts:
        budget.charge(1)
        var common = _gcd(_fiber(parts[0], real, budget), _fiber(parts[1], real, budget), budget)
        var imaginary_parts = _real_roots(common, budget)
        for imag in imaginary_parts:
            _emit(sink, checker, Complex(real.copy(), imag.copy()), multiplicity, budget)


def _valid(value: Complex, mut budget: Budget) raises PolyError -> Bool:
    var re = value.real()
    var im = value.imag()
    if not check_real_algebraic(re.witness(), budget) or not check_real_algebraic(im.witness(), budget):
        return False
    var normalized_re = _plain(_copy_witness(re.witness(), budget), budget)
    var normalized_im = _plain(_copy_witness(im.witness(), budget), budget)
    return normalized_re == re and normalized_im == im and normalized_re.is_zero() == re.is_zero() and normalized_im.is_zero() == im.is_zero()


def _check(p: DensePoly[Complex], proof: ComplexRootsCertificate, mut budget: Budget) raises PolyError -> Bool:
    p.check_invariants(budget)
    if p.is_zero() or len(proof.roots) != len(proof.multiplicities) or len(proof.roots) > p.degree():
        return False
    var count = 0
    for i in range(len(proof.roots)):
        budget.charge(1)
        var m = proof.multiplicities[i]
        if m <= 0 or m > p.degree() - count:
            return False
        count += m
        if not _valid(proof.roots[i], budget):
            return False
        if i and _compare(proof.roots[i - 1], proof.roots[i], budget) >= 0:
            return False
    if count != p.degree():
        return False
    for coefficient_index in range(len(p.coefficients)):
        ref c = p.coefficients[coefficient_index]
        if not _valid(c, budget):
            return False
    _allocation[Complex](1, budget)
    var coefficients: List[Complex] = [p.coefficients[p.degree()].copy()]
    var product = DensePoly[Complex].from_coefficients(0, NoContext(), coefficients^, budget)
    for i in range(len(proof.roots)):
        _allocation[Complex](2, budget)
        var entries: List[Complex] = [proof.roots[i].neg(budget), Complex.one(NoContext())]
        var factor = DensePoly[Complex].from_coefficients(0, NoContext(), entries^, budget)
        for _ in range(proof.multiplicities[i]):
            budget.charge(1)
            product = product.mul(factor, budget)
    if len(product.coefficients) != len(p.coefficients):
        return False
    for i in range(len(p.coefficients)):
        budget.charge(1)
        if not product.coefficients[i].sub(p.coefficients[i], budget).is_zero():
            return False
    return True


def _solve_to[S: ExactComplexRootSink](p: DensePoly[Complex], mut sink: S, mut budget: Budget) raises PolyError -> ExactComplexRootReceipt:
    if p.is_zero():
        raise PolyError(PolyErrorKind.DOMAIN, "the zero polynomial has no finite root list")
    var checker = _checker(p, budget)
    var rational = True
    for coefficient_index in range(len(p.coefficients)):
        ref c = p.coefficients[coefficient_index]
        budget.charge(1)
        rational = rational and c.real().witness().rational and c.imag().witness().rational
    if p.degree() > 0:
        if rational:
            _allocation[QQi](len(p.coefficients), budget)
            var coefficients = List[QQi](capacity=len(p.coefficients))
            for coefficient_index in range(len(p.coefficients)):
                ref c = p.coefficients[coefficient_index]
                _charge_size(c.size(), budget)
                coefficients.append(QQi(ExactComplex(c.real().witness().low.value, c.imag().witness().low.value)))
            var q = DensePoly[QQi].from_coefficients(0, NoContext(), coefficients^, budget)
            var blocks = squarefree(q, budget)
            for block in blocks.factors:
                _allocation[QQ](len(block[0].coefficients), budget)
                _allocation[QQ](len(block[0].coefficients), budget)
                var reals = List[QQ](capacity=len(block[0].coefficients))
                var imags = List[QQ](capacity=len(block[0].coefficients))
                for coefficient_index in range(len(block[0].coefficients)):
                    ref c = block[0].coefficients[coefficient_index]
                    _charge_size(c.size(), budget)
                    reals.append(QQ(c.value.real()))
                    imags.append(QQ(c.value.imag()))
                _block(reals, imags, NoContext(), block[1], sink, checker, budget)
        else:
            var blocks = squarefree(p, budget)
            for block in blocks.factors:
                _allocation[Real](len(block[0].coefficients), budget)
                _allocation[Real](len(block[0].coefficients), budget)
                var reals = List[Real](capacity=len(block[0].coefficients))
                var imags = List[Real](capacity=len(block[0].coefficients))
                for coefficient_index in range(len(block[0].coefficients)):
                    ref c = block[0].coefficients[coefficient_index]
                    _charge_size(c.size(), budget)
                    reals.append(c.real())
                    imags.append(c.imag())
                _block(reals, imags, NoContext(), block[1], sink, checker, budget)
    var receipt = ExactComplexRootReceipt(checker._items, checker._degree)
    if not checker.finish(receipt, budget):
        raise PolyError(PolyErrorKind.DOMAIN, "complex roots failed complete reconstruction")
    return receipt^


@fieldwise_init
struct _Collector(ExactComplexRootSink):
    var proof: ComplexRootsCertificate

    def write(mut self, var event: ExactComplexRootEvent, mut budget: Budget) raises PolyError:
        _append(self.proof, event.value.copy(), event.multiplicity, budget)


def _solve(p: DensePoly[Complex], mut budget: Budget) raises PolyError -> ComplexRoots:
    var sink = _Collector(ComplexRootsCertificate(List[Complex](), List[Int]()))
    _ = _solve_to(p, sink, budget)
    _sort(sink.proof, budget)
    _allocation[_Data](1, budget)
    var source = _dense_copy(p, 0, budget)
    var proof = ComplexRootsCertificate(List[Complex](), List[Int]())
    swap(proof, sink.proof)
    return ComplexRoots(ArcPointer(_Data(source^, proof^)))




@fieldwise_init
struct ExactComplexRootEvent(Copyable, Movable):
    """One distinct exact root and its positive multiplicity."""
    var value: Complex
    """The exact complex root, with component witnesses."""
    var multiplicity: Int
    """The root's algebraic multiplicity."""


@fieldwise_init
struct ExactComplexRootReceipt(Copyable, Movable):
    """Completion metadata; accept only after checking every root event."""
    var items: Int
    """The number of distinct roots emitted."""
    var degree: Int
    """The sum of the emitted multiplicities."""


trait ExactComplexRootSink(Movable):
    """Consume exact root events under the caller's shared work allowance."""
    def write(mut self, var event: ExactComplexRootEvent, mut budget: Budget) raises PolyError:
        """Consume one owned root witness.

        Args:
            event: The exact root and its positive multiplicity.
            budget: The shared production and consumption allowance.

        Raises:
            PolyError: Consumer failures, including LIMIT or CAPACITY.
        """
        ...


@fieldwise_init
struct ExactComplexRootChecker(ExactComplexRootSink, Movable):
    """Check roots incrementally through exact factor reconstruction.

    Completed root witnesses are discarded. The checker retains the original
    polynomial and a growing product, so its storage still depends on degree.
    """
    var _source: DensePoly[Complex]
    var _product: DensePoly[Complex]
    var _items: Int
    var _degree: Int
    var _failed: Bool

    def accept(mut self, var event: ExactComplexRootEvent, mut budget: Budget) raises PolyError -> Bool:
        """Validate and consume one distinct root; rejection invalidates the checker.

        Args:
            event: The exact root and its positive multiplicity.
            budget: The shared validation allowance.

        Returns:
            Whether the event is valid and the partial reconstruction remains sound.

        Raises:
            PolyError: LIMIT or CAPACITY on resource limits; failure invalidates the checker.
        """
        if self._failed:
            return False
        self._failed = True
        budget.charge(1)
        if event.multiplicity <= 0 or event.multiplicity > self._source.degree() - self._degree:
            return False
        if not _valid(event.value, budget):
            return False
        # The partial product vanishes exactly at previously accepted roots.
        var previous = Complex.zero(NoContext())
        for i in range(len(self._product.coefficients) - 1, -1, -1):
            previous = previous.mul(event.value, budget).add(self._product.coefficients[i], budget)
        if previous.is_zero():
            return False
        _allocation[Complex](2, budget)
        var entries: List[Complex] = [event.value.neg(budget), Complex.one(NoContext())]
        var factor = DensePoly[Complex].from_coefficients(0, NoContext(), entries^, budget)
        var expanded = _dense_copy(self._product, 0, budget)
        for _ in range(event.multiplicity):
            budget.charge(1)
            expanded = expanded.mul(factor, budget)
        # Exact division binds each partial product to the original polynomial.
        if not divmod(self._source, expanded, budget).remainder.is_zero():
            return False
        self._product = expanded^
        self._degree += event.multiplicity
        self._items += 1
        self._failed = False
        return True

    def write(mut self, var event: ExactComplexRootEvent, mut budget: Budget) raises PolyError:
        """Check an event directly as a producer sink.

        Args:
            event: The next owned root witness.
            budget: The shared allowance.

        Raises:
            PolyError: DOMAIN on rejected evidence, or LIMIT/CAPACITY on resource limits.
        """
        if not self.accept(event^, budget):
            raise PolyError(PolyErrorKind.DOMAIN, "invalid exact complex root event")

    def finish(self, receipt: ExactComplexRootReceipt, mut budget: Budget) raises PolyError -> Bool:
        """Verify complete reconstruction and the final event counts.

        Args:
            receipt: Metadata returned after successful production.
            budget: The shared final checking allowance.

        Returns:
            Whether every root is accounted for with its exact multiplicity.

        Raises:
            PolyError: LIMIT or CAPACITY on resource limits.
        """
        budget.charge(1)
        if self._failed or receipt.items != self._items or receipt.degree != self._degree or self._degree != self._source.degree():
            return False
        if len(self._source.coefficients) != len(self._product.coefficients):
            return False
        for i in range(len(self._source.coefficients)):
            if not self._source.coefficients[i].sub(self._product.coefficients[i], budget).is_zero():
                return False
        return True


def _checker(p: DensePoly[Complex], mut budget: Budget) raises PolyError -> ExactComplexRootChecker:
    p.check_invariants(budget)
    if p.is_zero():
        raise PolyError(PolyErrorKind.DOMAIN, "the zero polynomial has no finite root list")
    for coefficient_index in range(len(p.coefficients)):
        ref c = p.coefficients[coefficient_index]
        if not _valid(c, budget):
            raise PolyError(PolyErrorKind.DOMAIN, "invalid algebraic coefficient witness")
    _allocation[Complex](1, budget)
    var lead: List[Complex] = [p.coefficients[p.degree()].copy()]
    var product = DensePoly[Complex].from_coefficients(0, NoContext(), lead^, budget)
    return ExactComplexRootChecker(_dense_copy(p, 0, budget), product^, 0, 0, False)


def _emit[S: ExactComplexRootSink](mut sink: S, mut checker: ExactComplexRootChecker, var value: Complex, multiplicity: Int, mut budget: Budget) raises PolyError:
    var event = ExactComplexRootEvent(value^, multiplicity)
    checker.write(event.copy(), budget)
    sink.write(event^, budget)
