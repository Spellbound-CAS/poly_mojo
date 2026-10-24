"""Construct CAD data and accept the complete claim before returning it.

The compact backend keeps field-relative sections and no intermediate proof
histories. Explicit McCallum requests and restricted coordinate prefixes use
the compatible legacy backend, with its complete checker before return.
"""
from std.memory import ArcPointer
from .coefficients import QQ
from .sparse import SparsePoly
from .algebraic_types import RealAlgebraic
from .errors import Budget, PolyError, PolyErrorKind, _allocation, _sum
from .cad_types import CADCertificate, CADCoveringCertificate
from .cad_fiber import CADFiber
from .cad_projection import cad_projection
from .cad_final_types import CADFinalCell
from ._cad_tree_final import _FinalCAD, _verified_final_cad, _check_final_cad
from ._cad_final_output import _materialize_final_cells, _final_cell_count, _final_witness_point
from ._cad_final_codec import _encode_final, _decode_final
from ._cad_common import _coordinates


@fieldwise_init
struct _CADFinalData(Copyable, Movable):
    var compact: Optional[_FinalCAD]
    var legacy: Optional[CADCertificate]
    var covering: Optional[CADCoveringCertificate]
    var query: Bool


@fieldwise_init
struct CADFinalCertificate(Copyable, Movable):
    """Immutable retained evidence for a complete decomposition or conjunction.

    Inspect geometry through cells(), and metadata through variables(),
    factor_counts(), cell_count() and stack_count(). Compact evidence keeps
    relative roots until cells() requests absolute coordinates. A legacy
    payload supports explicit McCallum, fixed prefixes and old archives.
    Independent acceptance always binds the original inputs, variable order,
    fixed prefix and query through check_cad_final.
    """
    var _data: ArcPointer[_CADFinalData]

    def is_compact(self) -> Bool:
        """Report whether the retained representation uses relative sections.

        Returns:
            True for compact final evidence, False for compatible legacy evidence.
        """
        return Bool(self._data[].compact)

    def is_covering(self) -> Bool:
        """Report whether the claim includes a sign conjunction.

        Returns:
            True for a full-CAD conjunction claim.
        """
        return self._data[].query

    def variables(self, mut budget: Budget) raises PolyError -> List[Int64]:
        """Copy the cylindrical variable order.

        Args:
            budget: Shared storage and work allowance.

        Returns:
            Variables in base-to-top order.

        Raises:
            PolyError: DOMAIN for missing evidence, LIMIT/CAPACITY on resource limits.
        """
        ref data = self._data[]
        if data.compact:
            _allocation[Int64](len(data.compact.value().projection.variables), budget)
            return data.compact.value().projection.variables.copy()
        if data.legacy:
            _allocation[Int64](len(data.legacy.value().projection.variables), budget)
            return data.legacy.value().projection.variables.copy()
        if data.covering:
            _allocation[Int64](len(data.covering.value().decomposition.projection.variables), budget)
            return data.covering.value().decomposition.projection.variables.copy()
        raise PolyError(PolyErrorKind.DOMAIN, "missing CAD final evidence")

    def factor_counts(self, mut budget: Budget) raises PolyError -> List[Int]:
        """Count the retained projection factors at each level.

        Args:
            budget: Shared storage and work allowance.

        Returns:
            Factor counts in cylindrical order.

        Raises:
            PolyError: DOMAIN for missing evidence, LIMIT/CAPACITY on resource limits.
        """
        var counts = List[Int]()
        ref data = self._data[]
        if data.compact:
            for level in data.compact.value().projection.levels:
                _allocation[Int](1, budget)
                counts.append(len(level))
        elif data.legacy:
            for level in data.legacy.value().projection.levels:
                _allocation[Int](1, budget)
                counts.append(len(level))
        elif data.covering:
            for level in data.covering.value().decomposition.projection.levels:
                _allocation[Int](1, budget)
                counts.append(len(level))
        else:
            raise PolyError(PolyErrorKind.DOMAIN, "missing CAD final evidence")
        return counts^

    def cell_count(self, mut budget: Budget) raises PolyError -> Int:
        """Count all cells, including the root, without materializing samples.

        Args:
            budget: Shared traversal allowance.

        Returns:
            The number of retained cells at all levels.

        Raises:
            PolyError: DOMAIN for missing evidence, LIMIT/CAPACITY on resource limits.
        """
        budget.charge(1)
        ref data = self._data[]
        if data.compact:
            return _final_cell_count(data.compact.value(), budget)
        if data.legacy:
            return len(data.legacy.value().cells)
        if data.covering:
            return len(data.covering.value().decomposition.cells)
        raise PolyError(PolyErrorKind.DOMAIN, "missing CAD final evidence")

    def stack_count(self) -> Int:
        """Count retained stacks without constructing output coordinates.

        Returns:
            The number of stacks, or zero for an empty invalid representation.
        """
        ref data = self._data[]
        if data.compact:
            return len(data.compact.value().stacks)
        if data.legacy:
            return len(data.legacy.value().stacks)
        if data.covering:
            return len(data.covering.value().decomposition.stacks)
        return 0

    def satisfiable(self) -> Bool:
        """Read the retained conjunction truth value.

        Returns:
            True when the retained covering contains a witness; False for a decomposition-only claim.
        """
        ref data = self._data[]
        if not data.query:
            return False
        if data.compact:
            return data.compact.value().satisfiable
        return data.covering and data.covering.value().witness_cell >= 0

    def cells(self, mut budget: Budget) raises PolyError -> List[CADFinalCell]:
        """Materialize absolute coordinates once per retained cell.

        Args:
            budget: Shared allowance including absolute root construction.

        Returns:
            Breadth-first geometry, beginning with a root placeholder. Its signs belong to the bound original inputs.

        Raises:
            PolyError: DOMAIN for missing evidence, LIMIT/CAPACITY on resource limits.
        """
        ref data = self._data[]
        if data.compact:
            return _materialize_final_cells(data.compact.value(), budget)
        if data.legacy:
            return _legacy_cells(data.legacy.value(), budget)
        if data.covering:
            return _legacy_cells(data.covering.value().decomposition, budget)
        raise PolyError(PolyErrorKind.DOMAIN, "missing CAD final evidence")

    def legacy_decomposition(self) -> ref[ImmOrigin(origin_of(self._data[].legacy))] Optional[CADCertificate]:
        """Borrow a retained legacy decomposition for compatible archive transport.

        Returns:
            The optional legacy payload; absent for compact evidence and coverings.
        """
        return self._data[].legacy

    def sample(self, mut budget: Budget) raises PolyError -> Optional[List[RealAlgebraic]]:
        """Materialize a satisfying path without converting the other cells.

        Args:
            budget: Shared coordinate construction and storage allowance.

        Returns:
            Absolute coordinates in cylindrical order, or None without a satisfying conjunction.

        Raises:
            PolyError: DOMAIN for malformed evidence, LIMIT/CAPACITY on resource limits.
        """
        if not self.satisfiable():
            return None
        ref data = self._data[]
        if data.compact:
            return _final_witness_point(data.compact.value(), budget)
        ref proof = data.covering.value()
        var path = _coordinates(proof.decomposition.cells, proof.witness_cell, proof.decomposition.projection.variables, budget)
        _allocation[RealAlgebraic](len(path[0]), budget)
        var point = List[RealAlgebraic](capacity=len(path[0]))
        for entry in path[0]:
            point.append(entry[1].copy())
        return Optional(point^)

    def legacy_covering(self) -> ref[ImmOrigin(origin_of(self._data[].covering))] Optional[CADCoveringCertificate]:
        """Borrow a retained legacy covering for compatible archive transport.

        Returns:
            The optional legacy payload; absent for compact evidence and decompositions.
        """
        return self._data[].covering

    def compact_bytes(self, mut budget: Budget) raises PolyError -> List[UInt8]:
        """Encode compact final data in the version-one QCF payload.

        Args:
            budget: Shared encoding and output storage allowance.

        Returns:
            Exact final-data bytes, excluding producer search counters. The enclosing archive must bind and verify the claim.

        Raises:
            PolyError: DOMAIN for a legacy representation, LIMIT/CAPACITY on resource limits.
        """
        if not self._data[].compact:
            raise PolyError(PolyErrorKind.DOMAIN, "legacy CAD has no compact payload")
        return _encode_final(self._data[].compact.value(), budget)

    @staticmethod
    def from_compact_bytes(data: Span[UInt8, _], inputs: List[SparsePoly[QQ]], variables: List[Int64], masks: List[UInt8], query: Bool, mut budget: Budget) raises PolyError -> Self:
        """Decode and independently accept a compact full-space claim.

        Args:
            data: Exact QCF payload with no trailing bytes.
            inputs: Original polynomials supplied by the caller.
            variables: Expected cylindrical order.
            masks: Expected conjunction masks, or empty for a decomposition.
            query: Whether the caller requests a conjunction claim.
            budget: Shared decoding and complete verification allowance.

        Returns:
            Accepted compact evidence; no producer assertions are trusted.

        Raises:
            PolyError: ENCODING for rejected evidence, LIMIT/CAPACITY on resource limits.
        """
        var proof: _FinalCAD
        try:
            proof = _decode_final(data, budget)
        except e:
            if e.kind == PolyErrorKind.DOMAIN:
                raise PolyError(PolyErrorKind.ENCODING, e.message)
            raise e^
        _allocation[_CADFinalData](1, budget)
        var result = Self(ArcPointer(_CADFinalData(proof^, None, None, query)))
        if not check_cad_final(inputs, variables, masks, result, budget, query=query):
            raise PolyError(PolyErrorKind.ENCODING, "compact CAD failed final verification")
        return result^

    @staticmethod
    def from_legacy(var proof: CADCertificate, inputs: List[SparsePoly[QQ]], fixed: List[RealAlgebraic], mut budget: Budget) raises PolyError -> Self:
        """Accept and retain an old decomposition without regenerating it.

        Args:
            proof: Owned historical decomposition evidence.
            inputs: Original polynomial family.
            fixed: Expected fixed coordinate prefix.
            budget: Shared complete verification and storage allowance.

        Returns:
            Accepted legacy evidence in the new result interface.

        Raises:
            PolyError: ENCODING for rejected evidence, LIMIT/CAPACITY on resource limits.
        """
        var fiber = CADFiber(fixed, budget)
        if not fiber.check(inputs, proof, budget):
            raise PolyError(PolyErrorKind.ENCODING, "legacy CAD failed final verification")
        _allocation[_CADFinalData](1, budget)
        return Self(ArcPointer(_CADFinalData(None, proof^, None, False)))

    @staticmethod
    def from_legacy_covering(var proof: CADCoveringCertificate, inputs: List[SparsePoly[QQ]], masks: List[UInt8], fixed: List[RealAlgebraic], mut budget: Budget) raises PolyError -> Self:
        """Accept and retain an old full covering without regenerating it.

        Args:
            proof: Owned historical covering evidence.
            inputs: Original polynomial family.
            masks: Expected conjunction masks.
            fixed: Expected fixed coordinate prefix.
            budget: Shared complete verification and storage allowance.

        Returns:
            Accepted legacy covering in the new result interface.

        Raises:
            PolyError: ENCODING for rejected evidence, LIMIT/CAPACITY on resource limits.
        """
        var fiber = CADFiber(fixed, budget)
        if not fiber.check(inputs, masks, proof, budget):
            raise PolyError(PolyErrorKind.ENCODING, "legacy covering failed final verification")
        _allocation[_CADFinalData](1, budget)
        return Self(ArcPointer(_CADFinalData(None, None, proof^, True)))


def _legacy_cells(proof: CADCertificate, mut budget: Budget) raises PolyError -> List[CADFinalCell]:
    _allocation[CADFinalCell](len(proof.cells), budget)
    var cells = List[CADFinalCell](capacity=len(proof.cells))
    for cell in proof.cells:
        _allocation[Int](len(cell.signs), budget)
        var signs = List[Int](capacity=len(cell.signs))
        for sign in cell.signs:
            signs.append(sign.relation)
        cells.append(CADFinalCell(cell.parent, cell.level, cell.position, cell.sample.copy(), signs^))
    return cells^


def _final_masks(masks: List[UInt8], mut budget: Budget) raises PolyError -> List[Int]:
    _allocation[Int](len(masks), budget)
    var output = List[Int](capacity=len(masks))
    for mask in masks:
        output.append(Int(mask))
    return output^


def check_cad_final(inputs: List[SparsePoly[QQ]], variables: List[Int64], masks: List[UInt8], certificate: CADFinalCertificate, mut budget: Budget, *, fixed: List[RealAlgebraic] = List[RealAlgebraic](), query: Bool = False, vectorized: Bool = True) raises PolyError -> Bool:
    """Independently verify a complete retained CAD against the requested scope.

    Args:
        inputs: Original polynomial family.
        variables: Expected cylindrical variable order.
        masks: Conjunction masks, or empty for decomposition only.
        certificate: Retained compact or compatible legacy evidence.
        budget: Shared complete verification allowance.
        fixed: Expected section prefix; empty means all real space.
        query: Whether the claim includes the conjunction's truth.
        vectorized: Use checked exact rational SIMD batches; False selects the scalar reference.

    Returns:
        Whether the whole certificate proves this exact requested claim.

    Raises:
        PolyError: LIMIT/CAPACITY on resource limits; neither permits acceptance.
    """
    try:
        return _check_cad_final(inputs, variables, masks, certificate, fixed, query, vectorized, budget)
    except e:
        if e.kind == PolyErrorKind.DOMAIN:
            return False
        raise e^


def _check_cad_final(inputs: List[SparsePoly[QQ]], variables: List[Int64], masks: List[UInt8], certificate: CADFinalCertificate, fixed: List[RealAlgebraic], query: Bool, vectorized: Bool, mut budget: Budget) raises PolyError -> Bool:
    ref data = certificate._data[]
    if Int(Bool(data.compact)) + Int(Bool(data.legacy)) + Int(Bool(data.covering)) != 1 or data.query != query:
        return False
    if (query and len(masks) != len(inputs)) or (not query and len(masks)):
        return False
    for mask in masks:
        budget.charge(1)
        if mask > 7:
            return False
    if certificate.variables(budget) != variables:
        return False
    if data.compact:
        if len(fixed):
            return False
        return _check_final_cad(inputs, variables, _final_masks(masks, budget), data.compact.value(), budget, vectorized=vectorized)
    var fiber = CADFiber(fixed, budget)
    if data.legacy:
        return not query and fiber.check(inputs, data.legacy.value(), budget)
    return query and fiber.check(inputs, masks, data.covering.value(), budget)


def cad_final(inputs: List[SparsePoly[QQ]], variables: List[Int64], mut budget: Budget, *, fixed: List[RealAlgebraic] = List[RealAlgebraic](), masks: List[UInt8] = List[UInt8](), query: Bool = False, lazard: Optional[Bool] = None) raises PolyError -> CADFinalCertificate:
    """Construct a full CAD and accept its final evidence before returning.

    Args:
        inputs: Original polynomial family.
        variables: All variables in base-to-top order.
        budget: Complete allowance, including construction retries and final verification.
        fixed: Section prefix; a nonempty prefix uses the checked legacy backend.
        masks: Allowed sign masks for a conjunction; empty for decomposition only.
        query: Request a full-decomposition conjunction result.
        lazard: None selects Revision 15 with interval lifting; True selects plain Lazard with interval lifting; False selects checked legacy McCallum.

    Returns:
        Only independently accepted evidence. Relative terminal sections remain unmaterialized on the compact backend.

    Raises:
        PolyError: DOMAIN for invalid input, UNSUPPORTED for legacy nullification, LIMIT/CAPACITY on resource limits.
    """
    if (query and len(masks) != len(inputs)) or (not query and len(masks)):
        raise PolyError(PolyErrorKind.DOMAIN, "CAD query mask count mismatch")
    for mask in masks:
        budget.charge(1)
        if mask > 7:
            raise PolyError(PolyErrorKind.DOMAIN, "invalid CAD sign mask")
    if len(fixed) == 0 and len(variables) and (not lazard or lazard.value()):
        var proof = _verified_final_cad(inputs, variables, _final_masks(masks, budget), budget, brown=not Bool(lazard), interval=True)
        _allocation[_CADFinalData](1, budget)
        return CADFinalCertificate(ArcPointer(_CADFinalData(proof^, None, None, query)))
    var fiber = CADFiber(fixed, budget)
    var use_lazard = lazard.value() if lazard else True
    if query:
        var proof = fiber.cover(inputs, masks, variables, budget, lazard=use_lazard)
        return CADFinalCertificate.from_legacy_covering(proof^, inputs, masks, fixed, budget)
    var proof = fiber.decompose(inputs, variables, budget, lazard=use_lazard)
    return CADFinalCertificate.from_legacy(proof^, inputs, fixed, budget)
