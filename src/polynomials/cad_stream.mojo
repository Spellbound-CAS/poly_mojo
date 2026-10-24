"""Depth-first CAD witnesses with independent incremental coverage checking."""
from .errors import Budget, PolyError, PolyErrorKind, _allocation, _sum, _product
from .coefficients import QQ, NoContext
from .sparse import SparsePoly
from .algebraic_types import RealAlgebraic, AlgebraicComparisonCertificate, _copy_real, _rational
from .algebraic_specialization import AlgebraicPolynomialEvaluation
from .cad_types import CADStackCertificate, CADSignWitness, CADCell
from .cad_projection import CADProjectionCertificate, cad_projection, check_cad_projection, _empty
from ._cad_lifting import _stack_cached, _check_stack, _sample_lifting
from ._cad_common import _check_signs, _signs_witnessed, _check_cell_signs
from ._algebraic import _NormalizationCache, _Frame
from ._algebraic_filter import _child_frame, _installed
from ._cad_fiber_common import _validate_fixed, _copy_fixed, _selection


@fieldwise_init
struct CADStreamEvent(Copyable, Movable):
    """A complete stack or terminal sign witness in depth-first traversal order."""
    var stack: Optional[CADStackCertificate]
    """Present for an internal node; its parent is the preceding path level."""
    var evaluations: List[AlgebraicPolynomialEvaluation]
    """Input evaluations at a terminal sample, empty for a stack event."""
    var signs: List[AlgebraicComparisonCertificate]
    """Terminal input signs, empty for a stack event."""
    var witnesses: List[CADSignWitness]
    """One sign proof per input at a terminal sample, or empty when every sign has an exact evaluation."""


@fieldwise_init
struct CADStreamReceipt(Copyable, Movable):
    """Record a completed CAD traversal requiring its separately checked witness stream."""
    var stacks: Int
    """The number of complete stack events."""
    var cells: Int
    """The number of terminal sign-cell events."""


trait CADSink:
    """Consume a CAD projection followed by depth-first local witnesses."""

    def begin(mut self, projection: CADProjectionCertificate, mut budget: Budget) raises PolyError:
        """Receive the shared projection evidence once.

        Args:
            projection: The borrowed complete projection certificate.
            budget: The shared transport allowance.

        Raises:
            PolyError: On transport or resource failure.
        """
        ...

    def write(mut self, var event: CADStreamEvent, mut budget: Budget) raises PolyError:
        """Consume the next owned depth-first event.

        Args:
            event: The complete stack or terminal sign evidence.
            budget: The shared transport allowance.

        Raises:
            PolyError: On transport or resource failure.
        """
        ...


struct _CADTraversal(Movable):
    var fixed: List[RealAlgebraic]
    var samples: List[List[RealAlgebraic]]
    var positions: List[Int]
    var point: List[Tuple[Int64, RealAlgebraic]]
    var dimension: Int
    var done: Bool

    def __init__(out self):
        self.fixed = List[RealAlgebraic]()
        self.samples = List[List[RealAlgebraic]]()
        self.positions = List[Int]()
        self.point = List[Tuple[Int64, RealAlgebraic]]()
        self.dimension = 0
        self.done = False

    def enter(mut self, stack: CADStackCertificate, variables: List[Int64], mut budget: Budget) raises PolyError:
        var selected = _selection(stack, self.fixed, budget)
        _allocation[RealAlgebraic](selected[1] - selected[0], budget)
        var owned = List[RealAlgebraic](capacity=selected[1] - selected[0])
        for i in range(selected[0], selected[1]):
            owned.append(_copy_real(stack.samples[i], budget))
        var first = _copy_real(owned[0], budget)
        _allocation[Tuple[Int64, RealAlgebraic]](1, budget)
        _allocation[Int](1, budget)
        _allocation[List[RealAlgebraic]](1, budget)
        self.point.append((variables[len(self.point)], first^))
        self.samples.append(owned^)
        self.positions.append(selected[0])
        if selected[0] % 2 == 0:
            self.dimension += 1

    def advance(mut self, variables: List[Int64], mut budget: Budget) raises PolyError:
        while len(self.point):
            budget.charge(1)
            var level = len(self.point) - 1
            _ = self.point.pop()
            if self.positions[level] % 2 == 0:
                self.dimension -= 1
            self.positions[level] += 1
            if level >= len(self.fixed) and self.positions[level] < len(self.samples[level]):
                var next = _copy_real(self.samples[level][self.positions[level]], budget)
                _allocation[Tuple[Int64, RealAlgebraic]](1, budget)
                self.point.append((variables[level], next^))
                if self.positions[level] % 2 == 0:
                    self.dimension += 1
                return
            _ = self.positions.pop()
            _ = self.samples.pop()
        self.done = True


def _basis_only(projection: CADProjectionCertificate, mut budget: Budget) raises PolyError -> CADProjectionCertificate:
    var result = _empty(projection.variables, projection.lazard, budget)
    for i in range(len(projection.levels)):
        for p in projection.levels[i]:
            _allocation[SparsePoly[QQ]](1, budget)
            result.levels[i].append(p._copy(budget))
    return result^


def cad_to[S: CADSink](inputs: List[SparsePoly[QQ]], variables: List[Int64], mut sink: S, mut budget: Budget, lazard: Bool = False) raises PolyError -> CADStreamReceipt:
    """Produce a complete CAD as depth-first stack and terminal-cell witnesses.

    Parameters:
        S: The witness destination type.

    Args:
        inputs: The original rational polynomial family.
        variables: Variables in base-to-top cylindrical order.
        sink: The incremental destination, called first with projection evidence.
        budget: The shared computation and transport allowance.
        lazard: Whether to use Lazard projection and lifting.

    Returns:
        A completion receipt requiring the separately transported witnesses.

    Raises:
        PolyError: DOMAIN for invalid inputs, UNSUPPORTED for disallowed nullification, LIMIT or CAPACITY on work limits; failure can leave an incomplete destination.
    """
    return _cad_to(inputs, variables, List[RealAlgebraic](), sink, budget, lazard)


def _cad_to[S: CADSink](inputs: List[SparsePoly[QQ]], variables: List[Int64], fixed: List[RealAlgebraic], mut sink: S, mut budget: Budget, lazard: Bool) raises PolyError -> CADStreamReceipt:
    _validate_fixed(fixed, len(variables), budget)
    var projection = cad_projection(inputs, variables, budget, lazard=lazard)
    sink.begin(projection, budget)
    # The header has been transported. Lifting needs its basis, not its proofs.
    projection = _basis_only(projection, budget)
    var traversal = _CADTraversal()
    traversal.fixed = _copy_fixed(fixed, budget)
    var receipt = CADStreamReceipt(0, 0)
    # frames[k] writes the first k path coordinates in one generator; positions
    # records the stack sample each extension used, and liftings[k] the lifting
    # polynomial of each root of the level-k stack when frames[k] has a generator.
    var frames: List[_Frame] = [_Frame()]
    var positions = List[Int]()
    var liftings = List[List[Optional[SparsePoly[QQ]]]]()
    while not traversal.done:
        budget.charge(1)
        var level = len(traversal.point)
        var kept = 0
        while kept < len(positions) and kept < level and positions[kept] == traversal.positions[kept]:
            kept += 1
        while len(positions) > kept:
            _ = positions.pop()
            _ = frames.pop()
        for j in range(kept, level):
            budget.charge(1)
            var position = traversal.positions[j]
            var lifting: Optional[SparsePoly[QQ]] = None
            if position % 2 == 1 and j < len(liftings) and position // 2 < len(liftings[j]) and liftings[j][position // 2]:
                lifting = liftings[j][position // 2].value()._copy(budget)
            # As in the retained decomposition: a terminal cell whose sample
            # would be a second irrational coordinate is signed by its own
            # zero proofs, and needs no frame.
            var irrational = 0
            for e in range(j):
                if not traversal.point[e][1].rational:
                    irrational += 1
            if frames[j].generator and not traversal.point[j][1].rational and j + 1 == len(variables) and irrational < 2:
                var skipped = _Frame()
                skipped.valid = False
                frames.append(skipped^)
            else:
                frames.append(_child_frame(frames[j], variables[j], lifting, traversal.point[j][1], budget, terminal=j + 1 == len(variables)))
            positions.append(position)
        var cache = _NormalizationCache()
        cache.frame = _installed(frames[level], traversal.point, budget)
        if level < len(variables):
            var stack = _stack_cached(projection, level - 1, level, traversal.point, traversal.dimension, None, cache, budget)
            var roots = List[Optional[SparsePoly[QQ]]]()
            if frames[level].valid and frames[level].generator:
                for gap in range(len(stack.roots)):
                    roots.append(_sample_lifting(projection, stack, traversal.point, _sum(_product(2, gap), 1), budget))
            while len(liftings) > level:
                _ = liftings.pop()
            liftings.append(roots^)
            traversal.enter(stack, variables, budget)
            receipt.stacks = _sum(receipt.stacks, 1)
            sink.write(CADStreamEvent(stack^, List[AlgebraicPolynomialEvaluation](), List[AlgebraicComparisonCertificate](), List[CADSignWitness]()), budget)
        else:
            var evaluations = List[AlgebraicPolynomialEvaluation]()
            var signs = List[AlgebraicComparisonCertificate]()
            var witnesses = List[CADSignWitness]()
            _signs_witnessed(inputs, traversal.point, evaluations, signs, witnesses, cache, budget)
            receipt.cells = _sum(receipt.cells, 1)
            sink.write(CADStreamEvent(None, evaluations^, signs^, witnesses^), budget)
            traversal.advance(variables, budget)
    return receipt^


@fieldwise_init
struct CADStreamChecker(CADSink, Movable):
    """Check CAD coverage and signs while retaining only the active sample path."""
    var _inputs: List[SparsePoly[QQ]]
    var _projection: Optional[CADProjectionCertificate]
    var _traversal: _CADTraversal
    var _stacks: Int
    var _cells: Int
    var _failed: Bool

    @staticmethod
    def prepare(inputs: List[SparsePoly[QQ]], mut budget: Budget) raises PolyError -> Self:
        """Prepare checking for the original polynomial family.

        Args:
            inputs: The original rational polynomials.
            budget: The shared preparation allowance.

        Returns:
            A checker awaiting one projection header and its traversal.

        Raises:
            PolyError: DOMAIN for malformed input, LIMIT or CAPACITY on work limits.
        """
        _allocation[SparsePoly[QQ]](len(inputs), budget)
        var owned = List[SparsePoly[QQ]](capacity=len(inputs))
        for p in inputs:
            p.check_invariants(budget)
            owned.append(p._copy(budget))
        return Self(owned^, None, _CADTraversal(), 0, 0, False)

    def begin(mut self, projection: CADProjectionCertificate, mut budget: Budget) raises PolyError:
        """Validate shared projection evidence and retain only its checked basis.

        Args:
            projection: The complete projection header.
            budget: The shared checking allowance.

        Raises:
            PolyError: DOMAIN for invalid or repeated headers, LIMIT or CAPACITY on work limits; failure invalidates the checker.
        """
        if self._failed or self._projection:
            self._failed = True
            raise PolyError(PolyErrorKind.DOMAIN, "repeated CAD stream header")
        self._failed = True
        _validate_fixed(self._traversal.fixed, len(projection.variables), budget)
        if not check_cad_projection(self._inputs, projection, budget):
            raise PolyError(PolyErrorKind.DOMAIN, "invalid CAD stream projection")
        self._projection = _basis_only(projection, budget)
        self._failed = False

    def write(mut self, var event: CADStreamEvent, mut budget: Budget) raises PolyError:
        """Check a producer event directly as a streaming destination.

        Args:
            event: The next owned local certificate.
            budget: The shared checking allowance.

        Raises:
            PolyError: DOMAIN for rejected evidence, LIMIT or CAPACITY on work limits.
        """
        if not self.accept(event^, budget):
            raise PolyError(PolyErrorKind.DOMAIN, "invalid CAD stream event")

    def accept(mut self, var event: CADStreamEvent, mut budget: Budget) raises PolyError -> Bool:
        """Check and discard one complete stack or terminal sign witness.

        Args:
            event: The next owned depth-first event.
            budget: The shared checking allowance.

        Returns:
            Whether the event is valid; rejection permanently invalidates the checker.

        Raises:
            PolyError: LIMIT or CAPACITY on work limits; failure invalidates the checker.
        """
        if self._failed:
            return False
        self._failed = True
        if not self._projection or self._traversal.done:
            return False
        try:
            ref projection = self._projection.value()
            var level = len(self._traversal.point)
            if level < len(projection.variables):
                if not event.stack or len(event.evaluations) or len(event.signs) or len(event.witnesses) or not _check_stack(projection, level - 1, level, self._traversal.point, self._traversal.dimension, event.stack.value(), budget):
                    return False
                self._traversal.enter(event.stack.value(), projection.variables, budget)
                self._stacks = _sum(self._stacks, 1)
            else:
                if event.stack:
                    return False
                if len(event.witnesses) == 0:
                    if not _check_signs(self._inputs, self._traversal.point, event.evaluations, event.signs, budget):
                        return False
                else:
                    # The cell's level and sample are not read by the sign check.
                    var sample = self._traversal.point[level - 1][1].copy() if level else _rational(QQ.zero(NoContext()), budget)
                    var cell = CADCell(-1, level - 1, 0, sample^, event.evaluations.copy(), event.signs.copy(), event.witnesses.copy())
                    if not _check_cell_signs(self._inputs, self._traversal.point, cell, budget):
                        return False
                self._traversal.advance(projection.variables, budget)
                self._cells = _sum(self._cells, 1)
        except e:
            if e.kind == PolyErrorKind.DOMAIN:
                return False
            raise e^
        self._failed = False
        return True

    def finish(self, receipt: CADStreamReceipt, mut budget: Budget) raises PolyError -> Bool:
        """Verify exhaustion of every stack and the terminal receipt counts.

        Args:
            receipt: The producer's completed traversal receipt.
            budget: The shared checking allowance.

        Returns:
            True only when every branch and terminal sign required by the prepared scope was checked.

        Raises:
            PolyError: LIMIT on exhaustion.
        """
        budget.charge(1)
        return not self._failed and Bool(self._projection) and self._traversal.done and receipt.stacks == self._stacks and receipt.cells == self._cells
