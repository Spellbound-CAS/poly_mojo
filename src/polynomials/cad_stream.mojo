"""Depth-first CAD witnesses with independent incremental coverage checking."""
from .errors import Budget, PolyError, PolyErrorKind, _allocation, _sum
from .coefficients import QQ
from .sparse import SparsePoly
from .algebraic_types import RealAlgebraic, AlgebraicComparisonCertificate, _copy_real
from .algebraic_specialization import AlgebraicPolynomialEvaluation
from .cad_types import CADStackCertificate
from .cad_projection import CADProjectionCertificate, cad_projection, check_cad_projection, _empty
from ._cad_lifting import _stack, _check_stack
from ._cad_common import _signs, _check_signs
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
    while not traversal.done:
        budget.charge(1)
        var level = len(traversal.point)
        if level < len(variables):
            var stack = _stack(projection, level - 1, level, traversal.point, traversal.dimension, None, budget)
            traversal.enter(stack, variables, budget)
            receipt.stacks = _sum(receipt.stacks, 1)
            sink.write(CADStreamEvent(stack^, List[AlgebraicPolynomialEvaluation](), List[AlgebraicComparisonCertificate]()), budget)
        else:
            var evaluations = List[AlgebraicPolynomialEvaluation]()
            var signs = List[AlgebraicComparisonCertificate]()
            _signs(inputs, traversal.point, evaluations, signs, budget)
            receipt.cells = _sum(receipt.cells, 1)
            sink.write(CADStreamEvent(None, evaluations^, signs^), budget)
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
                if not event.stack or len(event.evaluations) or len(event.signs) or not _check_stack(projection, level - 1, level, self._traversal.point, self._traversal.dimension, event.stack.value(), budget):
                    return False
                self._traversal.enter(event.stack.value(), projection.variables, budget)
                self._stacks = _sum(self._stacks, 1)
            else:
                if event.stack or not _check_signs(self._inputs, self._traversal.point, event.evaluations, event.signs, budget):
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
