"""CAD on specified algebraic sections of the leading coordinate axes.

The ordinary CAD entry points cover all coordinates. CADFiber instead binds
an exact section prefix and proves coverage only over the remaining axes.
Certificate layouts are shared; verification must use the same prefix.
"""
from .errors import Budget, PolyError, PolyErrorKind
from .coefficients import QQ
from .sparse import SparsePoly
from .algebraic_types import RealAlgebraic
from .cad_types import CADCertificate, CADCoveringCertificate
from .cad import _cad, _check_cad
from .cad_covering import _cad_covering, _check_cad_covering
from .cylindrical_covering import CylindricalCoveringCertificate, _cylindrical_covering, _check
from .cad_coverings_search import CADCoveringsCertificate, _coverings, _check_coverings
from .cad_stream import CADSink, CADStreamReceipt, CADStreamChecker, _cad_to
from ._cad_common import _masks
from ._cad_fiber_common import _validate_fixed, _copy_fixed


struct CADFiber(Copyable, Movable):
    """Bind a zero-dimensional section prefix for CAD generation and checking.

    Each coordinate must be a section of its projection stack. Include a
    defining polynomial for each value among the inputs. Equivalent witnesses
    may use different defining polynomials or isolating intervals. Coordinates
    follow the leading variables in base-to-top order; an empty prefix covers
    all real space.
    """
    var _fixed: List[RealAlgebraic]

    def __init__(out self, fixed: List[RealAlgebraic], mut budget: Budget) raises PolyError:
        """Validate and copy the requested algebraic coordinate prefix.

        Args:
            fixed: Exact coordinates for the leading cylindrical variables.
            budget: The shared validation and storage allowance.

        Raises:
            PolyError: DOMAIN for invalid witnesses, LIMIT or CAPACITY on work limits.
        """
        _validate_fixed(fixed, len(fixed), budget)
        self._fixed = _copy_fixed(fixed, budget)

    def decompose(self, inputs: List[SparsePoly[QQ]], variables: List[Int64], mut budget: Budget, lazard: Bool = False) raises PolyError -> CADCertificate:
        """Decompose all remaining coordinates over the fixed section prefix.

        Args:
            inputs: Rational polynomials, including defining polynomials for the prefix.
            variables: All variables in base-to-top order, with fixed coordinates first.
            budget: The shared work and storage allowance.
            lazard: Whether to use Lazard projection and lifting.

        Returns:
            Evidence covering the specified fiber; check it with this CADFiber.

        Raises:
            PolyError: DOMAIN for invalid inputs or a prefix outside projection sections, UNSUPPORTED for McCallum nullification, LIMIT or CAPACITY on work limits.
        """
        return _cad(inputs, variables, self._fixed, budget, lazard)

    def cover(self, inputs: List[SparsePoly[QQ]], masks: List[UInt8], variables: List[Int64], mut budget: Budget, lazard: Bool = False) raises PolyError -> CADCoveringCertificate:
        """Decide a conjunction with a full decomposition of the requested fiber.

        Args:
            inputs: Rational polynomials, including defining polynomials for the prefix.
            masks: Allowed sign bits per input: 1 for negative, 2 for zero, 4 for positive.
            variables: All variables in base-to-top order, with fixed coordinates first.
            budget: The shared work and storage allowance.
            lazard: Whether to use Lazard projection and lifting.

        Returns:
            Evidence covering the specified fiber; check it with this CADFiber.

        Raises:
            PolyError: DOMAIN for invalid inputs or a prefix outside projection sections, UNSUPPORTED for McCallum nullification, LIMIT or CAPACITY on work limits.
        """
        return _cad_covering(inputs, masks, variables, self._fixed, budget, lazard)

    def solve(self, inputs: List[SparsePoly[QQ]], masks: List[UInt8], variables: List[Int64], mut budget: Budget, lazard: Bool = False) raises PolyError -> CylindricalCoveringCertificate:
        """Search only the requested fiber for a witness or conflict covering.

        Args:
            inputs: Rational polynomials, including defining polynomials for the prefix.
            masks: Allowed sign bits per input: 1 for negative, 2 for zero, 4 for positive.
            variables: All variables in base-to-top order, with fixed coordinates first.
            budget: The shared work and storage allowance.
            lazard: Whether to use Lazard projection and lifting.

        Returns:
            Evidence covering the specified fiber; check it with this CADFiber.

        Raises:
            PolyError: DOMAIN for invalid inputs or a prefix outside projection sections, UNSUPPORTED for McCallum nullification, LIMIT or CAPACITY on work limits.
        """
        return _cylindrical_covering(inputs, masks, variables, self._fixed, budget, lazard)

    def search(self, inputs: List[SparsePoly[QQ]], masks: List[UInt8], variables: List[Int64], mut budget: Budget, lazard: Bool = False) raises PolyError -> CADCoveringsCertificate:
        """Decide the conjunction on the fiber by a conflict-driven search over cylindrical coverings.

        Args:
            inputs: Rational polynomials over all variables.
            masks: Allowed sign bits per input: 1 for negative, 2 for zero, 4 for positive.
            variables: All variables in base-to-top order, with fixed coordinates first.
            budget: The shared work and storage allowance.
            lazard: Whether to use Lazard characterizations; False selects McCallum's.

        Returns:
            A witness on the fiber or a covering of it; check it with this CADFiber.

        Raises:
            PolyError: DOMAIN for invalid inputs, UNSUPPORTED for McCallum nullification, LIMIT or CAPACITY on work limits.
        """
        return _coverings(inputs, masks, variables, self._fixed, budget, lazard)

    def check(self, inputs: List[SparsePoly[QQ]], masks: List[UInt8], certificate: CADCoveringsCertificate, mut budget: Budget) raises PolyError -> Bool:
        """Verify a coverings witness or covering on exactly the bound prefix.

        Args:
            inputs: The original rational polynomial family.
            masks: The original allowed-sign masks.
            certificate: The proposed witness or covering.
            budget: The shared verification allowance.

        Returns:
            Whether the certificate proves its verdict on the fiber.

        Raises:
            PolyError: LIMIT or CAPACITY on work limits, DOMAIN for invalid masks.
        """
        _masks(inputs, masks, budget)
        try:
            return _check_coverings(inputs, masks, certificate, self._fixed, budget)
        except e:
            if e.kind == PolyErrorKind.DOMAIN:
                return False
            raise e^

    def check(self, inputs: List[SparsePoly[QQ]], certificate: CADCertificate, mut budget: Budget) raises PolyError -> Bool:
        """Verify the geometric claim on exactly the bound section prefix.

        Args:
            inputs: The original rational polynomial family.
            certificate: The proposed decomposition or feasibility evidence.
            budget: The shared verification allowance.

        Returns:
            Whether local witnesses, prefix binding, and required coverage are valid.

        Raises:
            PolyError: LIMIT or CAPACITY on work limits.
        """
        try:
            return _check_cad(inputs, certificate, self._fixed, budget)
        except e:
            if e.kind == PolyErrorKind.DOMAIN:
                return False
            raise e^

    def check(self, inputs: List[SparsePoly[QQ]], masks: List[UInt8], certificate: CADCoveringCertificate, mut budget: Budget) raises PolyError -> Bool:
        """Verify the geometric claim on exactly the bound section prefix.

        Args:
            inputs: The original rational polynomial family.
            masks: The original allowed-sign masks.
            certificate: The proposed decomposition or feasibility evidence.
            budget: The shared verification allowance.

        Returns:
            Whether local witnesses, prefix binding, and required coverage are valid.

        Raises:
            PolyError: LIMIT or CAPACITY on work limits, DOMAIN for invalid masks.
        """
        _masks(inputs, masks, budget)
        try:
            return _check_cad_covering(inputs, masks, certificate, self._fixed, budget)
        except e:
            if e.kind == PolyErrorKind.DOMAIN:
                return False
            raise e^

    def check(self, inputs: List[SparsePoly[QQ]], masks: List[UInt8], certificate: CylindricalCoveringCertificate, mut budget: Budget) raises PolyError -> Bool:
        """Verify the geometric claim on exactly the bound section prefix.

        Args:
            inputs: The original rational polynomial family.
            masks: The original allowed-sign masks.
            certificate: The proposed decomposition or feasibility evidence.
            budget: The shared verification allowance.

        Returns:
            Whether local witnesses, prefix binding, and required coverage are valid.

        Raises:
            PolyError: LIMIT or CAPACITY on work limits, DOMAIN for invalid masks.
        """
        _masks(inputs, masks, budget)
        try:
            return _check(inputs, masks, certificate, self._fixed, budget)
        except e:
            if e.kind == PolyErrorKind.DOMAIN:
                return False
            raise e^

    def decompose_to[S: CADSink](self, inputs: List[SparsePoly[QQ]], variables: List[Int64], mut sink: S, mut budget: Budget, lazard: Bool = False) raises PolyError -> CADStreamReceipt:
        """Stream complete local stacks and terminal signs for the requested fiber.

        Parameters:
            S: The witness destination type.

        Args:
            inputs: Rational polynomials including the prefix's defining polynomials.
            variables: All variables in base-to-top order, fixed coordinates first.
            sink: The consumer of the projection header and depth-first events.
            budget: The shared computation and transport allowance.
            lazard: Whether to use Lazard projection and lifting.

        Returns:
            Counts of visited stacks and terminal cells, requiring a matching fiber checker.

        Raises:
            PolyError: DOMAIN for invalid inputs or a nonsection prefix, UNSUPPORTED for McCallum nullification, LIMIT or CAPACITY on work limits; failure may leave an incomplete destination.
        """
        return _cad_to(inputs, variables, self._fixed, sink, budget, lazard)

    def checker(self, inputs: List[SparsePoly[QQ]], mut budget: Budget) raises PolyError -> CADStreamChecker:
        """Prepare incremental verification bound to the requested section prefix.

        Args:
            inputs: The original rational polynomial family.
            budget: The shared preparation allowance.

        Returns:
            A checker requiring every branch after the fixed prefix, and no other prefix branch.

        Raises:
            PolyError: DOMAIN for invalid input, LIMIT or CAPACITY on work limits.
        """
        var result = CADStreamChecker.prepare(inputs, budget)
        result._traversal.fixed = _copy_fixed(self._fixed, budget)
        return result^
