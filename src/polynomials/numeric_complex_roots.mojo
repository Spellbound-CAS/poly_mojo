"""Validated complex root covers using exact winding numbers and ball traces."""
from apn_mojo import Rational
from .coefficients import QQ, CoefficientRing
from .dense import DensePoly
from .root_stream import RootStreamReceipt
from .errors import Budget, PolyError, _allocation, _sum, _product
from ._root_ops import _compare
from .numeric_domains import NumericDomain, ComplexBallDomain
from .numeric_dense import NumericDensePoly
from .numeric_bounds import NumericBox, numeric_bounds, _box_point, _box_zero, _box_sub, _box_center, _box_valid, _box_charge
from .numeric_evaluate import EnclosureEvaluation, _evaluate_program, _check_program, _enclosure_context, _exact_dense_boxes, _numeric_dense_boxes
from ._numeric_program import _dense_program
from .numeric_options import NumericOptions, NumericOutcome, _accurate, _next_precision
from ._numeric_roots import _root_degree, _root_bound, _exact_box_evaluate, _translated_boxes, _disjoint, _require_success


@fieldwise_init
struct ComplexBoundarySegment(Copyable, Movable):
    """One oriented rectangle-boundary segment and its zero-free image evidence."""
    var start: NumericBox
    """Exact singleton starting point."""
    var end: NumericBox
    """Exact singleton ending point."""
    var evaluation: EnclosureEvaluation[ComplexBallDomain]
    """Trace for the translated polynomial over the entire segment."""


@fieldwise_init
struct ComplexRootRegion(Copyable, Movable):
    """A rectangle with a uniform root count including multiplicities."""
    var rectangle: NumericBox
    """The rational rectangle, whose boundary is certified root-free."""
    var count: Int
    """The positive interior count, including multiplicity."""
    var boundary: List[ComplexBoundarySegment]
    """A counterclockwise closed boundary cover with image witnesses."""


@fieldwise_init
struct ComplexRootResult(Copyable, Movable):
    """A disjoint complex root cover with explicit completeness and accuracy."""
    var outcome: NumericOutcome
    """SUCCESS only when all requested claims are established."""
    var complete: Bool
    """Checked region counts sum to the uniform polynomial degree."""
    var accurate: Bool
    """Every region meets the requested absolute/relative widths."""
    var isolated: Bool
    """Every region has count one; repeated roots remain clusters."""
    var regions: List[ComplexRootRegion]
    """Disjoint positive-count regions with independently checkable witnesses."""


trait ComplexRootSink:
    """Consume completed complex root witnesses without retaining their prefix."""

    def write(mut self, var region: ComplexRootRegion, mut budget: Budget) raises PolyError:
        """Consume one completed complex root witness.

        Args:
            region: The next owned witness in traversal order.
            budget: The shared transport allowance.

        Raises:
            PolyError: On transport or resource failure.
        """
        ...


@fieldwise_init
struct _ComplexCollector(ComplexRootSink):
    var regions: List[ComplexRootRegion]

    def write(mut self, var region: ComplexRootRegion, mut budget: Budget) raises PolyError:
        _allocation[ComplexRootRegion](1, budget)
        self.regions.append(region^)


def _segment_box(a: NumericBox, b: NumericBox, mut budget: Budget) raises PolyError -> NumericBox:
    _box_charge(a, budget)
    _box_charge(b, budget)
    var forward_x = _compare(a.real_low, b.real_low, budget) <= 0
    var forward_y = _compare(a.imag_low, b.imag_low, budget) <= 0
    return NumericBox(a.real_low.copy() if forward_x else b.real_low.copy(), b.real_low.copy() if forward_x else a.real_low.copy(), a.imag_low.copy() if forward_y else b.imag_low.copy(), b.imag_low.copy() if forward_y else a.imag_low.copy())


def _corners(rectangle: NumericBox, mut budget: Budget) raises PolyError -> List[NumericBox]:
    _allocation[NumericBox](4, budget)
    return [NumericBox.point(rectangle.real_low, rectangle.imag_low, budget), NumericBox.point(rectangle.real_high, rectangle.imag_low, budget), NumericBox.point(rectangle.real_high, rectangle.imag_high, budget), NumericBox.point(rectangle.real_low, rectangle.imag_high, budget)]


def _representative(coefficients: List[NumericBox], mut budget: Budget) raises PolyError -> List[NumericBox]:
    _allocation[NumericBox](len(coefficients), budget)
    var result = List[NumericBox]()
    for box in coefficients:
        result.append(_box_center(box, budget))
    return result^


def _crossing(a: NumericBox, b: NumericBox, mut budget: Budget) raises PolyError -> Int:
    var ay = a.imag_low.value.sign()
    var by = b.imag_low.value.sign()
    if (ay <= 0 and by > 0) or (by <= 0 and ay > 0):
        var cross = a.real_low.mul(b.imag_low, budget).sub(b.real_low.mul(a.imag_low, budget), budget)
        if ay <= 0 and by > 0 and cross.value.sign() > 0:
            return 1
        if by <= 0 and ay > 0 and cross.value.sign() < 0:
            return -1
    return 0


@fieldwise_init
struct _ContourEdge(Movable):
    var start: Int
    var end: Int


def _contour(coefficients: List[NumericBox], representative: List[NumericBox], rectangle: NumericBox, options: NumericOptions, allowance: Int, mut subdivisions: Int, mut budget: Budget) raises PolyError -> Optional[ComplexRootRegion]:
    # The local arena owns endpoints and their exact images. Subdivision shares
    # indices, so each new midpoint is evaluated once and no borrowed span moves.
    var points = _corners(rectangle, budget)
    _allocation[NumericBox](4, budget)
    var images = List[NumericBox](capacity=4)
    for point in points:
        images.append(_exact_box_evaluate(representative, point, budget))
        if _box_zero(images[len(images) - 1], budget):
            return None
    _allocation[_ContourEdge](4, budget)
    var pending = List[_ContourEdge](capacity=4)
    for i in range(3, -1, -1):
        pending.append(_ContourEdge(i, (i + 1) % 4))
    var boundary = List[ComplexBoundarySegment]()
    var used = 0
    var count = 0
    while len(pending):
        budget.charge(1)
        var edge = pending.pop()
        var segment_box = _segment_box(points[edge.start], points[edge.end], budget)
        var center = _box_center(segment_box, budget)
        var translated = _translated_boxes(coefficients, center, budget)
        var segment = _box_sub(segment_box, center, budget)
        var program = _dense_program(translated, segment, budget)
        var bits = options.initial_precision
        var proof = _evaluate_program[ComplexBallDomain](program, _enclosure_context[ComplexBallDomain](bits, options.max_precision), budget)
        var image = numeric_bounds[ComplexBallDomain](proof.value(budget), budget)
        # Avoid futile precision retries when exact interval dependency already
        # includes zero. Geometry, not extra bits, must resolve that case.
        if _box_zero(image, budget) and not _box_zero(_exact_box_evaluate(translated, segment, budget), budget):
            while bits < options.max_precision and _box_zero(image, budget):
                bits = _next_precision(bits, options.max_precision)
                proof = _evaluate_program[ComplexBallDomain](program, _enclosure_context[ComplexBallDomain](bits, options.max_precision), budget)
                image = numeric_bounds[ComplexBallDomain](proof.value(budget), budget)
        if not _box_zero(image, budget):
            count += _crossing(images[edge.start], images[edge.end], budget)
            _allocation[ComplexBoundarySegment](1, budget)
            _box_charge(points[edge.start], budget)
            _box_charge(points[edge.end], budget)
            boundary.append(ComplexBoundarySegment(points[edge.start].copy(), points[edge.end].copy(), proof^))
            continue
        if used == allowance or subdivisions == options.max_subdivisions:
            return None
        used += 1
        subdivisions += 1
        var value = _exact_box_evaluate(representative, center, budget)
        if _box_zero(value, budget):
            return None
        var middle = len(points)
        _allocation[NumericBox](2, budget)
        points.append(center^)
        images.append(value^)
        _allocation[_ContourEdge](2, budget)
        pending.append(_ContourEdge(middle, edge.end))
        pending.append(_ContourEdge(edge.start, middle))
    _box_charge(rectangle, budget)
    return ComplexRootRegion(rectangle.copy(), count, boundary^)


def _oriented_edge(a: NumericBox, b: NumericBox, rectangle: NumericBox, mut budget: Budget) raises PolyError -> Int:
    if not _box_point(a, budget) or not _box_point(b, budget) or not rectangle.contains(a, budget) or not rectangle.contains(b, budget):
        return -1
    if a.imag_low == rectangle.imag_low and b.imag_low == rectangle.imag_low and _compare(a.real_low, b.real_low, budget) < 0:
        return 0
    if a.real_low == rectangle.real_high and b.real_low == rectangle.real_high and _compare(a.imag_low, b.imag_low, budget) < 0:
        return 1
    if a.imag_low == rectangle.imag_high and b.imag_low == rectangle.imag_high and _compare(a.real_low, b.real_low, budget) > 0:
        return 2
    if a.real_low == rectangle.real_low and b.real_low == rectangle.real_low and _compare(a.imag_low, b.imag_low, budget) > 0:
        return 3
    return -1


def _check_region(coefficients: List[NumericBox], representative: List[NumericBox], region: ComplexRootRegion, mut budget: Budget) raises PolyError -> Bool:
    if not _box_valid(region.rectangle, budget) or _compare(region.rectangle.real_low, region.rectangle.real_high, budget) >= 0 or _compare(region.rectangle.imag_low, region.rectangle.imag_high, budget) >= 0 or region.count < 0 or len(region.boundary) < 4:
        return False
    var first = NumericBox.point(region.rectangle.real_low, region.rectangle.imag_low, budget)
    var previous = first.copy()
    var edge = 0
    var count = 0
    for segment in region.boundary:
        if segment.start != previous:
            return False
        var next_edge = _oriented_edge(segment.start, segment.end, region.rectangle, budget)
        if next_edge < edge or next_edge > edge + 1:
            return False
        edge = next_edge
        var segment_box = _segment_box(segment.start, segment.end, budget)
        var center = _box_center(segment_box, budget)
        var translated = _translated_boxes(coefficients, center, budget)
        var argument = _box_sub(segment_box, center, budget)
        if not _check_program[ComplexBallDomain](_dense_program(translated, argument, budget), segment.evaluation, budget):
            return False
        var enclosure = numeric_bounds[ComplexBallDomain](segment.evaluation.value(budget), budget)
        if _box_zero(enclosure, budget):
            return False
        var a = _exact_box_evaluate(representative, segment.start, budget)
        var b = _exact_box_evaluate(representative, segment.end, budget)
        if not enclosure.contains(a, budget) or not enclosure.contains(b, budget):
            return False
        count += _crossing(a, b, budget)
        previous = segment.end.copy()
    return edge == 3 and previous == first and count == region.count


def _complex_roots(coefficients: List[NumericBox], options: NumericOptions, mut budget: Budget, complete_only: Bool) raises PolyError -> ComplexRootResult:
    options._validate()
    var degree = _root_degree(coefficients, budget)
    if not degree:
        _require_success(NumericOutcome.UNCERTAIN_DEGREE, complete_only)
        return ComplexRootResult(NumericOutcome.UNCERTAIN_DEGREE, False, False, False, List[ComplexRootRegion]())
    if degree.value() == 0:
        return ComplexRootResult(NumericOutcome.SUCCESS, True, True, True, List[ComplexRootRegion]())
    var bound = _root_bound(coefficients, degree.value(), budget)
    var rectangle = NumericBox(bound.neg(budget), bound.copy(), bound.neg(budget), bound.copy())
    var representative = _representative(coefficients, budget)
    var subdivisions = 0
    var first = _contour(coefficients, representative, rectangle, options, options.max_subdivisions, subdivisions, budget)
    if not first or first.value().count != degree.value():
        _require_success(NumericOutcome.UNRESOLVED_BOUNDARY, complete_only)
        return ComplexRootResult(NumericOutcome.UNRESOLVED_BOUNDARY, False, False, False, List[ComplexRootRegion]())
    _allocation[ComplexRootRegion](1, budget)
    var pending: List[ComplexRootRegion] = [first.take()]
    return _refine_complex_regions(coefficients, representative, pending^, degree.value(), options, subdivisions, budget, complete_only)


def _refine_complex_regions(coefficients: List[NumericBox], representative: List[NumericBox], var pending: List[ComplexRootRegion], degree: Int, options: NumericOptions, mut subdivisions: Int, mut budget: Budget, complete_only: Bool) raises PolyError -> ComplexRootResult:
    var sink = _ComplexCollector(List[ComplexRootRegion]())
    var receipt = _refine_complex_regions_to(coefficients, representative, pending^, degree, options, subdivisions, budget, complete_only, sink)
    var completed = List[ComplexRootRegion]()
    swap(completed, sink.regions)
    return ComplexRootResult(receipt.outcome, receipt.complete, receipt.accurate, receipt.isolated, completed^)


def _refine_complex_regions_to[S: ComplexRootSink](coefficients: List[NumericBox], representative: List[NumericBox], var pending: List[ComplexRootRegion], degree: Int, options: NumericOptions, mut subdivisions: Int, mut budget: Budget, complete_only: Bool, mut sink: S) raises PolyError -> RootStreamReceipt:
    var items = 0
    var accurate = True
    var isolated = True
    var reason = NumericOutcome.SUCCESS
    while len(pending):
        budget.charge(1)
        var region = pending.pop()
        var enough = _accurate(region.rectangle, options, budget) and (not options.require_isolation or region.count == 1)
        if enough or subdivisions == options.max_subdivisions:
            accurate = accurate and _accurate(region.rectangle, options, budget)
            isolated = isolated and region.count == 1
            if not enough:
                reason = NumericOutcome.SUBDIVISION_LIMIT
            items = _sum(items, 1)
            sink.write(region^, budget)
            continue
        var split_ok = False
        var width = region.rectangle.real_high.sub(region.rectangle.real_low, budget)
        var height = region.rectangle.imag_high.sub(region.rectangle.imag_low, budget)
        var vertical = _compare(width, height, budget) >= 0
        # Three asymmetric rational cuts avoid common roots on coordinate axes.
        for numerator in range(2, 5):
            if subdivisions == options.max_subdivisions:
                break
            subdivisions += 1
            var fraction = QQ(Rational(numerator)).mul(QQ(Rational(7)).inverse(budget), budget)
            var left = region.rectangle.copy()
            var right = region.rectangle.copy()
            if vertical:
                var split = left.real_low.add(width.mul(fraction, budget), budget)
                left.real_high = split.copy()
                right.real_low = split^
            else:
                var split = left.imag_low.add(height.mul(fraction, budget), budget)
                left.imag_high = split.copy()
                right.imag_low = split^
            # Bounded failed contours leave work for alternative split lines.
            var allowance = min(options.max_subdivisions - subdivisions, _sum(32, _product(degree, 8)))
            var a = _contour(coefficients, representative, left, options, allowance, subdivisions, budget)
            if not a:
                continue
            allowance = min(options.max_subdivisions - subdivisions, _sum(32, _product(degree, 8)))
            var b = _contour(coefficients, representative, right, options, allowance, subdivisions, budget)
            if not b or a.value().count < 0 or b.value().count < 0 or a.value().count + b.value().count != region.count:
                continue
            _allocation[ComplexRootRegion](2, budget)
            if b.value().count:
                pending.append(b.take())
            if a.value().count:
                pending.append(a.take())
            split_ok = True
            break
        if not split_ok:
            accurate = accurate and _accurate(region.rectangle, options, budget)
            isolated = isolated and region.count == 1
            reason = NumericOutcome.SUBDIVISION_LIMIT if subdivisions == options.max_subdivisions else NumericOutcome.UNRESOLVED_BOUNDARY
            items = _sum(items, 1)
            sink.write(region^, budget)
    _require_success(reason, complete_only)
    return RootStreamReceipt(reason, True, accurate, isolated, items)


def _check_complex(coefficients: List[NumericBox], result: ComplexRootResult, options: NumericOptions, mut budget: Budget) raises PolyError -> Bool:
    options._validate()
    var degree = _root_degree(coefficients, budget)
    if not degree or not result.complete:
        return False
    var bound = _root_bound(coefficients, degree.value(), budget)
    var global_box = NumericBox(bound.neg(budget), bound.copy(), bound.neg(budget), bound.copy())
    var representative = _representative(coefficients, budget)
    var total = 0
    var accurate = True
    var isolated = True
    for i in range(len(result.regions)):
        ref region = result.regions[i]
        if region.count <= 0 or not global_box.contains(region.rectangle, budget) or not _check_region(coefficients, representative, region, budget):
            return False
        for j in range(i):
            if not _disjoint(region.rectangle, result.regions[j].rectangle, budget, True):
                return False
        total = _sum(total, region.count)
        accurate = accurate and _accurate(region.rectangle, options, budget)
        isolated = isolated and region.count == 1
    var success = accurate and (isolated or not options.require_isolation)
    return total == degree.value() and result.accurate == accurate and result.isolated == isolated and ((result.outcome == NumericOutcome.SUCCESS and success) or ((result.outcome == NumericOutcome.SUBDIVISION_LIMIT or result.outcome == NumericOutcome.UNRESOLVED_BOUNDARY) and not success))


def complex_roots_exact[C: CoefficientRing](poly: DensePoly[C], options: NumericOptions, mut budget: Budget, *, complete_only: Bool = False) raises PolyError -> ComplexRootResult:
    """Produce a validated complex cover of an exact ZZ/QQ/QQi polynomial's roots.

    Parameters:
        C: ZZ, QQ, or QQi.

    Args:
        poly: The original exact polynomial.
        options: Precision, accuracy, isolation and subdivision settings.
        budget: The shared work allowance.
        complete_only: Raise unless every requested guarantee succeeds.

    Returns:
        Root rectangles with multiplicity counts and winding witnesses.

    Raises:
        PolyError: DOMAIN for zero or invalid input, UNSUPPORTED for other domains, LIMIT or CAPACITY on work limits.
    """
    return _complex_roots(_exact_dense_boxes(poly, budget), options, budget, complete_only)


def complex_roots_numeric[D: NumericDomain](poly: NumericDensePoly[D], options: NumericOptions, mut budget: Budget, *, complete_only: Bool = False) raises PolyError -> ComplexRootResult:
    """Validate complex roots uniformly over a numerical coefficient family.

    Parameters:
        D: Any of the four numerical coefficient domains.

    Args:
        poly: The point or family polynomial.
        options: Precision, accuracy, isolation and subdivision settings.
        budget: The shared work allowance.
        complete_only: Raise unless every requested guarantee succeeds.

    Returns:
        A certified complete cover or an explicit inconclusive outcome.

    Raises:
        PolyError: DOMAIN for zero or invalid input, LIMIT or CAPACITY on work limits.
    """
    return _complex_roots(_numeric_dense_boxes(poly, budget), options, budget, complete_only)


def check_complex_roots_exact[C: CoefficientRing](poly: DensePoly[C], result: ComplexRootResult, options: NumericOptions, mut budget: Budget) raises PolyError -> Bool:
    """Check complex winding evidence and completeness without root search.

    Parameters:
        C: ZZ, QQ, or QQi.

    Args:
        poly: The original exact polynomial.
        result: The proposed complete root cover.
        options: The original requested guarantees.
        budget: The shared work allowance.

    Returns:
        Whether all counts, boundaries, disjointness and outcome flags hold.

    Raises:
        PolyError: DOMAIN for invalid input, LIMIT or CAPACITY on work limits.
    """
    return _check_complex(_exact_dense_boxes(poly, budget), result, options, budget)


def check_complex_roots_numeric[D: NumericDomain](poly: NumericDensePoly[D], result: ComplexRootResult, options: NumericOptions, mut budget: Budget) raises PolyError -> Bool:
    """Verify a uniform complex root cover for the original coefficient family.

    Parameters:
        D: The original numerical coefficient domain.

    Args:
        poly: The original point or family polynomial.
        result: The proposed complete root cover.
        options: The original requested guarantees.
        budget: The shared work allowance.

    Returns:
        Whether every root-count and global completeness obligation is valid.

    Raises:
        PolyError: DOMAIN for invalid input, LIMIT or CAPACITY on work limits.
    """
    return _check_complex(_numeric_dense_boxes(poly, budget), result, options, budget)


def _complex_roots_to[S: ComplexRootSink](coefficients: List[NumericBox], options: NumericOptions, mut budget: Budget, complete_only: Bool, mut sink: S) raises PolyError -> RootStreamReceipt:
    options._validate()
    var degree = _root_degree(coefficients, budget)
    if not degree:
        _require_success(NumericOutcome.UNCERTAIN_DEGREE, complete_only)
        return RootStreamReceipt(NumericOutcome.UNCERTAIN_DEGREE, False, False, False, 0)
    if degree.value() == 0:
        return RootStreamReceipt(NumericOutcome.SUCCESS, True, True, True, 0)
    var bound = _root_bound(coefficients, degree.value(), budget)
    var rectangle = NumericBox(bound.neg(budget), bound.copy(), bound.neg(budget), bound.copy())
    var representative = _representative(coefficients, budget)
    var subdivisions = 0
    var first = _contour(coefficients, representative, rectangle, options, options.max_subdivisions, subdivisions, budget)
    if not first or first.value().count != degree.value():
        _require_success(NumericOutcome.UNRESOLVED_BOUNDARY, complete_only)
        return RootStreamReceipt(NumericOutcome.UNRESOLVED_BOUNDARY, False, False, False, 0)
    _allocation[ComplexRootRegion](1, budget)
    var pending: List[ComplexRootRegion] = [first.take()]
    return _refine_complex_regions_to(coefficients, representative, pending^, degree.value(), options, subdivisions, budget, complete_only, sink)


def complex_roots_exact_to[C: CoefficientRing, S: ComplexRootSink](poly: DensePoly[C], options: NumericOptions, mut sink: S, mut budget: Budget, *, complete_only: Bool = False) raises PolyError -> RootStreamReceipt:
    """Produce complex root witnesses directly into an incremental consumer.

    Parameters:
        C: The source coefficient domain.
        S: The witness destination type.

    Args:
        poly: The original polynomial or coefficient family.
        options: The explicit accuracy and subdivision settings.
        sink: The consumer receiving completed witnesses.
        budget: The shared computation and transport allowance.
        complete_only: Raise unless all requested guarantees hold.

    Returns:
        A completion receipt; the receipt alone is not a root certificate.

    Raises:
        PolyError: DOMAIN for invalid input, UNSUPPORTED for unsupported domains, LIMIT or CAPACITY on work limits; failure can leave an incomplete destination.
    """
    return _complex_roots_to(_exact_dense_boxes(poly, budget), options, budget, complete_only, sink)


def complex_roots_numeric_to[D: NumericDomain, S: ComplexRootSink](poly: NumericDensePoly[D], options: NumericOptions, mut sink: S, mut budget: Budget, *, complete_only: Bool = False) raises PolyError -> RootStreamReceipt:
    """Produce complex root witnesses directly into an incremental consumer.

    Parameters:
        D: The source coefficient domain.
        S: The witness destination type.

    Args:
        poly: The original polynomial or coefficient family.
        options: The explicit accuracy and subdivision settings.
        sink: The consumer receiving completed witnesses.
        budget: The shared computation and transport allowance.
        complete_only: Raise unless all requested guarantees hold.

    Returns:
        A completion receipt; the receipt alone is not a root certificate.

    Raises:
        PolyError: DOMAIN for invalid input, UNSUPPORTED for unsupported domains, LIMIT or CAPACITY on work limits; failure can leave an incomplete destination.
    """
    return _complex_roots_to(_numeric_dense_boxes(poly, budget), options, budget, complete_only, sink)
