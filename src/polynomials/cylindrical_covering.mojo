"""Demand-driven cylindrical coverings with independently checked overlapping conflicts.

Source-guided by SMT-RAT fc28c1d, coveringng Algorithm.h and Sampling.h.
The global projection closure supplies delineability; lifting and conflict
search visit only uncovered intervals. No complete CAD is constructed.
"""
# Copyright (c) 2015-2022 Florian Corzilius, Gereon Kremer, Jasper Nalbach,
# Sebastian Junges, Stefan Schupp, Erika Abraham.
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
# THE SOFTWARE.
from .errors import Budget, PolyError, PolyErrorKind, _allocation, _product, _sum
from .coefficients import QQ, NoContext
from .sparse import SparsePoly
from .certificates import _equal
from .algebraic_types import RealAlgebraic, AlgebraicComparisonCertificate, _rational, _copy_real, _same_root, _check_comparison
from .algebraic import algebraic_sign
from .algebraic_specialization import AlgebraicPolynomialEvaluation, _evaluate_sparse_cached, check_algebraic_evaluate_sparse
from .cad_projection import CADProjectionCertificate, cad_projection, check_cad_projection, _level
from .cad_types import CADCell, CADStackCertificate
from ._cad_common import _cell, _coordinates, _signs, _signs_cached, _check_signs, _masks, _violation
from ._cad_lifting import _stack_cached, _check_stack, _position, _sample_lifting
from ._algebraic import _NormalizationCache, _Frame as _PointFrame
from ._algebraic_filter import _child_frame, _installed
from ._cad_fiber_common import _validate_fixed, _selection


@fieldwise_init
struct CADConflict(Copyable, Movable, Writable):
    """A certified violated constraint at an exact partial assignment."""
    var constraint: Int
    """The index of the violated input polynomial."""
    var evaluation: AlgebraicPolynomialEvaluation
    """Exact polynomial evaluation at the interval sample or root assignment."""
    var sign: AlgebraicComparisonCertificate
    """The value's sign, excluded by the constraint's allowed-sign mask."""


@fieldwise_init
struct CADConflictInterval(Copyable, Movable, Writable):
    """A conflicting interval spanning consecutive sections and sectors of a stack."""
    var parent: Int
    """The base cell whose stack contains the interval."""
    var lower: Int
    """The first included stack position, with zero denoting the lower unbounded sector."""
    var upper: Int
    """The last included stack position; section parity determines endpoint inclusion."""
    var position: Int
    """The stack sample at which the conflict is certified."""
    var proof: CADConflict
    """The violated constraint and its evaluation and sign proofs."""


@fieldwise_init
struct CylindricalCoveringCertificate(Copyable, Movable, Writable):
    """A selective feasibility witness or overlapping cylindrical conflict covering."""
    var projection: CADProjectionCertificate
    """The projection closure proving delineability and sign invariance."""
    var cells: List[CADCell]
    """Only visited base cells, with parents preceding children and root at index zero."""
    var stacks: List[CADStackCertificate]
    """Complete local root delineations for visited nonterminal cells."""
    var conflicts: List[CADConflictInterval]
    """Possibly overlapping intervals excluded by individual constraints."""
    var root_conflict: Optional[CADConflict]
    """An optional constant contradiction excluding all coordinates without lifting."""
    var witness_cell: Int
    """A fully assigned satisfying cell, or minus one for a complete conflict covering."""


@fieldwise_init
struct _Frame(Movable):
    var parent: Int
    var stack: Int
    var covered: List[Bool]


def _conflict(inputs: List[SparsePoly[QQ]], masks: List[UInt8], levels: List[Int], level: Int, point: List[Tuple[Int64, RealAlgebraic]], mut cache: _NormalizationCache, mut budget: Budget) raises PolyError -> Optional[CADConflict]:
    for i in range(len(inputs)):
        budget.charge(1)
        if levels[i] != level or masks[i] == 7:
            continue
        var evaluation = _evaluate_sparse_cached(inputs[i], point, cache, budget)
        var sign = algebraic_sign(evaluation.value, budget)
        if masks[i] & (UInt8(1) << UInt8(sign.relation + 1)) == 0:
            return CADConflict(i, evaluation^, sign^)
    return None


def _check_conflict(inputs: List[SparsePoly[QQ]], masks: List[UInt8], variables: List[Int64], point: List[Tuple[Int64, RealAlgebraic]], proof: CADConflict, mut budget: Budget) raises PolyError -> Bool:
    var i = proof.constraint
    if i < 0 or i >= len(inputs) or _level(inputs[i], variables, budget) >= len(point):
        return False
    var zero = _rational(QQ.zero(NoContext()), budget)
    return check_algebraic_evaluate_sparse(inputs[i], point, proof.evaluation, budget) and _check_comparison(proof.evaluation.value, zero, proof.sign, budget) and masks[i] & (UInt8(1) << UInt8(proof.sign.relation + 1)) == 0


def _range(projection: CADProjectionCertificate, stack: CADStackCertificate, constraint: Int, position: Int, mask: UInt8, mut budget: Budget) raises PolyError -> Tuple[Int, Int]:
    var lower = 0
    var upper = len(stack.samples) - 1
    if mask == 0:
        return (lower, upper)
    # Input factorizations occupy the first entries in the projection schedule.
    # Between consecutive roots of these factors the input sign cannot change.
    # Other polynomials' boundaries can therefore be crossed by this conflict.
    for factor in projection.factorizations[constraint].factors:
        budget.charge(1)
        if _level(factor[0], projection.variables, budget) != stack.level:
            continue
        var index = -1
        for j in range(len(projection.levels[stack.level])):
            budget.charge(1)
            if _equal(factor[0], projection.levels[stack.level][j], budget):
                index = j
                break
        if index < 0:
            raise PolyError(PolyErrorKind.DOMAIN, "missing conflict projection factor")
        for root in stack.mapping[index]:
            budget.charge(1)
            var boundary = _sum(_product(2, root), 1)
            if boundary == position:
                return (position, position)
            if boundary < position:
                lower = max(lower, boundary + 1)
            else:
                upper = min(upper, boundary - 1)
    return (lower, upper)


def _mark(mut covered: List[Bool], lower: Int, upper: Int, mut budget: Budget) raises PolyError:
    for i in range(lower, upper + 1):
        budget.charge(1)
        covered[i] = True


def _uncovered(stack: CADStackCertificate, covered: List[Bool], mut budget: Budget) raises PolyError -> Int:
    var zero = _position(stack.roots, _rational(QQ.zero(NoContext()), budget), budget)
    if not covered[zero]:
        return zero
    if not covered[0]:
        return 0
    if not covered[len(covered) - 1]:
        return len(covered) - 1
    # Prefer rational sector samples before algebraic boundary samples.
    for parity in range(2):
        for i in range(parity, len(covered), 2):
            budget.charge(1)
            if not covered[i]:
                return i
    return -1


def _push(mut result: CylindricalCoveringCertificate, parent: Int, mut frames: List[_Frame], fixed: List[RealAlgebraic], point_frame: _PointFrame, mut liftings: List[List[Optional[SparsePoly[QQ]]]], mut budget: Budget) raises PolyError:
    var point = _coordinates(result.cells, parent, result.projection.variables, budget)
    var level = result.cells[parent].level + 1
    var cache = _NormalizationCache()
    var framed = point_frame.valid and Bool(point_frame.generator)
    cache.frame = _installed(point_frame, point[0], budget)
    var stack = _stack_cached(result.projection, parent, level, point[0], point[1], _rational(QQ.zero(NoContext()), budget), cache, budget)
    # The lifting polynomial of each root, for the frames of section children.
    var roots = List[Optional[SparsePoly[QQ]]]()
    if framed:
        for gap in range(len(stack.roots)):
            roots.append(_sample_lifting(result.projection, stack, point[0], _sum(_product(2, gap), 1), budget))
    liftings.append(roots^)
    _allocation[Bool](len(stack.samples), budget)
    _allocation[_Frame](1, budget)
    _allocation[CADStackCertificate](1, budget)
    var selected = _selection(stack, fixed, budget)
    var covered = List[Bool](length=len(stack.samples), fill=True)
    for i in range(selected[0], selected[1]):
        budget.charge(1)
        covered[i] = False
    frames.append(_Frame(parent, len(result.stacks), covered^))
    result.stacks.append(stack^)


def cylindrical_covering(inputs: List[SparsePoly[QQ]], masks: List[UInt8], variables: List[Int64], mut budget: Budget, lazard: Bool = False) raises PolyError -> CylindricalCoveringCertificate:
    """Decide a sign conjunction by sampling uncovered intervals and pruning cylinders.

    Args:
        inputs: The rational constraint polynomials.
        masks: Allowed sign bits per polynomial: one for negative, two for zero, four for positive.
        variables: Distinct variables in base-to-top cylindrical order.
        budget: The shared work allowance.
        lazard: Whether to use Lazard projection and vanishing-order evaluation.

    Returns:
        A satisfying assignment or independently checkable overlapping conflict coverage.

    Raises:
        PolyError: DOMAIN for invalid inputs, UNSUPPORTED for McCallum nullification, LIMIT or CAPACITY on work limits.
    """
    return _cylindrical_covering(inputs, masks, variables, List[RealAlgebraic](), budget, lazard)


def _cylindrical_covering(inputs: List[SparsePoly[QQ]], masks: List[UInt8], variables: List[Int64], fixed: List[RealAlgebraic], mut budget: Budget, lazard: Bool) raises PolyError -> CylindricalCoveringCertificate:
    _masks(inputs, masks, budget)
    _validate_fixed(fixed, len(variables), budget)
    var projection = cad_projection(inputs, variables, budget, lazard=lazard)
    _allocation[CADCell](1, budget)
    var result = CylindricalCoveringCertificate(projection^, [_cell(-1, -1, 0, _rational(QQ.zero(NoContext()), budget), budget)], List[CADStackCertificate](), List[CADConflictInterval](), None, -1)
    _allocation[Int](len(inputs), budget)
    var levels = List[Int]()
    for f in inputs:
        budget.charge(1)
        levels.append(_level(f, variables, budget))
    var empty = List[Tuple[Int64, RealAlgebraic]]()
    var plain = _NormalizationCache()
    result.root_conflict = _conflict(inputs, masks, levels, -1, empty, plain, budget)
    if result.root_conflict:
        return result^
    if len(variables) == 0:
        _signs(inputs, empty, result.cells[0].evaluations, result.cells[0].signs, budget)
        result.witness_cell = 0
        return result^
    var frames = List[_Frame]()
    # Each cell's point frame, and each stack's root lifting polynomials.
    var point_frames: List[_PointFrame] = [_PointFrame()]
    var liftings = List[List[Optional[SparsePoly[QQ]]]]()
    _push(result, 0, frames, fixed, point_frames[0], liftings, budget)
    while len(frames):
        budget.charge(1)
        var frame = len(frames) - 1
        var parent = frames[frame].parent
        var stack_index = frames[frame].stack
        var position = _uncovered(result.stacks[stack_index], frames[frame].covered, budget)
        if position == -1:
            _ = frames.pop()
            if len(frames):
                frames[len(frames) - 1].covered[result.cells[parent].position] = True
            continue
        var point = _coordinates(result.cells, parent, variables, budget)
        var level = result.cells[parent].level + 1
        _allocation[Tuple[Int64, RealAlgebraic]](1, budget)
        point[0].append((variables[level], _copy_real(result.stacks[stack_index].samples[position], budget)))
        var lifting: Optional[SparsePoly[QQ]] = None
        if position % 2 == 1 and position // 2 < len(liftings[stack_index]) and liftings[stack_index][position // 2]:
            lifting = liftings[stack_index][position // 2].value()._copy(budget)
        var child_frame = _child_frame(point_frames[parent], variables[level], lifting, result.stacks[stack_index].samples[position], budget)
        var cache = _NormalizationCache()
        cache.frame = _installed(child_frame, point[0], budget)
        var conflict = _conflict(inputs, masks, levels, level, point[0], cache, budget)
        if conflict:
            var proof = conflict.take()
            var bounds = _range(result.projection, result.stacks[stack_index], proof.constraint, position, masks[proof.constraint], budget)
            if level < len(fixed):
                bounds = (position, position)
            _allocation[CADConflictInterval](1, budget)
            result.conflicts.append(CADConflictInterval(parent, bounds[0], bounds[1], position, proof^))
            _mark(frames[frame].covered, bounds[0], bounds[1], budget)
            continue
        var child = _cell(parent, level, position, result.stacks[stack_index].samples[position], budget)
        if level == len(variables) - 1:
            _signs_cached(inputs, point[0], child.evaluations, child.signs, cache, budget)
            result.witness_cell = len(result.cells)
        _allocation[CADCell](1, budget)
        result.cells.append(child^)
        point_frames.append(child_frame^)
        if result.witness_cell != -1:
            return result^
        _push(result, len(result.cells) - 1, frames, fixed, point_frames[len(point_frames) - 1], liftings, budget)
    return result^


def _check(inputs: List[SparsePoly[QQ]], masks: List[UInt8], certificate: CylindricalCoveringCertificate, fixed: List[RealAlgebraic], mut budget: Budget) raises PolyError -> Bool:
    if not check_cad_projection(inputs, certificate.projection, budget) or len(certificate.cells) == 0:
        return False
    ref variables = certificate.projection.variables
    _validate_fixed(fixed, len(variables), budget)
    ref cells = certificate.cells
    var size = len(cells)
    ref root = cells[0]
    if root.parent != -1 or root.level != -1 or root.position != 0 or not _same_root(root.sample, _rational(QQ.zero(NoContext()), budget), budget):
        return False
    if certificate.root_conflict:
        return certificate.witness_cell == -1 and size == 1 and len(certificate.stacks) == 0 and len(certificate.conflicts) == 0 and len(root.evaluations) == 0 and len(root.signs) == 0 and _check_conflict(inputs, masks, variables, List[Tuple[Int64, RealAlgebraic]](), certificate.root_conflict.value(), budget)
    var witness = certificate.witness_cell
    if witness < -1 or witness >= size:
        return False
    _allocation[Int](size, budget)
    _allocation[List[Bool]](_product(2, size), budget)
    var stack_indices = List[Int](length=size, fill=-1)
    var covered = List[List[Bool]](length=size, fill=List[Bool]())
    var occupied = List[List[Bool]](length=size, fill=List[Bool]())
    for i in range(len(certificate.stacks)):
        budget.charge(1)
        ref stack = certificate.stacks[i]
        var parent = stack.parent
        if parent < 0 or parent >= size or stack_indices[parent] != -1 or cells[parent].level >= len(variables) - 1:
            return False
        var point = _coordinates(cells, parent, variables, budget)
        if not _check_stack(certificate.projection, parent, cells[parent].level + 1, point[0], point[1], stack, budget):
            return False
        stack_indices[parent] = i
        _allocation[Bool](_product(2, len(stack.samples)), budget)
        var selected = _selection(stack, fixed, budget)
        covered[parent] = List[Bool](length=len(stack.samples), fill=True)
        for j in range(selected[0], selected[1]):
            budget.charge(1)
            covered[parent][j] = False
        occupied[parent] = List[Bool](length=len(stack.samples), fill=False)
    for i in range(size):
        budget.charge(1)
        ref cell = cells[i]
        if i:
            if cell.parent < 0 or cell.parent >= i or cell.level != cells[cell.parent].level + 1 or cell.level >= len(variables):
                return False
            var stack_index = stack_indices[cell.parent]
            if stack_index < 0:
                return False
            ref stack = certificate.stacks[stack_index]
            if cell.position < 0 or cell.position >= len(stack.samples) or occupied[cell.parent][cell.position] or not _same_root(cell.sample, stack.samples[cell.position], budget):
                return False
            if cell.level < len(fixed):
                # The child refines its stack sample; bind that sample to the
                # requested value with the same checked selection as traversal.
                var selected = _selection(stack, fixed, budget)
                if cell.position != selected[0]:
                    return False
            occupied[cell.parent][cell.position] = True
        if i == witness:
            if cell.level != len(variables) - 1:
                return False
            var point = _coordinates(cells, i, variables, budget)
            if not _check_signs(inputs, point[0], cell.evaluations, cell.signs, budget) or _violation(cell, masks, budget) != -1:
                return False
        elif len(cell.evaluations) or len(cell.signs) or stack_indices[i] < 0:
            return False
    for conflict in certificate.conflicts:
        budget.charge(1)
        var parent = conflict.parent
        if parent < 0 or parent >= size or stack_indices[parent] < 0:
            return False
        ref stack = certificate.stacks[stack_indices[parent]]
        if conflict.lower < 0 or conflict.upper >= len(stack.samples) or conflict.position < conflict.lower or conflict.position > conflict.upper:
            return False
        if stack.level < len(fixed):
            var selected = _selection(stack, fixed, budget)
            if conflict.lower != selected[0] or conflict.upper != selected[0]:
                return False
        var point = _coordinates(cells, parent, variables, budget)
        _allocation[Tuple[Int64, RealAlgebraic]](1, budget)
        point[0].append((variables[stack.level], _copy_real(stack.samples[conflict.position], budget)))
        if not _check_conflict(inputs, masks, variables, point[0], conflict.proof, budget):
            return False
        var bounds = _range(certificate.projection, stack, conflict.proof.constraint, conflict.position, masks[conflict.proof.constraint], budget)
        if conflict.lower < bounds[0] or conflict.upper > bounds[1]:
            return False
        _mark(covered[parent], conflict.lower, conflict.upper, budget)
    if witness != -1:
        return True
    # Children precede their parents in this traversal. Each completed child
    # excludes its entire base cell; intervals and children must cover the line.
    for i in range(size - 1, -1, -1):
        budget.charge(1)
        if stack_indices[i] < 0:
            return False
        for value in covered[i]:
            budget.charge(1)
            if not value:
                return False
        if i:
            covered[cells[i].parent][cells[i].position] = True
    return True


def check_cylindrical_covering(inputs: List[SparsePoly[QQ]], masks: List[UInt8], certificate: CylindricalCoveringCertificate, mut budget: Budget) raises PolyError -> Bool:
    """Verify selective conflict intervals and cylindrical coverage without searching.

    Args:
        inputs: The original polynomial constraints.
        masks: The original allowed-sign masks.
        certificate: The proposed selective feasibility or infeasibility certificate.
        budget: The shared work allowance.

    Returns:
        Whether exact evaluation, projection, interval expansion and coverage justify the result.

    Raises:
        PolyError: DOMAIN for invalid masks, LIMIT or CAPACITY on work limits.
    """
    _masks(inputs, masks, budget)
    try:
        return _check(inputs, masks, certificate, List[RealAlgebraic](), budget)
    except e:
        if e.kind == PolyErrorKind.DOMAIN:
            return False
        raise e^
