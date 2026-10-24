"""Absolute geometry requested after acceptance of a compact CAD.

This adapter creates coordinate values, not intermediate lifting proofs. Its
cost belongs to the caller that requests absolute cells; construction and
final verification do not call it.
"""
from .coefficients import QQ, NoContext
from .algebraic_types import RealAlgebraic, _rational
from .errors import Budget, PolyError, PolyErrorKind, _allocation, _sum, _product
from ._cad_tree_final import _FinalCAD, _terminal_signs
from ._cad_field_extension import _absolute_section
from .cad_final_types import CADFinalCell as _MaterializedFinalCell


def _final_cell_count(result: _FinalCAD, mut budget: Budget) raises PolyError -> Int:
    var count = 1
    for node in result.stacks:
        budget.charge(1)
        count = _sum(count, _sum(_product(2, len(node.stack.roots)), 1))
    return count


def _materialize_final_cells(result: _FinalCAD, mut budget: Budget) raises PolyError -> List[_MaterializedFinalCell]:
    """Materialize every accepted tree coordinate once in breadth-first order."""
    var count = _final_cell_count(result, budget)
    _allocation[_MaterializedFinalCell](count, budget)
    _allocation[Int](len(result.stacks), budget)
    var cells = List[_MaterializedFinalCell](capacity=count)
    cells.append(_MaterializedFinalCell(-1, -1, 0, _rational(QQ.zero(NoContext()), budget), List[Int]()))
    var starts = List[Int](capacity=len(result.stacks))
    for index in range(len(result.stacks)):
        budget.charge(1)
        ref node = result.stacks[index]
        var parent = 0
        if index:
            if node.parent < 0 or node.parent >= index:
                raise PolyError(PolyErrorKind.DOMAIN, "invalid compact CAD parent")
            parent = _sum(starts[node.parent], node.slot)
        starts.append(len(cells))
        # Absolute root selection can refine the embedding. Keep those
        # refinements local to this output operation.
        var field = node.field.copy()
        var size = _sum(_product(2, len(node.stack.roots)), 1)
        for slot in range(size):
            budget.charge(1)
            var sample: RealAlgebraic
            if slot % 2:
                sample = _absolute_section(node.stack, slot // 2, field, budget)
            else:
                sample = _rational(node.stack.sectors[slot // 2], budget)
            var signs = List[Int]()
            if node.level == len(result.projection.variables) - 1:
                signs = _terminal_signs(node, len(result.inputs), slot, budget)
            cells.append(_MaterializedFinalCell(parent, node.level, slot, sample^, signs^))
    return cells^


def _final_witness_point(result: _FinalCAD, mut budget: Budget) raises PolyError -> Optional[List[RealAlgebraic]]:
    """Materialize only the first satisfying path, in cylindrical order."""
    for index in range(len(result.stacks)):
        ref terminal = result.stacks[index]
        if terminal.level != len(result.projection.variables) - 1:
            continue
        for position in range(len(terminal.truth)):
            budget.charge(1)
            if not terminal.truth[position]:
                continue
            _allocation[Tuple[Int, Int]](len(result.projection.variables), budget)
            var path = List[Tuple[Int, Int]]()
            var cursor = index
            var slot = position
            while cursor >= 0:
                budget.charge(1)
                ref node = result.stacks[cursor]
                path.append((cursor, slot))
                slot = node.slot
                cursor = node.parent
            _allocation[RealAlgebraic](len(path), budget)
            var point = List[RealAlgebraic](capacity=len(path))
            for k in range(len(path) - 1, -1, -1):
                ref node = result.stacks[path[k][0]]
                slot = path[k][1]
                if slot % 2:
                    var field = node.field.copy()
                    point.append(_absolute_section(node.stack, slot // 2, field, budget))
                else:
                    point.append(_rational(node.stack.sectors[slot // 2], budget))
            return Optional(point^)
    return None
