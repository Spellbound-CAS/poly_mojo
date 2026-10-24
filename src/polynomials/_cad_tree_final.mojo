"""Private multivariate CAD with relative sections and a final acceptance gate.

All real space is represented by a complete stack tree. A terminal section
is just a relative root; only a section with children extends its sample field.
Original input signs are checked separately from Lazard lifting polynomials.
The initial private query language is a conjunction of input sign masks.
"""
from .coefficients import QQ, NoContext
from .dense import DensePoly
from .sparse import SparsePoly
from .errors import Budget, PolyError, PolyErrorKind, _allocation, _sum, _product
from ._algebraic import _Frame
from ._algebraic_filter import _constant
from ._cad_sample_field import _SampleField, _sample_field, _check_sample_field, _specialize_field, _lazard_field
from ._cad_root_cache import _RootCountCache
from ._cad_field_stack import _check_stack_in_field_cached, _FieldStack, _construct_field_stack, _check_stack_in_field
from ._cad_field_extension import _check_extension_data, _extend_section_owned, _check_extension, _FieldExtension, _FieldExtensionCache
from ._cad_projection_final import _FinalProjection, _project_values, _check_projection_values
from ._root_ops import _compare
from ._cad_field_verify import _CheckedFieldModuli


@fieldwise_init
struct _TreeStack(Copyable, Movable):
    var parent: Int
    var slot: Int
    var level: Int
    var field: _SampleField
    var parent_generator: DensePoly[QQ]
    var orders: List[List[Int]]
    var stack: _FieldStack
    var children: List[Int]
    var truth: List[Bool]


@fieldwise_init
struct _FinalCAD(Copyable, Movable):
    var inputs: List[SparsePoly[QQ]]
    var masks: List[Int]
    var projection: _FinalProjection
    var stacks: List[_TreeStack]
    var satisfiable: Bool
    var materializations: Int
    var retries: Int
    var field_map_hits: Int


def _query(masks: List[Int], signs: List[Int]) -> Bool:
    for i in range(len(masks)):
        if not (masks[i] & (1 << (signs[i] + 1))):
            return False
    return True


def _query_tail(masks: List[Int], signs: List[Int8], input_count: Int) -> Bool:
    # Shape and sign encoding are checked by the caller. Read the original
    # input suffix directly instead of allocating a new list for every cell.
    for i in range(len(masks)):
        if not (masks[i] & (1 << (Int(signs[len(signs) - input_count + i]) + 1))):
            return False
    return True


def _query_valid(inputs: List[SparsePoly[QQ]], masks: List[Int], mut budget: Budget) raises PolyError -> Bool:
    if len(masks) != 0 and len(masks) != len(inputs):
        return False
    for mask in masks:
        budget.charge(1)
        if mask < 0 or mask > 7:
            return False
    return True


def _tree_inputs(inputs: List[SparsePoly[QQ]], projection: _FinalProjection, level: Int, field: _SampleField, mut budget: Budget) raises PolyError -> Tuple[List[List[Int]], List[List[DensePoly[QQ]]]]:
    var orders = List[List[Int]]()
    var specialized = List[List[DensePoly[QQ]]]()
    var variable = projection.variables[level]
    for f in projection.levels[level]:
        var lifted = _lazard_field(f, variable, field, budget)
        _allocation[List[Int]](1, budget)
        _allocation[List[DensePoly[QQ]]](1, budget)
        # The tuple is temporary; transfer its buffers instead of cloning them.
        var lifted_orders = List[Int]()
        var lifted_coefficients = List[DensePoly[QQ]]()
        swap(lifted_orders, lifted[0])
        swap(lifted_coefficients, lifted[1])
        orders.append(lifted_orders^)
        specialized.append(lifted_coefficients^)
    if level == len(projection.variables) - 1:
        # The original polynomial may vanish identically although its Lazard
        # derivative does not. Its signs come from direct specialization.
        for f in inputs:
            _allocation[List[DensePoly[QQ]]](1, budget)
            specialized.append(_specialize_field(f, variable, field, budget))
    return (orders^, specialized^)


def _terminal_signs(stack: _TreeStack, input_count: Int, position: Int, mut budget: Budget) raises PolyError -> List[Int]:
    ref row = stack.stack.signs[position]
    if len(row) < input_count:
        raise PolyError(PolyErrorKind.DOMAIN, "missing original input signs")
    _allocation[Int](input_count, budget)
    var result = List[Int]()
    for i in range(len(row) - input_count, len(row)):
        result.append(Int(row[i]))
    return result^


def _new_tree_stack(inputs: List[SparsePoly[QQ]], projection: _FinalProjection, masks: List[Int], parent: Int, slot: Int, level: Int, var field: _SampleField, var parent_image: DensePoly[QQ], mut budget: Budget, interval: Bool = False) raises PolyError -> _TreeStack:
    var specialized = _tree_inputs(inputs, projection, level, field, budget)
    var stack = _construct_field_stack(specialized[1], field, budget, interval=interval)
    var orders = List[List[Int]]()
    swap(orders, specialized[0])
    var result = _TreeStack(parent, slot, level, field^, parent_image^, orders^, stack^, List[Int](), List[Bool]())
    if level == len(projection.variables) - 1:
        _allocation[Bool](len(result.stack.signs), budget)
        for position in range(len(result.stack.signs)):
            result.truth.append(_query_tail(masks, result.stack.signs[position], len(inputs)))
    return result^


def _construct_final_cad(inputs: List[SparsePoly[QQ]], variables: List[Int64], masks: List[Int], mut budget: Budget, *, brown: Bool = True, interval: Bool = False) raises PolyError -> _FinalCAD:
    if len(variables) == 0 or not _query_valid(inputs, masks, budget):
        raise PolyError(PolyErrorKind.DOMAIN, "private CAD requires variables and valid sign masks")
    var projection = _project_values(inputs, variables, budget, brown=brown)
    var start = _sample_field(_Frame(), budget)
    _allocation[_TreeStack](1, budget)
    var stacks: List[_TreeStack] = [_new_tree_stack(inputs, projection, masks, -1, -1, 0, start^, DensePoly[QQ](0, NoContext()), budget, interval)]
    var cursor = 0
    var materializations = 0
    var satisfiable = False
    var field_maps = _FieldExtensionCache()
    while cursor < len(stacks):
        budget.charge(1)
        var level = stacks[cursor].level
        if level == len(variables) - 1:
            for truth in stacks[cursor].truth:
                satisfiable = satisfiable or truth
            cursor += 1
            continue
        var size = _sum(_product(2, len(stacks[cursor].stack.roots)), 1)
        for slot in range(size):
            var field = stacks[cursor].field.copy()
            var image = DensePoly[QQ](0, NoContext())
            if slot % 2:
                var root = slot // 2
                ref section = stacks[cursor].stack.roots[root]
                if len(stacks[cursor].stack.basis[section.basis]) > 2 and section.low != section.high:
                    materializations = _sum(materializations, 1)
                var extension = _extend_section_owned(stacks[cursor].stack, root, variables[level], field^, field_maps, budget)
                field = _SampleField(DensePoly[QQ](0, NoContext()), None, List[Int64](), List[DensePoly[QQ]]())
                swap(field, extension.child)
                swap(image, extension.parent_generator)
            else:
                field.variables.append(variables[level])
                field.coordinates.append(_constant(stacks[cursor].stack.sectors[slot // 2], budget))
            var next = _new_tree_stack(inputs, projection, masks, cursor, slot, level + 1, field^, image^, budget, interval)
            _allocation[Int](1, budget)
            stacks[cursor].children.append(len(stacks))
            _allocation[_TreeStack](1, budget)
            stacks.append(next^)
        cursor += 1
    _allocation[SparsePoly[QQ]](len(inputs), budget)
    _allocation[Int](len(masks), budget)
    return _FinalCAD(inputs.copy(), masks.copy(), projection^, stacks^, satisfiable, materializations, 0, field_maps.hits)


def _check_final_cad_with_cache(inputs: List[SparsePoly[QQ]], variables: List[Int64], masks: List[Int], result: _FinalCAD, mut root_counts: _RootCountCache, mut budget: Budget, *, vectorized: Bool = True) raises PolyError -> Bool:
    # Internal entry: root_counts must be created by this verification run.
    try:
        if not _query_valid(inputs, masks, budget) or result.inputs != inputs or result.masks != masks or result.projection.variables != variables or len(variables) == 0:
            return False
        if not _check_projection_values(inputs, result.projection, budget) or len(result.stacks) == 0:
            return False
        _allocation[Bool](len(result.stacks), budget)
        var reached = List[Bool](length=len(result.stacks), fill=False)
        reached[0] = True
        var satisfiable = False
        var moduli = _CheckedFieldModuli()
        for i in range(len(result.stacks)):
            budget.charge(1)
            ref node = result.stacks[i]
            if not reached[i] or node.level < 0 or node.level >= len(variables):
                return False
            var checked_field = node.field.copy()
            var modulus_checked = not node.field.embedding or moduli.contains(node.field.modulus, budget)
            if i == 0:
                if node.parent != -1 or node.slot != -1 or node.level != 0 or len(node.field.variables) or node.field.embedding or not _check_sample_field(node.field, budget):
                    return False
            else:
                if node.parent < 0 or node.parent >= i:
                    return False
                ref parent = result.stacks[node.parent]
                if node.level != parent.level + 1 or node.slot < 0 or node.slot >= len(parent.children) or parent.children[node.slot] != i:
                    return False
                if node.slot % 2:
                    ref root = parent.stack.roots[node.slot // 2]
                    if not _check_extension_data(parent.field, parent.stack.basis[root.basis], root, variables[parent.level], node.field, node.parent_generator, budget, modulus_checked=modulus_checked):
                        return False
                else:
                    # Compare the shared prefix directly; do not copy the
                    # entire parent field just to append one rational value.
                    if node.field.modulus != parent.field.modulus or len(node.field.variables) != len(parent.field.variables) + 1 or len(node.field.coordinates) != len(parent.field.coordinates) + 1:
                        return False
                    for k in range(len(parent.field.variables)):
                        if node.field.variables[k] != parent.field.variables[k] or node.field.coordinates[k] != parent.field.coordinates[k]:
                            return False
                    if node.field.variables[len(parent.field.variables)] != variables[parent.level] or node.field.coordinates[len(parent.field.coordinates)] != _constant(parent.stack.sectors[node.slot // 2], budget):
                        return False
                    if Bool(node.field.embedding) != Bool(parent.field.embedding):
                        return False
                    if node.field.embedding:
                        ref original = parent.field.embedding.value()
                        ref narrowed = node.field.embedding.value()
                        if narrowed.polynomial != original.polynomial or not _check_sample_field(node.field, budget, modulus_checked=modulus_checked):
                            return False
                        if _compare(narrowed.low, original.low, budget) < 0 or _compare(narrowed.high, original.high, budget) > 0:
                            return False
            if len(node.field.variables) != node.level:
                return False
            for k in range(node.level):
                if node.field.variables[k] != variables[k]:
                    return False
            if not modulus_checked:
                moduli.add(node.field.modulus, budget)
            var specialized = _tree_inputs(inputs, result.projection, node.level, checked_field, budget)
            if node.orders != specialized[0] or not _check_stack_in_field_cached(specialized[1], node.stack, checked_field, root_counts, budget, vectorized=vectorized):
                return False
            if node.level == len(variables) - 1:
                if len(node.children) or len(node.truth) != len(node.stack.signs):
                    return False
                for position in range(len(node.stack.signs)):
                    var truth = _query_tail(masks, node.stack.signs[position], len(inputs))
                    if node.truth[position] != truth:
                        return False
                    satisfiable = satisfiable or truth
            else:
                if len(node.children) != _sum(_product(2, len(node.stack.roots)), 1) or len(node.truth):
                    return False
                for child in node.children:
                    if child <= i or child >= len(result.stacks) or reached[child]:
                        return False
                    reached[child] = True
        return result.satisfiable == satisfiable and root_counts.flush(budget, vectorized)
    except e:
        if e.kind == PolyErrorKind.DOMAIN:
            return False
        raise e^


def _check_final_cad(inputs: List[SparsePoly[QQ]], variables: List[Int64], masks: List[Int], result: _FinalCAD, mut budget: Budget, *, vectorized: Bool = True) raises PolyError -> Bool:
    var root_counts = _RootCountCache()
    return _check_final_cad_with_cache(inputs, variables, masks, result, root_counts, budget, vectorized=vectorized)


def _verified_final_cad(inputs: List[SparsePoly[QQ]], variables: List[Int64], masks: List[Int], mut budget: Budget, *, brown: Bool = True, interval: Bool = False) raises PolyError -> _FinalCAD:
    var retries = 0
    for attempt in range(3):
        if attempt == 1 and not interval or attempt == 2 and not brown:
            continue
        try:
            var result = _construct_final_cad(inputs, variables, masks, budget, brown=brown and attempt != 2, interval=interval and attempt == 0)
            result.retries = retries
            if _check_final_cad(inputs, variables, masks, result, budget):
                return result^
        except e:
            if e.kind != PolyErrorKind.DOMAIN:
                raise e^
        # At most an exact retry, then a plain-Lazard exact retry. All work
        # shares the caller's budget. Imported data is only checked/rejected.
        retries += 1
    raise PolyError(PolyErrorKind.DOMAIN, "final CAD verification failed")
