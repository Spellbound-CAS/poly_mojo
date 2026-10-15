"""Reusable scalar slots for a fixed numerical instruction schedule."""
from .errors import Budget, PolyError, PolyErrorKind, _allocation
from ._numeric_program import _NumericProgram


@fieldwise_init
struct _NumericStorage(Movable):
    var slots: List[Int]
    var count: Int


@always_inline
def _instruction_source[reuse: Bool](program: _NumericProgram, i: Int) -> Int:
    comptime if reuse:
        if len(program.aliases) and program.aliases[i] >= 0:
            return program.aliases[i]
    return i


def _storage_plan[reuse: Bool = False](program: _NumericProgram, mut budget: Budget) raises PolyError -> _NumericStorage:
    var n = len(program.instructions)
    if program.result < 0 or program.result >= n:
        raise PolyError(PolyErrorKind.DOMAIN, "invalid numerical program result")
    comptime if reuse:
        if len(program.aliases) and len(program.aliases) != n:
            raise PolyError(PolyErrorKind.DOMAIN, "invalid numerical aliases")
    _allocation[Int](n, budget)
    var last_use = List[Int](length=n, fill=-1)
    for i in range(n):
        budget.charge(1)
        var source = _instruction_source[reuse](program, i)
        if source != i:
            if source < 0 or source >= i or _instruction_source[reuse](program, source) != source:
                raise PolyError(PolyErrorKind.DOMAIN, "invalid numerical alias source")
            last_use[source] = i
            continue
        ref instruction = program.instructions[i]
        if instruction.operation:
            if instruction.left < 0 or instruction.right < 0 or instruction.left >= i or instruction.right >= i:
                raise PolyError(PolyErrorKind.DOMAIN, "invalid numerical program operand")
            last_use[_instruction_source[reuse](program, instruction.left)] = i
            last_use[_instruction_source[reuse](program, instruction.right)] = i
    last_use[_instruction_source[reuse](program, program.result)] = n
    _allocation[Int](n, budget)
    _allocation[Int](n, budget)
    var slots = List[Int](capacity=n)
    var free = List[Int](capacity=n)
    var count = 0
    for i in range(n):
        budget.charge(1)
        var source = _instruction_source[reuse](program, i)
        if source != i:
            slots.append(slots[source])
            if last_use[source] == i:
                free.append(slots[source])
            continue
        var slot: Int
        if len(free):
            slot = free.pop()
        else:
            slot = count
            count += 1
        slots.append(slot)
        ref instruction = program.instructions[i]
        if instruction.operation:
            var left = _instruction_source[reuse](program, instruction.left)
            var right = _instruction_source[reuse](program, instruction.right)
            if last_use[left] == i:
                free.append(slots[left])
            if right != left and last_use[right] == i:
                free.append(slots[right])
        if last_use[i] < i:
            free.append(slot)
    return _NumericStorage(slots^, count)
