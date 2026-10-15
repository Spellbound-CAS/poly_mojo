"""Deterministic budgeted lexicographic subsets without machine-word masks."""
from .errors import Budget, PolyError, PolyErrorKind, _allocation, _charge_size, _sum


struct _Subsets(Movable):
    var total: Int
    var indices: List[Int]

    def __init__(out self, total: Int, size: Int, mut budget: Budget) raises PolyError:
        if total < 0 or size < 0 or size > total:
            raise PolyError(PolyErrorKind.DOMAIN, "invalid subset dimensions")
        self.total = total
        _allocation[Int](size, budget)
        self.indices = List[Int](capacity=size)
        for i in range(size):
            budget.charge(1)
            self.indices.append(i)

    def advance(mut self, mut budget: Budget) raises PolyError -> Bool:
        var size = len(self.indices)
        var position = size - 1
        while position >= 0:
            budget.charge(1)
            if self.indices[position] < self.total - size + position:
                break
            position -= 1
        if position < 0:
            return False
        # Charge the suffix before changing it, leaving the cursor intact on
        # exhaustion as well as preserving the borrowed active-factor list.
        _charge_size(size - position, budget)
        self.indices[position] += 1
        for i in range(position + 1, size):
            self.indices[i] = self.indices[i - 1] + 1
        return True


struct _Recombination(Movable):
    # Only selection bookkeeping is shared. Callers retain their distinct
    # reconstruction, complementary-degree and exact-divisibility rules.
    var active: List[Int]
    var selection: _Subsets

    def __init__(out self, count: Int, mut budget: Budget) raises PolyError:
        _allocation[Int](count, budget)
        self.active = List[Int](capacity=count)
        for i in range(count):
            budget.charge(1)
            self.active.append(i)
        self.selection = _Subsets(count, 1 if count > 1 else 0, budget)

    def ready(self) -> Bool:
        return len(self.selection.indices) > 0

    def advance(mut self, mut budget: Budget) raises PolyError:
        if self.selection.advance(budget):
            return
        var size = _sum(len(self.selection.indices), 1)
        if size > len(self.active) // 2:
            size = 0
        var next_selection = _Subsets(len(self.active), size, budget)
        self.selection = next_selection^

    def accept(mut self, mut budget: Budget) raises PolyError:
        # Remove the selected indices even when the caller reconstructed the
        # complementary polynomial. Commit paired state only after all charges.
        var size = len(self.selection.indices)
        if size == 0:
            raise PolyError(PolyErrorKind.DOMAIN, "no recombination selection")
        var remaining = len(self.active) - size
        _allocation[Int](remaining, budget)
        var active = List[Int](capacity=remaining)
        var cursor = 0
        for i in range(len(self.active)):
            budget.charge(1)
            if cursor < size and self.selection.indices[cursor] == i:
                cursor += 1
            else:
                active.append(self.active[i])
        var selection = _Subsets(remaining, 1 if remaining > 1 else 0, budget)
        self.active = active^
        self.selection = selection^
