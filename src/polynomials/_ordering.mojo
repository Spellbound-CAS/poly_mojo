"""Stable index sorting and heap repair without moving owned term buffers."""
from .errors import Budget, PolyError, _allocation, _product
from .monomial import MonomialOrder


trait _IndexedOrder:
    comptime Item: Copyable & Deinitable

    @staticmethod
    def compare(items: List[Self.Item], a: Int, b: Int, order: MonomialOrder, mut budget: Budget) raises PolyError -> Int:
        ...


def _stable_indices[A: _IndexedOrder](items: List[A.Item], order: MonomialOrder, mut budget: Budget) raises PolyError -> List[Int]:
    var n = len(items)
    _allocation[Int](_product(n, 2), budget)
    var indices = List[Int](capacity=n)
    var scratch = List[Int](length=n, fill=0)
    for i in range(n):
        budget.charge(1)
        indices.append(i)
    var width = 1
    while width < n:
        var start = 0
        while start < n:
            var middle = start + min(width, n - start)
            var end = middle + min(width, n - middle)
            var left = start
            var right = middle
            for k in range(start, end):
                budget.charge(1)
                if right == end or (left < middle and A.compare(items, indices[left], indices[right], order, budget) <= 0):
                    scratch[k] = indices[left]
                    left += 1
                else:
                    scratch[k] = indices[right]
                    right += 1
            start = end
        swap(indices, scratch)
        width = n if width > n // 2 else width * 2
    return indices^


def _sift_min[A: _IndexedOrder](mut heap: List[Int], items: List[A.Item], start: Int, order: MonomialOrder, mut budget: Budget) raises PolyError:
    var value = heap[start]
    var position = start
    while position < len(heap) // 2:
        budget.charge(1)
        var child = 2 * position + 1
        if child + 1 < len(heap) and A.compare(items, heap[child + 1], heap[child], order, budget) < 0:
            child += 1
        heap[position] = heap[child]
        position = child
    while position > start:
        budget.charge(1)
        var parent = (position - 1) // 2
        if A.compare(items, value, heap[parent], order, budget) >= 0:
            break
        heap[position] = heap[parent]
        position = parent
    heap[position] = value
