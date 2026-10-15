# output: 6
from polynomials import CoefficientBuffer
from polynomials.coefficient_buffer import _view
from std.iter import iter, next

def main() raises:
    var values: List[Int] = [1, 2, 3]
    var buffer = CoefficientBuffer(values^)
    var cursor = iter(buffer)
    var total = next(cursor)
    buffer.append(9)
    total += next(cursor)
    total += next(cursor)
    var view = _view(buffer)
    assert view[0] == 1
    buffer.clear()
    assert len(buffer) == 0
    assert total == 6
    print(total)
