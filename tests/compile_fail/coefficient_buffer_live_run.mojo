# expect: invalidated interior reference
from polynomials import CoefficientBuffer
from polynomials.coefficient_buffer import _view

def main():
    var values: List[Int] = [1, 2, 3]
    var buffer = CoefficientBuffer(values^)
    var run = _view(buffer)._run(0)
    buffer.clear()
    print(run[0])
