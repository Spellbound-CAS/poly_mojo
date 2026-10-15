# expect: invalidated interior reference
from polynomials import CoefficientBuffer
from polynomials.coefficient_buffer import _view

def main():
    var values: List[Int] = [1, 2, 3]
    var buffer = CoefficientBuffer(values^)
    var view = _view(buffer)
    buffer.append(9)
    print(view[0])
