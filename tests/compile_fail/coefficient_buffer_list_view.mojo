# expect: invalidated interior reference
from polynomials.coefficient_buffer import _view

def main():
    var values: List[Int] = [1, 2, 3]
    var view = _view(values)
    values.append(9)
    print(view[0])
