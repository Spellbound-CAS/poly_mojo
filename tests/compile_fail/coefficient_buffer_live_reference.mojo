# expect: invalidated interior reference
from polynomials import CoefficientBuffer

def main():
    var values: List[Int] = [1, 2, 3]
    var buffer = CoefficientBuffer(values^)
    ref entry = buffer[0]
    _ = buffer.pop()
    print(entry)
