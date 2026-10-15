# expect: CAD point has the wrong number of coordinates
from poly import cad

def main() raises:
    var inputs = [[((1, 0), 1)]]
    _ = cad.cell_at(inputs, (0,))
