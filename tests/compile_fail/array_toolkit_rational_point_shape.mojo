# expect: rational evaluation point has the wrong number of coordinates
from poly import rational

def main() raises:
    var p = rational.from_terms([((1, 0), 1)], [((0, 0), 1)])
    _ = rational.evaluate(p, (1,))
