# expect: wrong number of coordinates
import poly


def main() raises:
    _ = poly.sparse.evaluate_many([((1, 0), 2)], [(1, 2, 3)])
