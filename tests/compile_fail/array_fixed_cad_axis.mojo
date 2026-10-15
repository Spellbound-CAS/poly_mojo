# expect: fixed CAD axes must be Int
from poly import cad

def main() raises:
    _ = cad.decompose([[(1, 1)]], fixed=[(True, 1)])
