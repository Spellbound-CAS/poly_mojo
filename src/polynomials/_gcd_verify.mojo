"""Exact verification of modular gcd candidates."""
from .errors import Budget, PolyError
from .coefficients import ZZ
from .sparse import SparsePoly
from .division import divide_exact


def _divides_both(f: SparsePoly[ZZ], h: SparsePoly[ZZ], g: SparsePoly[ZZ], mut budget: Budget) raises PolyError -> Bool:
    """Whether the nonzero candidate g divides both f and h exactly over the integers."""
    return Bool(divide_exact(f, g, budget)) and Bool(divide_exact(h, g, budget))
