# expect: cannot call function that may raise 'Error' in context that supports an error type of 'PolyError'
"""An apn_mojo error cannot escape a typed polynomial boundary."""
from apn_mojo import Rational
from polynomials import PolyError

def invalid(a: Rational, b: Rational) raises PolyError -> Rational:
    return a + b

def main():
    pass
