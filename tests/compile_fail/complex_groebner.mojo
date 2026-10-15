# expect: CoefficientField
"""Rounded complex domains cannot instantiate exact Gröbner bases."""
from polynomials import ComplexDomain, GroebnerBasis


def reject(value: GroebnerBasis[ComplexDomain]):
    _ = value


def main():
    pass
