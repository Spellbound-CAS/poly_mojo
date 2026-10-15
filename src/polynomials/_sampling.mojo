"""Prepared monomial metadata for repeated finite-field specializations."""
from .errors import Budget, PolyError, PolyErrorKind, _allocation, _sum
from .coefficients import ModInt, Modulus
from .sparse import SparsePoly, Term
from .dense import DensePoly
from ._powers import _Powers, _powers


def _point_powers(point: List[UInt64], ctx: Modulus, mut budget: Budget) raises PolyError -> List[_Powers[ModInt]]:
    _allocation[_Powers[ModInt]](len(point), budget)
    var result = List[_Powers[ModInt]](capacity=len(point))
    for coordinate in point:
        budget.charge(1)
        result.append(_powers(ModInt(coordinate, ctx), budget))
    return result^


struct _SamplingPlan(Movable):
    var main: Int64
    var degree: Int
    var exponents: List[Int]
    var factors: List[List[Tuple[Int, UInt32]]]

    def __init__(out self, terms: List[Term[ModInt]], main: Int64, variables: List[Int64], mut budget: Budget) raises PolyError:
        self.main = main
        self.degree = -1
        _allocation[Int](len(terms), budget)
        _allocation[List[Tuple[Int, UInt32]]](len(terms), budget)
        self.exponents = List[Int](capacity=len(terms))
        self.factors = List[List[Tuple[Int, UInt32]]](capacity=len(terms))
        for term in terms:
            budget.charge(1)
            _allocation[Tuple[Int, UInt32]](len(term.monomial.pairs), budget)
            var factors = List[Tuple[Int, UInt32]](capacity=len(term.monomial.pairs))
            var exponent = 0
            for pair in term.monomial.pairs:
                budget.charge(1)
                if pair[0] == main:
                    exponent = Int(UInt64(pair[1]))
                    continue
                var lo = 0
                var hi = len(variables)
                while lo < hi:
                    budget.charge(1)
                    var mid = lo + (hi - lo) // 2
                    if variables[mid] < pair[0]:
                        lo = mid + 1
                    else:
                        hi = mid
                if lo == len(variables) or variables[lo] != pair[0]:
                    raise PolyError(PolyErrorKind.DOMAIN, "sampling variables do not cover the monomial")
                factors.append((lo, pair[1]))
            self.exponents.append(exponent)
            self.degree = max(self.degree, exponent)
            self.factors.append(factors^)

    def sample(self, index: Int, mut powers: List[_Powers[ModInt]], ctx: Modulus, mut budget: Budget) raises PolyError -> ModInt:
        var result = ModInt.one(ctx)
        for factor in self.factors[index]:
            budget.charge(1)
            result = result.mul(powers[factor[0]].power(UInt64(factor[1]), budget), budget)
        return result^

    def evaluate(self, f: SparsePoly[ModInt], mut powers: List[_Powers[ModInt]], mut budget: Budget) raises PolyError -> DensePoly[ModInt]:
        var size = _sum(self.degree, 1) if self.degree >= 0 else 0
        _allocation[ModInt](size, budget)
        var values = List[ModInt](length=size, fill=ModInt.zero(f.ring))
        for i in range(len(f.terms)):
            budget.charge(1)
            var k = self.exponents[i]
            values[k] = values[k].add(f.terms[i].coefficient.mul(self.sample(i, powers, f.ring, budget), budget), budget)
        return DensePoly[ModInt].from_coefficients(self.main, f.ring, values^, budget)
