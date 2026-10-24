"""Reuse the exact batched endpoint kernel for rational sample fields."""
from .coefficients import QQ, NoContext
from .dense import DensePoly
from .errors import Budget, PolyError, PolyErrorKind, _allocation
from ._cad_final import _FinalStack, _FinalRoot, _check_roots
from ._cad_final_signs import _SignStats
from ._algebraic import _norm_key
from .certificates import _dense_equal
from ._dense_ops import _copy as _copy_poly
from std.collections import Dict


struct _CheckedFieldModuli(Movable):
    """Verifier-owned irreducibility facts; collisions use exact comparison."""
    var heads: Dict[Int, Int]
    var next: List[Int]
    var moduli: List[DensePoly[QQ]]

    def __init__(out self):
        self.heads = Dict[Int, Int]()
        self.next = List[Int]()
        self.moduli = List[DensePoly[QQ]]()

    def key(self, modulus: DensePoly[QQ], mut budget: Budget) raises PolyError -> Int:
        modulus.check_invariants(budget)
        try:
            return _norm_key(modulus)
        except e:
            raise PolyError(PolyErrorKind.CAPACITY, String(e))

    def contains(self, modulus: DensePoly[QQ], mut budget: Budget) raises PolyError -> Bool:
        var key = self.key(modulus, budget)
        var index = self.heads.get(key, -1)
        while index >= 0:
            budget.charge(1)
            if _dense_equal(modulus, self.moduli[index], budget):
                return True
            index = self.next[index]
        return False

    def add(mut self, modulus: DensePoly[QQ], mut budget: Budget) raises PolyError:
        var key = self.key(modulus, budget)
        _allocation[Int](1, budget)
        _allocation[DensePoly[QQ]](1, budget)
        self.next.append(self.heads.get(key, -1))
        self.heads[key] = len(self.moduli)
        self.moduli.append(_copy_poly(modulus, budget))


def _rational_root_check(basis: List[List[DensePoly[QQ]]], roots: List[_FinalRoot], mut stats: _SignStats, mut budget: Budget, vectorized: Bool) raises PolyError -> Bool:
    """Check root completeness after the caller checks basis identities/geometry.

    Each coefficient must be a rational constant. This only repacks the final
    basis and brackets; no producer facts or arithmetic histories are trusted.
    """
    _allocation[DensePoly[QQ]](len(basis), budget)
    var polynomials = List[DensePoly[QQ]]()
    for member in basis:
        _allocation[QQ](len(member), budget)
        var coefficients = List[QQ]()
        for c in member:
            c.check_invariants(budget)
            if c.variable != 0 or c.degree() > 0:
                return False
            coefficients.append(c.coefficients[0].copy() if not c.is_zero() else QQ.zero(NoContext()))
        polynomials.append(DensePoly[QQ].from_coefficients(0, NoContext(), coefficients^, budget))
    _allocation[_FinalRoot](len(roots), budget)
    var data = _FinalStack(polynomials^, List[QQ](), List[List[Int]](), roots.copy(), List[QQ](), List[List[Int]]())
    return _check_roots(data, stats, budget, vectorized)
