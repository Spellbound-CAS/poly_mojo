"""Iterative simultaneous Hensel lifting by base-prime correction digits (Z-4)."""
from apn_mojo import Integer
from ._integer_ops import _primitive, _divide_scalar
from .coefficients import ZZ, ModInt, NoContext, _work
from .dense import DensePoly, _slice, _trim
from .coefficient_buffer import _view
from .errors import Budget, PolyError, PolyErrorKind, _allocation
from ._coefficient_ops import _quotient
from ._dense_ops import _constant, _remainder, _exact
from ._euclidean import _extended
from ._images import _mod_image
from .primes import _symmetric


def _lift_image(f: DensePoly[ModInt], mut budget: Budget) raises PolyError -> DensePoly[ZZ]:
    _allocation[ZZ](len(f.coefficients), budget)
    var values = List[ZZ](capacity=len(f.coefficients))
    for coefficient_index in range(len(f.coefficients)):
        ref c = f.coefficients[coefficient_index]
        budget.charge(1)
        values.append(ZZ(Integer(Int(c.value))))
    return DensePoly[ZZ]._from_canonical(f.variable, NoContext(), values^)


def _truncate(f: DensePoly[ZZ], modulus: ZZ, mut budget: Budget) raises PolyError -> DensePoly[ZZ]:
    _allocation[ZZ](len(f.coefficients), budget)
    var values = List[ZZ](capacity=len(f.coefficients))
    for coefficient_index in range(len(f.coefficients)):
        ref c = f.coefficients[coefficient_index]
        budget.charge(1)
        values.append(ZZ(_symmetric(c.value, modulus.value, budget)))
    return DensePoly[ZZ].from_coefficients(f.variable, NoContext(), values^, budget)


def _trial_divide(f: DensePoly[ZZ], g: DensePoly[ZZ], mut budget: Budget) raises PolyError -> Optional[DensePoly[ZZ]]:
    f._same(g)
    if g.is_zero():
        raise PolyError(PolyErrorKind.DOMAIN, "zero integer polynomial divisor")
    var n = max(f.degree() - g.degree() + 1, 0)
    _allocation[ZZ](n, budget)
    var quotient = List[ZZ](length=n, fill=ZZ.zero(f.ring))
    var remainder = _slice(_view(f.coefficients), 0, len(f.coefficients), budget)
    while len(remainder) >= len(g.coefficients):
        budget.charge(1)
        var shift = len(remainder) - len(g.coefficients)
        var scalar = _quotient(remainder[len(remainder) - 1], g.coefficients[g.degree()], budget)
        if not scalar:
            return None
        var q = scalar.take()
        for i in range(len(g.coefficients)):
            budget.charge(1)
            remainder[shift + i] = remainder[shift + i].sub(q.mul(g.coefficients[i], budget), budget)
        quotient[shift] = q^
        _trim(remainder, budget)
    if len(remainder):
        return None
    return DensePoly[ZZ].from_coefficients(f.variable, f.ring, quotient^, budget)


@fieldwise_init
struct _HenselState(Movable):
    var lifted: List[DensePoly[ZZ]]
    var modulus: ZZ
    var inverses: List[DensePoly[ModInt]]


def _hensel_start(factors: List[DensePoly[ModInt]], mut budget: Budget) raises PolyError -> _HenselState:
    if len(factors) == 0:
        raise PolyError(PolyErrorKind.DOMAIN, "Hensel lifting requires modular factors")
    _allocation[DensePoly[ZZ]](len(factors), budget)
    var lifted = List[DensePoly[ZZ]](capacity=len(factors))
    for g in factors:
        budget.charge(1)
        lifted.append(_lift_image(g, budget))
    return _HenselState(lifted^, ZZ(Integer(Int(factors[0].ring.value))), List[DensePoly[ModInt]]())


def _hensel_extend(f: DensePoly[ZZ], factors: List[DensePoly[ModInt]], bound: ZZ, mut state: _HenselState, mut budget: Budget) raises PolyError:
    # The original leading coefficient is carried separately. For E=(f-lc*ΠG)/m,
    # delta_i = E * (lc*Π_{j!=i} g_j)^(-1) mod g_i. The Chinese remainder
    # identity then gives f=lc*Π(G_i+m*delta_i) mod m*p. Each correction has
    # degree below g_i, preserving monicity and all factor degrees.
    var ctx = factors[0].ring
    var p = ZZ(Integer(Int(ctx.value)))
    _work(state.modulus.size(), bound.size(), False, budget)
    if state.modulus.value > bound.value:
        return
    if len(state.inverses) == 0:
        var image = _mod_image(f, ctx, budget)
        _allocation[DensePoly[ModInt]](len(factors), budget)
        for g in factors:
            budget.charge(1)
            var complement = _exact(image, g, budget)
            var egcd = _extended(complement, g, budget)
            state.inverses.append(_remainder(egcd[1], g, budget))
    while True:
        budget.charge(1)
        _work(state.modulus.size(), bound.size(), False, budget)
        if state.modulus.value > bound.value:
            break
        var product = _constant(f, f.coefficients[f.degree()], budget)
        for g in state.lifted:
            budget.charge(1)
            product = product.mul(g, budget)
        var error = _mod_image(_divide_scalar(f.sub(product, budget), state.modulus, budget), ctx, budget)
        var next_modulus = state.modulus.mul(p, budget)
        for i in range(len(state.lifted)):
            budget.charge(1)
            var correction = _remainder(error.mul(state.inverses[i], budget), factors[i], budget)
            state.lifted[i] = _truncate(state.lifted[i].add(_lift_image(correction, budget).scale(state.modulus, budget), budget), next_modulus, budget)
        state.modulus = next_modulus^


def _hensel(f: DensePoly[ZZ], factors: List[DensePoly[ModInt]], bound: ZZ, mut budget: Budget) raises PolyError -> Tuple[List[DensePoly[ZZ]], ZZ]:
    var state = _hensel_start(factors, budget)
    _hensel_extend(f, factors, bound, state, budget)
    var lifted = List[DensePoly[ZZ]]()
    var modulus = ZZ.zero(NoContext())
    swap(lifted, state.lifted)
    swap(modulus, state.modulus)
    return (lifted^, modulus^)
