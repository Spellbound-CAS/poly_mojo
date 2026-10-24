"""Cheap reductions ahead of multivariate integer gcds.

Deflation: when every exponent of a variable in both inputs is a multiple
of J, substituting x^J -> x is an injective ring map that preserves gcds
over a field of characteristic zero (it maps coprime pairs to coprime
pairs), so gcd(f, h) is the inflated gcd of the deflated pair.

Image test: if lc_x(f) and lc_x(h) survive specialization of the other
variables at a point modulo p, the true gcd's image has its full degree in
x and divides both images. A constant image gcd therefore proves that the
gcd has degree zero in x, leaving only the gcd of the contents.
"""
from apn_mojo import Integer
from .errors import Budget, PolyError, _allocation, _charge_size
from .coefficients import CoefficientRing, ZZ
from .sparse import SparsePoly, Term
from .monomial import Monomial
from .primes import prime, _residue
from ._nmod_euclid import _Field, _gcd_words


def _integer_gcd(a: Int, b: Int) -> Int:
    var x = a
    var y = b
    while y:
        var t = x % y
        x = y
        y = t
    return x


def _exponent_gcds[C: CoefficientRing](p: SparsePoly[C], variables: List[Int64], mut out: List[Int], mut ones: Int, mut budget: Budget) raises PolyError:
    """Fold p's exponents into the gcds in out; `ones` counts the gcds equal to one, which no exponent changes."""
    for t in p.terms:
        if ones == len(out):
            return
        budget.charge(1)
        var position = 0
        for pair in t.monomial.pairs:
            while variables[position] != pair[0]:
                position += 1
            if out[position] != 1:
                out[position] = _integer_gcd(out[position], Int(pair[1]))
                if out[position] == 1:
                    ones += 1


def _deflation[C: CoefficientRing](f: SparsePoly[C], h: SparsePoly[C], variables: List[Int64], mut budget: Budget) raises PolyError -> List[Int]:
    """The gcd of all exponents of each variable in both inputs (one for an absent variable)."""
    _allocation[Int](len(variables), budget)
    var out = List[Int](length=len(variables), fill=0)
    var ones = 0
    _exponent_gcds(f, variables, out, ones, budget)
    _exponent_gcds(h, variables, out, ones, budget)
    for i in range(len(out)):
        if out[i] == 0:
            out[i] = 1
    return out^


def _scaled[C: CoefficientRing](p: SparsePoly[C], variables: List[Int64], exponents: List[Int], inflate: Bool, mut budget: Budget) raises PolyError -> SparsePoly[C]:
    """Divide (deflate) or multiply (inflate) each variable's exponents by its factor."""
    _allocation[Term[C]](len(p.terms), budget)
    var terms = List[Term[C]](capacity=len(p.terms))
    for t in p.terms:
        budget.charge(1)
        _charge_size(len(t.monomial.pairs), budget)
        var pairs = List[Tuple[Int64, UInt32]](capacity=len(t.monomial.pairs))
        var position = 0
        for pair in t.monomial.pairs:
            while variables[position] != pair[0]:
                position += 1
            var factor = UInt32(exponents[position])
            pairs.append((pair[0], pair[1] * factor if inflate else pair[1] // factor))
        terms.append(Term[C](Monomial(pairs^, budget), t.coefficient.copy()))
    return SparsePoly[C].from_terms(terms^, p.order, p.ring, budget)


def _deflates(exponents: List[Int]) -> Bool:
    for e in exponents:
        if e > 1:
            return True
    return False


def _degrees(p: SparsePoly[ZZ], variables: List[Int64], mut out: List[Int], mut budget: Budget) raises PolyError:
    """Raise out[i] to the degree of p in variables[i]."""
    for t in p.terms:
        budget.charge(1)
        var position = 0
        for pair in t.monomial.pairs:
            while variables[position] != pair[0]:
                position += 1
            out[position] = max(out[position], Int(pair[1]))


def _residues(p: SparsePoly[ZZ], field: _Field, mut budget: Budget) raises PolyError -> List[UInt64]:
    _allocation[UInt64](len(p.terms), budget)
    var out = List[UInt64](capacity=len(p.terms))
    for t in p.terms:
        budget.charge(1)
        out.append(field.enter(_residue(t.coefficient.value, field.p, budget)))
    return out^


def _univariate_image(p: SparsePoly[ZZ], residues: List[UInt64], degree: Int, at: Int, variables: List[Int64], powers: List[List[UInt64]], field: _Field, mut budget: Budget) raises PolyError -> List[UInt64]:
    """p modulo the field at the powers' points for every variable but variables[at]."""
    _allocation[UInt64](degree + 1, budget)
    var out = List[UInt64](length=degree + 1, fill=0)
    for i in range(len(p.terms)):
        budget.charge(1)
        ref pairs = p.terms[i].monomial.pairs
        _charge_size(len(pairs), budget)
        var value = residues[i]
        var e = 0
        var position = 0
        for pair in pairs:
            while variables[position] != pair[0]:
                position += 1
            if position == at:
                e = Int(pair[1])
            else:
                value = field.mul(value, powers[position][Int(pair[1])])
        out[e] = field.add(out[e], value)
    return out^


def _image_coprime(f: SparsePoly[ZZ], h: SparsePoly[ZZ], variable: Int64, variables: List[Int64], tops_f: List[Int], tops_h: List[Int], mut budget: Budget) raises PolyError -> Bool:
    """Whether one specialization modulo a prime proves deg_x gcd(f, h) = 0 (False is inconclusive)."""
    var field = _Field.of(prime(0, budget))
    var at = 0
    while variables[at] != variable:
        at += 1
    var df = tops_f[at]
    var dh = tops_h[at]
    var residues_f = _residues(f, field, budget)
    var residues_h = _residues(h, field, budget)
    for attempt in range(3):
        budget.charge(1)
        # Deterministic points 2, 3, … shifted per attempt; the main variable's row is unused.
        var powers = List[List[UInt64]]()
        for position in range(len(variables)):
            var top = max(tops_f[position], tops_h[position])
            _allocation[UInt64](top + 1, budget)
            var row = List[UInt64](capacity=top + 1)
            var point = field.enter(UInt64(2 + position + attempt * len(variables)))
            var value = field.enter(1)
            for _ in range(top + 1):
                row.append(value)
                value = field.mul(value, point)
            powers.append(row^)
        var a = _univariate_image(f, residues_f, df, at, variables, powers, field, budget)
        var b = _univariate_image(h, residues_h, dh, at, variables, powers, field, budget)
        if a[df] == 0 or b[dh] == 0:
            continue
        return len(_gcd_words(a^, b^, field, budget)) == 1
    return False
