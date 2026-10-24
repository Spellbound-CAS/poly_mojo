"""The étale certificate for trace-derived rational bases, and every fallback to the S-pairs.

A traced basis over QQ lies in the ideal by its exact replay; `_etale` then
makes it a Gröbner basis when the generators are n polynomials in n
variables whose image modulo a word prime is a reduced algebra of the
basis's dimension (see `_groebner_check._etale`). Every unmet condition
must decline, and the S-pairs must still prove the basis: a non-square
system, a positive-dimensional ideal, a prime whose image is not reduced
though the ideal is radical over QQ, an ideal that is not radical, and
primes whose image has other leading monomials. On systems where it
applies, the certified basis must equal the default mode's, which the
S-pairs prove.
"""
from apn_mojo import Rational
from polynomials import Budget, QQ, NoContext, SparsePoly, Term, Monomial, MonomialOrder, groebner
from polynomials._groebner_check import _ETALE_START, _check_basis, _etale, _proves_reduced
from polynomials._packed_groebner import _Packed, _pack
from support.check import Checker
from support.corpus import parse_polynomial
from support.rng import SplitMix64


def term(n: Int, exponents: List[Int], coefficient: Int, order: MonomialOrder, mut work: Budget) raises -> SparsePoly[QQ]:
    var pairs = List[Tuple[Int64, UInt32]]()
    for i in range(n):
        if exponents[i]:
            pairs.append((Int64(i), UInt32(exponents[i])))
    var t = Term[QQ](Monomial(pairs^, work), QQ(Rational(coefficient)))
    return SparsePoly[QQ].from_terms([t^], order, NoContext(), work)


def unit(n: Int, i: Int) -> List[Int]:
    var e = List[Int](length=n, fill=0)
    e[i] = 1
    return e^


def plus(a: List[Int], b: List[Int]) -> List[Int]:
    var e = a.copy()
    for i in range(len(b)):
        e[i] += b[i]
    return e^


def katsura(n: Int, order: MonomialOrder, mut work: Budget) raises -> List[SparsePoly[QQ]]:
    var v = n + 1
    var polys = List[SparsePoly[QQ]]()
    var linear = term(v, List[Int](length=v, fill=0), -1, order, work)
    linear = linear.add(term(v, unit(v, 0), 1, order, work), work)
    for i in range(1, v):
        linear = linear.add(term(v, unit(v, i), 2, order, work), work)
    polys.append(linear^)
    for m in range(n):
        var p = term(v, unit(v, m), -1, order, work)
        for i in range(-n, n + 1):
            var j = m - i
            if abs(i) <= n and abs(j) <= n:
                p = p.add(term(v, plus(unit(v, abs(i)), unit(v, abs(j))), 1, order, work), work)
        polys.append(p^)
    return polys^


def quadratics(n: Int, seed: UInt64, order: MonomialOrder, mut work: Budget) raises -> List[SparsePoly[QQ]]:
    var rng = SplitMix64(seed)
    var polys = List[SparsePoly[QQ]]()
    for _ in range(n):
        var p = term(n, List[Int](length=n, fill=0), rng.below(19) - 9, order, work)
        for i in range(n):
            for j in range(i, n):
                p = p.add(term(n, plus(unit(n, i), unit(n, j)), rng.below(19) - 9, order, work), work)
        polys.append(p^)
    return polys^


def grevlex(n: Int, mut work: Budget) raises -> MonomialOrder:
    var variables = List[Int64]()
    for i in range(n):
        variables.append(Int64(i))
    return MonomialOrder(MonomialOrder.GREVLEX, variables^, work)


def parsed(texts: List[String], order: MonomialOrder, mut work: Budget) raises -> List[SparsePoly[QQ]]:
    var polys = List[SparsePoly[QQ]]()
    for text in texts:
        polys.append(parse_polynomial(text, order, work))
    return polys^


def etale(generators: List[SparsePoly[QQ]], basis: List[SparsePoly[QQ]], order: MonomialOrder, start: UInt64, mut work: Budget) raises -> Bool:
    var leads = List[_Packed[8]]()
    var terms = 0
    for g in basis:
        var packed = _pack[QQ, 8](g, order, work)
        leads.append(packed[len(packed) - 1].monomial)
        terms += len(packed)
    return _etale[8](generators, leads, terms, order, work, start)


def check_case(generators: List[SparsePoly[QQ]], order: MonomialOrder, expected: Bool, name: String, mut c: Checker, mut work: Budget, start: UInt64 = _ETALE_START) raises:
    """The traced basis: `_etale`'s verdict as expected, proved either way, equal to the default mode's basis, and accepted by the S-pairs."""
    var traced = groebner(generators, order, work, trace=True)
    var default = groebner(generators, order, work)
    c.check(traced.basis == default.basis, String(name, ": traced basis equals the default mode's"))
    c.check(etale(generators, traced.basis, order, start, work) == expected, String(name, ": étale verdict"))
    c.check(_proves_reduced(generators, traced.basis, order, work, start), String(name, ": proved, by étale or the S-pairs"))
    var pairs = _check_basis(generators, traced.basis, order, work)
    c.check(Bool(pairs) and pairs.value(), String(name, ": the S-pairs prove it too"))


def main() raises:
    var work = Budget()
    var c = Checker("etale")
    # Zero-dimensional square systems whose image is reduced.
    for n in range(2, 6):
        var order = grevlex(n + 1, work)
        check_case(katsura(n, order, work), order, True, String("katsura-", n), c, work)
    for seed in range(3):
        var order = grevlex(3, work)
        check_case(quadratics(3, UInt64(seed + 1), order, work), order, True, String("random quadratics ", seed), c, work)
    var o2 = grevlex(2, work)
    # Not square: a redundant third generator.
    var three = katsura(1, o2, work)
    three.append(three[0].add(three[1], work))
    check_case(three, o2, False, "three generators in two variables", c, work)
    # Positive-dimensional: the line x0 = 0, so no pure power of x1 leads.
    check_case(parsed(["x0^2*x1", "x0*x1^2"], o2, work), o2, False, "positive-dimensional", c, work)
    # Radical over QQ, roots 1 +- sqrt 7, but (x0 - 1)^2 modulo 7, which divides
    # no coefficient: from 9 the first prime is 7, and the minimal polynomial
    # is not squarefree there.
    var sqrt7 = parsed(["x0^2 - 2*x0 - 6", "x1 - x0"], o2, work)
    check_case(sqrt7.copy(), o2, False, "radical ideal, image not reduced modulo 7", c, work, start=9)
    check_case(sqrt7^, o2, True, "radical ideal at the word primes", c, work)
    # Not radical: a double root.
    check_case(parsed(["x0^2 - 2*x0 + 1", "x1"], o2, work), o2, False, "double root", c, work)
    # x0 + x1 and x0 + 36 x1 have basis {x0, x1}, but modulo 7 and 5, both
    # tried from 9, the second generator is the first: other leading monomials.
    var collapse = parsed(["x0 + x1", "x0 + 36*x1"], o2, work)
    check_case(collapse.copy(), o2, False, "leading monomials differ modulo both primes", c, work, start=9)
    # x0 + 8 x1 differs modulo 7 only: the second prime, 5, decides.
    check_case(parsed(["x0 + x1", "x0 + 8*x1"], o2, work), o2, True, "leading monomials differ modulo the first prime", c, work, start=9)
    check_case(collapse^, o2, True, "leading monomials agree at the word primes", c, work)
    c.finish()
