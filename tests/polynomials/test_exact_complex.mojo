"""Gaussian rational field laws and exact complex Gröbner certificates."""
from apn_mojo import ExactComplex, Integer, Rational
from polynomials import Budget, PolyError, PolyErrorKind, QQi, NoContext, MonomialOrder, SparsePoly, DensePoly, groebner, reduce, membership, MembershipCertificate, check_membership, check_list_division, divmod, check_division, evaluate_dense
from support.check import Checker
from support.groebner import verify_basis
from support.rng import SplitMix64


def scalar_op(op: Int, a: QQi, b: QQi, mut work: Budget) raises PolyError -> QQi:
    if op == 0:
        return a.add(b, work)
    if op == 1:
        return a.sub(b, work)
    if op == 2:
        return a.mul(b, work)
    if op == 3:
        return a.neg(work)
    return a.inverse(work)


def coefficients(mut c: Checker) raises:
    var work = Budget()
    var zero = QQi.zero(NoContext())
    var one = QQi.one(NoContext())
    var i = QQi(ExactComplex(0, 1))
    var a = QQi(ExactComplex(Rational(1, 2), Rational(2, 3)))
    var b = QQi(ExactComplex(Rational(3, 5), Rational(-4, 7)))
    c.check(zero.is_zero() and not i.is_zero(), "imaginary coefficients are nonzero")
    c.check(a.context() == NoContext() and zero.context() == one.context(), "unique field context")
    c.check(a.size() >= 4 and String(a).byte_length() > 0, "storage size and formatting")
    c.check(a.add(b, work) == QQi(ExactComplex(Rational(11, 10), Rational(2, 21))), "hand-computed sum")
    c.check(a.sub(b, work) == QQi(ExactComplex(Rational(-1, 10), Rational(26, 21))), "hand-computed difference")
    c.check(a.mul(b, work) == QQi(ExactComplex(Rational(143, 210), Rational(4, 35))), "hand-computed product")
    c.check(a.inverse(work) == QQi(ExactComplex(Rational(18, 25), Rational(-24, 25))), "hand-computed reciprocal")
    c.check(a.neg(work) == QQi(ExactComplex(Rational(-1, 2), Rational(-2, 3))), "hand-computed negative")
    c.check(i.mul(i, work) == one.neg(work) and i.inverse(work) == i.neg(work), "imaginary unit identities")
    c.check(QQi(ExactComplex(Rational(6, 8), Rational(-2, 6))) == QQi(ExactComplex(Rational(3, 4), Rational(-1, 3))), "canonical rational components")

    var rng = SplitMix64(29047)
    for _ in range(40):
        var x = QQi(ExactComplex(Rational(rng.below(201) - 100, 17), Rational(rng.below(201) - 100, 23)))
        var y = QQi(ExactComplex(Rational(rng.below(201) - 100, 31), Rational(rng.below(201) - 100, 7)))
        var z = QQi(ExactComplex(Rational(rng.below(201) - 100, 11), Rational(rng.below(201) - 100, 13)))
        c.check(x.add(zero, work) == x and x.mul(one, work) == x, "field identities")
        c.check(x.sub(x, work).is_zero() and x.add(x.neg(work), work).is_zero(), "additive inverses")
        c.check(x.add(y, work) == y.add(x, work) and x.mul(y, work) == y.mul(x, work), "commutativity")
        c.check(x.add(y, work).add(z, work) == x.add(y.add(z, work), work), "addition associativity")
        c.check(x.mul(y, work).mul(z, work) == x.mul(y.mul(z, work), work), "multiplication associativity")
        c.check(x.mul(y.add(z, work), work) == x.mul(y, work).add(x.mul(z, work), work), "distributivity")
        if not x.is_zero():
            c.check(x.mul(x.inverse(work), work) == one, "multiplicative inverse")

    var big = Rational(Integer(1) << 2048, Integer(3))
    var large = QQi(ExactComplex(big, Rational(7, 11)))
    c.check(large.size() > a.size() and large.mul(large.inverse(work), work) == one, "arbitrary precision inverse")
    for op in range(5):
        var measured = Budget()
        var expected = scalar_op(op, a, b, measured)
        var cost = UInt64.MAX - measured.remaining()
        var exact = Budget(cost)
        c.check(scalar_op(op, a, b, exact) == expected and exact.remaining() == 0, "exact work allowance")
        var allowances: List[UInt64] = [0, cost - 1]
        for allowance in allowances:
            var tiny = Budget(allowance)
            var original_a = a.copy()
            var original_b = b.copy()
            try:
                _ = scalar_op(op, a, b, tiny)
                c.check(False, "coefficient operation ignored budget")
            except e:
                c.check(e.kind == PolyErrorKind.LIMIT, "typed coefficient budget failure")
            c.check(a == original_a and b == original_b, "failure preserves borrowed coefficients")
    var empty = Budget(0)
    try:
        _ = zero.inverse(empty)
        c.check(False, "zero inverse accepted")
    except e:
        c.check(e.kind == PolyErrorKind.DOMAIN and empty.remaining() == 0, "zero inverse DOMAIN before arithmetic")


def bases(mut c: Checker) raises:
    var work = Budget()
    var i = QQi(ExactComplex(0, 1))
    for kind in range(1, 4):
        var order = MonomialOrder(UInt8(kind), [0, 1], work)
        var x = SparsePoly[QQi].variable(0, order, NoContext(), work)
        var y = SparsePoly[QQi].variable(1, order, NoContext(), work)
        var one = SparsePoly[QQi].constant(QQi.one(NoContext()), order, work)
        var zero = SparsePoly[QQi].zero(order, NoContext(), work)
        var imag = SparsePoly[QQi].constant(i, order, work)
        var f = x.pow(2, work).sub(y.scale(i, work), work)
        var g = x.mul(y, work).sub(one, work)
        # x^2 = i y and xy = 1 imply x = i y^2, y^3 = -i, x^3 = i.
        var expected = List[SparsePoly[QQi]]()
        if kind == 1:
            expected.append(x.sub(y.pow(2, work).scale(i, work), work))
            expected.append(y.pow(3, work).add(imag, work))
        else:
            expected.append(f.copy())
            expected.append(g.copy())
            expected.append(y.pow(2, work).add(x.scale(i, work), work))
        var generators: List[SparsePoly[QQi]] = [f.scale(QQi(ExactComplex(2, 3)), work), g.scale(QQi(ExactComplex(Rational(-1, 2), Rational(1, 3))), work), zero.copy(), f.copy()]
        var basis = groebner(generators, order, work, track_cofactors=True)
        c.check(basis.basis == expected, "hand-derived complex reduced basis")
        verify_basis(basis, generators, c, work)
        var repeated = groebner(generators, order, work, track_cofactors=True)
        c.check(repeated.basis == basis.basis and repeated.cofactors.value() == basis.cofactors.value(), "deterministic basis and witnesses")
        var reverse: List[SparsePoly[QQi]] = [f.copy(), zero.copy(), generators[1].copy(), generators[0].copy()]
        var reversed_basis = groebner(reverse, order, work)
        c.check(reversed_basis.basis == expected and not reversed_basis.cofactors, "generator order and cofactor option")
        for k in range(len(basis.basis)):
            var proof = MembershipCertificate[QQi](basis.cofactors.value()[k].copy())
            c.check(check_membership(basis.basis[k], generators, proof, work), "complex basis cofactor identity")
        var cube = x.pow(3, work)
        var normal = reduce(cube, basis, work)
        c.check(normal.remainder == imag and check_list_division(cube, basis.basis, normal, work), "exact complex normal form and division certificate")
        var target = cube.sub(imag, work)
        var member = membership(target, generators, order, work)
        c.check(Bool(member) and check_membership(target, generators, member.value(), work), "complex ideal membership")
        var corrupt = member.take()
        corrupt.cofactors[0] = corrupt.cofactors[0].add(imag, work)
        c.check(not check_membership(target, generators, corrupt, work), "imaginary cofactor corruption rejected")
        _ = corrupt.cofactors.pop()
        c.check(not check_membership(target, generators, corrupt, work), "truncated complex witness rejected")
        c.check(not membership(one, generators, order, work), "one outside proper complex ideal")
        var impossible: List[SparsePoly[QQi]] = [x.sub(imag, work), x.add(imag, work)]
        var contradiction = membership(one, impossible, order, work)
        c.check(Bool(contradiction) and check_membership(one, impossible, contradiction.value(), work), "complex infeasibility certificate")
        c.check(groebner(List[SparsePoly[QQi]](), order, work).basis == List[SparsePoly[QQi]](), "empty complex basis")
        var zeros: List[SparsePoly[QQi]] = [zero.copy()]
        c.check(groebner(zeros, order, work).basis == List[SparsePoly[QQi]](), "all-zero complex basis")
        var zero_proof = membership(zero, List[SparsePoly[QQi]](), order, work)
        c.check(Bool(zero_proof) and check_membership(zero, List[SparsePoly[QQi]](), zero_proof.value(), work), "empty complex membership certificate")
        var before = generators.copy()
        var tiny = Budget(100)
        try:
            _ = groebner(generators, order, tiny, track_cofactors=True)
            c.check(False, "complex basis ignored budget")
        except e:
            c.check(e.kind == PolyErrorKind.LIMIT and generators == before, "basis LIMIT preserves generators")
        try:
            _ = groebner(generators, MonomialOrder(), work)
            c.check(False, "complex basis accepted KERNEL order")
        except e:
            c.check(e.kind == PolyErrorKind.DOMAIN, "complex basis rejects KERNEL order")


def dense(mut c: Checker) raises:
    var work = Budget()
    var i = QQi(ExactComplex(0, 1))
    var one = QQi.one(NoContext())
    var left = DensePoly[QQi].from_coefficients(0, NoContext(), [i.neg(work), one.copy()], work)
    var right = DensePoly[QQi].from_coefficients(0, NoContext(), [i.copy(), one.copy()], work)
    var product = left.mul(right, work)
    var expected = DensePoly[QQi].from_coefficients(0, NoContext(), [one.copy(), QQi.zero(NoContext()), one.copy()], work)
    c.check(product == expected, "(x-i)(x+i) = x^2+1")
    product.check_invariants(work)
    var division = divmod(product, left, work)
    c.check(division.quotient == right and division.remainder.is_zero() and check_division(product, left, division, work), "dense complex exact division")
    c.check(evaluate_dense(product, i, work).is_zero(), "exact complex root evaluation")
    var order = MonomialOrder(MonomialOrder.LEX, [0], work)
    c.check(DensePoly[QQi].from_sparse(product.to_sparse(order, work), work) == expected, "dense and sparse complex round trip")


def main() raises:
    var c = Checker("exact_complex")
    coefficients(c)
    bases(c)
    dense(c)
    c.finish()
