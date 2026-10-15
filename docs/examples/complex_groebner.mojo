"""Exact complex coefficients, Gröbner reduction, and membership witnesses."""
from apn_mojo import ExactComplex
from polynomials import Budget, QQi, NoContext, SparsePoly, MonomialOrder, groebner, reduce, membership, check_membership


def main() raises:
    var work = Budget()
    var order = MonomialOrder(MonomialOrder.LEX, [0, 1], work)
    var x = SparsePoly[QQi].variable(0, order, NoContext(), work)
    var y = SparsePoly[QQi].variable(1, order, NoContext(), work)
    var i = QQi(ExactComplex(0, 1))
    var imag = SparsePoly[QQi].constant(i, order, work)
    var one = SparsePoly[QQi].constant(QQi.one(NoContext()), order, work)
    var generators: List[SparsePoly[QQi]] = [x.pow(2, work).sub(y.scale(i, work), work), x.mul(y, work).sub(one, work)]
    var basis = groebner(generators, order, work)
    var expected: List[SparsePoly[QQi]] = [x.sub(y.pow(2, work).scale(i, work), work), y.pow(3, work).add(imag, work)]
    if basis.basis != expected:
        raise Error("unexpected exact complex basis")
    print("Exact complex basis elements:", len(basis.basis))
    var cube = x.pow(3, work)
    if reduce(cube, basis, work).remainder != imag:
        raise Error("x^3 should reduce to i")
    print("x^3 reduces exactly to i")
    var target = cube.sub(imag, work)
    var proof = membership(target, generators, order, work)
    if not proof or not check_membership(target, generators, proof.value(), work):
        raise Error("complex membership witness did not check")
    print("Complex membership certificate checked")
