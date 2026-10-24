"""Compare private final stacks with the existing certified CAD and exact samples."""
from apn_mojo import Rational
from polynomials import Budget, DensePoly, QQ, NoContext, SparsePoly, MonomialOrder, cad, check_cad
from polynomials._cad_final import _verified_stack, _check_final_stack
from polynomials._cad_final_signs import _SignStats
from polynomials._root_ops import _sign_at
from support.check import Checker


def poly(values: List[Int], mut work: Budget) raises -> DensePoly[QQ]:
    var coefficients = List[QQ]()
    for value in values:
        coefficients.append(QQ(Rational(value)))
    return DensePoly[QQ].from_coefficients(0, NoContext(), coefficients^, work)


def reference(mut c: Checker, mut work: Budget) raises:
    var first = poly([-1, 0, 1], work)
    var second = poly([-2, 0, 1], work)
    var repeated = poly([1, -2, 1], work)
    var inputs: List[DensePoly[QQ]] = [first.copy(), second.copy(), repeated.copy(), first.mul(second, work)]
    var proposed = _verified_stack(inputs, work)
    var sparse = List[SparsePoly[QQ]]()
    for f in inputs:
        sparse.append(f.to_sparse(MonomialOrder(), work))
    for mode in range(2):
        var baseline = cad(sparse, [Int64(0)], work, lazard=Bool(mode))
        c.check(check_cad(sparse, baseline, work), "baseline CAD certificate checks")
        c.check(len(baseline.cells) == len(proposed.signs) + 1, "same complete ordered cell decomposition")
        for i in range(len(proposed.signs)):
            ref cell = baseline.cells[i + 1]
            c.check(cell.position == i, "same section/sector ordering")
            for j in range(len(inputs)):
                c.check(cell.signs[j].relation == proposed.signs[i][j], "same input sign on every cell")


def generated(mut c: Checker, mut work: Budget) raises:
    # Fixed-seed families: shared rational factors, repeated factors,
    # no-real-root factors, signed scalars, and coefficient denominators.
    var seed = UInt64(0x516CAD)
    for family in range(32):
        seed = seed * UInt64(6364136223846793005) + 1
        var a = Int((seed >> 32) % 9) - 4
        seed = seed * UInt64(6364136223846793005) + 1
        var b = Int((seed >> 32) % 9) - 4
        var left = poly([-a, 1], work)
        var right = poly([-b, 1], work)
        var positive = poly([family + 1, 0, 1], work)
        var repeated = left.mul(left, work).mul(right, work)
        var signed = right.mul(positive, work).scale(QQ(Rational(-3, family + 2)), work)
        var inputs: List[DensePoly[QQ]] = [left.copy(), right.copy(), repeated.copy(), signed.copy()]
        var stack = _verified_stack(inputs, work)
        var scalar = _SignStats()
        c.check(_check_final_stack(inputs, stack, scalar, work, vectorized=False), "generated scalar/SIMD agreement")
        c.check(len(stack.roots) == (1 if a == b else 2), "known rational root count")
        for i in range(len(stack.sectors)):
            for j in range(len(inputs)):
                c.check(stack.signs[2 * i][j] == _sign_at(inputs[j], stack.sectors[i], work), "generated direct sector oracle")
        for i in range(len(stack.roots)):
            var x = min(a, b) if i == 0 else max(a, b)
            var sample = QQ(Rational(x))
            for j in range(len(inputs)):
                c.check(stack.signs[2 * i + 1][j] == _sign_at(inputs[j], sample, work), "generated direct exact-section oracle")
        var bad = stack.copy()
        bad.signs[0][0] = -bad.signs[0][0]
        var vector = _SignStats()
        c.check(not _check_final_stack(inputs, bad, vector, work), "generated forged result rejected")

        # No separately supplied factors may be needed to discover distinct
        # multiplicities within a single input polynomial.
        inputs = [repeated.mul(right, work).mul(right, work).mul(positive, work)]
        stack = _verified_stack(inputs, work)
        c.check(len(stack.roots) == (1 if a == b else 2), "mixed multiplicities in one input")
        for i in range(len(stack.sectors)):
            c.check(stack.signs[2 * i][0] == _sign_at(inputs[0], stack.sectors[i], work), "single-input sector oracle")
        for i in range(len(stack.roots)):
            c.check(stack.signs[2 * i + 1][0] == 0, "single-input section oracle")


def main() raises:
    var c = Checker("final CAD differential")
    var work = Budget()
    reference(c, work)
    generated(c, work)
    c.finish()
