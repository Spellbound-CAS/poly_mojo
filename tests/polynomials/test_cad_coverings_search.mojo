"""Coverings search: witnesses, conflicts, degenerate inputs, budgets, and agreement with the certified covering."""
from polynomials import Budget, QQ, SparsePoly, MonomialOrder, PolyErrorKind, PolyError
from polynomials import coverings_solve, check_coverings, CADCoveringsCertificate, CADFiber, CoveringBound
from polynomials.cylindrical_covering import cylindrical_covering, check_cylindrical_covering
from polynomials.algebraic_specialization import algebraic_evaluate_sparse
from polynomials.algebraic import algebraic_sign
from polynomials.algebraic_types import RealAlgebraic
from support.check import Checker
from apn_mojo import Rational
from polynomials import real_rational, real_root, DensePoly, NoContext
from support.corpus import parse_polynomial


def polys(texts: List[String], mut work: Budget) raises -> List[SparsePoly[QQ]]:
    var order = MonomialOrder()
    var result = List[SparsePoly[QQ]]()
    for t in texts:
        result.append(parse_polynomial(t, order, work))
    return result^


def axes(n: Int) -> List[Int64]:
    var result = List[Int64]()
    for i in range(n):
        result.append(Int64(i))
    return result^


def satisfies(inputs: List[SparsePoly[QQ]], masks: List[UInt8], point: List[RealAlgebraic], mut work: Budget) raises -> Bool:
    """Substitute the point exactly and check every sign against its mask."""
    var assignment = List[Tuple[Int64, RealAlgebraic]]()
    for i in range(len(point)):
        assignment.append((Int64(i), point[i].copy()))
    for i in range(len(inputs)):
        var value = algebraic_evaluate_sparse(inputs[i], assignment, work)
        var sign = algebraic_sign(value.value, work).relation
        if masks[i] & (UInt8(1) << UInt8(sign + 1)) == 0:
            return False
    return True


def verdict(texts: List[String], masks: List[UInt8], n: Int, mut c: Checker, what: String, expected: Bool, mut work: Budget) raises:
    var inputs = polys(texts, work)
    var outcome = coverings_solve(inputs, masks, axes(n), work)
    c.check(outcome.satisfiable() == expected, what)
    c.check(check_coverings(inputs, masks, outcome, work), what + ": the certificate checks")
    if outcome.satisfiable():
        c.check(len(outcome.witness) == n and satisfies(inputs, masks, outcome.witness, work), what + ": the point satisfies every constraint")
    var mccallum = coverings_solve(inputs, masks, axes(n), work, lazard=False)
    c.check(mccallum.satisfiable() == expected and check_coverings(inputs, masks, mccallum, work), what + ": McCallum characterizations agree and check")


struct Lcg(Movable):
    var state: UInt64

    def __init__(out self, seed: UInt64):
        self.state = seed

    def next(mut self, bound: Int) -> Int:
        self.state = self.state * 6364136223846793005 + 1442695040888963407
        return Int((self.state >> 33) % UInt64(bound))


def random_poly(mut rng: Lcg, n: Int, degree: Int) -> String:
    """Up to four random terms plus a constant, with distinct monomials and nonzero coefficients."""
    var monomials = List[String]()
    var coefficients = List[Int]()
    var terms = 1 + rng.next(4)
    for t in range(terms + 1):
        var c = rng.next(7) - 3
        var monomial = String()
        if t < terms:
            var remaining = 1 + rng.next(degree)
            for v in range(n):
                var e = rng.next(remaining + 1) if v + 1 < n else remaining
                remaining -= e
                if e:
                    monomial.write("*" if monomial.byte_length() else "", "x", v, "^", e)
        var found = -1
        for i in range(len(monomials)):
            if monomials[i] == monomial:
                found = i
        if found >= 0:
            coefficients[found] += c
        else:
            monomials.append(monomial)
            coefficients.append(c)
    var text = String()
    for i in range(len(monomials)):
        if coefficients[i] == 0:
            continue
        if text.byte_length():
            text.write(" + ")
        if monomials[i].byte_length() == 0:
            text.write(coefficients[i])
        else:
            text.write(coefficients[i], "*", monomials[i])
    if text.byte_length() == 0:
        return String("x0 - 1")
    return text^


def main() raises:
    var work = Budget()
    var c = Checker("coverings search")
    verdict(["x0^2 + x1^2 - 1"], [1], 2, c, "disk interior", True, work)
    verdict(["x0^2 + x1^2 + 1"], [3], 2, c, "positive definite is never nonpositive", False, work)
    verdict(["x0^2 + x1^2 - 1", "x0 - 1"], [2, 4], 2, c, "circle beyond its tangent", False, work)
    verdict(["x0^2 + x1^2 - 1", "x0 - 1"], [2, 6], 2, c, "circle touches its tangent", True, work)
    verdict(["x0^2 - 2"], [2], 1, c, "irrational section", True, work)
    verdict(["x0^2 - 2", "x1^2 - x0", "x0*x1 - 1"], [2, 2, 4], 2, c, "tower of sections", True, work)
    verdict(["x0^2 - 2", "x1^2 - 3", "x0*x1 + x1 - 5"], [2, 2, 4], 2, c, "two irrational coordinates, no product large enough", False, work)
    verdict(["x0 - 1", "x0", "x0^2 - 4"], [1, 4, 6], 1, c, "overlapping conflicts", False, work)
    verdict(["x0^2 + 1", "x1^2 + x2^2 - 1"], [1, 2], 3, c, "lower-dimensional conflict", False, work)
    verdict(["x0", "x0*x1"], [2, 4], 2, c, "Lazard nullification over a section", False, work)
    verdict(["x0", "x0*x1"], [2, 2], 2, c, "nullified input is zero above the section", True, work)
    verdict(["x0^2 - 2", "x0^2*x1 - 2*x1"], [2, 4], 2, c, "nullification over an irrational section", False, work)
    verdict(["x0^2 - 2", "x0^2*x1 - 2*x1"], [2, 2], 2, c, "nullified at an irrational section, zero above it", True, work)
    verdict(["x0^2 - 2", "x0^2*x1^2 - 2*x1^2 + x1 - 1"], [2, 1], 2, c, "leading coefficient vanishing at an irrational section", True, work)
    verdict(["x0^3 - 2", "x1^2 - x0", "x1^3 - 2*x1"], [2, 2, 2], 2, c, "no common root over a cubic field", False, work)
    verdict(["x0*x1 - 1", "x0^2 + x1^2 - 2"], [2, 2], 2, c, "hyperbola meets circle at two points", True, work)
    verdict(["x0*x1 - 1", "x0^2 + x1^2 - 1"], [2, 1], 2, c, "hyperbola outside the unit disk", False, work)
    verdict(["x0^4*x1^2 + x0^2*x1^4 - 3*x0^2*x1^2 + 1"], [1], 2, c, "Motzkin is nonnegative", False, work)
    verdict(["x0^4*x1^2 + x0^2*x1^4 - 3*x0^2*x1^2 + 1"], [3], 2, c, "Motzkin vanishes at (1, 1)", True, work)
    verdict(["x0^2 + x1^2 + x2^2 - 1", "x0 + x1 + x2 - 3"], [2, 2], 3, c, "sphere misses a far plane", False, work)
    verdict(["x0^2 + x1^2 + x2^2 - 1", "x0 + x1 + x2 - 1"], [2, 2], 3, c, "sphere meets a near plane", True, work)
    # Degenerate inputs.
    verdict(["1"], [2], 2, c, "constant contradiction", False, work)
    verdict(["1"], [4], 0, c, "constant truth in zero variables", True, work)
    verdict(["x0"], [0], 1, c, "an empty mask excludes the line", False, work)
    verdict(["x0"], [7], 1, c, "a full mask allows the line", True, work)
    verdict(["0", "x0 - 3"], [2, 2], 1, c, "a zero input is zero", True, work)
    verdict(["0"], [5], 1, c, "a zero input is never nonzero", False, work)
    verdict(["x1"], [4], 2, c, "an unused base variable", True, work)
    var empty = coverings_solve(List[SparsePoly[QQ]](), List[UInt8](), List[Int64](), work)
    c.check(empty.satisfiable() and len(empty.witness) == 0, "empty conjunction in zero variables")
    var line = polys(["x0"], work)
    var plane = polys(["x0*x1"], work)
    var rejected = False
    try:
        _ = coverings_solve(line, [UInt8(8)], axes(1), work)
    except e:
        rejected = e.kind == PolyErrorKind.DOMAIN
    c.check(rejected, "invalid mask rejected")
    rejected = False
    try:
        _ = coverings_solve(plane, [UInt8(4)], axes(1), work)
    except e:
        rejected = e.kind == PolyErrorKind.DOMAIN
    c.check(rejected, "variables must cover the inputs")
    # Budgets: exhaustion raises and never reads as a verdict; charges are deterministic.
    var system = polys(["x0^2 + x1^2 + x2^2 - 1", "x0 + x1 + x2 - 3"], work)
    var masks: List[UInt8] = [2, 2]
    var exhausted = False
    try:
        var tiny = Budget(200)
        _ = coverings_solve(system, masks, axes(3), tiny)
    except e:
        exhausted = e.kind == PolyErrorKind.LIMIT
    c.check(exhausted, "a small budget raises LIMIT")
    var first = Budget()
    _ = coverings_solve(system, masks, axes(3), first)
    var second = Budget()
    _ = coverings_solve(system, masks, axes(3), second)
    c.check(first.remaining() == second.remaining() and first.remaining() < UInt64.MAX, "budget charges repeat exactly")
    var used = UInt64.MAX - first.remaining()
    var exact = Budget(used)
    var outcome = coverings_solve(system, masks, axes(3), exact)
    c.check(not outcome.satisfiable() and exact.remaining() == 0, "the measured charge suffices exactly")
    exhausted = False
    try:
        var short = Budget(used - 1)
        _ = coverings_solve(system, masks, axes(3), short)
    except e:
        exhausted = e.kind == PolyErrorKind.LIMIT
    c.check(exhausted, "one step less raises LIMIT")
    # Corrupted coverings are rejected.
    var sphere = polys(["x0^2 + x1^2 + x2^2 - 1", "x0 + x1 + x2 - 3"], work)
    var sphere_masks: List[UInt8] = [2, 2]
    var proof = coverings_solve(sphere, sphere_masks, axes(3), work)
    c.check(not proof.satisfiable() and len(proof.nodes) > 1 and check_coverings(sphere, sphere_masks, proof, work), "sphere and far plane: covering checks")
    var bad = proof.copy()
    _ = bad.nodes[0].intervals.pop()
    c.check(not check_coverings(sphere, sphere_masks, bad, work), "dropped interval leaves a gap")
    bad = proof.copy()
    _ = bad.nodes[0].intervals.pop(0)
    c.check(not check_coverings(sphere, sphere_masks, bad, work), "dropped first interval leaves minus infinity uncovered")
    bad = proof.copy()
    bad.projections.clear()
    c.check(not check_coverings(sphere, sphere_masks, bad, work), "missing projection factorizations rejected")
    bad = proof.copy()
    for k in range(len(bad.nodes[0].intervals)):
        bad.nodes[0].intervals[k].polynomials.clear()
    c.check(not check_coverings(sphere, sphere_masks, bad, work), "emptied characterizations rejected")
    bad = proof.copy()
    for k in range(len(bad.nodes[0].intervals)):
        if bad.nodes[0].intervals[k].lower:
            bad.nodes[0].intervals[k].lower.value().root += 1
            break
    c.check(not check_coverings(sphere, sphere_masks, bad, work), "shifted bound root rejected")
    bad = proof.copy()
    for k in range(len(bad.nodes)):
        for j in range(len(bad.nodes[k].intervals)):
            if bad.nodes[k].intervals[j].constraint >= 0:
                bad.nodes[k].intervals[j].constraint = 1 - bad.nodes[k].intervals[j].constraint
    c.check(not check_coverings(sphere, sphere_masks, bad, work), "swapped conflict constraints rejected")
    bad = proof.copy()
    bad.nodes.clear()
    c.check(not check_coverings(sphere, sphere_masks, bad, work), "an empty tree is not a witness")
    var disk = polys(["x0^2 + x1^2 - 1"], work)
    var disk_masks: List[UInt8] = [1]
    var inside = coverings_solve(disk, disk_masks, axes(2), work)
    c.check(inside.satisfiable() and check_coverings(disk, disk_masks, inside, work), "disk witness checks")
    var moved = inside.copy()
    moved.witness[0] = inside.witness[0].copy()
    moved.signs[0].relation = 1
    c.check(not check_coverings(disk, disk_masks, moved, work), "false witness sign rejected")
    c.check(not check_coverings(disk, [UInt8(4)], inside, work), "witness checked against the masks")
    # A fixed prefix: the circle above x0 = 1/2 and above x0 = sqrt(2).
    var fiber = CADFiber([real_rational(QQ(Rational(1, 2)), work)], work)
    var above = fiber.search(disk, disk_masks, axes(2), work, lazard=True)
    c.check(above.satisfiable() and fiber.check(disk, disk_masks, above, work), "fiber above 1/2 meets the disk")
    var root2 = real_root(DensePoly[QQ].from_coefficients(0, NoContext(), [QQ(Rational(-2)), QQ(Rational(0)), QQ(Rational(1))], work), QQ(Rational(1)), QQ(Rational(2)), work)
    var outside = CADFiber([root2.copy()], work)
    var missed = outside.search(disk, disk_masks, axes(2), work, lazard=True)
    c.check(not missed.satisfiable() and outside.check(disk, disk_masks, missed, work), "fiber above sqrt(2) misses the disk")
    c.check(not fiber.check(disk, disk_masks, missed, work), "a fiber certificate does not check on another fiber")
    # Agreement with the certified covering on deterministic random systems.
    var rng = Lcg(20261021)
    var agreed = 0
    var satisfiable = 0
    var choices: List[UInt8] = [1, 2, 2, 3, 4, 5, 6]
    var trials = 120
    for trial in range(trials):
        var n = 2 if trial < 80 else 3
        var count = 1 + rng.next(4)
        var texts = List[String]()
        var trial_masks = List[UInt8]()
        for _ in range(count):
            texts.append(random_poly(rng, n, 2 if n == 3 else 3))
            trial_masks.append(choices[rng.next(len(choices))])
        var inputs = polys(texts, work)
        var mine = coverings_solve(inputs, trial_masks, axes(n), work)
        var reference = cylindrical_covering(inputs, trial_masks, axes(n), work, lazard=True)
        var same = mine.satisfiable() == (reference.witness_cell >= 0)
        if not same:
            var shown = String()
            for t in texts:
                shown.write(t, "; ")
            print("disagreement:", shown, "masks", trial_masks[0])
        c.check(same, "verdict agrees with the certified covering")
        c.check(check_cylindrical_covering(inputs, trial_masks, reference, work), "reference certificate checks")
        c.check(check_coverings(inputs, trial_masks, mine, work), "random certificate checks")
        if mine.satisfiable():
            satisfiable += 1
            c.check(satisfies(inputs, trial_masks, mine.witness, work), "random witness satisfies every constraint")
        if same:
            agreed += 1
    print("random systems:", satisfiable, "satisfiable of", agreed, "agreeing")
    c.check(satisfiable >= 10 and trials - satisfiable >= 10, "random systems mix both verdicts")
    c.finish()
