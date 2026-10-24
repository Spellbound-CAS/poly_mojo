"""Exact LLL invariants and independently planted lattice recombination cases."""
from apn_mojo import Integer, Rational
from polynomials import Budget, PolyErrorKind, ZZ, QQ, NoContext, DensePoly, Modulus, factor_mod_p, factor_integer, check_factorization
from polynomials._lattice import _lll, _lll_knapsack, _lll_knapsack_integer, _q, _orthogonal, _exact_norms, _modular_norms
from polynomials._lattice_factor import _van_hoeij, _cld_bound, _modular_partition, _cld_values, _lattice_precision, _trace_precision, _cld_precision, _uses_traces
from polynomials._dense_ops import _derivative
from polynomials._hensel import _hensel, _hensel_start, _hensel_extend, _divmod_monic, _truncate
from polynomials._gcd import _mod_image
from polynomials._dense_ops import _monic
from polynomials.factor import _bound, _large_recombine
from support.corpus import parse_polynomial, as_zz
from polynomials import MonomialOrder
from support.check import Checker


def polynomial(a: Int, b: Int, d: Int, mut budget: Budget) raises -> DensePoly[ZZ]:
    return DensePoly[ZZ].from_coefficients(0, NoContext(), [ZZ(Integer(a)), ZZ(Integer(b)), ZZ(Integer(d))], budget)


def main() raises:
    var work = Budget()
    var c = Checker("lattice")
    var input: List[List[ZZ]] = [[ZZ(Integer(1)), ZZ(Integer(1)), ZZ(Integer(1))], [ZZ(Integer(-1)), ZZ(Integer(0)), ZZ(Integer(2))], [ZZ(Integer(3)), ZZ(Integer(5)), ZZ(Integer(6))]]
    var tracked_work = Budget()
    var result = _lll(input, tracked_work)
    var untracked_work = Budget()
    var untracked = _lll[False](input, untracked_work)
    c.check(len(untracked.transform) == 0, "untracked LLL omits transformation storage")
    c.check(untracked.basis == result.basis and untracked.mu == result.mu and untracked.norms == result.norms, "transformation tracking does not affect exact reduction")
    c.check(untracked_work.remaining() > tracked_work.remaining(), "untracked reduction avoids transformation work")
    var saved = input.copy()
    for limit in [UInt64(0), UInt64(50), UInt64(500)]:
        var tiny = Budget(limit)
        try:
            var limited = _lll[False](input, tiny)
            c.check(limited.basis == result.basis, "limited LLL returns the complete basis")
        except e:
            c.check(e.kind == PolyErrorKind.LIMIT, "untracked LLL fails closed on exhaustion")
        c.check(input == saved, "LLL exhaustion preserves its input matrix")
    var dependent: List[List[ZZ]] = [[ZZ(Integer(1)), ZZ(Integer(2))], [ZZ(Integer(2)), ZZ(Integer(4))]]
    try:
        _ = _lll[False](dependent, work)
        c.check(False, "dependent untracked lattice accepted")
    except e:
        c.check(e.kind == PolyErrorKind.DOMAIN, "untracked LLL rejects a rank-deficient matrix")
    # Floating steering, exact answers: the knapsack reduction returns a basis
    # whose norms equal an independent exact Gram-Schmidt pass.
    var knapsack: List[List[ZZ]] = [[ZZ(Integer(1)), ZZ(Integer(0)), ZZ(Integer(0)), ZZ(Integer(123456789))], [ZZ(Integer(0)), ZZ(Integer(1)), ZZ(Integer(0)), ZZ(Integer(987654321))], [ZZ(Integer(0)), ZZ(Integer(0)), ZZ(Integer(1)), ZZ(Integer(555555555))], [ZZ(Integer(0)), ZZ(Integer(0)), ZZ(Integer(0)), ZZ(Integer(1000000007))]]
    var steered = _lll_knapsack(knapsack, work)
    var exact_check = _orthogonal(steered.basis, work)
    var same_norms = True
    var volume = QQ.one(NoContext())
    for i in range(4):
        same_norms = same_norms and steered.norms[i] == exact_check.norms[i]
        volume = volume.mul(steered.norms[i], work)
    c.check(same_norms, "knapsack LLL norms are exact")
    c.check(volume == _q(ZZ(Integer(1000000007).copy()), work).mul(_q(ZZ(Integer(1000000007)), work), work), "knapsack LLL keeps the lattice determinant")
    var recomputed = _orthogonal(result.basis, work)
    for i in range(3):
        c.check(result.norms[i] == recomputed.norms[i], "incremental exact Gram-Schmidt norm")
        for j in range(i):
            c.check(result.mu[i][j] == recomputed.mu[i][j], "incremental exact Gram-Schmidt coefficient")
    for i in range(3):
        for j in range(3):
            var total = ZZ.zero(NoContext())
            for k in range(3):
                total = total.add(result.transform[i][k].mul(input[k][j], work), work)
            c.check(total == result.basis[i][j], "B=U*A")
    ref u = result.transform
    var det = u[0][0].mul(u[1][1].mul(u[2][2], work).sub(u[1][2].mul(u[2][1], work), work), work).sub(u[0][1].mul(u[1][0].mul(u[2][2], work).sub(u[1][2].mul(u[2][0], work), work), work), work).add(u[0][2].mul(u[1][0].mul(u[2][1], work).sub(u[1][1].mul(u[2][0], work), work), work), work)
    c.check(det.value == Integer(1) or det.value == Integer(-1), "unimodular row transformation")
    var half = QQ(Rational(1, 2))
    var delta = QQ(Rational(3, 4))
    for i in range(3):
        for j in range(i):
            ref mu = result.mu[i][j]
            c.check(mu.value <= half.value and mu.value >= half.neg(work).value, "size-reduced Gram-Schmidt coefficient")
        if i:
            var lower = delta.sub(result.mu[i][i - 1].mul(result.mu[i][i - 1], work), work).mul(result.norms[i - 1], work)
            c.check(result.norms[i].value >= lower.value, "exact Lovasz condition")
    var ctx = Modulus(UInt64(97), work)
    for nonmonic in range(2):
        var left = polynomial(-2 + nonmonic, 0, 1 + nonmonic, work)
        var right = polynomial(-3, 0, 1, work)
        var f = left.mul(right, work)
        var modular = factor_mod_p(_monic(_mod_image(f, ctx, work), work), work, seed=0)
        var state = _hensel_start(modular, work)
        _hensel_extend(f, modular, ZZ(Integer(10000)), state, work)
        var old_modulus = state.modulus.copy()
        var cached_work = Budget()
        _hensel_extend(f, modular, ZZ(Integer(100000000)), state, cached_work)
        var word_work = Budget()
        var fresh = _hensel(f, modular, ZZ(Integer(100000000)), word_work)
        c.check(state.modulus == fresh[1] and state.modulus.value > old_modulus.value, "resumed lifting reaches fresh precision")
        var same = len(state.lifted) == len(fresh[0])
        for i in range(len(state.lifted)):
            same = same and state.lifted[i] == fresh[0][i]
        c.check(same, "resumed and fresh Hensel factors agree")
        # Against restarting the same big-integer tree lift; _hensel itself
        # may take the word-size path.
        var fresh_work = Budget()
        var restarted = _hensel_start(modular, fresh_work)
        _hensel_extend(f, modular, ZZ(Integer(100000000)), restarted, fresh_work)
        c.check(cached_work.remaining() > fresh_work.remaining(), "resumption saves correction and inverse work")
        var factors = _van_hoeij(f, modular, _bound(f, work), work)
        var product = polynomial(1, 0, 0, work)
        for factor in factors:
            product = product.mul(factor, work)
        c.check(product == f, "lattice factors reconstruct the input")
        c.check(len(factors) == 2, "expected irreducible factor count")
        if not nonmonic:
            c.check((factors[0] == left and factors[1] == right) or (factors[1] == left and factors[0] == right), "planted irrational quadratics recovered")
    # The CLD bound dominates every coefficient of g' (f/g) for true factors g.
    var g1 = DensePoly[ZZ].from_coefficients(0, NoContext(), [ZZ(Integer(-7)), ZZ(Integer(3)), ZZ(Integer(0)), ZZ(Integer(11))], work)
    var g2 = DensePoly[ZZ].from_coefficients(0, NoContext(), [ZZ(Integer(5)), ZZ(Integer(-2)), ZZ(Integer(9))], work)
    var g3 = DensePoly[ZZ].from_coefficients(0, NoContext(), [ZZ(Integer(1000)), ZZ(Integer(1)), ZZ(Integer(-40)), ZZ(Integer(0)), ZZ(Integer(13))], work)
    var planted = g1.mul(g2, work).mul(g3, work)
    var dominated = True
    for which in range(3):
        var g = g1.copy() if which == 0 else (g2.copy() if which == 1 else g3.copy())
        var h = g2.mul(g3, work) if which == 0 else (g1.mul(g3, work) if which == 1 else g1.mul(g2, work))
        var cld = _derivative(g, work).mul(h, work)
        for j in range(planted.degree()):
            var value = cld.coefficients[j].value.copy() if j < len(cld.coefficients) else Integer(0)
            if value.sign() < 0:
                value = -value
            dominated = dominated and not (value > _cld_bound(planted, j, work).value)
    c.check(dominated, "CLD bound dominates logarithmic derivative coefficients")
    # Series coefficients of g' (f/g) from both ends equal the division's,
    # for lifted factors of a nonmonic f.
    var planted_images = factor_mod_p(_monic(_mod_image(planted, Modulus(10007, work), work), work), work, seed=0)
    var planted_state = _hensel_start(planted_images, work)
    _hensel_extend(planted, planted_images, ZZ(Integer(1) << 200), planted_state, work)
    var every = List[Int]()
    for j in range(planted.degree()):
        every.append(j)
    var series_agree = len(planted_images) >= 2
    for g in planted_state.lifted:
        var quotient = _divmod_monic(planted, g, planted_state.modulus, work)[0].copy()
        var cld = _truncate(_derivative(g, work).mul(quotient, work), planted_state.modulus, work)
        var values = _cld_values(planted, g, every, planted_state.modulus, work)
        series_agree = series_agree and Bool(values)
        if values:
            for j in range(len(every)):
                var expected = cld.coefficients[j].copy() if j < len(cld.coefficients) else ZZ.zero(NoContext())
                series_agree = series_agree and values.value()[j] == expected
    c.check(series_agree, "CLD series coefficients equal the division's")
    # The first lattice precision is the smaller of the trace and CLD ones,
    # and it names the columns the lattice takes.
    var lifting_prime = ZZ(Integer(10007))
    var modular_count = len(planted_images)
    var trace_first = _trace_precision(planted, modular_count, lifting_prime, work)
    var cld_first = _cld_precision(planted, modular_count, lifting_prime, work)
    var chosen = _lattice_precision(planted, modular_count, lifting_prime, work).value().copy()
    var traces = _uses_traces(planted, modular_count, lifting_prime, work)
    var smaller = cld_first.copy() if trace_first.value > cld_first.value else trace_first.copy()
    c.check(chosen == smaller and traces == (chosen == trace_first and not (trace_first.value > cld_first.value)), "lattice precision takes the cheaper columns")
    var split = polynomial(1, 0, 0, work)
    for i in range(1, 21):
        split = split.mul(polynomial(-i, 1, 0, work), work)
    var public_result = factor_integer(split, work)
    c.check(len(public_result.factors) == 20 and check_factorization(split, public_result, work), "large public dispatch peels easy singleton groups")
    var cyclotomic = DensePoly[ZZ].from_sparse(as_zz(parse_polynomial("x0^20 - x0^18 + x0^16 - x0^14 + x0^12 - x0^10 + x0^8 - x0^6 + x0^4 - x0^2 + 1", MonomialOrder(), work), work), work)
    var completely_split = factor_mod_p(_monic(_mod_image(cyclotomic, Modulus(89, work), work), work), work, seed=0)
    var hard = _large_recombine(cyclotomic, completely_split, work)
    c.check(len(completely_split) == 20 and len(hard) == 1 and hard[0] == cyclotomic, "Phi44 irreducibility survives twenty-way modular splitting")
    # The Int128 knapsack reduction takes the same steps as the Integer one.
    var wide_lattice = List[List[ZZ]]()
    var lattice_seed = Integer(987654321)
    var lattice_top = Integer(1) << 80
    for i in range(30):
        var row = List[ZZ]()
        for j in range(30):
            row.append(ZZ(Integer(1) if i == j else Integer(0)))
        lattice_seed = (lattice_seed * Integer(6364136223846793005) + Integer(1442695040888963407)) % lattice_top
        row.append(ZZ(lattice_seed - (lattice_top >> 1)))
        wide_lattice.append(row^)
    var lattice_last = List[ZZ]()
    for _ in range(30):
        lattice_last.append(ZZ(Integer(0)))
    lattice_last.append(ZZ(lattice_top))
    wide_lattice.append(lattice_last^)
    var narrow = _lll_knapsack(wide_lattice, work)
    var wide = _lll_knapsack_integer(wide_lattice, work)
    var agree = len(narrow.basis) == len(wide.basis)
    for i in range(len(narrow.basis)):
        for j in range(len(narrow.basis[i])):
            agree = agree and narrow.basis[i][j] == wide.basis[i][j]
        agree = agree and narrow.norms[i].value == wide.norms[i].value
    c.check(agree, "Int128 knapsack reduction equals the Integer one")
    # Multimodular norms equal the fraction-free ones, on the reduced basis
    # and on the original one with its 2^80 entries.
    var modular_agree = True
    for lattice in [narrow.basis.copy(), wide_lattice.copy()]:
        var exact = _exact_norms(lattice, work)
        var modular = _modular_norms(lattice, work)
        modular_agree = modular_agree and Bool(modular) and len(modular.value()) == len(exact)
        if modular:
            for i in range(len(exact)):
                modular_agree = modular_agree and modular.value()[i] == exact[i]
    c.check(modular_agree, "multimodular Gram-Schmidt norms equal the exact recurrence")
    # The modular partition: groups come back in pivot order, a
    # non-partition gives none, and dependent rows reduce away.
    var groups: List[List[ZZ]] = [[ZZ(Integer(1)), ZZ(Integer(1)), ZZ(Integer(0))], [ZZ(Integer(0)), ZZ(Integer(0)), ZZ(Integer(1))]]
    var mixed: List[List[ZZ]] = [[ZZ(Integer(1)), ZZ(Integer(2)), ZZ(Integer(0))], [ZZ(Integer(0)), ZZ(Integer(1)), ZZ(Integer(1))]]
    var deficient: List[List[ZZ]] = [[ZZ(Integer(1)), ZZ(Integer(2)), ZZ(Integer(3))], [ZZ(Integer(2)), ZZ(Integer(4)), ZZ(Integer(6))], [ZZ(Integer(0)), ZZ(Integer(1)), ZZ(Integer(1))], [ZZ(Integer(1)), ZZ(Integer(0)), ZZ(Integer(5))]]
    var found = _modular_partition(groups, 3, work)
    var singletons = _modular_partition(deficient, 3, work)
    var one = QQ(Rational(1))
    var zero = QQ(Rational(0))
    var exact = len(found) == 2 and found[0][0] == one and found[0][1] == one and found[0][2] == zero and found[1][0] == zero and found[1][1] == zero and found[1][2] == one
    exact = exact and len(singletons) == 3 and singletons[0][0] == one and singletons[1][1] == one and singletons[2][2] == one and singletons[0][1] == zero
    c.check(exact and len(_modular_partition(mixed, 3, work)) == 0, "modular partition groups, rejection and dependent rows")
    c.finish()
