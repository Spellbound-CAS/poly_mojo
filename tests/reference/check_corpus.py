"""Independent checks of the reference corpus, tests/reference/corpus.txt.

The corpus's expected values come from SymPy (tests/reference/gen_corpus.py).
This script re-checks them with plain Python integers and fractions, sharing
no code with the package or with the generator, so that a case enters the
package's tests only if two independent computations agree on it. It also
checks the format (canonical terms, terms in decreasing order, the
normalizations each operation documents) and the claims behind particular
cases (for instance that the unlucky prime really is unlucky).

What is proved here, by operation:

| op | checked |
|---|---|
| compare | the three orders of Definition R-2, implemented here |
| divmod, pseudo_divmod | the identities and the degree bound (Theorems V-1, V-2) |
| gcd | normalization (Definition G-1), both cofactor identities, and coprimality of the cofactors (Euclid on specializations, so maximality is proved); over ZZ the content |
| primes | each is prime (deterministic Miller-Rabin), decreasing, and none is skipped |
| resultant | the Sylvester determinant (Definition E-1, rows x^(n-1) f first) and each principal subresultant coefficient as its minor; multivariate cases at specializations |
| normalize_fraction | equality of fractions, the monic denominator (Definition F-1), and coprimality |
| squarefree | the product identity, monic and squarefree factors, pairwise coprime (Definition Z-1) |
| factor_mod_p | the product identity, monic factors, irreducibility by distinct-degree factorization, the sort order |
| factor, factor_integer | the product identity, normalization, distinct squarefree pairwise coprime factors, the sort order; irreducibility where degree patterns modulo small primes prove it (counted, not required) |
| groebner, membership | monic, reduced, sorted bases; each generator reduces to 0; Buchberger's criterion on every pair (Theorem B-2); the remainder and the membership answer. Not checked: that the basis lies in the generators' ideal |
| sturm, count_roots, isolate | Definition T-1 recomputed; Sturm counts (Theorem T-2) for the counts, for each isolating interval, and for all real roots within Cauchy's bound |

Run from the repository root:

    python3 tests/reference/check_corpus.py [-v]
"""
import re
import sys
from fractions import Fraction
from math import gcd
from pathlib import Path

CORPUS = Path(__file__).with_name('corpus.txt')
VERBOSE = '-v' in sys.argv[1:]


class CorpusError(Exception):
    pass


def require(condition, case, message):
    if not condition:
        raise CorpusError(f'{case.name}: {message}')


# ---------------------------------------------------------------- orders

def order_key(kind):
    """The sort key of exponent vectors (aligned with the order's variable list, largest first)."""
    if kind == 'lex':
        return lambda e: e
    if kind == 'grlex':
        return lambda e: (sum(e), e)
    if kind == 'grevlex':
        return lambda e: (sum(e), tuple(-x for x in reversed(e)))
    raise ValueError(kind)


# ---------------------------------------------------------------- cases and parsing

class Case:
    def __init__(self, name):
        self.name = name
        self.op = None
        self.ring = None          # 'QQ', 'ZZ', or a prime p
        self.vars = None          # variable indices, largest first
        self.order = 'lex'
        self.ins = []             # (key, text)
        self.outs = []

    @property
    def modulus(self):
        return self.ring if isinstance(self.ring, int) else 0

    def key(self):
        return order_key(self.order)

    def values(self, which, key):
        return [text for k, text in (self.ins if which == 'in' else self.outs) if k == key]

    def one(self, which, key):
        found = self.values(which, key)
        require(len(found) == 1, self, f'expected one "{which} {key}", found {len(found)}')
        return found[0]


NUMBER = r'\d+(?:/\d+)?'
MONOMIAL = r'x\d+(?:\^\d+)?(?:\*x\d+(?:\^\d+)?)*'
TERM = re.compile(rf'(?:({NUMBER})|(?:({NUMBER})\*)?({MONOMIAL}))')
FACTOR = re.compile(r'x(\d+)(?:\^(\d+))?')


def parse_number(text, case):
    value = Fraction(text)
    require(str(value) == text.lstrip('+'), case, f'non-canonical number {text!r}')
    return value


def parse_poly(text, case):
    """A dict from exponent tuples to coefficients; checks canonical form and term order."""
    if text == '0':
        return {}
    parts = re.split(r' ([+-]) ', text)
    terms, signs = parts[0::2], ['+'] + parts[1::2]
    if terms[0].startswith('-'):
        signs[0], terms[0] = '-', terms[0][1:]
    poly, printed = {}, []
    for sign, term in zip(signs, terms):
        match = TERM.fullmatch(term)
        require(match is not None, case, f'bad term {term!r} in {text!r}')
        if match.group(1) is not None:
            coeff, monomial = match.group(1), ''
        else:
            coeff, monomial = match.group(2) or '1', match.group(3)
        c = parse_number(coeff, case)
        require(c != 0, case, f'zero coefficient in {text!r}')
        require(not (match.group(2) and c == 1), case, f'coefficient 1 written out in {text!r}')
        if sign == '-':
            c = -c
        if case.modulus:
            require(c.denominator == 1 and 0 < c < case.modulus, case, f'coefficient outside [1, p) in {text!r}')
            c = int(c)
        elif case.ring == 'ZZ':
            require(c.denominator == 1, case, f'non-integer coefficient in {text!r}')
        exponents = [0] * len(case.vars)
        last = -1
        for factor in monomial.split('*') if monomial else []:
            fm = FACTOR.fullmatch(factor)
            index, e = int(fm.group(1)), int(fm.group(2) or 1)
            require(fm.group(2) is None or e >= 2, case, f'exponent written as 0 or 1 in {text!r}')
            require(index > last, case, f'variables not increasing in {term!r}')
            require(index in case.vars, case, f'variable x{index} not declared')
            last = index
            exponents[case.vars.index(index)] = e
        e = tuple(exponents)
        require(e not in poly, case, f'repeated monomial in {text!r}')
        poly[e] = c
        printed.append(e)
    key = case.key()
    require(all(key(a) > key(b) for a, b in zip(printed, printed[1:])), case,
            f'terms not in decreasing {case.order} order: {text!r}')
    return poly


def read_corpus(path):
    cases, case = [], None
    for number, raw in enumerate(path.read_text().splitlines(), 1):
        line = raw.rstrip('\n')
        if not line or line.startswith('#'):
            continue
        word, _, rest = line.partition(' ')
        if word == 'case':
            if case is not None:
                raise CorpusError(f'line {number}: case {rest} starts before "end"')
            case = Case(rest)
            continue
        if case is None:
            raise CorpusError(f'line {number}: {word!r} outside a case')
        if word == 'op':
            case.op = rest
        elif word == 'ring':
            case.ring = int(rest[3:]) if rest.startswith('GF ') else rest
        elif word == 'var':
            case.vars = [int(rest)]
        elif word == 'order':
            kind, *indices = rest.split(' ')
            case.order, case.vars = kind, [int(i) for i in indices]
        elif word in ('in', 'out'):
            key, _, value = rest.partition(' ')
            (case.ins if word == 'in' else case.outs).append((key, value))
        elif word == 'end':
            cases.append(case)
            case = None
        else:
            raise CorpusError(f'line {number}: unknown line {line!r}')
    if case is not None:
        raise CorpusError(f'case {case.name} has no "end"')
    names = [c.name for c in cases]
    duplicates = {n for n in names if names.count(n) > 1}
    if duplicates:
        raise CorpusError(f'repeated case names: {sorted(duplicates)}')
    return cases


# ---------------------------------------------------------------- sparse arithmetic

def norm(c, m):
    return c % m if m else c


def inverse(c, m):
    return pow(int(c), -1, m) if m else 1 / Fraction(c)


def p_add(a, b, m=0, sign=1):
    out = dict(a)
    for e, c in b.items():
        v = norm(out.get(e, 0) + sign * c, m)
        if v:
            out[e] = v
        else:
            out.pop(e, None)
    return out


def p_sub(a, b, m=0):
    return p_add(a, b, m, -1)


def p_scale(a, c, m=0):
    return {e: v for e, v in ((e, norm(x * c, m)) for e, x in a.items()) if v}


def p_mul(a, b, m=0):
    out = {}
    for e1, c1 in a.items():
        for e2, c2 in b.items():
            e = tuple(x + y for x, y in zip(e1, e2))
            out[e] = norm(out.get(e, 0) + c1 * c2, m)
    return {e: c for e, c in out.items() if c}


def p_const(c, nvars):
    return {tuple([0] * nvars): c} if c else {}


def p_shift_mul(a, e, c, m=0):
    """a times the term c * X^e."""
    return {tuple(x + y for x, y in zip(k, e)): norm(v * c, m) for k, v in a.items() if norm(v * c, m)}


def leading(p, key):
    e = max(p, key=key)
    return e, p[e]


def divides(a, b):
    return all(x <= y for x, y in zip(a, b))


def reduce_full(f, basis, key, m=0):
    """The remainder of f on division by basis (Theorem V-4): no term divisible by a leading monomial."""
    leads = [leading(g, key) for g in basis]
    p, remainder = dict(f), {}
    while p:
        e, c = leading(p, key)
        for g, (eg, cg) in zip(basis, leads):
            if divides(eg, e):
                shift = tuple(x - y for x, y in zip(e, eg))
                p = p_sub(p, p_shift_mul(g, shift, c * inverse(cg, m), m), m)
                break
        else:
            remainder[e] = c
            del p[e]
    return remainder


def s_polynomial(f, g, key, m=0):
    ef, cf = leading(f, key)
    eg, cg = leading(g, key)
    lcm = tuple(max(x, y) for x, y in zip(ef, eg))
    return p_sub(p_shift_mul(f, tuple(x - y for x, y in zip(lcm, ef)), inverse(cf, m), m),
                 p_shift_mul(g, tuple(x - y for x, y in zip(lcm, eg)), inverse(cg, m), m), m)


def evaluate_at(p, position, value):
    """Substitute a rational for the variable at `position`."""
    out = {}
    for e, c in p.items():
        k = e[:position] + (0,) + e[position + 1:]
        out[k] = out.get(k, 0) + c * Fraction(value) ** e[position]
    return {k: c for k, c in out.items() if c}


def specialize(p, keep, point):
    """Substitute point[i] for every variable but position `keep`: a dense univariate polynomial."""
    coeffs = {}
    for e, c in p.items():
        v = Fraction(c)
        for i, x in enumerate(e):
            if i != keep:
                v *= Fraction(point[i]) ** x
        coeffs[e[keep]] = coeffs.get(e[keep], 0) + v
    return trim([coeffs.get(i, Fraction(0)) for i in range(max(coeffs, default=-1) + 1)])


def degree_in(p, position):
    return max((e[position] for e in p), default=-1)


# ---------------------------------------------------------------- dense univariate arithmetic (coefficients ascending)

def trim(a):
    a = list(a)
    while a and a[-1] == 0:
        a.pop()
    return a


def to_dense(p):
    if not p:
        return []
    d = max(e[0] for e in p)
    return trim([Fraction(p.get((i,), 0)) for i in range(d + 1)])


def to_sparse(a):
    return {(i,): c for i, c in enumerate(a) if c}


def d_deg(a):
    return len(a) - 1


def d_add(a, b, sign=1):
    n = max(len(a), len(b))
    return trim([(a[i] if i < len(a) else 0) + sign * (b[i] if i < len(b) else 0) for i in range(n)])


def d_mul(a, b):
    if not a or not b:
        return []
    out = [Fraction(0)] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        for j, y in enumerate(b):
            out[i + j] += x * y
    return trim(out)


def d_divmod(a, b):
    a, q = list(a), [Fraction(0)] * max(len(a) - len(b) + 1, 1)
    while len(a) >= len(b) and a:
        shift, c = len(a) - len(b), a[-1] / b[-1]
        q[shift] = c
        for i, y in enumerate(b):
            a[i + shift] -= c * y
        a = trim(a)
    return trim(q), a


def d_monic(a):
    return [c / a[-1] for c in a] if a else []


def d_gcd(a, b):
    a, b = trim(a), trim(b)
    while b:
        a, b = b, d_divmod(a, b)[1]
    return d_monic(a)


def d_diff(a):
    return trim([i * a[i] for i in range(1, len(a))])


def d_eval(a, x):
    v = Fraction(0)
    for c in reversed(a):
        v = v * x + c
    return v


def content(p):
    """The positive gcd of an integer polynomial's coefficients."""
    g = 0
    for c in p.values():
        g = gcd(g, int(c))
    return g


# ---------------------------------------------------------------- dense arithmetic modulo p (coefficients ascending, ints)

def m_trim(a):
    a = list(a)
    while a and a[-1] == 0:
        a.pop()
    return a


def m_mul(a, b, p):
    if not a or not b:
        return []
    out = [0] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        if x:
            for j, y in enumerate(b):
                out[i + j] = (out[i + j] + x * y) % p
    return m_trim(out)


def m_divmod(a, b, p):
    a, inv = list(a), pow(b[-1], -1, p)
    q = [0] * max(len(a) - len(b) + 1, 1)
    while len(a) >= len(b) and a:
        shift, c = len(a) - len(b), a[-1] * inv % p
        q[shift] = c
        for i, y in enumerate(b):
            a[i + shift] = (a[i + shift] - c * y) % p
        a = m_trim(a)
    return m_trim(q), a


def m_gcd(a, b, p):
    a, b = m_trim(a), m_trim(b)
    while b:
        a, b = b, m_divmod(a, b, p)[1]
    if not a:
        return a
    inv = pow(a[-1], -1, p)
    return [c * inv % p for c in a]


def m_powmod(base, n, f, p):
    result, base = [1], m_divmod(base, f, p)[1]
    while n:
        if n & 1:
            result = m_divmod(m_mul(result, base, p), f, p)[1]
        base = m_divmod(m_mul(base, base, p), f, p)[1]
        n >>= 1
    return result


def m_ddf(f, p):
    """Distinct-degree factorization of a monic squarefree f modulo p: the degrees of its irreducible factors."""
    degrees, h, d = [], [0, 1], 0
    while len(f) - 1 >= 2 * (d + 1):
        d += 1
        h = m_powmod(h, p, f, p)
        diff = m_trim([(h[i] if i < len(h) else 0) - (1 if i == 1 else 0) for i in range(max(len(h), 2))])
        diff = [c % p for c in diff]
        g = m_gcd(f, diff, p)
        if len(g) > 1:
            degrees += [d] * ((len(g) - 1) // d)
            f = m_divmod(f, g, p)[0]
            h = m_divmod(h, f, p)[1]
    if len(f) > 1:
        degrees.append(len(f) - 1)
    return degrees


# ---------------------------------------------------------------- primes

def is_prime(n):
    """Deterministic Miller-Rabin; these bases decide every n < 3.3 * 10^24."""
    if n < 2:
        return False
    small = (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37)
    for q in small:
        if n % q == 0:
            return n == q
    d, s = n - 1, 0
    while d % 2 == 0:
        d, s = d // 2, s + 1
    for a in small:
        x = pow(a, d, n)
        if x in (1, n - 1):
            continue
        for _ in range(s - 1):
            x = x * x % n
            if x == n - 1:
                break
        else:
            return False
    return True


SMALL_PRIMES = [q for q in range(3, 400) if is_prime(q)]


# ---------------------------------------------------------------- determinants and resultants

def determinant(rows):
    """Exact Gaussian elimination over the rationals."""
    a = [[Fraction(x) for x in row] for row in rows]
    n, det = len(a), Fraction(1)
    for col in range(n):
        pivot = next((r for r in range(col, n) if a[r][col] != 0), None)
        if pivot is None:
            return Fraction(0)
        if pivot != col:
            a[col], a[pivot] = a[pivot], a[col]
            det = -det
        det *= a[col][col]
        for r in range(col + 1, n):
            factor = a[r][col] / a[col][col]
            if factor:
                for c in range(col, n):
                    a[r][c] -= factor * a[col][c]
    return det


def subresultant_matrix(f, g, j):
    """Rows x^(n-j-1) f, ..., f, x^(m-j-1) g, ..., g; columns x^(m+n-j-1) down to x^j (j = 0: Sylvester)."""
    m, n = d_deg(f), d_deg(g)
    top = m + n - j - 1
    rows = []
    for poly, count in ((f, n - j), (g, m - j)):
        for shift in range(count - 1, -1, -1):
            rows.append([poly[k - shift] if 0 <= k - shift < len(poly) else 0 for k in range(top, j - 1, -1)])
    return rows


def resultant(f, g):
    return determinant(subresultant_matrix(f, g, 0))


def psc(f, g, j):
    return determinant(subresultant_matrix(f, g, j))


# ---------------------------------------------------------------- Sturm sequences

def sturm_sequence(f):
    """Definition T-1: f, f', then the negated remainders, without normalization."""
    seq = [f, d_diff(f)]
    while True:
        r = d_divmod(seq[-2], seq[-1])[1]
        if not r:
            return seq
        seq.append([-c for c in r])


def sign_changes(seq, x):
    signs = [s for s in ((d_eval(p, x) > 0) - (d_eval(p, x) < 0) for p in seq) if s]
    return sum(1 for a, b in zip(signs, signs[1:]) if a != b)


def sturm_count(f, a, b):
    seq = sturm_sequence(f)
    return sign_changes(seq, a) - sign_changes(seq, b)


def cauchy_bound(f):
    return 1 + max(abs(c / f[-1]) for c in f[:-1]) if len(f) > 1 else Fraction(1)


# ---------------------------------------------------------------- coprimality by specialization

POINTS = [(3, -2, 5, 7, -11, 13), (1, 4, -3, 2, 9, -5), (-7, 6, 1, -4, 3, 8), (2, -9, 4, 11, -1, 6)]


def coprime(a, b):
    """True when a and b have no common factor of positive degree (a, b not both zero).

    A common factor of positive degree in the variable at position k survives every
    specialization of the other variables that keeps both degrees in that variable,
    so one coprime specialization per variable proves coprimality.
    """
    if not a or not b:
        other = b or a
        return all(sum(e) == 0 for e in other)
    nvars = len(next(iter(a)))
    for k in range(nvars):
        da, db = degree_in(a, k), degree_in(b, k)
        if da <= 0 or db <= 0:
            continue
        for point in POINTS:
            sa, sb = specialize(a, k, point), specialize(b, k, point)
            if d_deg(sa) == da and d_deg(sb) == db:
                if len(d_gcd(sa, sb)) == 1:
                    break
        else:
            return False
    return True


# ---------------------------------------------------------------- checks by operation

def check_compare(case):
    a, b = parse_poly(case.one('in', 'a'), case), parse_poly(case.one('in', 'b'), case)
    require(len(a) == 1 and len(b) == 1 and set(a.values()) == {1} == set(b.values()), case, 'operands must be monomials')
    key = case.key()
    ka, kb = key(next(iter(a))), key(next(iter(b)))
    expected = (ka > kb) - (ka < kb)
    require(int(case.one('out', 'compare')) == expected, case, f'compare should be {expected}')


def check_divmod(case):
    f, g = to_dense(parse_poly(case.one('in', 'f'), case)), to_dense(parse_poly(case.one('in', 'g'), case))
    q, r = to_dense(parse_poly(case.one('out', 'quotient'), case)), to_dense(parse_poly(case.one('out', 'remainder'), case))
    require(d_add(d_mul(q, g), r) == f, case, 'f != q g + r')
    require(d_deg(r) < d_deg(g), case, 'deg r >= deg g')


def check_pseudo_divmod(case):
    f, g = to_dense(parse_poly(case.one('in', 'f'), case)), to_dense(parse_poly(case.one('in', 'g'), case))
    q, r = to_dense(parse_poly(case.one('out', 'quotient'), case)), to_dense(parse_poly(case.one('out', 'remainder'), case))
    delta = max(d_deg(f) - d_deg(g) + 1, 0)
    require(int(case.one('out', 'delta')) == delta, case, 'delta')
    require(d_add(d_mul(q, g), r) == [c * g[-1] ** delta for c in f], case, 'lc(g)^delta f != q g + r')
    require(d_deg(r) < d_deg(g), case, 'deg r >= deg g')


def check_gcd(case):
    f, h = parse_poly(case.one('in', 'f'), case), parse_poly(case.one('in', 'h'), case)
    g = parse_poly(case.one('out', 'gcd'), case)
    key = case.key()
    if not f and not h:
        require(not g and not case.values('out', 'cofactor_f'), case, 'gcd(0, 0) = 0 without cofactors')
        return
    require(bool(g), case, 'nonzero inputs have a nonzero gcd')
    lc = leading(g, key)[1]
    if case.ring == 'QQ':
        require(lc == 1, case, 'the gcd over QQ is monic in the order')
    else:
        require(lc > 0, case, 'the gcd over ZZ has a positive leading coefficient')
        require(content(g) == gcd(content(f), content(h)), case, 'the gcd over ZZ carries the gcd of the contents')
    cf, ch = parse_poly(case.one('out', 'cofactor_f'), case), parse_poly(case.one('out', 'cofactor_h'), case)
    require(p_mul(g, cf) == f and p_mul(g, ch) == h, case, 'cofactor identities')
    require(coprime(cf, ch), case, 'the cofactors are not coprime: the gcd is not maximal')
    if case.ring == 'ZZ' and cf and ch:
        require(gcd(content(cf), content(ch)) == 1, case, 'the cofactors share a content')
    if case.name == 'gcd.uni.unlucky_first_prime':
        p0 = prime_sequence_head()
        require(m_trim([int(c) % p0 for c in to_dense(cf)]) == m_trim([int(c) % p0 for c in to_dense(ch)]),
                case, 'the cofactors should agree modulo the first prime')


_PRIMES = []


def prime_sequence_head():
    return _PRIMES[0]


def check_primes(case):
    primes = [int(t) for t in case.values('out', 'prime')]
    require(primes, case, 'no primes')
    require(all(is_prime(p) for p in primes), case, 'a listed number is not prime')
    require(all(a > b for a, b in zip(primes, primes[1:])), case, 'not decreasing')
    require(primes[0] < 2 ** 62, case, 'the first prime is not below 2^62')
    bounds = [2 ** 62] + primes
    for high, low in zip(bounds, bounds[1:]):
        require(not any(is_prime(n) for n in range(low + 1, high)), case, f'a prime between {low} and {high} is missing')
    _PRIMES[:] = primes


def check_resultant(case):
    v = int(case.one('in', 'variable'))
    f, g = parse_poly(case.one('in', 'f'), case), parse_poly(case.one('in', 'g'), case)
    res = parse_poly(case.one('out', 'resultant'), case)
    position = case.vars.index(v)
    require(position not in {i for e in res for i, x in enumerate(e) if x}, case, 'the resultant contains the variable')
    if len(case.vars) == 1:
        fd, gd = to_dense(f), to_dense(g)
        expected = resultant(fd, gd)
        require(to_dense(res) == trim([expected]), case, f'resultant should be {expected}')
        pscs = {int(t.split()[0]): Fraction(t.split()[1]) for t in case.values('out', 'psc')}
        require(sorted(pscs) == list(range(min(d_deg(fd), d_deg(gd)))), case, 'psc indices')
        for j, value in pscs.items():
            require(psc(fd, gd, j) == value, case, f'psc {j} should be {psc(fd, gd, j)}')
        return
    df, dg = degree_in(f, position), degree_in(g, position)
    checked = 0
    for point in POINTS + [(5, 3, -2, 1, 1, 1), (-3, 7, 2, 1, 1, 1)]:
        sf, sg = specialize(f, position, point), specialize(g, position, point)
        if d_deg(sf) != df or d_deg(sg) != dg:
            continue
        value = sum((Fraction(c) * prod_point(e, point, position) for e, c in res.items()), Fraction(0))
        require(resultant(sf, sg) == value, case, f'resultant disagrees at {point}')
        checked += 1
    require(checked >= 3, case, 'too few specializations keep the degrees')


def prod_point(e, point, skip):
    v = Fraction(1)
    for i, x in enumerate(e):
        if i != skip:
            v *= Fraction(point[i]) ** x
    return v


def check_fraction(case):
    p, q = parse_poly(case.one('in', 'p'), case), parse_poly(case.one('in', 'q'), case)
    num, den = parse_poly(case.one('out', 'numerator'), case), parse_poly(case.one('out', 'denominator'), case)
    require(bool(den) and leading(den, case.key())[1] == 1, case, 'the denominator is monic in the order')
    require(p_mul(num, q) == p_mul(den, p), case, 'num/den != p/q')
    if not p:
        require(not num and den == p_const(Fraction(1), len(case.vars)), case, '0 normalizes to 0/1')
    else:
        require(coprime(num, den), case, 'numerator and denominator share a factor')


def check_squarefree(case):
    f = to_dense(parse_poly(case.one('in', 'f'), case))
    unit = parse_number(case.one('out', 'unit'), case)
    require(unit == f[-1], case, 'the unit is the leading coefficient')
    factors = []
    for text in case.values('out', 'factor'):
        e, poly = text.split(' ', 1)
        factors.append((int(e), to_dense(parse_poly(poly, case))))
    exps = [e for e, _ in factors]
    require(exps == sorted(set(exps)) and all(e >= 1 for e in exps), case, 'exponents distinct and increasing')
    product = [unit]
    for e, a in factors:
        require(d_deg(a) >= 1 and a[-1] == 1, case, 'factors are monic of positive degree')
        require(len(d_gcd(a, d_diff(a))) == 1, case, 'a factor is not squarefree')
        for _ in range(e):
            product = d_mul(product, a)
    require(product == f, case, 'f != c prod a_i^i')
    for i, (_, a) in enumerate(factors):
        for _, b in factors[i + 1:]:
            require(len(d_gcd(a, b)) == 1, case, 'factors are not pairwise coprime')


def check_factor_mod_p(case):
    p = case.modulus
    f = m_trim([int(c) for c in to_dense(parse_poly(case.one('in', 'f'), case))])
    require(f[-1] == 1, case, 'the input is monic')
    require(len(m_gcd(f, m_trim([i * f[i] % p for i in range(1, len(f))]), p)) == 1, case, 'the input is squarefree')
    factors = [m_trim([int(c) for c in to_dense(parse_poly(t, case))]) for t in case.values('out', 'factor')]
    product = [1]
    for a in factors:
        require(a[-1] == 1, case, 'factors are monic')
        require(m_ddf(a, p) == [len(a) - 1], case, f'factor {a} is reducible')
        product = m_mul(product, a, p)
    require(product == f, case, 'the product of the factors is not f')
    keys = [(len(a), list(reversed(a))) for a in factors]
    require(keys == sorted(keys) and len(set(map(str, keys))) == len(keys), case, 'factors not sorted by degree, then coefficients')


IRREDUCIBILITY = {'proved': 0, 'unproved': 0}


def irreducible_by_degrees(a):
    """Whether degree patterns modulo small primes prove the integer polynomial a irreducible over Q."""
    n = d_deg(a)
    if n <= 1:
        return True
    possible = set(range(n + 1))
    ints = [int(c) for c in a]
    for p in SMALL_PRIMES:
        if ints[-1] % p == 0:
            continue
        fp = m_trim([c % p for c in ints])
        inv = pow(fp[-1], -1, p)
        fp = [c * inv % p for c in fp]
        if len(m_gcd(fp, m_trim([i * fp[i] % p for i in range(1, len(fp))]), p)) != 1:
            continue
        sums = {0}
        for d in m_ddf(fp, p):
            sums |= {s + d for s in sums}
        possible &= sums
        if possible <= {0, n}:
            return True
    return False


def check_factor(case):
    f = to_dense(parse_poly(case.one('in', 'f'), case))
    unit = parse_number(case.one('out', 'unit'), case)
    factors = []
    for text in case.values('out', 'factor'):
        e, poly = text.split(' ', 1)
        factors.append((int(e), to_dense(parse_poly(poly, case))))
    product = [unit]
    for e, a in factors:
        require(e >= 1 and d_deg(a) >= 1, case, 'factors have positive degree and exponent')
        if case.ring == 'QQ':
            require(a[-1] == 1, case, 'factors over QQ are monic')
        else:
            require(a[-1] > 0 and content(to_sparse(a)) == 1, case, 'factors over ZZ are primitive with positive leading coefficient')
        require(len(d_gcd(a, d_diff(a))) == 1, case, 'a factor is not squarefree')
        for _ in range(e):
            product = d_mul(product, a)
    require(product == f, case, 'f != unit * prod f_i^e_i')
    for i, (_, a) in enumerate(factors):
        for _, b in factors[i + 1:]:
            require(len(d_gcd(a, b)) == 1, case, 'factors are not pairwise coprime')
    keys = [(d_deg(a), list(reversed(a))) for _, a in factors]
    require(keys == sorted(keys), case, 'factors not sorted by degree, then coefficients')
    for _, a in factors:
        lcm_den = 1
        for c in a:
            lcm_den = lcm_den * c.denominator // gcd(lcm_den, c.denominator)
        integral = [c * lcm_den for c in a]
        IRREDUCIBILITY['proved' if irreducible_by_degrees(integral) else 'unproved'] += 1


def check_basis(case, generators, basis):
    m, key = case.modulus, case.key()
    leads = [leading(g, key) for g in basis]
    require(all(c == 1 for _, c in leads), case, 'basis elements are monic')
    require(all(key(a[0]) > key(b[0]) for a, b in zip(leads, leads[1:])), case, 'basis not sorted by leading monomial, largest first')
    for i, g in enumerate(basis):
        for j, (e, _) in enumerate(leads):
            if i != j:
                require(not any(divides(e, t) for t in g), case, 'basis not reduced')
    for f in generators:
        require(not reduce_full(f, basis, key, m), case, 'a generator does not reduce to 0')
    for i in range(len(basis)):
        for j in range(i + 1, len(basis)):
            require(not reduce_full(s_polynomial(basis[i], basis[j], key, m), basis, key, m), case,
                    f'S-polynomial of elements {i} and {j} does not reduce to 0 (Theorem B-2)')


def check_groebner(case):
    generators = [parse_poly(t, case) for t in case.values('in', 'generator')]
    basis = [parse_poly(t, case) for t in case.values('out', 'basis')]
    if all(not g for g in generators):
        require(not basis, case, 'the zero ideal has the empty basis')
        return
    check_basis(case, generators, basis)


def check_membership(case):
    generators = [parse_poly(t, case) for t in case.values('in', 'generator')]
    basis = [parse_poly(t, case) for t in case.values('out', 'basis')]
    check_basis(case, generators, basis)
    f = parse_poly(case.one('in', 'f'), case)
    r = reduce_full(f, basis, case.key())
    require(r == parse_poly(case.one('out', 'remainder'), case), case, 'remainder')
    require(case.one('out', 'member') == ('true' if not r else 'false'), case, 'membership answer')


def check_sturm(case):
    f = to_dense(parse_poly(case.one('in', 'f'), case))
    listed = [to_dense(parse_poly(t, case)) for t in case.values('out', 'sturm')]
    require(listed == sturm_sequence(f), case, 'Sturm sequence differs from Definition T-1')


def check_count_roots(case):
    f = to_dense(parse_poly(case.one('in', 'f'), case))
    a, b = parse_number(case.one('in', 'a'), case), parse_number(case.one('in', 'b'), case)
    require(a < b and d_eval(f, a) != 0 and d_eval(f, b) != 0, case, 'needs a < b, f(a) f(b) != 0')
    require(int(case.one('out', 'count')) == sturm_count(f, a, b), case, 'count differs from the Sturm count')


def check_isolate(case):
    f = to_dense(parse_poly(case.one('in', 'f'), case))
    r = cauchy_bound(f)
    total = sturm_count(f, -r, r)
    require(int(case.one('out', 'count')) == total, case, f'count should be {total}')
    roots = case.values('out', 'root')
    require(len(roots) == total, case, 'one root line per distinct real root')
    previous = None
    for text in roots:
        kind, *values = text.split(' ')
        ends = [parse_number(v, case) for v in values]
        if kind == 'exact':
            require(len(ends) == 1 and d_eval(f, ends[0]) == 0, case, f'{text} is not a root')
            low = high = ends[0]
        else:
            require(kind == 'interval' and len(ends) == 2, case, f'bad root line {text!r}')
            low, high = ends
            require(low < high and d_eval(f, low) != 0 and d_eval(f, high) != 0, case, f'bad interval {text}')
            require(sturm_count(f, low, high) == 1, case, f'{text} does not isolate one root')
        require(previous is None or previous < low, case, 'roots not sorted and disjoint')
        previous = high


CHECKS = {
    'compare': check_compare, 'divmod': check_divmod, 'pseudo_divmod': check_pseudo_divmod,
    'gcd': check_gcd, 'primes': check_primes, 'resultant': check_resultant,
    'normalize_fraction': check_fraction, 'squarefree': check_squarefree,
    'factor_mod_p': check_factor_mod_p, 'factor': check_factor, 'factor_integer': check_factor,
    'groebner': check_groebner, 'membership': check_membership, 'sturm': check_sturm,
    'count_roots': check_count_roots, 'isolate': check_isolate,
}


def main():
    cases = read_corpus(CORPUS)
    # The prime sequence first: the unlucky-prime case refers to it.
    cases.sort(key=lambda c: c.op != 'primes')
    failures, counts = [], {}
    for case in cases:
        try:
            require(case.op in CHECKS, case, f'unknown op {case.op!r}')
            require(case.ring in ('QQ', 'ZZ') or isinstance(case.ring, int), case, 'unknown ring')
            if case.op != 'primes':
                require(case.vars, case, 'no "var" or "order" line')
            if isinstance(case.ring, int):
                require(is_prime(case.ring), case, 'GF modulus is not prime')
            CHECKS[case.op](case)
            counts[case.op] = counts.get(case.op, 0) + 1
            if VERBOSE:
                print(f'ok {case.name}')
        except CorpusError as error:
            failures.append(str(error))
    for op in sorted(counts):
        print(f'{op:20} {counts[op]:4} cases')
    print(f'irreducibility over QQ proved by degree patterns for {IRREDUCIBILITY["proved"]} of '
          f'{IRREDUCIBILITY["proved"] + IRREDUCIBILITY["unproved"]} factors (the others rely on the SymPy generator)')
    if failures:
        print(f'{len(failures)} FAILED:')
        for failure in failures:
            print('  ' + failure)
        sys.exit(1)
    print(f'all {len(cases)} cases agree')


if __name__ == '__main__':
    main()
