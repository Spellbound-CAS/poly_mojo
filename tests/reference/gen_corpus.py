"""Generate exact regression answers with SymPy 1.14.0, independently of Mojo.

The frozen input-only cases preserve the original corpus coverage without a
vendor-specific random generator. Generation never reads expected answers
from corpus.txt or imports check_corpus.py. The latter validates its output
using Python's standard library. SymPy is an optional regeneration dependency, not a package
or test-suite dependency. All arithmetic and root isolation here are exact.

Run from the repository root with Python 3.11+ and sympy==1.14.0 installed:

    python tests/reference/gen_corpus.py --output /tmp/corpus.txt
    python tests/reference/gen_corpus.py --check

The CLI computes all cases before writing, so a failed calculation cannot
leave a partially regenerated corpus. Without either option it writes stdout.
"""
import argparse
from dataclasses import dataclass
from pathlib import Path
import re
import sys

import sympy as sp
from sympy.polys.matrices import DomainMatrix
from sympy.polys.orderings import monomial_key


SYMPY_VERSION = '1.14.0'
DIRECTORY = Path(__file__).resolve().parent


@dataclass
class Case:
    name: str
    op: str
    ring: str
    order: str
    variables: tuple
    lines: list
    inputs: dict

    @property
    def modulus(self):
        return int(self.ring.split()[1]) if self.ring.startswith('GF ') else None

    @property
    def domain(self):
        return sp.GF(self.modulus) if self.modulus else {'QQ': sp.QQ, 'ZZ': sp.ZZ}[self.ring]

    def one(self, key):
        values = self.inputs[key]
        if len(values) != 1:
            raise ValueError(f'{self.name}: expected one input {key}')
        return values[0]

    def poly(self, text):
        """Parse the corpus's restricted polynomial grammar without eval."""
        terms = {}
        for term in text.replace(' - ', ' + -').split(' + '):
            coefficient = sp.S.One
            if term.startswith('-'):
                coefficient = -coefficient
                term = term[1:]
            exponents = [0] * len(self.variables)
            for factor in term.split('*'):
                if re.fullmatch(r'\d+(?:/\d+)?', factor):
                    coefficient *= sp.Rational(factor)
                else:
                    match = re.fullmatch(r'x(\d+)(?:\^(\d+))?', factor)
                    if match is None:
                        raise ValueError(f'{self.name}: invalid polynomial factor {factor!r}')
                    variable = sp.Symbol('x' + match[1])
                    exponents[self.variables.index(variable)] += int(match[2] or 1)
            key = tuple(exponents)
            terms[key] = terms.get(key, sp.S.Zero) + coefficient
        return sp.Poly.from_dict(terms, self.variables, domain=self.domain)

    def input_poly(self, key):
        return self.poly(self.one(key))

    def format(self, value):
        p = sp.Poly(value, *self.variables, domain=self.domain)
        if p.is_zero:
            return '0'
        parts = []
        for exponents, coefficient in p.terms(order=self.order):
            if self.modulus:
                coefficient = int(coefficient) % self.modulus
            factors = []
            for variable, exponent in sorted(zip(self.variables, exponents), key=lambda item: int(str(item[0])[1:])):
                if exponent:
                    factors.append(str(variable) + (f'^{exponent}' if exponent != 1 else ''))
            if not factors or abs(coefficient) != 1:
                factors.insert(0, str(abs(coefficient)))
            body = '*'.join(factors)
            sign = (' - ' if parts else '-') if coefficient < 0 else (' + ' if parts else '')
            parts.append(sign + body)
        return ''.join(parts)


def read_inputs(path):
    """Read input-only case blocks, rejecting output records and malformed headers."""
    cases, names, lines = [], set(), None
    for number, line in enumerate(path.read_text().splitlines(), 1):
        if not line or line.startswith('#'):
            continue
        if line.startswith('case '):
            if lines is not None:
                raise ValueError(f'line {number}: nested case')
            lines = [line]
        elif line == 'end':
            if lines is None:
                raise ValueError(f'line {number}: end outside case')
            headers, inputs = {}, {}
            for item in lines[1:]:
                key, value = item.split(' ', 1)
                if key == 'in':
                    field, text = value.split(' ', 1)
                    inputs.setdefault(field, []).append(text)
                elif key in ('op', 'ring', 'var', 'order') and key not in headers:
                    headers[key] = value
                else:
                    raise ValueError(f'line {number}: invalid or duplicate field {key!r}')
            if not {'op', 'ring'} <= headers.keys():
                raise ValueError(f'line {number}: missing operation or ring')
            if ('var' in headers) + ('order' in headers) != (headers['op'] != 'primes'):
                raise ValueError(f'line {number}: expected one variable/order declaration except for primes')
            order, *indices = headers.get('order', 'lex ' + headers.get('var', '')).split()
            if order not in ('lex', 'grlex', 'grevlex') or len(indices) != len(set(indices)):
                raise ValueError(f'line {number}: invalid order or repeated variable')
            name = lines[0][5:]
            if name in names:
                raise ValueError(f'line {number}: duplicate case {name}')
            names.add(name)
            case = Case(name, headers['op'], headers['ring'], order,
                        tuple(sp.Symbol('x' + str(int(i))) for i in indices), lines, inputs)
            if case.modulus and not sp.isprime(case.modulus):
                raise ValueError(f'{name}: modulus is not prime')
            _ = case.domain
            cases.append(case)
            lines = None
        elif lines is None:
            raise ValueError(f'line {number}: field outside case')
        else:
            lines.append(line)
    if lines is not None or not cases:
        raise ValueError('unterminated or empty input corpus')
    return cases


def principal_subresultant(f, g, variable, index):
    """Compute Definition E-1's Sylvester minor with an exact domain determinant."""
    f, g = sp.Poly(f.as_expr(), variable), sp.Poly(g.as_expr(), variable)
    m, n = f.degree(), g.degree()
    columns = range(m + n - index - 1, index - 1, -1)
    rows = [[poly.nth(column - shift) if column >= shift else 0 for column in columns]
            for poly, count in ((f, n - index), (g, m - index))
            for shift in range(count - 1, -1, -1)]
    matrix = DomainMatrix.from_Matrix(sp.Matrix(rows))
    return matrix.domain.to_sympy(matrix.det())


def factor_key(case, factor):
    coefficients = factor.all_coeffs()
    if case.modulus:
        coefficients = [int(c) % case.modulus for c in coefficients]
    return factor.degree(), coefficients


def answers(case):
    """Compute output records without consulting any stored expected answers."""
    op, fmt = case.op, case.format
    if op == 'compare':
        key = monomial_key(case.order)
        a = key(case.input_poly('a').monoms()[0])
        b = key(case.input_poly('b').monoms()[0])
        return [f'compare {int(a > b) - int(a < b)}']
    if op in ('divmod', 'pseudo_divmod'):
        f, g = case.input_poly('f'), case.input_poly('g')
        delta = max(f.degree() - g.degree() + 1, 0)
        if op == 'pseudo_divmod':
            f = f.mul_ground(g.LC() ** delta)
        q, r = f.div(g)
        return ([f'delta {delta}'] if op == 'pseudo_divmod' else []) + [f'quotient {fmt(q)}', f'remainder {fmt(r)}']
    if op == 'gcd':
        f, h = case.input_poly('f'), case.input_poly('h')
        g = f.gcd(h)
        if not g.is_zero:
            lc = g.terms(order=case.order)[0][1]
            g = g.mul_ground(1 / lc if case.ring == 'QQ' else sp.sign(lc))
        result = [f'gcd {fmt(g)}']
        if not g.is_zero:
            result += [f'cofactor_f {fmt(f.exquo(g))}', f'cofactor_h {fmt(h.exquo(g))}']
        return result
    if op == 'primes':
        p, result = 2**62, []
        for _ in range(100):
            p = sp.prevprime(p)
            result.append(f'prime {p}')
        return result
    if op == 'resultant':
        f, g = case.input_poly('f'), case.input_poly('g')
        variable = sp.Symbol('x' + case.one('variable'))
        result = [f'resultant {fmt(principal_subresultant(f, g, variable, 0))}']
        if len(case.variables) == 1:
            result += [f'psc {j} {principal_subresultant(f, g, variable, j)}'
                       for j in range(min(f.degree(), g.degree()))]
        return result
    if op == 'normalize_fraction':
        p, q = case.input_poly('p'), case.input_poly('q')
        common = p.gcd(q)
        numerator, denominator = p.exquo(common), q.exquo(common)
        scale = 1 / denominator.terms(order=case.order)[0][1]
        return [f'numerator {fmt(numerator.mul_ground(scale))}', f'denominator {fmt(denominator.mul_ground(scale))}']
    if op in ('squarefree', 'factor', 'factor_integer', 'factor_mod_p'):
        f = case.input_poly('f')
        unit, factors = f.sqf_list() if op == 'squarefree' else f.factor_list()
        normalized = []
        for factor, exponent in factors:
            if case.ring == 'QQ':
                unit *= factor.LC() ** exponent
                factor = factor.monic()
            normalized.append((factor, exponent))
        normalized.sort(key=lambda item: item[1] if op == 'squarefree' else factor_key(case, item[0]))
        if op == 'factor_mod_p':
            if unit != 1 or any(exponent != 1 for _, exponent in normalized):
                raise ValueError(f'{case.name}: expected monic squarefree input')
            return [f'factor {fmt(factor)}' for factor, _ in normalized]
        return [f'unit {unit}'] + [f'factor {exponent} {fmt(factor)}' for factor, exponent in normalized]
    if op in ('groebner', 'membership'):
        generators = [case.poly(text).as_expr() for text in case.inputs.get('generator', [])]
        basis = sp.groebner(generators, *case.variables, order=case.order, domain=case.domain)
        result = [f'basis {fmt(poly)}' for poly in basis.polys]
        if op == 'membership':
            remainder = basis.reduce(case.input_poly('f').as_expr())[1]
            result += [f'remainder {fmt(remainder)}', f'member {str(remainder == 0).lower()}']
        return result
    if op == 'sturm':
        # Definition T-1 retains repeated factors; SymPy.sturm squarefrees first.
        f = case.input_poly('f')
        sequence = [f, f.diff()]
        while not sequence[-1].is_zero:
            remainder = -sequence[-2].rem(sequence[-1])
            if remainder.is_zero:
                break
            sequence.append(remainder)
        return [f'sturm {fmt(poly)}' for poly in sequence]
    if op == 'count_roots':
        f = case.input_poly('f')
        a, b = sp.Rational(case.one('a')), sp.Rational(case.one('b'))
        if a >= b or f.eval(a) == 0 or f.eval(b) == 0:
            raise ValueError(f'{case.name}: endpoints must be ordered nonroots')
        return [f'count {f.count_roots(a, b)}']
    if op == 'isolate':
        f = case.input_poly('f').sqf_part()
        # Separate rational roots explicitly; intervals() can leave nonintegral
        # rational roots in open intervals until a later refinement.
        rational = sorted(f.ground_roots())
        for root in rational:
            f = f.exquo(sp.Poly(case.variables[0] - root, *case.variables, domain=sp.QQ))
        roots = [(root, root) for root in rational]
        roots += [ends for ends, _ in f.intervals(eps=sp.Rational(1, 2**50))]
        roots.sort()
        return [f'count {len(roots)}'] + [f'root exact {a}' if a == b else f'root interval {a} {b}' for a, b in roots]
    raise ValueError(f'{case.name}: unknown operation {op}')


def generate(path):
    if sp.__version__ != SYMPY_VERSION:
        raise RuntimeError(f'regeneration requires sympy=={SYMPY_VERSION}, found {sp.__version__}')
    blocks = [f'# Reference corpus for the polynomials package. Generated by tests/reference/gen_corpus.py\n'
              f'# with SymPy {SYMPY_VERSION}; checked by tests/reference/check_corpus.py.\n'
              '# Format: tests/reference/README.md. Do not edit by hand.']
    for case in read_inputs(path):
        try:
            output = answers(case)
        except Exception as error:
            raise ValueError(f'{case.name}: {error}') from error
        blocks.append('\n'.join(case.lines + ['out ' + line for line in output] + ['end']))
    return '\n\n'.join(blocks) + '\n'


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--inputs', type=Path, default=DIRECTORY / 'corpus_inputs.txt')
    mode = parser.add_mutually_exclusive_group()
    mode.add_argument('--output', type=Path)
    mode.add_argument('--check', action='store_true', help='compare regenerated answers with corpus.txt')
    args = parser.parse_args()
    text = generate(args.inputs)
    if args.check:
        if text != (DIRECTORY / 'corpus.txt').read_text():
            parser.exit(1, 'corpus.txt differs from regenerated answers\n')
        count = sum(line.startswith('case ') for line in text.splitlines())
        print(f'all {count} cases reproduce exactly')
    elif args.output:
        args.output.write_text(text)
    else:
        sys.stdout.write(text)


if __name__ == '__main__':
    main()
