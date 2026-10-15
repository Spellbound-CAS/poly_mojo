"""Regeneration regressions; run with Python and sympy==1.14.0 installed."""
from pathlib import Path
from contextlib import redirect_stderr
from io import StringIO
import tempfile
import unittest
from unittest.mock import patch

import sympy as sp

from gen_corpus import answers, generate, main, read_inputs


class GeneratorTests(unittest.TestCase):
    def case(self, body):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / 'inputs.txt'
            path.write_text('case example\n' + body + '\nend\n')
            return read_inputs(path)[0]

    def test_rejects_expected_answers_in_inputs(self):
        with self.assertRaisesRegex(ValueError, 'invalid or duplicate field'):
            self.case('op gcd\nring QQ\nvar 0\nin f x0\nin h x0\nout gcd x0')

    def test_rejects_ambiguous_headers(self):
        for header in ('var 0\norder lex 0', 'var 0\nvar 1', 'order lex 0 0'):
            with self.subTest(header=header), self.assertRaises(ValueError):
                self.case('op gcd\nring QQ\n' + header)

    def test_rejects_composite_modulus(self):
        with self.assertRaisesRegex(ValueError, 'not prime'):
            self.case('op factor_mod_p\nring GF 15\nvar 0\nin f x0')

    def test_generation_does_not_read_expected_output(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / 'inputs.txt'
            path.write_text('case example\nop divmod\nring QQ\nvar 0\nin f x0^2 - 1\nin g x0 - 1\nend\n')
            text = generate(path)
            self.assertIn('out quotient x0 + 1\nout remainder 0\n', text)
            (path.parent / 'corpus.txt').write_text('deliberately corrupt expected output')
            self.assertEqual(generate(path), text)

    def test_wrong_sympy_version_fails_closed(self):
        with patch('sympy.__version__', '0.0'), self.assertRaisesRegex(RuntimeError, 'requires sympy'):
            generate(Path('does-not-exist'))

    def test_check_rejects_stale_corpus(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            (root / 'corpus_inputs.txt').write_text('case example\nop gcd\nring QQ\nvar 0\nin f x0\nin h x0\nend\n')
            (root / 'corpus.txt').write_text('stale output')
            with patch('gen_corpus.DIRECTORY', root), patch('sys.argv', ['gen_corpus.py', '--check']):
                with redirect_stderr(StringIO()), self.assertRaises(SystemExit) as caught:
                    main()
            self.assertEqual(caught.exception.code, 1)
            self.assertEqual((root / 'corpus.txt').read_text(), 'stale output')

    def test_failed_generation_preserves_output_file(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            (root / 'corpus_inputs.txt').write_text('case example\nop unsupported\nring QQ\nvar 0\nend\n')
            output = root / 'corpus.txt'
            output.write_text('existing output')
            with patch('gen_corpus.DIRECTORY', root), patch('sys.argv', ['gen_corpus.py', '--output', str(output)]):
                with self.assertRaisesRegex(ValueError, 'unknown operation'):
                    main()
            self.assertEqual(output.read_text(), 'existing output')

    def test_reversed_variable_order_normalization(self):
        case = self.case('op normalize_fraction\nring QQ\norder lex 1 0\nin p x0 + 1\nin q 5*x1 - 2*x0')
        self.assertEqual(answers(case), ['numerator 1/5*x0 + 1/5', 'denominator x1 - 2/5*x0'])

    def test_resultant_sign_when_degrees_are_swapped(self):
        for f, g, expected in [('x0^3', 'x0 - 1', '-1'), ('x0 - 1', 'x0^3', '1')]:
            case = self.case(f'op resultant\nring ZZ\nvar 0\nin variable 0\nin f {f}\nin g {g}')
            self.assertEqual(answers(case), [f'resultant {expected}', f'psc 0 {expected}'])

    def test_pseudo_division_scales_even_when_degree_drops(self):
        case = self.case('op pseudo_divmod\nring ZZ\nvar 0\nin f x0^2 - 1\nin g 2*x0')
        self.assertEqual(answers(case), ['delta 2', 'quotient 2*x0', 'remainder -4'])

    def test_field_factors_use_unsigned_coefficients_and_sort_order(self):
        case = self.case('op factor_mod_p\nring GF 5\nvar 0\nin f x0^2 + 4')
        self.assertEqual(answers(case), ['factor x0 + 1', 'factor x0 + 4'])

    def test_rational_factor_unit_and_multiplicity(self):
        case = self.case('op factor\nring QQ\nvar 0\nin f -6*x0^2 + 4*x0 - 2/3')
        self.assertEqual(answers(case), ['unit -6', 'factor 2 x0 - 1/3'])

    def test_repeated_factor_sturm_sequence(self):
        case = self.case('op sturm\nring QQ\nvar 0\nin f x0^3 - 3*x0 + 2')
        self.assertEqual(answers(case), ['sturm x0^3 - 3*x0 + 2', 'sturm 3*x0^2 - 3', 'sturm 2*x0 - 2'])

    def test_exact_nonintegral_root_and_irrational_intervals(self):
        case = self.case('op isolate\nring QQ\nvar 0\nin f x0^3 - 1/3*x0^2 - 2*x0 + 2/3')
        result = answers(case)
        self.assertEqual(result[0], 'count 3')
        self.assertEqual(result[2], 'root exact 1/3')
        bounds = []
        for record in (result[1], result[3]):
            _, _, a, b = record.split()
            a, b = sp.Rational(a), sp.Rational(b)
            self.assertLess(a, b)
            self.assertLess(b - a, sp.Rational(1, 2**40))
            self.assertLess((a*a - 2) * (b*b - 2), 0)
            bounds.append((a, b))
        self.assertLess(bounds[0][1], sp.Rational(1, 3))
        self.assertLess(sp.Rational(1, 3), bounds[1][0])


if __name__ == '__main__':
    unittest.main()
