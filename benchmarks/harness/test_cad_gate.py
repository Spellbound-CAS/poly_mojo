"""The CAD interval gate accepts equal or narrower answers and rejects coarser ones."""
from pathlib import Path
import tempfile
import unittest

from harness.cad_gate import compare

BASE = '2\n1 1,-1 Q:1/2 A:-2,0,1:1:3/2\n0 0,0 Q:1 A:-3,0,1:3/2:2\n'


def write(directory, name, text):
    path = Path(directory) / name
    path.write_text(text)
    return path


class Gate(unittest.TestCase):
    def test_narrower_moved_samples_pass(self):
        with tempfile.TemporaryDirectory() as d:
            new = '2\n1 1,-1 Q:1 A:-2,0,1:5/4:3/2\n0 0,0 Q:1 A:-3,0,1:7/4:15/8\n'
            self.assertEqual(compare(write(d, 'a', BASE), write(d, 'b', new)), [])

    def test_wider_interval_fails(self):
        with tempfile.TemporaryDirectory() as d:
            new = '2\n1 1,-1 Q:1/2 A:-2,0,1:0:2\n0 0,0 Q:1 A:-3,0,1:3/2:2\n'
            self.assertTrue(any('wider' in f for f in compare(write(d, 'a', BASE), write(d, 'b', new))))

    def test_changed_signs_or_cells_fail(self):
        with tempfile.TemporaryDirectory() as d:
            signs = '2\n1 1,1 Q:1/2 A:-2,0,1:1:3/2\n0 0,0 Q:1 A:-3,0,1:3/2:2\n'
            self.assertTrue(compare(write(d, 'a', BASE), write(d, 'b', signs)))
            fewer = '1\n1 1,-1 Q:1/2 A:-2,0,1:1:3/2\n'
            self.assertTrue(compare(write(d, 'a', BASE), write(d, 'c', fewer)))
            interval = '2\n1 1,-1 A:-2,0,1:1:3/2 A:-2,0,1:1:3/2\n0 0,0 Q:1 A:-3,0,1:3/2:2\n'
            self.assertTrue(compare(write(d, 'a', BASE), write(d, 'd', interval)))

    def test_exact_replacing_interval_passes(self):
        with tempfile.TemporaryDirectory() as d:
            exact = '2\n1 1,-1 Q:1/2 Q:5/4\n0 0,0 Q:1 A:-3,0,1:3/2:2\n'
            self.assertEqual(compare(write(d, 'a', BASE), write(d, 'b', exact)), [])


if __name__ == '__main__':
    unittest.main()
