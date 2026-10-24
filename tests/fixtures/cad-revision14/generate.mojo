from poly import cad, io, algebraic, dense, sparse
from std.testing import assert_true

def save[T: Copyable & Deinitable](name: String, value: T) raises:
    with open(String("tests/fixtures/cad-revision14/", name, ".txt"), "w") as f:
        f.write(io.dumps(value))

def main() raises:
    var quadratic = cad.decompose([[(0, -2), (2, 1)]], lazard=False)
    assert_true(not quadratic.evidence().projection.lazard)
    save("quadratic", quadratic)
    save("unsatisfiable", cad.cover([cad.positive([(1, 1)]), cad.nonpositive([(1, 1)])]))
    save("fixed", cad.decompose([[((1, 0), 1), ((0, 1), -1)]], fixed=[(1, 2)]))
    var a = dense.root_between([-2, 0, 1], 1, 2)
    var p = sparse.asarray([(0, algebraic.negative(a)), (1, algebraic.asreal(1))])
    save("algebraic", cad.decompose([p.copy()]))
    print("legacy CAD archive fixtures accepted and written")
