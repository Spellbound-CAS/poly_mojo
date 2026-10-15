from poly import cad, dense, algebraic


def main() raises:
    var a = dense.root_between([-2, 0, 1], 1, 2)
    var p = [(0, algebraic.negative(a)), (1, algebraic.asreal(1))]
    var solution = cad.solve([cad.zero(p), cad.positive([(1, 1)])])
    print("satisfiable:", solution.satisfiable())
    print("checked:", cad.check(solution))
    print("exact sample:", cad.sample(solution).value()[0] == a)
    var decomposition = cad.decompose([p.copy()])
    print("cells on the line:", len(cad.cells(decomposition)))
    var checker = cad.checker([p.copy()])
    var receipt = cad.decompose_to([p.copy()], checker)
    from polynomials import Budget
    var budget = Budget()
    print("stream checked:", checker.finish(receipt, budget))
    print("stream stacks and cells:", receipt.stacks, receipt.cells)
