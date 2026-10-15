from poly import dense, algebraic, evaluation
from apn_mojo import ArithmeticContext, FloatFormat, BallContext
from polynomials import NumericOptions


def main() raises:
    var root = dense.root_between([-2, 0, 1], 1, 2)
    var exact = dense.asarray([root.copy(), algebraic.asreal(1)])
    var rounded = dense.asarray(exact, context=ArithmeticContext(format=FloatFormat(precision=80)))
    print("precision:", rounded[0].precision())
    var checked = dense.enclose(exact, 2, context=BallContext(precision=80))
    print("checked:", evaluation.check(checked))
    var adaptive = dense.adaptive_evaluate(exact, 2, options=NumericOptions(initial_precision=8, max_precision=128, absolute_bits=60))
    print("accurate:", adaptive.accurate(), evaluation.check(adaptive))
