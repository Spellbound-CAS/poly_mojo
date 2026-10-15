# expect: format
from poly import dense
from apn_mojo import ArithmeticContext, FloatFormat

def main() raises:
    _ = dense.format(dense.asarray([1, 2], context=ArithmeticContext(format=FloatFormat(precision=53))))
