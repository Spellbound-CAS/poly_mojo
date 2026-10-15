# expect: CAD requires exact real coefficients
from poly import cad, sparse
from apn_mojo import BallContext

def main() raises:
    var input = sparse.asarray([(1, 1)], context=BallContext(precision=80))
    _ = cad.project([input.copy()])
