# expect: Sturm queries require integer or rational points
from poly import dense, roots

def main() raises:
    var plan = dense.sturm_sequence([-2, 0, 1])
    _ = roots.sign_changes(plan, Float64(0))
