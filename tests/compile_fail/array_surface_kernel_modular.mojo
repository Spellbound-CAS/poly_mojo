# expect: to_kernel_bytes
from poly import sparse, modint

def main() raises:
    _ = sparse.to_kernel_bytes([(1, modint(2, 7))])
