# expect: poly.io supports prepared polynomial values, batches, and supported coefficient scalars
from poly import io

def main() raises:
    _ = io.dumps(List[Int]([1, 2]))
