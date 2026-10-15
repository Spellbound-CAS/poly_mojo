"""Functional polynomial arrays, explicit numerical contexts, and shared work limits."""
from . import dense, sparse
from ._common import Values, coefficients, modint
from ._conversion import to_dense, to_sparse
from . import algebraic
from . import roots
from . import cad
from . import series, basis
from . import rational
from . import batch

from . import evaluation

from . import exact

from . import io
