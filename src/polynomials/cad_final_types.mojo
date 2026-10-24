"""Materialized geometry of an independently verified compact CAD."""
from .algebraic_types import RealAlgebraic


@fieldwise_init
struct CADFinalCell(Copyable, Movable):
    """An absolute cell coordinate and labels, without intermediate proofs.

    The enclosing CADFinalCertificate proves geometry and completeness. This
    output value is not a replacement for that certificate.
    """
    var parent: Int
    """Parent index in the materialized breadth-first list; -1 at its root."""
    var level: Int
    """Zero-based coordinate level; -1 for the root placeholder."""
    var position: Int
    """Alternating stack position: sectors are even and sections are odd."""
    var sample: RealAlgebraic
    """The absolute last coordinate, computed only on an output request."""
    var signs: List[Int]
    """Original-input signs on terminal cells; empty on other levels."""
