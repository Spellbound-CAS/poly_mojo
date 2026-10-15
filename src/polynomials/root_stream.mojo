"""Completion receipts for separately transported root-count witnesses."""
from .numeric_options import NumericOutcome


@fieldwise_init
struct RootStreamReceipt(Copyable, Movable):
    """Describe a completed root stream whose local witnesses must also be checked."""
    var outcome: NumericOutcome
    """SUCCESS only when the requested guarantees hold."""
    var complete: Bool
    """Whether all roots have certified counts and coverage."""
    var accurate: Bool
    """Whether all positive-count pieces meet the requested accuracy."""
    var isolated: Bool
    """Whether each positive-count piece contains a single root."""
    var items: Int
    """The number of cells or regions in the transported stream."""
