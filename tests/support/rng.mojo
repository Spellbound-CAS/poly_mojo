"""SplitMix64, so generated tests are deterministic and failing seeds can be replayed."""


struct SplitMix64(Movable):
    """Each draw adds 0x9e3779b97f4a7c15 to the state and returns the state, finalized."""

    var state: UInt64
    """The generator's state."""

    def __init__(out self, seed: UInt64):
        """Seed a generator.

        Args:
            seed: The initial state.
        """
        self.state = seed

    def next(mut self) -> UInt64:
        """Draw 64 bits.

        Returns:
            The next value.
        """
        self.state += 0x9E3779B97F4A7C15
        var z = self.state
        z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9
        z = (z ^ (z >> 27)) * 0x94D049BB133111EB
        return z ^ (z >> 31)

    def below(mut self, n: Int) -> Int:
        """Draw an integer in `[0, n)`.

        Args:
            n: The bound, positive.

        Returns:
            The value.
        """
        return Int(self.next() % UInt64(n))
