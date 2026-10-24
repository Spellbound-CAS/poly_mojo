"""Sparse polynomials stored as packed exponent keys and a coefficient list.

A monomial over variables 0..n-1 in LEX order becomes one word by mixed-radix
(Kronecker) substitution with radix_i = (largest exponent of variable i) + 1,
variable 0 most significant, so integer order is LEX order. Keys increase with
the terms. Products, powers and exact quotients know their exponent bounds
exactly (the leading forms in each variable multiply without cancellation over
an integral domain), so their kernels emit keys already in the result's radix
and no Monomial is built until the terms are materialized.
"""
from .coefficients import CoefficientRing
from .errors import Budget, PolyError, _allocation, _charge_size, _product
from .monomial import Monomial, MonomialOrder
from .sparse import SparsePoly, Term

comptime _KEY_LIMIT = UInt64(1) << 62
"""Keys stay below this, as in the packed product kernels."""


struct _PackedPoly[C: CoefficientRing](Copyable, Movable, Sized):
    """Increasing mixed-radix keys over variables 0..n-1 with their nonzero coefficients."""
    var order: MonomialOrder
    var ring: Self.C.Context
    var radix: List[UInt64]
    """radix_i = 1 + the largest exponent of variable i (1 for an absent variable)."""
    var keys: List[UInt64]
    var values: List[Self.C]

    def __init__(out self, var order: MonomialOrder, var ring: Self.C.Context, var radix: List[UInt64], var keys: List[UInt64], var values: List[Self.C]):
        self.order = order^
        self.ring = ring^
        self.radix = radix^
        self.keys = keys^
        self.values = values^

    def __len__(self) -> Int:
        return len(self.keys)

    def weights(self) -> List[UInt64]:
        """weight_i = radix_(i+1) * ... * radix_(n-1): the place value of variable i."""
        return _weights(self.radix)

    @staticmethod
    def from_sparse(p: SparsePoly[Self.C], mut budget: Budget) raises PolyError -> Optional[Self]:
        """The packed form of p, or None unless p is in LEX over variables 0..n-1 with keys below 2^62."""
        if p.order.kind != MonomialOrder.LEX or not p.order._identity():
            return None
        var n = len(p.order.variables)
        var top = List[UInt64](length=n, fill=0)
        for term in p.terms:
            if not p.order._within(term.monomial):
                return None
            for pair in term.monomial.pairs:
                top[Int(pair[0])] = max(top[Int(pair[0])], UInt64(pair[1]))
        var radix = _radix_of(top)
        if not radix:
            return None
        var weights = _weights(radix.value())
        _charge_size(len(p.terms), budget)
        _allocation[UInt64](len(p.terms), budget)
        _allocation[Self.C](len(p.terms), budget)
        var keys = List[UInt64](capacity=len(p.terms))
        var values = List[Self.C](capacity=len(p.terms))
        for term in p.terms:
            var key = UInt64(0)
            for pair in term.monomial.pairs:
                key += UInt64(pair[1]) * weights[Int(pair[0])]
            keys.append(key)
            values.append(term.coefficient.copy())
        return Self(p.order.copy(), p.ring.copy(), radix.take(), keys^, values^)

    def to_sparse(self, mut budget: Budget) raises PolyError -> SparsePoly[Self.C]:
        """The canonical sparse form."""
        var weights = self.weights()
        var n = len(self.radix)
        var result = SparsePoly[Self.C].zero(self.order, self.ring, budget)
        _charge_size(_product(len(self.keys), n + 1), budget)
        _allocation[Term[Self.C]](len(self.keys), budget)
        result.terms.reserve(len(self.keys))
        for t in range(len(self.keys)):
            result.terms.append(Term[Self.C](_monomial(self.keys[t], weights, budget), self.values[t].copy()))
        return result^

    def repacked(self, radix: List[UInt64]) -> List[UInt64]:
        """The keys under another radix with every radix_i at least this one's."""
        var old = self.weights()
        var new = _weights(radix)
        var n = len(self.radix)
        var out = List[UInt64](capacity=len(self.keys))
        for key in self.keys:
            var rest = key
            var k = UInt64(0)
            for i in range(n):
                var e = rest // old[i]
                rest -= e * old[i]
                k += e * new[i]
            out.append(k)
        return out^

    def __eq__(self, other: Self) -> Bool:
        if self.order != other.order or self.radix != other.radix or self.keys != other.keys:
            return False
        for i in range(len(self.values)):
            if not (self.values[i] == other.values[i]):
                return False
        return True


def _weights(radix: List[UInt64]) -> List[UInt64]:
    var n = len(radix)
    var weight = List[UInt64](length=n, fill=1)
    var i = n - 2
    while i >= 0:
        weight[i] = weight[i + 1] * radix[i + 1]
        i -= 1
    return weight^


def _radix_of(top: List[UInt64]) -> Optional[List[UInt64]]:
    """radix_i = top_i + 1, or None when the key range reaches 2^62."""
    var radix = List[UInt64](capacity=len(top))
    var total = UInt128(1)
    for e in top:
        total *= UInt128(e + 1)
        if total >= UInt128(_KEY_LIMIT):
            return None
        radix.append(e + 1)
    return radix^


def _monomial(key: UInt64, weights: List[UInt64], mut budget: Budget) raises PolyError -> Monomial:
    var n = len(weights)
    var count = 0
    var rest = key
    for i in range(n):
        var e = rest // weights[i]
        rest -= e * weights[i]
        count += Int(e != 0)
    _allocation[Tuple[Int64, UInt32]](count, budget)
    var pairs = List[Tuple[Int64, UInt32]](capacity=count)
    rest = key
    for i in range(n):
        var e = rest // weights[i]
        rest -= e * weights[i]
        if e:
            pairs.append((Int64(i), UInt32(e)))
    var m = Monomial()
    m.pairs = pairs^
    return m^
