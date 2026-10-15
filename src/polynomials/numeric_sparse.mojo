"""Sparse numerical polynomials with stable accumulation and certain-zero trimming."""
from std.sys import size_of
from .errors import Budget, PolyError, PolyErrorKind, _allocation, _sum, _product, _charge_size
from .monomial import Monomial, MonomialOrder
from .numeric_domains import NumericDomain, ComplexBallDomain, _numeric_cost, _numeric_zero, _numeric_nonzero, _numeric_integer, _numeric_apply, _numeric_power, _context_cost
from .numeric_dense import NumericDensePoly
from ._ordering import _IndexedOrder, _stable_indices, _sift_min
from ._assignment import _AssignmentPlan
from ._numeric_powers import _NumericPowerCache


@fieldwise_init
struct NumericTerm[D: NumericDomain](Copyable, Movable):
    """A monomial with an apn_mojo numerical coefficient.

    Parameters:
        D: The numerical domain adapter.
    """
    var monomial: Monomial
    """The canonical sparse monomial."""
    var coefficient: Self.D.Scalar
    """The finite point value or bounded coefficient enclosure."""


struct _NumericTermOrder[D: NumericDomain](_IndexedOrder):
    comptime Item = NumericTerm[Self.D]

    @staticmethod
    @always_inline
    def compare(items: List[Self.Item], a: Int, b: Int, order: MonomialOrder, mut budget: Budget) raises PolyError -> Int:
        return order._compare(items[a].monomial, items[b].monomial, budget)


def _numeric_sorted[D: NumericDomain](terms: List[NumericTerm[D]], order: MonomialOrder, mut budget: Budget) raises PolyError -> List[Int]:
    for term in terms:
        order.validate(term.monomial, budget)
        _charge_size(_numeric_cost[D](term.coefficient), budget)
    return _stable_indices[_NumericTermOrder[D]](terms, order, budget)


def _canonical_terms[D: NumericDomain](terms: List[NumericTerm[D]], order: MonomialOrder, mut budget: Budget) raises PolyError -> Bool:
    var canonical = True
    for i in range(len(terms)):
        order.validate(terms[i].monomial, budget)
        _charge_size(_numeric_cost[D](terms[i].coefficient), budget)
        if _numeric_zero[D](terms[i].coefficient):
            canonical = False
        if i and order._compare(terms[i - 1].monomial, terms[i].monomial, budget) >= 0:
            canonical = False
    return canonical


struct NumericSparsePoly[D: NumericDomain](Copyable, Movable):
    """A sparse multivariate numerical polynomial with explicit monomial order.

    Parameters:
        D: The numerical domain adapter.
    """
    var order: MonomialOrder
    """The explicit term order."""
    var terms: List[NumericTerm[Self.D]]
    """Increasing distinct monomials with coefficients not certainly zero."""

    def __init__(out self, var order: MonomialOrder):
        """Construct zero by taking ownership of its monomial order.

        Args:
            order: The validated owned order.
        """
        self.order = order^
        self.terms = List[NumericTerm[Self.D]]()

    @staticmethod
    def from_terms(var terms: List[NumericTerm[Self.D]], order: MonomialOrder, context: Self.D.Context, mut budget: Budget) raises PolyError -> Self:
        """Sort stably, combine duplicate terms and remove only certain zeros.

        Args:
            terms: The owned terms in their original accumulation order.
            order: The output monomial order.
            context: Explicit settings for duplicate coefficient sums.
            budget: The shared work allowance.

        Returns:
            The canonical numerical polynomial.

        Raises:
            PolyError: DOMAIN for invalid values or monomials, LIMIT or CAPACITY on work limits.
        """
        _ = _context_cost[Self.D](context)
        _charge_size(len(order.variables), budget)
        var result = Self(order.copy())
        var indices = _numeric_sorted(terms, order, budget)
        _allocation[NumericTerm[Self.D]](len(terms), budget)
        for index in indices:
            var term = terms[index].copy()
            if _numeric_zero[Self.D](term.coefficient):
                continue
            var last = len(result.terms) - 1
            if last >= 0 and order._compare(result.terms[last].monomial, term.monomial, budget) == 0:
                result.terms[last].coefficient = _numeric_apply[Self.D](result.terms[last].coefficient, term.coefficient, 0, context, budget)
            else:
                # Keep an intermediate zero until its duplicate group ends:
                # adding the next coefficient to zero still rounds it.
                if last >= 0 and _numeric_zero[Self.D](result.terms[last].coefficient):
                    _ = result.terms.pop()
                result.terms.append(term^)
        if len(result.terms) and _numeric_zero[Self.D](result.terms[len(result.terms) - 1].coefficient):
            _ = result.terms.pop()
        return result^

    def _same(self, other: Self, mut budget: Budget) raises PolyError:
        _charge_size(_sum(len(self.order.variables), len(other.order.variables)), budget)
        if self.order != other.order:
            raise PolyError(PolyErrorKind.DOMAIN, "different numerical monomial orders")

    def degree_bound(self, mut budget: Budget) raises PolyError -> Int:
        """Return the total-degree upper bound, or minus one for zero.

        Args:
            budget: The shared work allowance.

        Returns:
            The largest total degree among stored monomials.

        Raises:
            PolyError: CAPACITY for an unrepresentable degree or LIMIT on exhaustion.
        """
        var result = -1
        for term in self.terms:
            var degree = term.monomial.degree(budget)
            if degree > UInt64(Int.MAX):
                raise PolyError(PolyErrorKind.CAPACITY, "numerical total degree exceeds capacity")
            result = max(result, Int(degree))
        return result

    def certified_degree(self, mut budget: Budget) raises PolyError -> Optional[Int]:
        """Determine the total degree shared by every represented polynomial.

        Args:
            budget: The shared work allowance.

        Returns:
            The degree, minus one for zero, or None for uncertain top terms.

        Raises:
            PolyError: DOMAIN for invalid values, LIMIT or CAPACITY on work limits.
        """
        var degree = self.degree_bound(budget)
        if degree < 0:
            return -1
        for term in self.terms:
            _charge_size(_numeric_cost[Self.D](term.coefficient), budget)
            if term.monomial.degree(budget) == UInt64(degree) and _numeric_nonzero[Self.D](term.coefficient):
                return degree
        return None

    def _combine(self, other: Self, subtract: Bool, context: Self.D.Context, mut budget: Budget) raises PolyError -> Self:
        self._same(other, budget)
        _ = _context_cost[Self.D](context)
        var left_canonical = _canonical_terms(self.terms, self.order, budget)
        var right_canonical = _canonical_terms(other.terms, self.order, budget)
        if left_canonical and right_canonical:
            return _merge_numeric(self, other, subtract, context, budget)
        # Public fields can have been edited. Preserve the normalizing fallback.
        _allocation[NumericTerm[Self.D]](_sum(len(self.terms), len(other.terms)), budget)
        for term in self.terms:
            _charge_size(_sum(len(term.monomial.pairs), _numeric_cost[Self.D](term.coefficient)), budget)
        var values = self.terms.copy()
        var zero = _numeric_integer[Self.D](0, context, budget)
        for term in other.terms:
            _charge_size(_sum(len(term.monomial.pairs), _numeric_cost[Self.D](term.coefficient)), budget)
            var coefficient = _numeric_apply[Self.D](zero, term.coefficient, 1, context, budget) if subtract else term.coefficient.copy()
            values.append(NumericTerm[Self.D](term.monomial.copy(), coefficient^))
        return Self.from_terms(values^, self.order, context, budget)

    def add(self, other: Self, context: Self.D.Context, mut budget: Budget) raises PolyError -> Self:
        """Add terms, accumulating equal monomials left before right.

        Args:
            other: The polynomial in the same monomial order.
            context: Explicit output settings.
            budget: The shared work allowance.

        Returns:
            The rounded or enclosing sum.

        Raises:
            PolyError: DOMAIN for incompatible values, LIMIT or CAPACITY on work limits.
        """
        return self._combine(other, False, context, budget)

    def sub(self, other: Self, context: Self.D.Context, mut budget: Budget) raises PolyError -> Self:
        """Subtract the right polynomial using explicit rounded coefficient negation.

        Args:
            other: The polynomial in the same monomial order.
            context: Explicit output settings.
            budget: The shared work allowance.

        Returns:
            The rounded or enclosing difference.

        Raises:
            PolyError: DOMAIN for incompatible values, LIMIT or CAPACITY on work limits.
        """
        return self._combine(other, True, context, budget)

    def neg(self, context: Self.D.Context, mut budget: Budget) raises PolyError -> Self:
        """Negate coefficients at the specified precision.

        Args:
            context: Explicit output settings.
            budget: The shared work allowance.

        Returns:
            The rounded or enclosing negation.

        Raises:
            PolyError: DOMAIN for invalid values, LIMIT or CAPACITY on work limits.
        """
        _charge_size(len(self.order.variables), budget)
        return Self(self.order.copy())._combine(self, True, context, budget)

    def mul(self, other: Self, context: Self.D.Context, mut budget: Budget) raises PolyError -> Self:
        """Multiply with stable left-major accumulation of equal products.

        Args:
            other: The polynomial in the same monomial order.
            context: Explicit output settings.
            budget: The shared work allowance.

        Returns:
            The rounded or enclosing product.

        Raises:
            PolyError: DOMAIN for incompatible values, LIMIT or CAPACITY on work limits.
        """
        self._same(other, budget)
        var count = _product(len(self.terms), len(other.terms))
        comptime stream_threshold = 4096 if Self.D == ComplexBallDomain else 256
        if count >= stream_threshold and self.order.kind != MonomialOrder.KERNEL:
            var left_canonical = _canonical_terms(self.terms, self.order, budget)
            var right_canonical = _canonical_terms(other.terms, self.order, budget)
            if left_canonical and right_canonical:
                return _numeric_heap_product(self, other, context, budget)
        _allocation[NumericTerm[Self.D]](count, budget)
        var values = List[NumericTerm[Self.D]](capacity=count)
        for a in self.terms:
            for b in other.terms:
                values.append(NumericTerm[Self.D](a.monomial.mul(b.monomial, budget), _numeric_apply[Self.D](a.coefficient, b.coefficient, 2, context, budget)))
        return Self.from_terms(values^, self.order, context, budget)

    def scale(self, scalar: Self.D.Scalar, context: Self.D.Context, mut budget: Budget, *, divide: Bool = False) raises PolyError -> Self:
        """Scale terms, optionally dividing by a certified nonzero scalar.

        Args:
            scalar: The finite scalar or bounded enclosure.
            context: Explicit output settings.
            budget: The shared work allowance.
            divide: Select division with a nonzero-divisor proof.

        Returns:
            The rounded or enclosing scaled polynomial.

        Raises:
            PolyError: DOMAIN for invalid or uncertain divisors, LIMIT or CAPACITY on work limits.
        """
        _charge_size(_numeric_cost[Self.D](scalar), budget)
        if divide and not _numeric_nonzero[Self.D](scalar):
            raise PolyError(PolyErrorKind.DOMAIN, "scalar divisor is not certified nonzero")
        _allocation[NumericTerm[Self.D]](len(self.terms), budget)
        var values = List[NumericTerm[Self.D]](capacity=len(self.terms))
        for term in self.terms:
            _charge_size(len(term.monomial.pairs), budget)
            values.append(NumericTerm[Self.D](term.monomial.copy(), _numeric_apply[Self.D](term.coefficient, scalar, 3 if divide else 2, context, budget)))
        return Self.from_terms(values^, self.order, context, budget)

    def pow(self, exponent: UInt64, context: Self.D.Context, mut budget: Budget) raises PolyError -> Self:
        """Compute a nonnegative power by deterministic binary powering.

        Args:
            exponent: The nonnegative power.
            context: Explicit output settings.
            budget: The shared work allowance.

        Returns:
            The rounded or enclosing power.

        Raises:
            PolyError: DOMAIN for invalid values, LIMIT or CAPACITY on work limits.
        """
        var one = NumericTerm[Self.D](Monomial(), _numeric_integer[Self.D](1, context, budget))
        var result = Self.from_terms([one^], self.order, context, budget)
        _allocation[NumericTerm[Self.D]](len(self.terms), budget)
        for term in self.terms:
            _charge_size(_sum(len(term.monomial.pairs), _numeric_cost[Self.D](term.coefficient)), budget)
        var base = self.copy()
        var n = exponent
        while n:
            budget.charge(1)
            if n & 1:
                result = result.mul(base, context, budget)
            n >>= 1
            if n:
                base = base.mul(base, context, budget)
        return result^

    def derivative(self, variable: Int64, context: Self.D.Context, mut budget: Budget) raises PolyError -> Self:
        """Differentiate with respect to one variable.

        Args:
            variable: The differentiation variable.
            context: Explicit output settings.
            budget: The shared work allowance.

        Returns:
            The rounded or enclosing partial derivative.

        Raises:
            PolyError: DOMAIN for invalid values, LIMIT or CAPACITY on work limits.
        """
        _allocation[NumericTerm[Self.D]](len(self.terms), budget)
        var values = List[NumericTerm[Self.D]]()
        for term in self.terms:
            var exponent = term.monomial.exponent(variable, budget)
            if not exponent:
                continue
            var pairs = List[Tuple[Int64, UInt32]]()
            _allocation[Tuple[Int64, UInt32]](len(term.monomial.pairs), budget)
            for pair in term.monomial.pairs:
                if pair[0] != variable:
                    pairs.append(pair)
                elif pair[1] > 1:
                    pairs.append((variable, pair[1] - 1))
            var factor = _numeric_integer[Self.D](Int(UInt64(exponent)), context, budget)
            values.append(NumericTerm[Self.D](Monomial(pairs^, budget), _numeric_apply[Self.D](term.coefficient, factor, 2, context, budget)))
        return Self.from_terms(values^, self.order, context, budget)

    def evaluate(self, assignments: List[Tuple[Int64, Self.D.Scalar]], context: Self.D.Context, mut budget: Budget) raises PolyError -> Self.D.Scalar:
        """Evaluate in canonical term order with cached variable powers.

        Args:
            assignments: Distinct variable/value pairs, allowing unused variables.
            context: Explicit output settings.
            budget: The shared work allowance.

        Returns:
            The rounded value or rigorous enclosure.

        Raises:
            PolyError: DOMAIN for missing or repeated variables or invalid values, LIMIT or CAPACITY on work limits.
        """
        if len(self.terms) >= 256 and len(assignments) <= len(self.terms) // 8:
            return self._evaluate[True](assignments, context, budget)
        return self._evaluate[False](assignments, context, budget)

    def _evaluate[shared: Bool](self, assignments: List[Tuple[Int64, Self.D.Scalar]], context: Self.D.Context, mut budget: Budget) raises PolyError -> Self.D.Scalar:
        var plan = _AssignmentPlan()
        var squares = _NumericPowerCache[Self.D]()
        var capacity = min(len(self.terms), 16)
        plan.reserve_small(len(assignments), capacity, budget)
        for i in range(len(assignments)):
            _charge_size(_numeric_cost[Self.D](assignments[i][1]), budget)
            plan.assign(assignments[i][0], i, budget)
        _allocation[Self.D.Scalar](capacity, budget)
        var powers = List[Self.D.Scalar](capacity=capacity)
        var result = _numeric_integer[Self.D](0, context, budget)
        for term in self.terms:
            _charge_size(_numeric_cost[Self.D](term.coefficient), budget)
            var value = term.coefficient.copy()
            for pair in term.monomial.pairs:
                var slot = plan.power(pair, budget)
                if slot == len(powers):
                    var assigned = plan.assignment(pair[0], budget)
                    if assigned < 0:
                        raise PolyError(PolyErrorKind.DOMAIN, "missing numerical assignment")
                    _allocation[Self.D.Scalar](1, budget)
                    comptime if shared:
                        powers.append(squares.power(assigned, assignments[assigned][1], UInt64(pair[1]), context, budget))
                    else:
                        powers.append(_numeric_power[Self.D](assignments[assigned][1], UInt64(pair[1]), context, budget))
                value = _numeric_apply[Self.D](value, powers[slot], 2, context, budget)
            result = _numeric_apply[Self.D](result, value, 0, context, budget)
        return result^

    @staticmethod
    def from_dense(poly: NumericDensePoly[Self.D], order: MonomialOrder, context: Self.D.Context, mut budget: Budget) raises PolyError -> Self:
        """Convert a dense univariate polynomial into sparse form.

        Args:
            poly: The dense polynomial.
            order: The output monomial order.
            context: Explicit arithmetic settings.
            budget: The shared work allowance.

        Returns:
            A value-preserving sparse representation.

        Raises:
            PolyError: DOMAIN for invalid orders, LIMIT or CAPACITY on work limits.
        """
        if UInt64(max(0, poly.degree_bound())) > UInt64(UInt32.MAX):
            raise PolyError(PolyErrorKind.CAPACITY, "dense degree exceeds sparse exponent capacity")
        _allocation[NumericTerm[Self.D]](len(poly.coefficients), budget)
        var values = List[NumericTerm[Self.D]]()
        for i in range(len(poly.coefficients)):
            var monomial = Monomial()
            if i:
                monomial = Monomial([(poly.variable, UInt32(i))], budget)
            _charge_size(_numeric_cost[Self.D](poly.coefficients[i]), budget)
            values.append(NumericTerm[Self.D](monomial^, poly.coefficients[i].copy()))
        return Self.from_terms(values^, order, context, budget)

    def to_dense(self, variable: Int64, context: Self.D.Context, mut budget: Budget) raises PolyError -> NumericDensePoly[Self.D]:
        """Convert a polynomial in one specified variable into dense form.

        Args:
            variable: The retained variable, including for constants and zero.
            context: Explicit settings for gap zeros.
            budget: The shared work allowance.

        Returns:
            A value-preserving dense representation.

        Raises:
            PolyError: DOMAIN for other variables, LIMIT or CAPACITY before dense allocation.
        """
        var count = 0
        for term in self.terms:
            term.monomial.check_invariants(budget)
            for pair in term.monomial.pairs:
                if pair[0] != variable:
                    raise PolyError(PolyErrorKind.DOMAIN, "multivariate polynomial cannot become dense")
            count = max(count, _sum(Int(UInt64(term.monomial.exponent(variable, budget))), 1))
        _allocation[Self.D.Scalar](count, budget)
        var values = List[Self.D.Scalar](length=count, fill=_numeric_integer[Self.D](0, context, budget))
        for term in self.terms:
            _charge_size(_numeric_cost[Self.D](term.coefficient), budget)
            values[Int(UInt64(term.monomial.exponent(variable, budget)))] = term.coefficient.copy()
        return NumericDensePoly[Self.D].from_coefficients(variable, values^, budget)


def _merge_numeric[D: NumericDomain](a: NumericSparsePoly[D], b: NumericSparsePoly[D], subtract: Bool, context: D.Context, mut budget: Budget) raises PolyError -> NumericSparsePoly[D]:
    _charge_size(len(a.order.variables), budget)
    var result = NumericSparsePoly[D](a.order.copy())
    _allocation[NumericTerm[D]](_sum(len(a.terms), len(b.terms)), budget)
    result.terms.reserve(_sum(len(a.terms), len(b.terms)))
    var zero = _numeric_integer[D](0, context, budget)
    var i = 0
    var j = 0
    while i < len(a.terms) or j < len(b.terms):
        budget.charge(1)
        var relation = -1 if j == len(b.terms) else (1 if i == len(a.terms) else a.order._compare(a.terms[i].monomial, b.terms[j].monomial, budget))
        if relation < 0:
            _charge_size(_sum(len(a.terms[i].monomial.pairs), _numeric_cost[D](a.terms[i].coefficient)), budget)
            result.terms.append(a.terms[i].copy())
            i += 1
        else:
            _charge_size(_sum(len(b.terms[j].monomial.pairs), _numeric_cost[D](b.terms[j].coefficient)), budget)
            var coefficient = _numeric_apply[D](zero, b.terms[j].coefficient, 1, context, budget) if subtract else b.terms[j].coefficient.copy()
            # from_terms ignores an individually rounded zero before combining.
            if relation == 0:
                if not _numeric_zero[D](coefficient):
                    coefficient = _numeric_apply[D](a.terms[i].coefficient, coefficient, 0, context, budget)
                else:
                    _charge_size(_numeric_cost[D](a.terms[i].coefficient), budget)
                    coefficient = a.terms[i].coefficient.copy()
                i += 1
            if not _numeric_zero[D](coefficient):
                result.terms.append(NumericTerm[D](b.terms[j].monomial.copy(), coefficient^))
            j += 1
    return result^


@fieldwise_init
struct _NumericProduct[D: NumericDomain](Copyable, Movable):
    var row: Int
    var column: Int
    var monomial: Monomial
    # One owning scalar slot avoids associated-field partial-move issues in
    # Mojo 1.1 and keeps comparisons independent of coefficient ownership.
    var coefficient: List[Self.D.Scalar]


struct _NumericProductOrder[D: NumericDomain](_IndexedOrder):
    comptime Item = _NumericProduct[Self.D]

    @staticmethod
    @always_inline
    def compare(items: List[Self.Item], a: Int, b: Int, order: MonomialOrder, mut budget: Budget) raises PolyError -> Int:
        var relation = order._compare(items[a].monomial, items[b].monomial, budget)
        if relation:
            return relation
        return Int(items[a].row > items[b].row) - Int(items[a].row < items[b].row)


def _numeric_heap_product[D: NumericDomain](a: NumericSparsePoly[D], b: NumericSparsePoly[D], context: D.Context, mut budget: Budget) raises PolyError -> NumericSparsePoly[D]:
    _ = _context_cost[D](context)
    _charge_size(len(a.order.variables), budget)
    var result = NumericSparsePoly[D](a.order.copy())
    _allocation[_NumericProduct[D]](len(a.terms), budget)
    _allocation[Int](len(a.terms), budget)
    var cursors = List[_NumericProduct[D]](capacity=len(a.terms))
    var heap = List[Int](capacity=len(a.terms))
    # Never exchange the operands: row order breaks monomial ties exactly as
    # the original left-major product list followed by a stable sort did.
    for i in range(len(a.terms)):
        var monomial = a.terms[i].monomial.mul(b.terms[0].monomial, budget)
        var coefficient = _numeric_apply[D](a.terms[i].coefficient, b.terms[0].coefficient, 2, context, budget)
        _allocation[D.Scalar](1, budget)
        cursors.append(_NumericProduct[D](i, 0, monomial^, [coefficient^]))
        heap.append(i)
    for i in range(len(heap) // 2 - 1, -1, -1):
        budget.charge(1)
        _sift_min[_NumericProductOrder[D]](heap, cursors, i, a.order, budget)
    while len(heap):
        budget.charge(1)
        var index = heap[0]
        ref cursor = cursors[index]
        var monomial = Monomial()
        swap(monomial, cursor.monomial)
        var term = NumericTerm[D](monomial^, cursor.coefficient.pop())
        if not _numeric_zero[D](term.coefficient):
            var last = len(result.terms) - 1
            if last >= 0 and a.order._compare(result.terms[last].monomial, term.monomial, budget) == 0:
                result.terms[last].coefficient = _numeric_apply[D](result.terms[last].coefficient, term.coefficient, 0, context, budget)
            else:
                # Retain a cancelled group until its final duplicate has been
                # accumulated: adding to rounded zero is still an operation.
                if last >= 0 and _numeric_zero[D](result.terms[last].coefficient):
                    _ = result.terms.pop()
                _ = _product(_sum(len(result.terms), 1), size_of[NumericTerm[D]]())
                _allocation[NumericTerm[D]](1, budget)
                result.terms.append(term^)
        cursor.column += 1
        if cursor.column < len(b.terms):
            ref left = a.terms[cursor.row]
            ref right = b.terms[cursor.column]
            cursor.monomial = left.monomial.mul(right.monomial, budget)
            cursor.coefficient.append(_numeric_apply[D](left.coefficient, right.coefficient, 2, context, budget))
        else:
            var last_index = heap.pop()
            if len(heap):
                heap[0] = last_index
        if len(heap):
            _sift_min[_NumericProductOrder[D]](heap, cursors, 0, a.order, budget)
    if len(result.terms) and _numeric_zero[D](result.terms[len(result.terms) - 1].coefficient):
        _ = result.terms.pop()
    return result^
