"""Coefficient-first Newton, Chebyshev, and Bernstein representations."""
from ._io_meta import _Record, _Basis
from std.iter import Iterable, Iterator, StopIteration
from std.builtin.rebind import downcast
from std.memory import ArcPointer
from polynomials import Budget, PolyError, PolyErrorKind, CoefficientField, DensePoly, Modulus, ModInt
from . import dense
from ._common import _Kind, _CoreKind, _Element, _Scalar, _Ring, _FieldKind, _Join, _WithModulus, _unwrap, _wrap, _convert, _allocation, _compact, _sum, _domain, _sequence, _sequence_coefficients, _data_context, _values_from_core, Values

comptime _Field[k: Int] = downcast[_Ring[k], CoefficientField]
from polynomials import BasisPolynomial


@fieldwise_init
struct Array[k: Int](Copyable, Movable, Iterable, Sized, _WithModulus, _Record):
    """Immutable basis coefficients with retained representation metadata.

    Parameters:
        k: The inferred exact field: QQ, QQi, GF(p), AA, or QQbar.
    """
    comptime _Format = _Basis
    comptime _Item = _Scalar[Self.k]
    comptime _Power = Int
    var _core: ArcPointer[BasisPolynomial[_Field[Self.k]]]
    comptime IteratorType[iterable_mut: Bool, //, iterable_origin: Origin[mut=iterable_mut]]: Iterator = _BasisIterator[Self.k, ImmOrigin(iterable_origin)]
    """A borrowing coefficient iterator.

    Parameters:
        iterable_mut: Whether the source origin permits mutation.
        iterable_origin: The lifetime of the source.
    """

    def __iter__(ref self) -> Self.IteratorType[origin_of(self)]:
        """Iterate in ascending coefficient order.

        Returns:
            A borrowing iterator over all retained coefficients.
        """
        return _BasisIterator[Self.k, ImmOrigin(origin_of(self))](Pointer(to=self).as_imm(), 0)

    def __len__(self) -> Int:
        """Return the coefficient extent, including trailing zeros.

        Returns:
            The number of retained coefficients.
        """
        return len(self._core[].coefficients)

    def __getitem__(self, index: Int) -> _Scalar[Self.k]:
        """Read a coefficient at a valid zero-based index.

        Args:
            index: An index from zero through len(self)-1.

        Returns:
            An independent exact scalar.
        """
        return _unwrap[Self.k](rebind[_Ring[Self.k]](self._core[].coefficients[index]))

    def modulus(self) -> UInt64:
        """Return the retained prime, or zero for QQ and QQi.

        Returns:
            The prime modulus or zero.
        """
        comptime if Self.k == 3:
            return rebind[Modulus](self._core[].ring).value
        return 0

    def domain(self) -> String:
        """Return the retained coefficient field.

        Returns:
            QQ, QQi, or GF(p).
        """
        return _domain[Self.k](self.modulus())

    def basis(self) -> String:
        """Return the name of the retained basis.

        Returns:
            The name newton, chebyshev, or bernstein.
        """
        if self._core[].kind == 0:
            return "newton"
        if self._core[].kind == 1:
            return "chebyshev"
        return "bernstein"


@fieldwise_init
struct _BasisIterator[k: Int, origin: ImmOrigin](Iterator):
    comptime Element = _Scalar[Self.k]
    var source: Pointer[Array[Self.k], Self.origin]
    var position: Int

    def __next__(mut self) raises StopIteration -> Self.Element:
        if self.position == len(self.source[]):
            raise StopIteration()
        var value = self.source[][self.position]
        self.position += 1
        return value^

    def bounds(self) -> Tuple[Int, Optional[Int]]:
        var remaining = len(self.source[]) - self.position
        return remaining, Optional(remaining)


def _kind(name: String) raises PolyError -> UInt8:
    if name == "newton":
        return 0
    if name == "chebyshev":
        return 1
    if name == "bernstein":
        return 2
    raise PolyError(PolyErrorKind.DOMAIN, "basis must be newton, chebyshev, or bernstein")


def _shape(kind: UInt8, count: Int, nodes: Int) raises PolyError:
    if (kind == 0 and nodes != max(0, count - 1)) or (kind != 0 and nodes):
        raise PolyError(PolyErrorKind.DOMAIN, "Newton coefficients require exactly one fewer nodes; other bases take no nodes")
    if kind == 2 and count == 0:
        raise PolyError(PolyErrorKind.DOMAIN, "Bernstein coefficients must retain at least degree zero")


def _wrap_basis[k: Int](var core: BasisPolynomial[_Field[k]], mut budget: Budget) raises PolyError -> Array[k]:
    _shape(core.kind, len(core.coefficients), len(core.nodes))
    _compact(core.coefficients, budget)
    _compact(core.nodes, budget)
    _allocation[BasisPolynomial[_Field[k]]](1, budget)
    return Array[k](ArcPointer(core^))


def _copy_values[k: Int, source: Int](values: List[_Field[source]], ctx: _Ring[k].Context, mut budget: Budget) raises PolyError -> List[_Field[k]]:
    _allocation[_Field[k]](len(values), budget)
    var result = List[_Field[k]](capacity=len(values))
    for c in values:
        result.append(rebind_var[_Field[k]](_wrap[k](_convert[k](_unwrap[source](rebind[_Ring[source]](c)), ctx, budget))))
    return result^


def _promote[k: Int, source: Int](values: Array[source], ctx: _Ring[k].Context, mut budget: Budget) raises PolyError -> Array[k]:
    comptime if k == source:
        budget.charge(1)
        return rebind[Array[k]](values).copy()
    else:
        ref p = values._core[]
        return _wrap_basis[k](BasisPolynomial[_Field[k]](0, p.kind, rebind[_Field[k].Context](ctx).copy(), _copy_values[k, source](p.coefficients, ctx, budget), _copy_values[k, source](p.nodes, ctx, budget)), budget)


def asarray[I: Iterable, N: Iterable](values: I, mut budget: Budget, *, basis: String, nodes: N) raises PolyError -> Array[_FieldKind[_Join[_Kind[_Element[I]], _Kind[_Element[N]]]]]:
    """Interpret the supplied values as basis coefficients, preserving every trailing zero.

    Parameters:
        I: The coefficient iterable type.
        N: The node iterable type.

    Args:
        values: Ascending coefficients in the named basis.
        budget: The shared work allowance, including preparation and result storage.
        basis: The name newton, chebyshev, or bernstein.
        nodes: Newton nodes; empty for other bases.

    Returns:
        An immutable basis representation; integer inputs promote to QQ.

    Raises:
        PolyError: DOMAIN for incompatible inputs or unavailable divisions, LIMIT for exhausted work, or CAPACITY for unrepresentable storage.
    """
    comptime k = _FieldKind[_Join[_Kind[_Element[I]], _Kind[_Element[N]]]]
    var kind = _kind(basis)
    var coefficients = _sequence(values, budget)
    var sites = _sequence(nodes, budget)
    var ctx = _data_context[k](coefficients.modulus(), sites.modulus(), budget)
    var c = rebind_var[List[_Field[k]]](_sequence_coefficients[k](coefficients, ctx, budget))
    var n = rebind_var[List[_Field[k]]](_sequence_coefficients[k](sites, ctx, budget))
    _shape(kind, len(c), len(n))
    return _wrap_basis[k](BasisPolynomial[_Field[k]](0, kind, rebind[_Field[k].Context](ctx).copy(), c^, n^), budget)


def asarray[I: Iterable, N: Iterable](values: I, *, max_steps: UInt64 = UInt64.MAX, basis: String, nodes: N) raises PolyError -> Array[_FieldKind[_Join[_Kind[_Element[I]], _Kind[_Element[N]]]]]:
    """Interpret the supplied values as basis coefficients, preserving every trailing zero.

    Parameters:
        I: The coefficient iterable type.
        N: The node iterable type.

    Args:
        values: Ascending coefficients in the named basis.
        max_steps: The complete call's work allowance; defaults to UInt64.MAX.
        basis: The name newton, chebyshev, or bernstein.
        nodes: Newton nodes; empty for other bases.

    Returns:
        An immutable basis representation; integer inputs promote to QQ.

    Raises:
        PolyError: DOMAIN for incompatible inputs or unavailable divisions, LIMIT for exhausted work, or CAPACITY for unrepresentable storage.
    """
    var budget = Budget(max_steps)
    return asarray(values, budget, basis=basis, nodes=nodes)


def asarray[I: Iterable](values: I, mut budget: Budget, *, basis: String) raises PolyError -> Array[_FieldKind[_Kind[_Element[I]]]]:
    """Interpret the supplied values as basis coefficients, preserving every trailing zero.

    Parameters:
        I: The coefficient iterable type.

    Args:
        values: The coefficient iterable.
        budget: The shared work allowance, including preparation and result storage.
        basis: The name newton, chebyshev, or bernstein.

    Returns:
        An immutable basis representation; integer inputs promote to QQ.

    Raises:
        PolyError: DOMAIN for incompatible inputs or unavailable divisions, LIMIT for exhausted work, or CAPACITY for unrepresentable storage.
    """
    comptime assert _Kind[_Element[I]] >= 0, "basis coefficients must be exact numbers"
    return rebind_var[Array[_FieldKind[_Kind[_Element[I]]]]](asarray(values, budget, basis=basis, nodes=List[Int]()))


def asarray[I: Iterable](values: I, *, max_steps: UInt64 = UInt64.MAX, basis: String) raises PolyError -> Array[_FieldKind[_Kind[_Element[I]]]]:
    """Interpret the supplied values as basis coefficients, preserving every trailing zero.

    Parameters:
        I: The coefficient iterable type.

    Args:
        values: The coefficient iterable.
        max_steps: The complete call's work allowance; defaults to UInt64.MAX.
        basis: The name newton, chebyshev, or bernstein.

    Returns:
        An immutable basis representation; integer inputs promote to QQ.

    Raises:
        PolyError: DOMAIN for incompatible inputs or unavailable divisions, LIMIT for exhausted work, or CAPACITY for unrepresentable storage.
    """
    var budget = Budget(max_steps)
    return asarray(values, budget, basis=basis)


def from_polynomial[I: Iterable, N: Iterable](values: I, mut budget: Budget, *, basis: String, nodes: N, degree: Int = -1) raises PolyError -> Array[_FieldKind[_Join[_Kind[_Element[I]], _Kind[_Element[N]]]]]:
    """Convert ascending monomial coefficients to the requested basis.

    Parameters:
        I: The coefficient iterable type.
        N: The node iterable type.

    Args:
        values: Ascending monomial coefficients.
        budget: The shared work allowance, including preparation and result storage.
        basis: The name newton, chebyshev, or bernstein.
        nodes: Newton nodes; empty for other bases.
        degree: The retained Bernstein degree; -1 chooses the polynomial degree, at least zero.

    Returns:
        An immutable basis representation; integer inputs promote to QQ.

    Raises:
        PolyError: DOMAIN for incompatible inputs or unavailable divisions, LIMIT for exhausted work, or CAPACITY for unrepresentable storage.
    """
    comptime k = _FieldKind[_Join[_Kind[_Element[I]], _Kind[_Element[N]]]]
    var kind = _kind(basis)
    var coefficients = _sequence(values, budget)
    var sites = _sequence(nodes, budget)
    var ctx = _data_context[k](coefficients.modulus(), sites.modulus(), budget)
    var c = rebind_var[List[_Field[k]]](_sequence_coefficients[k](coefficients, ctx, budget))
    var n = rebind_var[List[_Field[k]]](_sequence_coefficients[k](sites, ctx, budget))
    var p = DensePoly[_Field[k]].from_coefficients(0, rebind[_Field[k].Context](ctx).copy(), c^, budget)
    if kind == 0:
        if degree != -1:
            raise PolyError(PolyErrorKind.DOMAIN, "degree= is only for Bernstein conversion")
        return _wrap_basis[k](BasisPolynomial[_Field[k]].newton(p, n, budget), budget)
    if len(n):
        raise PolyError(PolyErrorKind.DOMAIN, "nodes= is only for Newton conversion")
    if kind == 1:
        if degree != -1:
            raise PolyError(PolyErrorKind.DOMAIN, "degree= is only for Bernstein conversion")
        return _wrap_basis[k](BasisPolynomial[_Field[k]].chebyshev(p, budget), budget)
    if degree < -1:
        raise PolyError(PolyErrorKind.DOMAIN, "invalid Bernstein degree")
    return _wrap_basis[k](BasisPolynomial[_Field[k]].bernstein(p, max(0, p.degree()) if degree == -1 else degree, budget), budget)


def from_polynomial[I: Iterable, N: Iterable](values: I, *, max_steps: UInt64 = UInt64.MAX, basis: String, nodes: N, degree: Int = -1) raises PolyError -> Array[_FieldKind[_Join[_Kind[_Element[I]], _Kind[_Element[N]]]]]:
    """Convert ascending monomial coefficients to the requested basis.

    Parameters:
        I: The coefficient iterable type.
        N: The node iterable type.

    Args:
        values: Ascending monomial coefficients.
        max_steps: The complete call's work allowance; defaults to UInt64.MAX.
        basis: The name newton, chebyshev, or bernstein.
        nodes: Newton nodes; empty for other bases.
        degree: The retained Bernstein degree; -1 chooses the polynomial degree, at least zero.

    Returns:
        An immutable basis representation; integer inputs promote to QQ.

    Raises:
        PolyError: DOMAIN for incompatible inputs or unavailable divisions, LIMIT for exhausted work, or CAPACITY for unrepresentable storage.
    """
    var budget = Budget(max_steps)
    return from_polynomial(values, budget, basis=basis, nodes=nodes, degree=degree)


def from_polynomial[I: Iterable](values: I, mut budget: Budget, *, basis: String, degree: Int = -1) raises PolyError -> Array[_FieldKind[_Kind[_Element[I]]]]:
    """Convert ascending monomial coefficients to the requested basis.

    Parameters:
        I: The coefficient iterable type.

    Args:
        values: The coefficient iterable.
        budget: The shared work allowance, including preparation and result storage.
        basis: The name newton, chebyshev, or bernstein.
        degree: The retained Bernstein degree; -1 chooses the polynomial degree, at least zero.

    Returns:
        An immutable basis representation; integer inputs promote to QQ.

    Raises:
        PolyError: DOMAIN for incompatible inputs or unavailable divisions, LIMIT for exhausted work, or CAPACITY for unrepresentable storage.
    """
    comptime assert _Kind[_Element[I]] >= 0, "basis coefficients must be exact numbers"
    return rebind_var[Array[_FieldKind[_Kind[_Element[I]]]]](from_polynomial(values, budget, basis=basis, nodes=List[Int](), degree=degree))


def from_polynomial[I: Iterable](values: I, *, max_steps: UInt64 = UInt64.MAX, basis: String, degree: Int = -1) raises PolyError -> Array[_FieldKind[_Kind[_Element[I]]]]:
    """Convert ascending monomial coefficients to the requested basis.

    Parameters:
        I: The coefficient iterable type.

    Args:
        values: The coefficient iterable.
        max_steps: The complete call's work allowance; defaults to UInt64.MAX.
        basis: The name newton, chebyshev, or bernstein.
        degree: The retained Bernstein degree; -1 chooses the polynomial degree, at least zero.

    Returns:
        An immutable basis representation; integer inputs promote to QQ.

    Raises:
        PolyError: DOMAIN for incompatible inputs or unavailable divisions, LIMIT for exhausted work, or CAPACITY for unrepresentable storage.
    """
    var budget = Budget(max_steps)
    return from_polynomial(values, budget, basis=basis, degree=degree)


def to_polynomial[k: Int](values: Array[k], mut budget: Budget) raises PolyError -> dense.Array[k]:
    """Expand a basis representation into ascending monomial coefficients.

    Parameters:
        k: The inferred coefficient field.

    Args:
        values: The prepared basis coefficients.
        budget: The shared work allowance, including preparation and result storage.

    Returns:
        A canonical dense polynomial in the same field.

    Raises:
        PolyError: DOMAIN for incompatible inputs or unavailable divisions, LIMIT for exhausted work, or CAPACITY for unrepresentable storage.
    """
    return dense.Array[k](rebind_var[DensePoly[_Ring[k]]](values._core[].to_polynomial(budget)), budget)


def to_polynomial[k: Int](values: Array[k], *, max_steps: UInt64 = UInt64.MAX) raises PolyError -> dense.Array[k]:
    """Expand a basis representation into ascending monomial coefficients.

    Parameters:
        k: The inferred coefficient field.

    Args:
        values: The prepared basis coefficients.
        max_steps: The complete call's work allowance; defaults to UInt64.MAX.

    Returns:
        A canonical dense polynomial in the same field.

    Raises:
        PolyError: DOMAIN for incompatible inputs or unavailable divisions, LIMIT for exhausted work, or CAPACITY for unrepresentable storage.
    """
    var budget = Budget(max_steps)
    return to_polynomial(values, budget)


def coefficients[k: Int](values: Array[k], mut budget: Budget) raises PolyError -> Values[_Scalar[k], k]:
    """Copy the retained coefficients without dropping trailing zeros.

    Parameters:
        k: The inferred coefficient field.

    Args:
        values: The prepared basis representation.
        budget: The shared work allowance, including preparation and result storage.

    Returns:
        An immutable sequence retaining its field metadata.

    Raises:
        PolyError: DOMAIN for incompatible inputs or unavailable divisions, LIMIT for exhausted work, or CAPACITY for unrepresentable storage.
    """
    var copied = _copy_values[k, k](values._core[].coefficients, rebind[_Ring[k].Context](values._core[].ring), budget)
    return _values_from_core[k](rebind_var[List[_Ring[k]]](copied^), rebind[_Ring[k].Context](values._core[].ring), budget)


def coefficients[k: Int](values: Array[k], *, max_steps: UInt64 = UInt64.MAX) raises PolyError -> Values[_Scalar[k], k]:
    """Copy the retained coefficients without dropping trailing zeros.

    Parameters:
        k: The inferred coefficient field.

    Args:
        values: The prepared basis representation.
        max_steps: The complete call's work allowance; defaults to UInt64.MAX.

    Returns:
        An immutable sequence retaining its field metadata.

    Raises:
        PolyError: DOMAIN for incompatible inputs or unavailable divisions, LIMIT for exhausted work, or CAPACITY for unrepresentable storage.
    """
    var budget = Budget(max_steps)
    return coefficients(values, budget)


def nodes[k: Int](values: Array[k], mut budget: Budget) raises PolyError -> Values[_Scalar[k], k]:
    """Copy the retained nodes without dropping trailing zeros.

    Parameters:
        k: The inferred coefficient field.

    Args:
        values: The prepared basis representation.
        budget: The shared work allowance, including preparation and result storage.

    Returns:
        An immutable sequence retaining its field metadata.

    Raises:
        PolyError: DOMAIN for incompatible inputs or unavailable divisions, LIMIT for exhausted work, or CAPACITY for unrepresentable storage.
    """
    var copied = _copy_values[k, k](values._core[].nodes, rebind[_Ring[k].Context](values._core[].ring), budget)
    return _values_from_core[k](rebind_var[List[_Ring[k]]](copied^), rebind[_Ring[k].Context](values._core[].ring), budget)


def nodes[k: Int](values: Array[k], *, max_steps: UInt64 = UInt64.MAX) raises PolyError -> Values[_Scalar[k], k]:
    """Copy the retained nodes without dropping trailing zeros.

    Parameters:
        k: The inferred coefficient field.

    Args:
        values: The prepared basis representation.
        max_steps: The complete call's work allowance; defaults to UInt64.MAX.

    Returns:
        An immutable sequence retaining its field metadata.

    Raises:
        PolyError: DOMAIN for incompatible inputs or unavailable divisions, LIMIT for exhausted work, or CAPACITY for unrepresentable storage.
    """
    var budget = Budget(max_steps)
    return nodes(values, budget)


def evaluate[k: Int, T: Copyable & Deinitable](values: Array[k], point: T, mut budget: Budget) raises PolyError -> _Scalar[_Join[k, _Kind[T]]]:
    """Evaluate directly in the retained basis, promoting the point and coefficients together.

    Parameters:
        k: The inferred coefficient field.
        T: The exact point type.

    Args:
        values: The prepared basis representation.
        point: An exact scalar in a compatible domain.
        budget: The shared work allowance, including preparation and result storage.

    Returns:
        The exact value from nested Newton, Clenshaw, or de Casteljau evaluation.

    Raises:
        PolyError: DOMAIN for incompatible inputs or unavailable divisions, LIMIT for exhausted work, or CAPACITY for unrepresentable storage.
    """
    comptime target = _Join[k, _Kind[T]]
    var modulus = UInt64(0)
    comptime if T == ModInt:
        modulus = rebind[ModInt](point).modulus
    var ctx = _data_context[target](values.modulus(), modulus, budget)
    var p = _promote[target](values, ctx, budget)
    var x = rebind_var[_Field[target]](_wrap[target](_convert[target](point, ctx, budget)))
    return _unwrap[target](rebind_var[_Ring[target]](p._core[].evaluate(x, budget)))


def evaluate[k: Int, T: Copyable & Deinitable](values: Array[k], point: T, *, max_steps: UInt64 = UInt64.MAX) raises PolyError -> _Scalar[_Join[k, _Kind[T]]]:
    """Evaluate directly in the retained basis, promoting the point and coefficients together.

    Parameters:
        k: The inferred coefficient field.
        T: The exact point type.

    Args:
        values: The prepared basis representation.
        point: An exact scalar in a compatible domain.
        max_steps: The complete call's work allowance; defaults to UInt64.MAX.

    Returns:
        The exact value from nested Newton, Clenshaw, or de Casteljau evaluation.

    Raises:
        PolyError: DOMAIN for incompatible inputs or unavailable divisions, LIMIT for exhausted work, or CAPACITY for unrepresentable storage.
    """
    var budget = Budget(max_steps)
    return evaluate(values, point, budget)


def subdivide[k: Int, T: Copyable & Deinitable](values: Array[k], point: T, mut budget: Budget) raises PolyError -> Tuple[Array[_Join[k, _Kind[T]]], Array[_Join[k, _Kind[T]]]]:
    """Subdivide a Bernstein representation at an exact parameter.

    Parameters:
        k: The inferred coefficient field.
        T: The exact parameter type.

    Args:
        values: Bernstein coefficients.
        point: The split parameter; values outside [0,1] give affine extrapolations.
        budget: The shared work allowance, including preparation and result storage.

    Returns:
        Representations of p(point*x) and p(point+(1-point)*x), retaining degree.

    Raises:
        PolyError: DOMAIN for incompatible inputs or unavailable divisions, LIMIT for exhausted work, or CAPACITY for unrepresentable storage.
    """
    comptime target = _Join[k, _Kind[T]]
    var modulus = UInt64(0)
    comptime if T == ModInt:
        modulus = rebind[ModInt](point).modulus
    var ctx = _data_context[target](values.modulus(), modulus, budget)
    var p = _promote[target](values, ctx, budget)
    var x = rebind_var[_Field[target]](_wrap[target](_convert[target](point, ctx, budget)))
    var split = p._core[].subdivide(x, budget)
    var left = BasisPolynomial[_Field[target]](0, 2, x.context(), List[_Field[target]](), List[_Field[target]]())
    var right = BasisPolynomial[_Field[target]](0, 2, x.context(), List[_Field[target]](), List[_Field[target]]())
    swap(left, split[0])
    swap(right, split[1])
    return _wrap_basis[target](left^, budget), _wrap_basis[target](right^, budget)


def subdivide[k: Int, T: Copyable & Deinitable](values: Array[k], point: T, *, max_steps: UInt64 = UInt64.MAX) raises PolyError -> Tuple[Array[_Join[k, _Kind[T]]], Array[_Join[k, _Kind[T]]]]:
    """Subdivide a Bernstein representation at an exact parameter.

    Parameters:
        k: The inferred coefficient field.
        T: The exact parameter type.

    Args:
        values: Bernstein coefficients.
        point: The split parameter; values outside [0,1] give affine extrapolations.
        max_steps: The complete call's work allowance; defaults to UInt64.MAX.

    Returns:
        Representations of p(point*x) and p(point+(1-point)*x), retaining degree.

    Raises:
        PolyError: DOMAIN for incompatible inputs or unavailable divisions, LIMIT for exhausted work, or CAPACITY for unrepresentable storage.
    """
    var budget = Budget(max_steps)
    return subdivide(values, point, budget)


def to_core[k: Int](values: Array[k], mut budget: Budget, *, variable: Int64 = 0) raises PolyError -> BasisPolynomial[_Field[k]]:
    """Export independent mutable core coefficients and nodes.

    Parameters:
        k: The inferred coefficient field.

    Args:
        values: The prepared representation.
        budget: The shared work allowance, including preparation and result storage.
        variable: The core variable ID assigned to the exported value.

    Returns:
        An independent core basis representation.

    Raises:
        PolyError: DOMAIN for incompatible inputs or unavailable divisions, LIMIT for exhausted work, or CAPACITY for unrepresentable storage.
    """
    ref p = values._core[]
    var ctx = rebind[_Ring[k].Context](p.ring).copy()
    return BasisPolynomial[_Field[k]](variable, p.kind, p.ring.copy(), _copy_values[k, k](p.coefficients, ctx, budget), _copy_values[k, k](p.nodes, ctx, budget))


def to_core[k: Int](values: Array[k], *, max_steps: UInt64 = UInt64.MAX, variable: Int64 = 0) raises PolyError -> BasisPolynomial[_Field[k]]:
    """Export independent mutable core coefficients and nodes.

    Parameters:
        k: The inferred coefficient field.

    Args:
        values: The prepared representation.
        max_steps: The complete call's work allowance; defaults to UInt64.MAX.
        variable: The core variable ID assigned to the exported value.

    Returns:
        An independent core basis representation.

    Raises:
        PolyError: DOMAIN for incompatible inputs or unavailable divisions, LIMIT for exhausted work, or CAPACITY for unrepresentable storage.
    """
    var budget = Budget(max_steps)
    return to_core(values, budget, variable=variable)


def from_core[C: CoefficientField](values: BasisPolynomial[C], mut budget: Budget) raises PolyError -> Array[_CoreKind[C]]:
    """Validate and copy a core basis representation, preserving coefficient extent.

    Parameters:
        C: The inferred supported exact coefficient field.

    Args:
        values: A mutable core basis representation.
        budget: The shared work allowance, including preparation and result storage.

    Returns:
        An independent representation using array axis zero.

    Raises:
        PolyError: DOMAIN for incompatible inputs or unavailable divisions, LIMIT for exhausted work, or CAPACITY for unrepresentable storage.
    """
    comptime k = _CoreKind[C]
    comptime assert k >= 1, "basis arrays support QQ, QQi, and ModInt"
    if values.kind > 2:
        raise PolyError(PolyErrorKind.DOMAIN, "invalid core basis tag")
    _shape(values.kind, len(values.coefficients), len(values.nodes))
    var modulus = UInt64(0)
    comptime if C == ModInt:
        modulus = rebind[Modulus](values.ring).value
    var ctx = _data_context[k](modulus, 0, budget)
    return _wrap_basis[k](BasisPolynomial[_Field[k]](0, values.kind, rebind[_Field[k].Context](ctx).copy(), _copy_values[k, k](rebind[List[_Field[k]]](values.coefficients), ctx, budget), _copy_values[k, k](rebind[List[_Field[k]]](values.nodes), ctx, budget)), budget)


def from_core[C: CoefficientField](values: BasisPolynomial[C], *, max_steps: UInt64 = UInt64.MAX) raises PolyError -> Array[_CoreKind[C]]:
    """Validate and copy a core basis representation, preserving coefficient extent.

    Parameters:
        C: The inferred supported exact coefficient field.

    Args:
        values: A mutable core basis representation.
        max_steps: The complete call's work allowance; defaults to UInt64.MAX.

    Returns:
        An independent representation using array axis zero.

    Raises:
        PolyError: DOMAIN for incompatible inputs or unavailable divisions, LIMIT for exhausted work, or CAPACITY for unrepresentable storage.
    """
    var budget = Budget(max_steps)
    return from_core(values, budget)


def _rebuild[k: Int](polynomial: dense.Array[k], model: Array[k], degree: Int, mut budget: Budget) raises PolyError -> Array[k]:
    ref p = rebind[DensePoly[_Field[k]]](polynomial._core[])
    if model._core[].kind == 0:
        return _wrap_basis[k](BasisPolynomial[_Field[k]].newton(p, model._core[].nodes, budget), budget)
    if model._core[].kind == 1:
        return _wrap_basis[k](BasisPolynomial[_Field[k]].chebyshev(p, budget), budget)
    return _wrap_basis[k](BasisPolynomial[_Field[k]].bernstein(p, degree, budget), budget)


def add[ka: Int, kb: Int](a: Array[ka], b: Array[kb], mut budget: Budget) raises PolyError -> Array[_Join[ka, kb]]:
    """Add representations of the same basis kind. The result retains the left Newton nodes.

    Parameters:
        ka: The left coefficient field.
        kb: The right coefficient field.

    Args:
        a: The left basis representation.
        b: The right representation with a compatible field and matching basis kind.
        budget: The shared work allowance, including preparation and result storage.

    Returns:
        A new basis array; insufficient retained Newton nodes raise DOMAIN.

    Raises:
        PolyError: DOMAIN for incompatible inputs or unavailable divisions, LIMIT for exhausted work, or CAPACITY for unrepresentable storage.
    """
    if a.basis() != b.basis():
        raise PolyError(PolyErrorKind.DOMAIN, "basis arithmetic requires matching basis kinds; convert explicitly")
    comptime k = _Join[ka, kb]
    var p = to_polynomial(a, budget)
    var q = to_polynomial(b, budget)
    var result = rebind_var[dense.Array[k]](dense.add(p, q, budget))
    var model = _promote[k](a, result._core[].ring, budget)
    return _rebuild[k](result, model, max(0, max(len(a), len(b)) - 1), budget)


def add[ka: Int, kb: Int](a: Array[ka], b: Array[kb], *, max_steps: UInt64 = UInt64.MAX) raises PolyError -> Array[_Join[ka, kb]]:
    """Add representations of the same basis kind. The result retains the left Newton nodes.

    Parameters:
        ka: The left coefficient field.
        kb: The right coefficient field.

    Args:
        a: The left basis representation.
        b: The right representation with a compatible field and matching basis kind.
        max_steps: The complete call's work allowance; defaults to UInt64.MAX.

    Returns:
        A new basis array; insufficient retained Newton nodes raise DOMAIN.

    Raises:
        PolyError: DOMAIN for incompatible inputs or unavailable divisions, LIMIT for exhausted work, or CAPACITY for unrepresentable storage.
    """
    var budget = Budget(max_steps)
    return add(a, b, budget)


def subtract[ka: Int, kb: Int](a: Array[ka], b: Array[kb], mut budget: Budget) raises PolyError -> Array[_Join[ka, kb]]:
    """Subtract representations of the same basis kind. The result retains the left Newton nodes.

    Parameters:
        ka: The left coefficient field.
        kb: The right coefficient field.

    Args:
        a: The left basis representation.
        b: The right representation with a compatible field and matching basis kind.
        budget: The shared work allowance, including preparation and result storage.

    Returns:
        A new basis array; insufficient retained Newton nodes raise DOMAIN.

    Raises:
        PolyError: DOMAIN for incompatible inputs or unavailable divisions, LIMIT for exhausted work, or CAPACITY for unrepresentable storage.
    """
    if a.basis() != b.basis():
        raise PolyError(PolyErrorKind.DOMAIN, "basis arithmetic requires matching basis kinds; convert explicitly")
    comptime k = _Join[ka, kb]
    var p = to_polynomial(a, budget)
    var q = to_polynomial(b, budget)
    var result = rebind_var[dense.Array[k]](dense.subtract(p, q, budget))
    var model = _promote[k](a, result._core[].ring, budget)
    return _rebuild[k](result, model, max(0, max(len(a), len(b)) - 1), budget)


def subtract[ka: Int, kb: Int](a: Array[ka], b: Array[kb], *, max_steps: UInt64 = UInt64.MAX) raises PolyError -> Array[_Join[ka, kb]]:
    """Subtract representations of the same basis kind. The result retains the left Newton nodes.

    Parameters:
        ka: The left coefficient field.
        kb: The right coefficient field.

    Args:
        a: The left basis representation.
        b: The right representation with a compatible field and matching basis kind.
        max_steps: The complete call's work allowance; defaults to UInt64.MAX.

    Returns:
        A new basis array; insufficient retained Newton nodes raise DOMAIN.

    Raises:
        PolyError: DOMAIN for incompatible inputs or unavailable divisions, LIMIT for exhausted work, or CAPACITY for unrepresentable storage.
    """
    var budget = Budget(max_steps)
    return subtract(a, b, budget)


def multiply[ka: Int, kb: Int](a: Array[ka], b: Array[kb], mut budget: Budget) raises PolyError -> Array[_Join[ka, kb]]:
    """Multiply representations of the same basis kind. The result retains the left Newton nodes.

    Parameters:
        ka: The left coefficient field.
        kb: The right coefficient field.

    Args:
        a: The left basis representation.
        b: The right representation with a compatible field and matching basis kind.
        budget: The shared work allowance, including preparation and result storage.

    Returns:
        A new basis array; insufficient retained Newton nodes raise DOMAIN.

    Raises:
        PolyError: DOMAIN for incompatible inputs or unavailable divisions, LIMIT for exhausted work, or CAPACITY for unrepresentable storage.
    """
    if a.basis() != b.basis():
        raise PolyError(PolyErrorKind.DOMAIN, "basis arithmetic requires matching basis kinds; convert explicitly")
    comptime k = _Join[ka, kb]
    var p = to_polynomial(a, budget)
    var q = to_polynomial(b, budget)
    var result = rebind_var[dense.Array[k]](dense.multiply(p, q, budget))
    var model = _promote[k](a, result._core[].ring, budget)
    return _rebuild[k](result, model, max(0, _sum(max(0, len(a) - 1), max(0, len(b) - 1))), budget)


def multiply[ka: Int, kb: Int](a: Array[ka], b: Array[kb], *, max_steps: UInt64 = UInt64.MAX) raises PolyError -> Array[_Join[ka, kb]]:
    """Multiply representations of the same basis kind. The result retains the left Newton nodes.

    Parameters:
        ka: The left coefficient field.
        kb: The right coefficient field.

    Args:
        a: The left basis representation.
        b: The right representation with a compatible field and matching basis kind.
        max_steps: The complete call's work allowance; defaults to UInt64.MAX.

    Returns:
        A new basis array; insufficient retained Newton nodes raise DOMAIN.

    Raises:
        PolyError: DOMAIN for incompatible inputs or unavailable divisions, LIMIT for exhausted work, or CAPACITY for unrepresentable storage.
    """
    var budget = Budget(max_steps)
    return multiply(a, b, budget)


def negative[k: Int](values: Array[k], mut budget: Budget) raises PolyError -> Array[k]:
    """Negate the represented polynomial and return coefficients in the same basis.

    Parameters:
        k: The inferred coefficient field.

    Args:
        values: The basis representation; integration uses a zero constant.
        budget: The shared work allowance, including preparation and result storage.

    Returns:
        A new representation; Newton nodes are retained and must cover the resulting degree.

    Raises:
        PolyError: DOMAIN for incompatible inputs or unavailable divisions, LIMIT for exhausted work, or CAPACITY for unrepresentable storage.
    """
    var polynomial = to_polynomial(values, budget)
    var result = dense.negative(polynomial, budget)
    return _rebuild[k](rebind[dense.Array[k]](result), values, max(0, len(values) - 1), budget)


def negative[k: Int](values: Array[k], *, max_steps: UInt64 = UInt64.MAX) raises PolyError -> Array[k]:
    """Negate the represented polynomial and return coefficients in the same basis.

    Parameters:
        k: The inferred coefficient field.

    Args:
        values: The basis representation; integration uses a zero constant.
        max_steps: The complete call's work allowance; defaults to UInt64.MAX.

    Returns:
        A new representation; Newton nodes are retained and must cover the resulting degree.

    Raises:
        PolyError: DOMAIN for incompatible inputs or unavailable divisions, LIMIT for exhausted work, or CAPACITY for unrepresentable storage.
    """
    var budget = Budget(max_steps)
    return negative(values, budget)


def derivative[k: Int](values: Array[k], mut budget: Budget) raises PolyError -> Array[k]:
    """Differentiate the represented polynomial and return coefficients in the same basis.

    Parameters:
        k: The inferred coefficient field.

    Args:
        values: The basis representation; integration uses a zero constant.
        budget: The shared work allowance, including preparation and result storage.

    Returns:
        A new representation; Newton nodes are retained and must cover the resulting degree.

    Raises:
        PolyError: DOMAIN for incompatible inputs or unavailable divisions, LIMIT for exhausted work, or CAPACITY for unrepresentable storage.
    """
    var polynomial = to_polynomial(values, budget)
    var result = dense.derivative(polynomial, budget)
    return _rebuild[k](rebind[dense.Array[k]](result), values, max(0, len(values) - 2), budget)


def derivative[k: Int](values: Array[k], *, max_steps: UInt64 = UInt64.MAX) raises PolyError -> Array[k]:
    """Differentiate the represented polynomial and return coefficients in the same basis.

    Parameters:
        k: The inferred coefficient field.

    Args:
        values: The basis representation; integration uses a zero constant.
        max_steps: The complete call's work allowance; defaults to UInt64.MAX.

    Returns:
        A new representation; Newton nodes are retained and must cover the resulting degree.

    Raises:
        PolyError: DOMAIN for incompatible inputs or unavailable divisions, LIMIT for exhausted work, or CAPACITY for unrepresentable storage.
    """
    var budget = Budget(max_steps)
    return derivative(values, budget)


def integral[k: Int](values: Array[k], mut budget: Budget) raises PolyError -> Array[k]:
    """Integrate the represented polynomial and return coefficients in the same basis.

    Parameters:
        k: The inferred coefficient field.

    Args:
        values: The basis representation; integration uses a zero constant.
        budget: The shared work allowance, including preparation and result storage.

    Returns:
        A new representation; Newton nodes are retained and must cover the resulting degree.

    Raises:
        PolyError: DOMAIN for incompatible inputs or unavailable divisions, LIMIT for exhausted work, or CAPACITY for unrepresentable storage.
    """
    var polynomial = to_polynomial(values, budget)
    var result = dense.integral(polynomial, budget)
    return _rebuild[k](rebind[dense.Array[k]](result), values, len(values), budget)


def integral[k: Int](values: Array[k], *, max_steps: UInt64 = UInt64.MAX) raises PolyError -> Array[k]:
    """Integrate the represented polynomial and return coefficients in the same basis.

    Parameters:
        k: The inferred coefficient field.

    Args:
        values: The basis representation; integration uses a zero constant.
        max_steps: The complete call's work allowance; defaults to UInt64.MAX.

    Returns:
        A new representation; Newton nodes are retained and must cover the resulting degree.

    Raises:
        PolyError: DOMAIN for incompatible inputs or unavailable divisions, LIMIT for exhausted work, or CAPACITY for unrepresentable storage.
    """
    var budget = Budget(max_steps)
    return integral(values, budget)


from . import _numeric as _num
from polynomials import NumericDomain, NumericBasisPolynomial


@fieldwise_init
struct NumericArray[D: NumericDomain](Copyable, Movable, Iterable, Sized, _Record):
    """Immutable numerical coefficients in a Newton, Chebyshev, or Bernstein basis.

    Parameters:
        D: The numerical coefficient domain.
    """
    comptime _Format = _Basis
    comptime _Item = Self.D.Scalar
    comptime _Power = Int
    var _core: ArcPointer[NumericBasisPolynomial[Self.D]]

    comptime IteratorType[iterable_mut: Bool, //, iterable_origin: Origin[mut=iterable_mut]]: Iterator = _NumericIterator[Self.D, ImmOrigin(iterable_origin)]
    """Borrow the representation's complete coefficient sequence.

    Parameters:
        iterable_mut: Whether the borrowing origin is mutable.
        iterable_origin: The source origin.
    """

    def __iter__(ref self) -> Self.IteratorType[origin_of(self)]:
        """Iterate over every basis coefficient, including trailing zeros.

        Returns:
            An iterator borrowing the immutable representation.
        """
        return _NumericIterator[Self.D, ImmOrigin(origin_of(self))](Pointer(to=self).as_imm(), 0)

    def __len__(self) -> Int:
        """Return the full coefficient extent, including zeros.

        Returns:
            The retained coefficient count.
        """
        return len(self._core[].coefficients)

    def __getitem__(self, index: Int) -> Self.D.Scalar:
        """Copy a coefficient at a valid zero-based index.

        Args:
            index: An index below len(self).

        Returns:
            An independent coefficient value.
        """
        debug_assert(index >= 0 and index < len(self), "coefficient index out of range")
        return self._core[].coefficients[index].copy()

    def basis(self) -> String:
        """Return the retained basis name.

        Returns:
            Newton, Chebyshev, or Bernstein in lowercase.
        """
        return "newton" if self._core[].kind == 0 else "chebyshev" if self._core[].kind == 1 else "bernstein"


@fieldwise_init
struct _NumericIterator[D: NumericDomain, origin: ImmOrigin](Iterator):
    comptime Element = Self.D.Scalar
    var source: Pointer[NumericArray[Self.D], Self.origin]
    var position: Int

    def __next__(mut self) raises StopIteration -> Self.D.Scalar:
        if self.position == len(self.source[]):
            raise StopIteration()
        var result = self.source[][self.position]
        self.position += 1
        return result^


def _numeric_wrap[D: NumericDomain](var core: NumericBasisPolynomial[D], context: D.Context, mut budget: Budget) raises PolyError -> NumericArray[D]:
    _num._context[D](context, budget)
    _compact(core.coefficients, budget)
    _compact(core.nodes, budget)
    _allocation[NumericBasisPolynomial[D]](1, budget)
    return NumericArray[D](ArcPointer(core^))


def asarray[I: Iterable, Context: Copyable & Deinitable, N: Iterable](values: I, mut budget: Budget, *, context: Context, basis: String, nodes: N) raises PolyError -> NumericArray[_num._Domain[Context, _num._CombinedScalar[_Element[I], _Element[N]]]]:
    """Construct a numerical basis representation with explicit precision.

    Parameters:
        I: Inferred from the input types.
        Context: Inferred from the input types.
        N: Inferred from the input types.

    Args:
        values: Coefficients, terms, or a prepared value.
        budget: The shared work allowance.
        context: Explicit output precision and rounding or enclosure settings.
        basis: One of newton, chebyshev, or bernstein.
        nodes: Newton nodes; ignored only when empty for other bases.

    Returns:
        Numerical basis coefficients with retained nodes or representation degree.

    Raises:
        PolyError: DOMAIN for invalid input, or LIMIT/CAPACITY on resource limits.
    """
    comptime D = _num._Domain[Context, _num._CombinedScalar[_Element[I], _Element[N]]]
    var sites = _num._sequence[D](nodes, rebind[D.Context](context), budget)
    var kind = _kind(basis)
    var c = _num._sequence[D](values, rebind[D.Context](context), budget)
    _shape(kind, len(c), len(sites))
    return _numeric_wrap(NumericBasisPolynomial[D](0, kind, c^, sites^), rebind[D.Context](context), budget)


def asarray[I: Iterable, Context: Copyable & Deinitable, N: Iterable](values: I, *, context: Context, basis: String, nodes: N, max_steps: UInt64 = UInt64.MAX) raises PolyError -> NumericArray[_num._Domain[Context, _num._CombinedScalar[_Element[I], _Element[N]]]]:
    """Construct a numerical basis representation with explicit precision.

    Parameters:
        I: Inferred from the input types.
        Context: Inferred from the input types.
        N: Inferred from the input types.

    Args:
        values: Coefficients, terms, or a prepared value.
        context: Explicit output precision and rounding or enclosure settings.
        basis: One of newton, chebyshev, or bernstein.
        nodes: Newton nodes; ignored only when empty for other bases.
        max_steps: The complete call's work allowance.

    Returns:
        Numerical basis coefficients with retained nodes or representation degree.

    Raises:
        PolyError: DOMAIN for invalid input, or LIMIT/CAPACITY on resource limits.
    """
    var budget = Budget(max_steps)
    return asarray(values, budget, context=context, basis=basis, nodes=nodes)


def asarray[I: Iterable, Context: Copyable & Deinitable](values: I, mut budget: Budget, *, context: Context, basis: String) raises PolyError -> NumericArray[_num._Domain[Context, _Element[I]]]:
    """Construct a numerical basis representation with explicit precision.

    Parameters:
        I: Inferred from the input types.
        Context: Inferred from the input types.

    Args:
        values: Coefficients, terms, or a prepared value.
        budget: The shared work allowance.
        context: Explicit output precision and rounding or enclosure settings.
        basis: One of newton, chebyshev, or bernstein.

    Returns:
        Numerical basis coefficients with retained nodes or representation degree.

    Raises:
        PolyError: DOMAIN for invalid input, or LIMIT/CAPACITY on resource limits.
    """
    return rebind_var[NumericArray[_num._Domain[Context, _Element[I]]]](asarray(values, budget, context=context, basis=basis, nodes=List[Int]()))


def asarray[I: Iterable, Context: Copyable & Deinitable](values: I, *, context: Context, basis: String, max_steps: UInt64 = UInt64.MAX) raises PolyError -> NumericArray[_num._Domain[Context, _Element[I]]]:
    """Construct a numerical basis representation with explicit precision.

    Parameters:
        I: Inferred from the input types.
        Context: Inferred from the input types.

    Args:
        values: Coefficients, terms, or a prepared value.
        context: Explicit output precision and rounding or enclosure settings.
        basis: One of newton, chebyshev, or bernstein.
        max_steps: The complete call's work allowance.

    Returns:
        Numerical basis coefficients with retained nodes or representation degree.

    Raises:
        PolyError: DOMAIN for invalid input, or LIMIT/CAPACITY on resource limits.
    """
    var budget = Budget(max_steps)
    return asarray(values, budget, context=context, basis=basis)


def from_polynomial[I: Iterable, Context: Copyable & Deinitable, N: Iterable](values: I, mut budget: Budget, *, context: Context, basis: String, nodes: N, degree: Int = -1) raises PolyError -> NumericArray[_num._Domain[Context, _num._CombinedScalar[_Element[I], _Element[N]]]]:
    """Convert numerical monomial coefficients to a named basis.

    Parameters:
        I: Inferred from the input types.
        Context: Inferred from the input types.
        N: Inferred from the input types.

    Args:
        values: Coefficients, terms, or a prepared value.
        budget: The shared work allowance.
        context: Explicit output precision and rounding or enclosure settings.
        basis: One of newton, chebyshev, or bernstein.
        nodes: Newton nodes; ignored only when empty for other bases.
        degree: Retained Bernstein degree, or -1 to infer it.

    Returns:
        Numerical basis coefficients with retained nodes or representation degree.

    Raises:
        PolyError: DOMAIN for invalid input, or LIMIT/CAPACITY on resource limits.
    """
    comptime D = _num._Domain[Context, _num._CombinedScalar[_Element[I], _Element[N]]]
    var sites = _num._sequence[D](nodes, rebind[D.Context](context), budget)
    var kind = _kind(basis)
    var p = _num._dense[D](values, rebind[D.Context](context), budget)
    if kind != 0 and len(sites):
        raise PolyError(PolyErrorKind.DOMAIN, "only Newton representations accept nodes")
    if kind == 0:
        if degree != -1:
            raise PolyError(PolyErrorKind.DOMAIN, "degree is only used for Bernstein representations")
        return _numeric_wrap(NumericBasisPolynomial[D].newton(p._core[], sites, rebind[D.Context](context), budget), rebind[D.Context](context), budget)
    if kind == 1:
        if degree != -1:
            raise PolyError(PolyErrorKind.DOMAIN, "degree is only used for Bernstein representations")
        return _numeric_wrap(NumericBasisPolynomial[D].chebyshev(p._core[], rebind[D.Context](context), budget), rebind[D.Context](context), budget)
    if degree < -1:
        raise PolyError(PolyErrorKind.DOMAIN, "invalid Bernstein degree")
    var n = max(0, len(p) - 1) if degree == -1 else degree
    return _numeric_wrap(NumericBasisPolynomial[D].bernstein(p._core[], n, rebind[D.Context](context), budget), rebind[D.Context](context), budget)


def from_polynomial[I: Iterable, Context: Copyable & Deinitable, N: Iterable](values: I, *, context: Context, basis: String, nodes: N, degree: Int = -1, max_steps: UInt64 = UInt64.MAX) raises PolyError -> NumericArray[_num._Domain[Context, _num._CombinedScalar[_Element[I], _Element[N]]]]:
    """Convert numerical monomial coefficients to a named basis.

    Parameters:
        I: Inferred from the input types.
        Context: Inferred from the input types.
        N: Inferred from the input types.

    Args:
        values: Coefficients, terms, or a prepared value.
        context: Explicit output precision and rounding or enclosure settings.
        basis: One of newton, chebyshev, or bernstein.
        nodes: Newton nodes; ignored only when empty for other bases.
        degree: Retained Bernstein degree, or -1 to infer it.
        max_steps: The complete call's work allowance.

    Returns:
        Numerical basis coefficients with retained nodes or representation degree.

    Raises:
        PolyError: DOMAIN for invalid input, or LIMIT/CAPACITY on resource limits.
    """
    var budget = Budget(max_steps)
    return from_polynomial(values, budget, context=context, basis=basis, nodes=nodes, degree=degree)


def from_polynomial[I: Iterable, Context: Copyable & Deinitable](values: I, mut budget: Budget, *, context: Context, basis: String, degree: Int = -1) raises PolyError -> NumericArray[_num._Domain[Context, _Element[I]]]:
    """Convert numerical monomial coefficients to a named basis.

    Parameters:
        I: Inferred from the input types.
        Context: Inferred from the input types.

    Args:
        values: Coefficients, terms, or a prepared value.
        budget: The shared work allowance.
        context: Explicit output precision and rounding or enclosure settings.
        basis: One of newton, chebyshev, or bernstein.
        degree: Retained Bernstein degree, or -1 to infer it.

    Returns:
        Numerical basis coefficients with retained nodes or representation degree.

    Raises:
        PolyError: DOMAIN for invalid input, or LIMIT/CAPACITY on resource limits.
    """
    return rebind_var[NumericArray[_num._Domain[Context, _Element[I]]]](from_polynomial(values, budget, context=context, basis=basis, nodes=List[Int](), degree=degree))


def from_polynomial[I: Iterable, Context: Copyable & Deinitable](values: I, *, context: Context, basis: String, degree: Int = -1, max_steps: UInt64 = UInt64.MAX) raises PolyError -> NumericArray[_num._Domain[Context, _Element[I]]]:
    """Convert numerical monomial coefficients to a named basis.

    Parameters:
        I: Inferred from the input types.
        Context: Inferred from the input types.

    Args:
        values: Coefficients, terms, or a prepared value.
        context: Explicit output precision and rounding or enclosure settings.
        basis: One of newton, chebyshev, or bernstein.
        degree: Retained Bernstein degree, or -1 to infer it.
        max_steps: The complete call's work allowance.

    Returns:
        Numerical basis coefficients with retained nodes or representation degree.

    Raises:
        PolyError: DOMAIN for invalid input, or LIMIT/CAPACITY on resource limits.
    """
    var budget = Budget(max_steps)
    return from_polynomial(values, budget, context=context, basis=basis, degree=degree)


def to_polynomial[D: NumericDomain](values: NumericArray[D], mut budget: Budget, *, context: D.Context) raises PolyError -> dense.NumericArray[D]:
    """Expand numerical basis coefficients into monomial coefficients.

    Parameters:
        D: The numerical coefficient domain.

    Args:
        values: Coefficients, terms, or a prepared value.
        budget: The shared work allowance.
        context: Explicit output precision and rounding or enclosure settings.

    Returns:
        An independent dense array of numerical monomial coefficients.

    Raises:
        PolyError: DOMAIN for invalid input, or LIMIT/CAPACITY on resource limits.
    """
    var result = values._core[].to_polynomial(context, budget)
    return dense.NumericArray[D](result^, budget)


def to_polynomial[D: NumericDomain](values: NumericArray[D], *, context: D.Context, max_steps: UInt64 = UInt64.MAX) raises PolyError -> dense.NumericArray[D]:
    """Expand numerical basis coefficients into monomial coefficients.

    Parameters:
        D: The numerical coefficient domain.

    Args:
        values: Coefficients, terms, or a prepared value.
        context: Explicit output precision and rounding or enclosure settings.
        max_steps: The complete call's work allowance.

    Returns:
        An independent dense array of numerical monomial coefficients.

    Raises:
        PolyError: DOMAIN for invalid input, or LIMIT/CAPACITY on resource limits.
    """
    var budget = Budget(max_steps)
    return to_polynomial(values, budget, context=context)


def coefficients[D: NumericDomain](values: NumericArray[D], mut budget: Budget) raises PolyError -> _num.Values[D]:
    """Copy the numerical coefficients without trimming zeros.

    Parameters:
        D: The numerical coefficient domain.

    Args:
        values: Coefficients, terms, or a prepared value.
        budget: The shared work allowance.

    Returns:
        A numerical sequence retaining every entry, including zeros.

    Raises:
        PolyError: DOMAIN for invalid input, or LIMIT/CAPACITY on resource limits.
    """
    ref stored = values._core[].coefficients
    _allocation[D.Scalar](len(stored), budget)
    for c in stored:
        _ = _num._bounds(c, budget)
    return _num.Values[D](stored.copy(), budget)


def coefficients[D: NumericDomain](values: NumericArray[D], *, max_steps: UInt64 = UInt64.MAX) raises PolyError -> _num.Values[D]:
    """Copy the numerical coefficients without trimming zeros.

    Parameters:
        D: The numerical coefficient domain.

    Args:
        values: Coefficients, terms, or a prepared value.
        max_steps: The complete call's work allowance.

    Returns:
        A numerical sequence retaining every entry, including zeros.

    Raises:
        PolyError: DOMAIN for invalid input, or LIMIT/CAPACITY on resource limits.
    """
    var budget = Budget(max_steps)
    return coefficients(values, budget)


def nodes[D: NumericDomain](values: NumericArray[D], mut budget: Budget) raises PolyError -> _num.Values[D]:
    """Copy the numerical nodes without trimming zeros.

    Parameters:
        D: The numerical coefficient domain.

    Args:
        values: Coefficients, terms, or a prepared value.
        budget: The shared work allowance.

    Returns:
        A numerical sequence retaining every entry, including zeros.

    Raises:
        PolyError: DOMAIN for invalid input, or LIMIT/CAPACITY on resource limits.
    """
    ref stored = values._core[].nodes
    _allocation[D.Scalar](len(stored), budget)
    for c in stored:
        _ = _num._bounds(c, budget)
    return _num.Values[D](stored.copy(), budget)


def nodes[D: NumericDomain](values: NumericArray[D], *, max_steps: UInt64 = UInt64.MAX) raises PolyError -> _num.Values[D]:
    """Copy the numerical nodes without trimming zeros.

    Parameters:
        D: The numerical coefficient domain.

    Args:
        values: Coefficients, terms, or a prepared value.
        max_steps: The complete call's work allowance.

    Returns:
        A numerical sequence retaining every entry, including zeros.

    Raises:
        PolyError: DOMAIN for invalid input, or LIMIT/CAPACITY on resource limits.
    """
    var budget = Budget(max_steps)
    return nodes(values, budget)


def evaluate[D: NumericDomain, X: Copyable & Deinitable](values: NumericArray[D], point: X, mut budget: Budget, *, context: D.Context) raises PolyError -> D.Scalar:
    """Evaluate using the recurrence for the retained numerical basis.

    Parameters:
        D: The numerical coefficient domain.
        X: Inferred from the input types.

    Args:
        values: Coefficients, terms, or a prepared value.
        point: The evaluation or expansion point.
        budget: The shared work allowance.
        context: Explicit output precision and rounding or enclosure settings.

    Returns:
        The rounded value or enclosing scalar.

    Raises:
        PolyError: DOMAIN for invalid input, or LIMIT/CAPACITY on resource limits.
    """
    return values._core[].evaluate(_num._value[D](point, context, budget), context, budget)


def evaluate[D: NumericDomain, X: Copyable & Deinitable](values: NumericArray[D], point: X, *, context: D.Context, max_steps: UInt64 = UInt64.MAX) raises PolyError -> D.Scalar:
    """Evaluate using the recurrence for the retained numerical basis.

    Parameters:
        D: The numerical coefficient domain.
        X: Inferred from the input types.

    Args:
        values: Coefficients, terms, or a prepared value.
        point: The evaluation or expansion point.
        context: Explicit output precision and rounding or enclosure settings.
        max_steps: The complete call's work allowance.

    Returns:
        The rounded value or enclosing scalar.

    Raises:
        PolyError: DOMAIN for invalid input, or LIMIT/CAPACITY on resource limits.
    """
    var budget = Budget(max_steps)
    return evaluate(values, point, budget, context=context)


def subdivide[D: NumericDomain, X: Copyable & Deinitable](values: NumericArray[D], point: X, mut budget: Budget, *, context: D.Context) raises PolyError -> Tuple[NumericArray[D], NumericArray[D]]:
    """Subdivide a numerical Bernstein representation at a parameter.

    Parameters:
        D: The numerical coefficient domain.
        X: Inferred from the input types.

    Args:
        values: Coefficients, terms, or a prepared value.
        point: The evaluation or expansion point.
        budget: The shared work allowance.
        context: Explicit output precision and rounding or enclosure settings.

    Returns:
        The two Bernstein pieces, each retaining the original degree.

    Raises:
        PolyError: DOMAIN for invalid input, or LIMIT/CAPACITY on resource limits.
    """
    var result = values._core[].subdivide(_num._value[D](point, context, budget), context, budget)
    var left = NumericBasisPolynomial[D](0, 2, List[D.Scalar](), List[D.Scalar]())
    var right = NumericBasisPolynomial[D](0, 2, List[D.Scalar](), List[D.Scalar]())
    swap(left, result[0])
    swap(right, result[1])
    return _numeric_wrap(left^, context, budget), _numeric_wrap(right^, context, budget)


def subdivide[D: NumericDomain, X: Copyable & Deinitable](values: NumericArray[D], point: X, *, context: D.Context, max_steps: UInt64 = UInt64.MAX) raises PolyError -> Tuple[NumericArray[D], NumericArray[D]]:
    """Subdivide a numerical Bernstein representation at a parameter.

    Parameters:
        D: The numerical coefficient domain.
        X: Inferred from the input types.

    Args:
        values: Coefficients, terms, or a prepared value.
        point: The evaluation or expansion point.
        context: Explicit output precision and rounding or enclosure settings.
        max_steps: The complete call's work allowance.

    Returns:
        The two Bernstein pieces, each retaining the original degree.

    Raises:
        PolyError: DOMAIN for invalid input, or LIMIT/CAPACITY on resource limits.
    """
    var budget = Budget(max_steps)
    return subdivide(values, point, budget, context=context)


def _numeric_rebuild[D: NumericDomain](p: dense.NumericArray[D], model: NumericArray[D], degree: Int, context: D.Context, mut budget: Budget) raises PolyError -> NumericArray[D]:
    if model._core[].kind == 0:
        return _numeric_wrap(NumericBasisPolynomial[D].newton(p._core[], model._core[].nodes, context, budget), context, budget)
    if model._core[].kind == 1:
        return _numeric_wrap(NumericBasisPolynomial[D].chebyshev(p._core[], context, budget), context, budget)
    return _numeric_wrap(NumericBasisPolynomial[D].bernstein(p._core[], max(0, degree), context, budget), context, budget)


def add[D: NumericDomain](a: NumericArray[D], b: NumericArray[D], mut budget: Budget, *, context: D.Context) raises PolyError -> NumericArray[D]:
    """Add numerical polynomials, retaining the left operand's basis.

    Newton results retain the left nodes and raise if those nodes cannot represent the result.

    Parameters:
        D: The numerical coefficient domain.

    Args:
        a: The first prepared value.
        b: The second prepared value.
        budget: The shared work allowance.
        context: Explicit output precision and rounding or enclosure settings.

    Returns:
        Numerical basis coefficients with retained nodes or representation degree.

    Raises:
        PolyError: DOMAIN for invalid input, or LIMIT/CAPACITY on resource limits.
    """
    var left = to_polynomial(a, budget, context=context)
    var right = to_polynomial(b, budget, context=context)
    var result = dense.add(left, right, budget, context=context)
    return _numeric_rebuild(result, a, max(len(a), len(b))-1, context, budget)


def add[D: NumericDomain](a: NumericArray[D], b: NumericArray[D], *, context: D.Context, max_steps: UInt64 = UInt64.MAX) raises PolyError -> NumericArray[D]:
    """Add numerical polynomials, retaining the left operand's basis.

    Newton results retain the left nodes and raise if those nodes cannot represent the result.

    Parameters:
        D: The numerical coefficient domain.

    Args:
        a: The first prepared value.
        b: The second prepared value.
        context: Explicit output precision and rounding or enclosure settings.
        max_steps: The complete call's work allowance.

    Returns:
        Numerical basis coefficients with retained nodes or representation degree.

    Raises:
        PolyError: DOMAIN for invalid input, or LIMIT/CAPACITY on resource limits.
    """
    var budget = Budget(max_steps)
    return add(a, b, budget, context=context)


def subtract[D: NumericDomain](a: NumericArray[D], b: NumericArray[D], mut budget: Budget, *, context: D.Context) raises PolyError -> NumericArray[D]:
    """Subtract numerical polynomials, retaining the left operand's basis.

    Newton results retain the left nodes and raise if those nodes cannot represent the result.

    Parameters:
        D: The numerical coefficient domain.

    Args:
        a: The first prepared value.
        b: The second prepared value.
        budget: The shared work allowance.
        context: Explicit output precision and rounding or enclosure settings.

    Returns:
        Numerical basis coefficients with retained nodes or representation degree.

    Raises:
        PolyError: DOMAIN for invalid input, or LIMIT/CAPACITY on resource limits.
    """
    var left = to_polynomial(a, budget, context=context)
    var right = to_polynomial(b, budget, context=context)
    var result = dense.subtract(left, right, budget, context=context)
    return _numeric_rebuild(result, a, max(len(a), len(b))-1, context, budget)


def subtract[D: NumericDomain](a: NumericArray[D], b: NumericArray[D], *, context: D.Context, max_steps: UInt64 = UInt64.MAX) raises PolyError -> NumericArray[D]:
    """Subtract numerical polynomials, retaining the left operand's basis.

    Newton results retain the left nodes and raise if those nodes cannot represent the result.

    Parameters:
        D: The numerical coefficient domain.

    Args:
        a: The first prepared value.
        b: The second prepared value.
        context: Explicit output precision and rounding or enclosure settings.
        max_steps: The complete call's work allowance.

    Returns:
        Numerical basis coefficients with retained nodes or representation degree.

    Raises:
        PolyError: DOMAIN for invalid input, or LIMIT/CAPACITY on resource limits.
    """
    var budget = Budget(max_steps)
    return subtract(a, b, budget, context=context)


def multiply[D: NumericDomain](a: NumericArray[D], b: NumericArray[D], mut budget: Budget, *, context: D.Context) raises PolyError -> NumericArray[D]:
    """Multiply numerical polynomials, retaining the left operand's basis.

    Newton results retain the left nodes and raise if those nodes cannot represent the result.

    Parameters:
        D: The numerical coefficient domain.

    Args:
        a: The first prepared value.
        b: The second prepared value.
        budget: The shared work allowance.
        context: Explicit output precision and rounding or enclosure settings.

    Returns:
        Numerical basis coefficients with retained nodes or representation degree.

    Raises:
        PolyError: DOMAIN for invalid input, or LIMIT/CAPACITY on resource limits.
    """
    var left = to_polynomial(a, budget, context=context)
    var right = to_polynomial(b, budget, context=context)
    var result = dense.multiply(left, right, budget, context=context)
    return _numeric_rebuild(result, a, _sum(max(0, len(a)-1), max(0, len(b)-1)), context, budget)


def multiply[D: NumericDomain](a: NumericArray[D], b: NumericArray[D], *, context: D.Context, max_steps: UInt64 = UInt64.MAX) raises PolyError -> NumericArray[D]:
    """Multiply numerical polynomials, retaining the left operand's basis.

    Newton results retain the left nodes and raise if those nodes cannot represent the result.

    Parameters:
        D: The numerical coefficient domain.

    Args:
        a: The first prepared value.
        b: The second prepared value.
        context: Explicit output precision and rounding or enclosure settings.
        max_steps: The complete call's work allowance.

    Returns:
        Numerical basis coefficients with retained nodes or representation degree.

    Raises:
        PolyError: DOMAIN for invalid input, or LIMIT/CAPACITY on resource limits.
    """
    var budget = Budget(max_steps)
    return multiply(a, b, budget, context=context)


def negative[D: NumericDomain](values: NumericArray[D], mut budget: Budget, *, context: D.Context) raises PolyError -> NumericArray[D]:
    """Negate the numerical polynomial, retaining its basis.

    Parameters:
        D: The numerical coefficient domain.

    Args:
        values: Coefficients, terms, or a prepared value.
        budget: The shared work allowance.
        context: Explicit output precision and rounding or enclosure settings.

    Returns:
        Numerical basis coefficients with retained nodes or representation degree.

    Raises:
        PolyError: DOMAIN for invalid input, or LIMIT/CAPACITY on resource limits.
    """
    var p = to_polynomial(values, budget, context=context)
    var result = dense.negative(p, budget, context=context)
    return _numeric_rebuild(result, values, len(values)-1, context, budget)


def negative[D: NumericDomain](values: NumericArray[D], *, context: D.Context, max_steps: UInt64 = UInt64.MAX) raises PolyError -> NumericArray[D]:
    """Negate the numerical polynomial, retaining its basis.

    Parameters:
        D: The numerical coefficient domain.

    Args:
        values: Coefficients, terms, or a prepared value.
        context: Explicit output precision and rounding or enclosure settings.
        max_steps: The complete call's work allowance.

    Returns:
        Numerical basis coefficients with retained nodes or representation degree.

    Raises:
        PolyError: DOMAIN for invalid input, or LIMIT/CAPACITY on resource limits.
    """
    var budget = Budget(max_steps)
    return negative(values, budget, context=context)


def derivative[D: NumericDomain](values: NumericArray[D], mut budget: Budget, *, context: D.Context) raises PolyError -> NumericArray[D]:
    """Return the numerical derivative in the retained basis.

    Parameters:
        D: The numerical coefficient domain.

    Args:
        values: Coefficients, terms, or a prepared value.
        budget: The shared work allowance.
        context: Explicit output precision and rounding or enclosure settings.

    Returns:
        Numerical basis coefficients with retained nodes or representation degree.

    Raises:
        PolyError: DOMAIN for invalid input, or LIMIT/CAPACITY on resource limits.
    """
    var p = to_polynomial(values, budget, context=context)
    var result = dense.derivative(p, budget, context=context)
    return _numeric_rebuild(result, values, max(0, len(values)-2), context, budget)


def derivative[D: NumericDomain](values: NumericArray[D], *, context: D.Context, max_steps: UInt64 = UInt64.MAX) raises PolyError -> NumericArray[D]:
    """Return the numerical derivative in the retained basis.

    Parameters:
        D: The numerical coefficient domain.

    Args:
        values: Coefficients, terms, or a prepared value.
        context: Explicit output precision and rounding or enclosure settings.
        max_steps: The complete call's work allowance.

    Returns:
        Numerical basis coefficients with retained nodes or representation degree.

    Raises:
        PolyError: DOMAIN for invalid input, or LIMIT/CAPACITY on resource limits.
    """
    var budget = Budget(max_steps)
    return derivative(values, budget, context=context)


def integral[D: NumericDomain](values: NumericArray[D], mut budget: Budget, *, context: D.Context) raises PolyError -> NumericArray[D]:
    """Return the numerical integral in the retained basis.

    Parameters:
        D: The numerical coefficient domain.

    Args:
        values: Coefficients, terms, or a prepared value.
        budget: The shared work allowance.
        context: Explicit output precision and rounding or enclosure settings.

    Returns:
        Numerical basis coefficients with retained nodes or representation degree.

    Raises:
        PolyError: DOMAIN for invalid input, or LIMIT/CAPACITY on resource limits.
    """
    var p = to_polynomial(values, budget, context=context)
    var result = dense.integral(p, budget, context=context)
    return _numeric_rebuild(result, values, _sum(len(values), 0), context, budget)


def integral[D: NumericDomain](values: NumericArray[D], *, context: D.Context, max_steps: UInt64 = UInt64.MAX) raises PolyError -> NumericArray[D]:
    """Return the numerical integral in the retained basis.

    Parameters:
        D: The numerical coefficient domain.

    Args:
        values: Coefficients, terms, or a prepared value.
        context: Explicit output precision and rounding or enclosure settings.
        max_steps: The complete call's work allowance.

    Returns:
        Numerical basis coefficients with retained nodes or representation degree.

    Raises:
        PolyError: DOMAIN for invalid input, or LIMIT/CAPACITY on resource limits.
    """
    var budget = Budget(max_steps)
    return integral(values, budget, context=context)


def to_core[D: NumericDomain](values: NumericArray[D], mut budget: Budget, *, variable: Int64 = 0) raises PolyError -> NumericBasisPolynomial[D]:
    """Copy a numerical basis into independent mutable core storage.

    Parameters:
        D: The numerical coefficient domain.

    Args:
        values: Coefficients, terms, or a prepared value.
        budget: The shared work allowance.
        variable: The core variable ID assigned to the exported value.

    Returns:
        An independent mutable core representation, including its metadata.

    Raises:
        PolyError: DOMAIN for invalid input, or LIMIT/CAPACITY on resource limits.
    """
    var c = coefficients(values, budget).to_list(budget)
    var n = nodes(values, budget).to_list(budget)
    return NumericBasisPolynomial[D](variable, values._core[].kind, c^, n^)


def to_core[D: NumericDomain](values: NumericArray[D], *, max_steps: UInt64 = UInt64.MAX, variable: Int64 = 0) raises PolyError -> NumericBasisPolynomial[D]:
    """Copy a numerical basis into independent mutable core storage.

    Parameters:
        D: The numerical coefficient domain.

    Args:
        values: Coefficients, terms, or a prepared value.
        max_steps: The complete call's work allowance.
        variable: The core variable ID assigned to the exported value.

    Returns:
        An independent mutable core representation, including its metadata.

    Raises:
        PolyError: DOMAIN for invalid input, or LIMIT/CAPACITY on resource limits.
    """
    var budget = Budget(max_steps)
    return to_core(values, budget, variable=variable)


def from_core[D: NumericDomain](values: NumericBasisPolynomial[D], mut budget: Budget, *, context: D.Context) raises PolyError -> NumericArray[D]:
    """Validate and copy a numerical core basis representation.

    Parameters:
        D: The numerical coefficient domain.

    Args:
        values: Coefficients, terms, or a prepared value.
        budget: The shared work allowance.
        context: Explicit output precision and rounding or enclosure settings.

    Returns:
        Numerical basis coefficients with retained nodes or representation degree.

    Raises:
        PolyError: DOMAIN for invalid input, or LIMIT/CAPACITY on resource limits.
    """
    if values.kind > 2:
        raise PolyError(PolyErrorKind.DOMAIN, "invalid numerical basis metadata")
    _shape(values.kind, len(values.coefficients), len(values.nodes))
    _num._context[D](context, budget)
    _allocation[D.Scalar](len(values.coefficients), budget)
    _allocation[D.Scalar](len(values.nodes), budget)
    for c in values.coefficients:
        _ = _num._bounds(c, budget)
    for c in values.nodes:
        _ = _num._bounds(c, budget)
    return _numeric_wrap(NumericBasisPolynomial[D](0, values.kind, values.coefficients.copy(), values.nodes.copy()), context, budget)


def from_core[D: NumericDomain](values: NumericBasisPolynomial[D], *, context: D.Context, max_steps: UInt64 = UInt64.MAX) raises PolyError -> NumericArray[D]:
    """Validate and copy a numerical core basis representation.

    Parameters:
        D: The numerical coefficient domain.

    Args:
        values: Coefficients, terms, or a prepared value.
        context: Explicit output precision and rounding or enclosure settings.
        max_steps: The complete call's work allowance.

    Returns:
        Numerical basis coefficients with retained nodes or representation degree.

    Raises:
        PolyError: DOMAIN for invalid input, or LIMIT/CAPACITY on resource limits.
    """
    var budget = Budget(max_steps)
    return from_core(values, budget, context=context)
