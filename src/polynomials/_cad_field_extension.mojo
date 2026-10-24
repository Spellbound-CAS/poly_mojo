"""Extend a sample field only when a relative section needs a further lift."""
from .coefficients import QQ, NoContext
from .dense import DensePoly
from .errors import Budget, PolyError, PolyErrorKind, _allocation, _sum
from .algebraic_types import RealAlgebraic, _rational
from ._algebraic import _Frame, _root
from ._algebraic_filter import _constant, _compose, _field_mul, _combination_norm, _combined_frame, _reduce, _gcd_over_field
from ._cad_sample_field import _SampleField, _sample_field, _check_sample_field, _ksign
from ._cad_final import _FinalRoot
from ._cad_field_stack import _FieldStack, _krefine
from ._root_ops import _sign_at
from .descartes import _unique_root


@fieldwise_init
struct _FieldExtension(Copyable, Movable):
    var child: _SampleField
    var parent_generator: DensePoly[QQ]


@fieldwise_init
struct _FieldMapEntry(Copyable, Movable):
    var parent_modulus: DensePoly[QQ]
    var parent_variables: List[Int64]
    var parent_coordinates: List[DensePoly[QQ]]
    var relative: List[DensePoly[QQ]]
    var extension: _FieldExtension


struct _FieldExtensionCache(Copyable, Movable):
    var entries: List[_FieldMapEntry]
    var minimum_degree: Int
    var replacement: Int
    var hits: Int

    def __init__(out self, minimum_degree: Int = 16):
        self.entries = List[_FieldMapEntry]()
        self.minimum_degree = minimum_degree
        self.replacement = 0
        self.hits = 0


def _sample_generated_extension(parent: _SampleField, rows: List[DensePoly[QQ]], sample: RealAlgebraic, variable: Int64, mut cache: _FieldExtensionCache, mut budget: Budget) raises PolyError -> Optional[_FieldExtension]:
    """Use beta itself when the relative equation determines gamma over QQ(beta)."""
    ref modulus = sample.polynomial
    if modulus.degree() < parent.modulus.degree():
        return None
    # The polynomial maps do not depend on the chosen real conjugate. The
    # new embedding does, and must replace the cached one. This producer-only
    # cache is discarded before return; final verification checks every map
    # and real embedding independently of its contents.
    var cacheable = modulus.degree() >= cache.minimum_degree
    if cacheable:
        for entry in cache.entries:
            budget.charge(1)
            ref child = entry.extension.child
            if child.modulus == modulus and child.variables[len(child.variables) - 1] == variable and entry.parent_modulus == parent.modulus and entry.parent_variables == parent.variables and entry.parent_coordinates == parent.coordinates and entry.relative == rows:
                cache.hits += 1
                var result = entry.extension.copy()
                result.child.embedding = sample.copy()
                return result^
    # Transpose P(gamma,beta) into QQ(beta)[gamma]. A linear gcd with the
    # parent's modulus gives gamma=c(beta), hence one field for all coordinates.
    # This is the sample-first branch of the historical _child_frame producer.
    var width = 0
    for row in rows:
        width = max(width, len(row.coefficients))
    _allocation[DensePoly[QQ]](width, budget)
    var specialized = List[DensePoly[QQ]]()
    for k in range(width):
        _allocation[QQ](len(rows), budget)
        var entries = List[QQ]()
        for row in rows:
            entries.append(row.coefficients[k].copy() if k < len(row.coefficients) else QQ.zero(NoContext()))
        specialized.append(_reduce(DensePoly[QQ].from_coefficients(0, NoContext(), entries^, budget), modulus, budget))
    _allocation[DensePoly[QQ]](len(parent.modulus.coefficients), budget)
    var definition = List[DensePoly[QQ]]()
    for c in parent.modulus.coefficients:
        definition.append(_constant(c, budget))
    # A high-degree primitive field can need more reconstructed coefficient
    # bits than the default allowance provides. Keep this bounded and charge
    # the caller's budget; exact division still verifies every candidate.
    var common = _gcd_over_field(definition, specialized, modulus, budget, prime_limit=4096 if modulus.degree() >= 64 else 512)
    if len(common) != 2:
        return None
    var image = common[0].neg(budget)
    var child = _SampleField(modulus.copy(), sample.copy(), parent.variables.copy(), List[DensePoly[QQ]]())
    _allocation[Int64](1, budget)
    child.variables.append(variable)
    _allocation[DensePoly[QQ]](_sum(len(parent.coordinates), 1), budget)
    for coordinate in parent.coordinates:
        child.coordinates.append(_compose(coordinate, image, modulus, budget))
    child.coordinates.append(DensePoly[QQ].from_coefficients(0, NoContext(), [QQ.zero(NoContext()), QQ.one(NoContext())], budget))
    var result = _FieldExtension(child^, image^)
    if cacheable:
        _allocation[_FieldMapEntry](1, budget)
        var entry = _FieldMapEntry(parent.modulus.copy(), parent.variables.copy(), parent.coordinates.copy(), rows.copy(), result.copy())
        # Retain only recent expensive maps, with bounded lookup and storage.
        if len(cache.entries) < 128:
            cache.entries.append(entry^)
        else:
            cache.entries[cache.replacement] = entry^
            cache.replacement = (cache.replacement + 1) % 128
    return result^


def _absolute_section(stack: _FieldStack, index: Int, mut parent: _SampleField, mut budget: Budget) raises PolyError -> RealAlgebraic:
    """One conversion at a nonterminal section, never at a coefficient."""
    var root = stack.roots[index].copy()
    if root.low == root.high:
        return _rational(root.low, budget)
    ref rows = stack.basis[root.basis]
    var norm: DensePoly[QQ]
    if parent.embedding:
        norm = _combination_norm(rows, parent.modulus, 0, budget)
    else:
        _allocation[QQ](len(rows), budget)
        var coefficients = List[QQ]()
        for r in rows:
            coefficients.append(r.coefficients[0].copy() if not r.is_zero() else QQ.zero(NoContext()))
        norm = DensePoly[QQ].from_coefficients(0, NoContext(), coefficients^, budget)
    while True:
        budget.charge(1)
        if root.low == root.high:
            return _rational(root.low, budget)
        if _sign_at(norm, root.low, budget) and _sign_at(norm, root.high, budget) and _unique_root(norm, root.low, root.high, budget):
            return _root(norm, root.low, root.high, budget)
        _krefine(root, stack.basis, parent, budget)


def _extend_section(stack: _FieldStack, index: Int, variable: Int64, mut parent: _SampleField, mut budget: Budget) raises PolyError -> _FieldExtension:
    var cache = _FieldExtensionCache()
    return _extend_section_cached(stack, index, variable, parent, cache, budget)


def _extend_section_owned(stack: _FieldStack, index: Int, variable: Int64, var parent: _SampleField, mut cache: _FieldExtensionCache, mut budget: Budget) raises PolyError -> _FieldExtension:
    """Reuse a private parent copy when a linear section stays in its field."""
    ref rows = stack.basis[stack.roots[index].basis]
    if len(rows) == 2:
        var image = DensePoly[QQ](0, NoContext())
        if parent.embedding:
            image = DensePoly[QQ].from_coefficients(0, NoContext(), [QQ.zero(NoContext()), QQ.one(NoContext())], budget)
        parent.variables.append(variable)
        parent.coordinates.append(rows[0].neg(budget))
        return _FieldExtension(parent^, image^)
    return _extend_section_cached(stack, index, variable, parent, cache, budget)


def _extend_section_cached(stack: _FieldStack, index: Int, variable: Int64, mut parent: _SampleField, mut cache: _FieldExtensionCache, mut budget: Budget) raises PolyError -> _FieldExtension:
    var z = DensePoly[QQ].from_coefficients(0, NoContext(), [QQ.zero(NoContext()), QQ.one(NoContext())], budget)
    ref rows = stack.basis[stack.roots[index].basis]
    if len(rows) == 2:
        # A monic linear section already belongs to the parent field.
        var child = parent.copy()
        child.variables.append(variable)
        child.coordinates.append(rows[0].neg(budget))
        return _FieldExtension(child^, z^ if parent.embedding else DensePoly[QQ](0, NoContext()))
    # Preserve the original embedding on same-field children even when root
    # selection refines the working parent. Defer copying the coordinate maps
    # until a branch actually needs them; new-field branches transport them.
    var original_embedding = parent.embedding.copy()
    var sample = _absolute_section(stack, index, parent, budget)
    if sample.rational:
        var child = parent.copy()
        child.embedding = original_embedding^
        child.variables.append(variable)
        child.coordinates.append(_constant(sample.low, budget))
        return _FieldExtension(child^, z^ if parent.embedding else DensePoly[QQ](0, NoContext()))
    if not parent.embedding:
        var child = _SampleField(sample.polynomial.copy(), sample.copy(), parent.variables.copy(), parent.coordinates.copy())
        child.variables.append(variable)
        child.coordinates.append(z^)
        return _FieldExtension(child^, DensePoly[QQ](0, NoContext()))
    var generated = _sample_generated_extension(parent, rows, sample, variable, cache, budget)
    if generated:
        return generated.take()
    var frame = _Frame()
    frame.generator = parent.embedding.copy()
    frame.variables = parent.variables.copy()
    frame.representations = parent.coordinates.copy()
    # Transport gamma as an extra temporary coordinate. The resulting image
    # is retained once per extension so the verifier can bind the two fields.
    frame.representations.append(z^)
    var t = 1
    while True:
        budget.charge(1)
        var combined = _combined_frame(frame, stack.basis[stack.roots[index].basis], sample, t, budget)
        if combined:
            var result = combined.take()
            var beta = result.representations.pop()
            var image = result.representations.pop()
            result.representations.append(beta^)
            result.variables.append(variable)
            return _FieldExtension(_sample_field(result, budget), image^)
        t = _sum(t, 1)


def _check_extension_data(parent: _SampleField, basis: List[DensePoly[QQ]], root: _FinalRoot, variable: Int64, child: _SampleField, image: DensePoly[QQ], mut budget: Budget, *, modulus_checked: Bool = False) raises PolyError -> Bool:
    """Bind a checked parent and checked relative root to one child field."""
    try:
        if not _check_sample_field(child, budget, modulus_checked=modulus_checked) or len(child.variables) != _sum(len(parent.variables), 1):
            return False
        if child.variables[len(parent.variables)] != variable:
            return False
        var checked = child.copy()
        image.check_invariants(budget)
        if image.variable != 0 or image.degree() >= child.modulus.degree():
            return False
        if not _compose(parent.modulus, image, child.modulus, budget).is_zero():
            return False
        if parent.embedding:
            ref gamma = parent.embedding.value()
            if _ksign(image.sub(_constant(gamma.low, budget), budget), checked, budget) <= 0 or _ksign(image.sub(_constant(gamma.high, budget), budget), checked, budget) >= 0:
                return False
        elif not image.is_zero():
            return False
        for i in range(len(parent.variables)):
            if child.variables[i] != parent.variables[i] or child.coordinates[i] != _compose(parent.coordinates[i], image, child.modulus, budget):
                return False
        ref beta = child.coordinates[len(parent.variables)]
        var value = DensePoly[QQ](0, NoContext())
        for k in range(len(basis) - 1, -1, -1):
            value = _field_mul(value, beta, child.modulus, budget).add(_compose(basis[k], image, child.modulus, budget), budget)
        if not value.is_zero():
            return False
        var lower = _ksign(beta.sub(_constant(root.low, budget), budget), checked, budget)
        var upper = _ksign(beta.sub(_constant(root.high, budget), budget), checked, budget)
        return lower == 0 if root.low == root.high else lower > 0 and upper < 0
    except e:
        if e.kind == PolyErrorKind.DOMAIN:
            return False
        raise e^


def _check_extension(parent: _SampleField, basis: List[DensePoly[QQ]], root: _FinalRoot, variable: Int64, extension: _FieldExtension, mut budget: Budget, *, modulus_checked: Bool = False) raises PolyError -> Bool:
    return _check_extension_data(parent, basis, root, variable, extension.child, extension.parent_generator, budget, modulus_checked=modulus_checked)
