"""A relative section becomes a checked child field only for another lift."""
from apn_mojo import Rational
from polynomials import Budget, QQ, NoContext, DensePoly, MonomialOrder, real_root
from polynomials._algebraic import _Frame
from polynomials._cad_sample_field import _sample_field, _specialize_field
from polynomials._cad_field_stack import _construct_field_stack, _check_field_stack
from polynomials._cad_field_extension import _FieldExtensionCache, _extend_section, _extend_section_cached, _extend_section_owned, _check_extension
from support.corpus import parse_polynomial
from support.check import Checker


def main() raises:
    var c = Checker("relative section child field")
    var work = Budget()
    var order = MonomialOrder()
    var frame = _Frame()
    var base = _sample_field(frame, work)
    var inputs: List[List[DensePoly[QQ]]] = [_specialize_field(parse_polynomial("x0^2 - 2", order, work), 0, base, work)]
    var stack = _construct_field_stack(inputs, base, work)
    c.check(_check_field_stack(inputs, stack, base, work), "rational-base relative stack")
    var child = _extend_section(stack, 1, 0, base, work)
    c.check(_check_extension(base, stack.basis[stack.roots[1].basis], stack.roots[1], 0, child, work), "positive sqrt(2) child embedding")
    var parent = child.child.copy()
    for expr in ["x1^2 - 3", "x1^2 - x0", "x1 - x0", "x1^2 - 1"]:
        inputs = [_specialize_field(parse_polynomial(expr, order, work), 1, parent, work)]
        stack = _construct_field_stack(inputs, parent, work)
        c.check(_check_field_stack(inputs, stack, parent, work), "relative stack over sqrt(2)")
        var index = len(stack.roots) - 1
        child = _extend_section(stack, index, 1, parent, work)
        c.check(_check_extension(parent, stack.basis[stack.roots[index].basis], stack.roots[index], 1, child, work), "section transports the entire parent field")
        var owned_cache = _FieldExtensionCache()
        var parent_coordinates = parent.coordinates.copy()
        var owned_child = _extend_section_owned(stack, index, 1, parent.copy(), owned_cache, work)
        c.check(_check_extension(parent, stack.basis[stack.roots[index].basis], stack.roots[index], 1, owned_child, work), "owned section extension preserves the checked real embedding")
        c.check(owned_child.child.modulus == child.child.modulus and owned_child.child.coordinates == child.child.coordinates and owned_child.parent_generator == child.parent_generator, "owned and borrowed extensions transport identical coordinates")
        c.check(parent.coordinates == parent_coordinates, "extending an owned sample copy preserves its parent")
        if expr == "x1^2 - x0":
            var simple = DensePoly[QQ].from_coefficients(0, NoContext(), [QQ(Rational(-2)), QQ.zero(NoContext()), QQ.zero(NoContext()), QQ.zero(NoContext()), QQ.one(NoContext())], work)
            c.check(child.child.modulus == simple, "fourth-root sample generates the child without adding sqrt(2)")
            c.check(child.parent_generator == DensePoly[QQ].from_coefficients(0, NoContext(), [QQ.zero(NoContext()), QQ.zero(NoContext()), QQ.one(NoContext())], work), "parent generator is exactly beta squared")
            var negative = _extend_section(stack, 0, 1, parent, work)
            c.check(negative.child.modulus == simple and _check_extension(parent, stack.basis[stack.roots[0].basis], stack.roots[0], 1, negative, work), "the negative embedding uses the same simple child modulus")
            var cache = _FieldExtensionCache(minimum_degree=0)
            var cached_positive = _extend_section_cached(stack, index, 1, parent, cache, work)
            var cached_negative = _extend_section_cached(stack, 0, 1, parent, cache, work)
            c.check(len(cache.entries) == 1 and cache.hits == 1 and cached_positive.parent_generator == cached_negative.parent_generator, "conjugate sections reuse their exact coordinate maps")
            c.check(_check_extension(parent, stack.basis[stack.roots[0].basis], stack.roots[0], 1, cached_negative, work) and not _check_extension(parent, stack.basis[stack.roots[index].basis], stack.roots[index], 1, cached_negative, work), "cache reuse retains the requested real embedding")
            cache.entries[0].extension.parent_generator = DensePoly[QQ](0, NoContext())
            var forged = _extend_section_cached(stack, 0, 1, parent, cache, work)
            c.check(not _check_extension(parent, stack.basis[stack.roots[0].basis], stack.roots[0], 1, forged, work), "producer cache contents cannot bypass final verification")
        var bad = child.copy()
        bad.parent_generator = DensePoly[QQ](0, NoContext())
        c.check(not _check_extension(parent, stack.basis[stack.roots[index].basis], stack.roots[index], 1, bad, work), "forged parent-generator image rejected")
        bad = child.copy()
        bad.child.coordinates[1] = DensePoly[QQ](0, NoContext())
        c.check(not _check_extension(parent, stack.basis[stack.roots[index].basis], stack.roots[index], 1, bad, work), "forged section image rejected")
        var third = child.child.copy()
        var next: List[List[DensePoly[QQ]]] = [_specialize_field(parse_polynomial("x2 - x0 - x1", order, work), 2, third, work)]
        var lifted = _construct_field_stack(next, third, work)
        c.check(len(lifted.roots) == 1 and _check_field_stack(next, lifted, third, work), "further lift uses the transported field")
    c.finish()
