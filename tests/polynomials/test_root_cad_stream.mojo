"""Complete and corrupted incremental root and cylindrical witnesses."""
from apn_mojo import Rational, BallContext
from polynomials import *
from support.check import Checker


@fieldwise_init
struct RealCollector(DescartesRootSink):
    var squarefree: DensePoly[ZZ]
    var cofactor: DensePoly[ZZ]
    var bound_exponent: Int
    var items: List[DescartesRootItem]

    def begin(mut self, squarefree: DensePoly[ZZ], cofactor: DensePoly[ZZ], bound_exponent: Int, mut budget: Budget) raises PolyError:
        budget.charge(1)
        self.squarefree = squarefree.copy()
        self.cofactor = cofactor.copy()
        self.bound_exponent = bound_exponent

    def write(mut self, var item: DescartesRootItem, mut budget: Budget) raises PolyError:
        budget.charge(1)
        self.items.append(item^)


@fieldwise_init
struct ComplexCollector(ComplexRootSink):
    var regions: List[ComplexRootRegion]

    def write(mut self, var region: ComplexRootRegion, mut budget: Budget) raises PolyError:
        budget.charge(1)
        self.regions.append(region^)


@fieldwise_init
struct CADCollector(CADSink):
    var projection: Optional[CADProjectionCertificate]
    var events: List[CADStreamEvent]

    def begin(mut self, projection: CADProjectionCertificate, mut budget: Budget) raises PolyError:
        budget.charge(1)
        self.projection = projection.copy()

    def write(mut self, var event: CADStreamEvent, mut budget: Budget) raises PolyError:
        budget.charge(1)
        self.events.append(event^)


def main() raises:
    var c = Checker("root and CAD streams")
    var work = Budget()
    var f = DensePoly[QQ].from_coefficients(0, NoContext(), [QQ(Rational(-2)), QQ(Rational(0)), QQ(Rational(1))], work)
    var options = NumericOptions(absolute_bits=2, max_subdivisions=1024)
    var checker = DescartesStreamChecker.from_dense(f, options, work)
    var receipt = real_roots_exact_to(f, options, checker, work)
    c.check(checker.finish(receipt, work), "direct real producer to independent checker")
    receipt.items += 1
    c.check(not checker.finish(receipt, work), "real receipt count binding")
    var collector = RealCollector(DensePoly[ZZ](0, NoContext()), DensePoly[ZZ](0, NoContext()), 0, List[DescartesRootItem]())
    receipt = real_roots_exact_to(f, options, collector, work)
    var materialized = real_roots_exact(f, options, work)
    c.check(len(materialized.certificate.leaves) == len(collector.items) and receipt.items == len(collector.items), "same real traversal")
    checker = DescartesStreamChecker.from_dense(f, options, work)
    checker.begin(collector.squarefree, collector.cofactor, collector.bound_exponent, work)
    c.check(not checker.finish(receipt, work), "empty prefix cannot prove completeness")
    for i in range(len(collector.items)):
        ref leaf = materialized.certificate.leaves[i]
        c.check(collector.items[i].leaf.low == leaf.low and collector.items[i].leaf.high == leaf.high and collector.items[i].leaf.variations == leaf.variations, "materialized/streamed real leaves agree")
        c.check(checker.accept(collector.items[i].copy(), work), "consume real leaf")
    c.check(checker.finish(receipt, work), "independent real replay completion")
    c.check(not checker.accept(collector.items[0].copy(), work) and not checker.finish(receipt, work), "extra real witness permanently rejected")
    checker = DescartesStreamChecker.from_dense(f, options, work)
    checker.begin(collector.squarefree, collector.cofactor, collector.bound_exponent, work)
    var corrupt = collector.items[0].copy()
    corrupt.leaf.variations += 1
    c.check(not checker.accept(corrupt^, work), "corrupt local real bound")
    checker = DescartesStreamChecker.from_dense(f, options, work)
    checker.begin(collector.squarefree, collector.cofactor, collector.bound_exponent, work)
    c.check(not checker.accept(collector.items[1].copy(), work), "a leaf out of order breaks coverage")
    checker = DescartesStreamChecker.from_dense(f, options, work)
    checker.begin(collector.squarefree, collector.cofactor, collector.bound_exponent, work)
    var tiny = Budget(0)
    try:
        _ = checker.accept(collector.items[0].copy(), tiny)
        c.check(False, "root stream budget ignored")
    except e:
        c.check(e.kind == PolyErrorKind.LIMIT and not checker.finish(receipt, work), "failed root stream cannot complete")
    var bad_squarefree = collector.squarefree.copy()
    bad_squarefree.coefficients[0] = ZZ(7)
    checker = DescartesStreamChecker.from_dense(f, options, work)
    try:
        checker.begin(bad_squarefree, collector.cofactor, collector.bound_exponent, work)
        c.check(False, "wrong squarefree header accepted")
    except e:
        c.check(e.kind == PolyErrorKind.DOMAIN, "independent squarefree header validation")
    checker = DescartesStreamChecker.from_dense(f, options, work)
    try:
        checker.begin(collector.squarefree, collector.cofactor, collector.bound_exponent - 1, work)
        c.check(False, "too small a bound accepted")
    except e:
        c.check(e.kind == PolyErrorKind.DOMAIN, "independent bound validation")
    var complex_checker = ComplexRootStreamChecker.from_dense(f, options, work)
    receipt = complex_roots_exact_to(f, options, complex_checker, work)
    c.check(complex_checker.finish(receipt, work), "direct complex producer to checker")
    var complex_collector = ComplexCollector(List[ComplexRootRegion]())
    receipt = complex_roots_exact_to(f, options, complex_collector, work)
    complex_checker = ComplexRootStreamChecker.from_dense(f, options, work)
    c.check(not complex_checker.finish(receipt, work), "missing complex regions")
    c.check(complex_checker.accept(complex_collector.regions[0].copy(), work), "complex local count")
    c.check(not complex_checker.accept(complex_collector.regions[0].copy(), work), "overlapping duplicate region rejected")
    var bad_region = complex_collector.regions[0].copy()
    bad_region.count += 1
    complex_checker = ComplexRootStreamChecker.from_dense(f, options, work)
    c.check(not complex_checker.accept(bad_region^, work), "corrupt winding count rejected")
    bad_region = complex_collector.regions[0].copy()
    _ = bad_region.boundary.pop()
    complex_checker = ComplexRootStreamChecker.from_dense(f, options, work)
    c.check(not complex_checker.accept(bad_region^, work), "truncated contour rejected")
    var constant = DensePoly[QQ].from_coefficients(0, NoContext(), [QQ(Rational(3))], work)
    complex_checker = ComplexRootStreamChecker.from_dense(constant, options, work)
    receipt = complex_roots_exact_to(constant, options, complex_checker, work)
    c.check(complex_checker.finish(receipt, work) and receipt.items == 0, "constant has empty complex root stream")
    var numeric = numeric_from_dense[BallDomain](f, BallContext(128), work)
    var numeric_checker = RealRootStreamChecker.from_numeric(numeric, options, work)
    receipt = real_roots_numeric_to(numeric, options, numeric_checker, work)
    c.check(numeric_checker.finish(receipt, work), "numerical real stream")
    complex_checker = ComplexRootStreamChecker.from_numeric(numeric, options, work)
    receipt = complex_roots_numeric_to(numeric, options, complex_checker, work)
    c.check(complex_checker.finish(receipt, work), "numerical complex stream")
    var x = SparsePoly[QQ].variable(0, MonomialOrder(), NoContext(), work)
    var inputs: List[SparsePoly[QQ]] = [x.copy()]
    for lazard in [False, True]:
        var cad_checker = CADStreamChecker.prepare(inputs, work)
        var cad_receipt = cad_to(inputs, [Int64(0)], cad_checker, work, lazard=lazard)
        c.check(cad_checker.finish(cad_receipt, work) and cad_receipt.stacks == 1 and cad_receipt.cells == 3, "three cylindrical cells of x=0")
        cad_receipt.cells += 1
        c.check(not cad_checker.finish(cad_receipt, work), "CAD receipt binding")
    var cad_collector = CADCollector(None, List[CADStreamEvent]())
    var cad_receipt = cad_to(inputs, [Int64(0)], cad_collector, work)
    var cad_checker = CADStreamChecker.prepare(inputs, work)
    cad_checker.begin(cad_collector.projection.value(), work)
    for i in range(len(cad_collector.events) - 1):
        c.check(cad_checker.accept(cad_collector.events[i].copy(), work), "CAD prefix accepted locally")
    c.check(not cad_checker.finish(cad_receipt, work), "CAD missing branch rejected globally")
    c.check(cad_checker.accept(cad_collector.events[len(cad_collector.events) - 1].copy(), work) and cad_checker.finish(cad_receipt, work), "CAD final branch completes coverage")
    c.check(not cad_checker.accept(cad_collector.events[0].copy(), work), "extra CAD event rejected")
    cad_checker = CADStreamChecker.prepare(inputs, work)
    cad_checker.begin(cad_collector.projection.value(), work)
    c.check(not cad_checker.accept(cad_collector.events[1].copy(), work), "out-of-order CAD sign event rejected")
    cad_checker = CADStreamChecker.prepare(inputs, work)
    cad_checker.begin(cad_collector.projection.value(), work)
    var bad_event = cad_collector.events[0].copy()
    _ = bad_event.stack.value().samples.pop()
    c.check(not cad_checker.accept(bad_event^, work), "incomplete stack rejected")
    cad_checker = CADStreamChecker.prepare(inputs, work)
    cad_checker.begin(cad_collector.projection.value(), work)
    c.check(cad_checker.accept(cad_collector.events[0].copy(), work), "valid stack before signs")
    bad_event = cad_collector.events[1].copy()
    bad_event.signs[0].relation = 0
    c.check(not cad_checker.accept(bad_event^, work), "corrupt terminal sign rejected")
    var exact = CADSignWitness(0, List[RealAlgebraic](), -1, List[DensePoly[QQ]]())
    cad_checker = CADStreamChecker.prepare(inputs, work)
    cad_checker.begin(cad_collector.projection.value(), work)
    bad_event = cad_collector.events[0].copy()
    bad_event.witnesses.append(exact.copy())
    c.check(not cad_checker.accept(bad_event^, work), "stack event with witnesses rejected")
    cad_checker = CADStreamChecker.prepare(inputs, work)
    cad_checker.begin(cad_collector.projection.value(), work)
    c.check(cad_checker.accept(cad_collector.events[0].copy(), work), "valid stack before witnesses")
    bad_event = cad_collector.events[1].copy()
    c.check(len(bad_event.witnesses) == len(inputs), "terminal event carries one witness per input")
    bad_event.witnesses.append(exact.copy())
    c.check(not cad_checker.accept(bad_event^, work), "extra terminal witness rejected")
    cad_checker = CADStreamChecker.prepare(inputs, work)
    cad_checker.begin(cad_collector.projection.value(), work)
    _ = cad_checker.accept(cad_collector.events[0].copy(), work)
    bad_event = cad_collector.events[1].copy()
    bad_event.witnesses[0] = exact.copy()
    if len(bad_event.evaluations) == 0:
        c.check(not cad_checker.accept(bad_event^, work), "exact witness without its evaluation rejected")
    cad_checker = CADStreamChecker.prepare(inputs, work)
    cad_checker.begin(cad_collector.projection.value(), work)
    tiny = Budget(0)
    try:
        _ = cad_checker.accept(cad_collector.events[0].copy(), tiny)
        c.check(False, "CAD stream budget ignored")
    except e:
        c.check(e.kind == PolyErrorKind.LIMIT and not cad_checker.finish(cad_receipt, work), "failed CAD stream cannot complete")
    inputs.append(SparsePoly[QQ].variable(1, MonomialOrder(), NoContext(), work))
    cad_checker = CADStreamChecker.prepare(inputs, work)
    cad_receipt = cad_to(inputs, [Int64(0), Int64(1)], cad_checker, work)
    c.check(cad_checker.finish(cad_receipt, work) and cad_receipt.stacks == 4 and cad_receipt.cells == 9, "two-dimensional stack backtracking")
    var empty = List[SparsePoly[QQ]]()
    cad_checker = CADStreamChecker.prepare(empty, work)
    cad_receipt = cad_to(empty, List[Int64](), cad_checker, work)
    c.check(cad_checker.finish(cad_receipt, work) and cad_receipt.cells == 1 and cad_receipt.stacks == 0, "zero-dimensional CAD stream")
    c.finish()
