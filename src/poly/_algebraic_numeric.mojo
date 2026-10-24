"""Exact-root refinement, certainly rounded conversion and retained enclosure witnesses."""
from apn_mojo import Integer, Rational, Float, Complex as RoundedComplex, ArithmeticContext, ComplexContext, BallContext, to_float_if_certain
from polynomials import Budget, PolyError, PolyErrorKind, QQ, CoefficientRing, DensePoly, SparsePoly, RealAlgebraic, check_real_algebraic, real_rational, evaluate_dense
from polynomials import NumericDomain, FloatDomain, BallDomain, ComplexBallDomain, NumericBox, numeric_value, numeric_bounds, NumericDensePoly, NumericSparsePoly, NumericTerm, MonomialOrder
from ._algebraic_value import Real, Complex, _copy_witness, _plain
from ._common import _allocation, _charge_size, _sum, _product


@fieldwise_init
struct _Witness(Copyable, Movable):
    var real: RealAlgebraic
    var imag: RealAlgebraic


def _initial(value: Real, mut budget: Budget) raises PolyError -> RealAlgebraic:
    if not check_real_algebraic(value.witness(), budget):
        raise PolyError(PolyErrorKind.DOMAIN, "invalid algebraic conversion witness")
    var normalized = _plain(_copy_witness(value.witness(), budget), budget)
    if value != normalized or value.is_zero() != normalized.is_zero():
        raise PolyError(PolyErrorKind.DOMAIN, "invalid algebraic conversion cache")
    # A nonminimal witness can describe a rational value. Recognize it before
    # deciding rounding traps; nearby irrational points would signal inexact.
    if normalized._canonical[].polynomial.degree() == 1:
        return real_rational(normalized._canonical[].polynomial.coefficients[0].neg(budget), budget)
    return _copy_witness(value.witness(), budget)


def _width(bits: Int, mut budget: Budget) raises PolyError -> QQ:
    if bits < 1:
        raise PolyError(PolyErrorKind.DOMAIN, "algebraic conversion needs positive precision")
    _charge_size(_sum(bits, 1), budget)
    try:
        return QQ(Rational(Integer(1), Integer(1) << bits))
    except e:
        raise PolyError(PolyErrorKind.CAPACITY, String(e))


def _refine(var witness: RealAlgebraic, width: QQ, mut budget: Budget) raises PolyError -> RealAlgebraic:
    if witness.rational:
        return witness^
    # A checked witness changes sign across its interval at its single root.
    var low_sign = evaluate_dense(witness.polynomial, witness.low, budget).value.sign()
    while witness.high.sub(witness.low, budget).sub(width, budget).value.sign() > 0:
        budget.charge(1)
        var middle = witness.low.add(witness.high, budget).mul(QQ(Rational(2)).inverse(budget), budget)
        var middle_sign = evaluate_dense(witness.polynomial, middle, budget).value.sign()
        if middle_sign == 0:
            return real_rational(middle, budget)
        if middle_sign == low_sign:
            witness.low = middle^
        else:
            witness.high = middle^
    if not check_real_algebraic(witness, budget):
        raise PolyError(PolyErrorKind.DOMAIN, "algebraic refinement failed verification")
    return witness^



def _magnitude(value: QQ, mut budget: Budget) raises PolyError -> QQ:
    if value.value.sign() < 0:
        return value.neg(budget)
    _charge_size(value.size(), budget)
    return value.copy()


def _localize(var witness: RealAlgebraic, mut budget: Budget) raises PolyError -> RealAlgebraic:
    # Locate a nonzero value to within a factor of two before spending the
    # requested mantissa bits. Exponent range does not consume max_precision.
    while not witness.rational:
        budget.charge(1)
        var low = _magnitude(witness.low, budget)
        var high = _magnitude(witness.high, budget)
        if low.sub(high, budget).value.sign() > 0:
            swap(low, high)
        if witness.low.value.sign() == witness.high.value.sign() and not low.is_zero() and high.sub(low.mul(QQ(Rational(2)), budget), budget).value.sign() <= 0:
            break
        var target = witness.high.sub(witness.low, budget).mul(QQ(Rational(2)).inverse(budget), budget)
        witness = _refine(witness^, target, budget)
    return witness^


def _relative_width(witness: RealAlgebraic, bits: Int, mut budget: Budget) raises PolyError -> QQ:
    var scale = _magnitude(witness.low, budget)
    var high = _magnitude(witness.high, budget)
    if high.sub(scale, budget).value.sign() > 0:
        scale = high^
    if scale.is_zero():
        scale = QQ(Rational(1))
    return scale.mul(_width(bits, budget), budget)


def _same(source: RealAlgebraic, refined: RealAlgebraic, mut budget: Budget) raises PolyError -> Bool:
    if not check_real_algebraic(source, budget) or not check_real_algebraic(refined, budget):
        return False
    if refined.low.sub(source.low, budget).value.sign() < 0 or source.high.sub(refined.high, budget).value.sign() < 0:
        return False
    if refined.rational:
        return evaluate_dense(source.polynomial, refined.low, budget).is_zero()
    return not source.rational and source.polynomial.sub(refined.polynomial, budget).is_zero()


def _box(proof: _Witness, mut budget: Budget) raises PolyError -> NumericBox:
    _charge_size(_sum(_sum(proof.real.low.size(), proof.real.high.size()), _sum(proof.imag.low.size(), proof.imag.high.size())), budget)
    return NumericBox(proof.real.low.copy(), proof.real.high.copy(), proof.imag.low.copy(), proof.imag.high.copy())


def _parts[C: CoefficientRing](value: C, mut budget: Budget) raises PolyError -> _Witness:
    comptime if C == Real:
        return _Witness(_initial(rebind[Real](value), budget), real_rational(QQ(Rational(0)), budget))
    else:
        return _Witness(_initial(rebind[Complex](value).real(), budget), _initial(rebind[Complex](value).imag(), budget))


def _check_parts[C: CoefficientRing, D: NumericDomain](source: C, value: D.Scalar, real: RealAlgebraic, imag: RealAlgebraic, mut budget: Budget) raises PolyError -> Bool:
    comptime if C == Real:
        if not _same(rebind[Real](source).witness(), real, budget) or not imag.rational or not imag.low.is_zero() or not check_real_algebraic(imag, budget):
            return False
    else:
        if not _same(rebind[Complex](source).real().witness(), real, budget) or not _same(rebind[Complex](source).imag().witness(), imag, budget):
            return False
    _charge_size(_sum(_sum(real.low.size(), real.high.size()), _sum(imag.low.size(), imag.high.size())), budget)
    var box = NumericBox(real.low.copy(), real.high.copy(), imag.low.copy(), imag.high.copy())
    return numeric_bounds[D](value, budget).contains(box, budget)


def _check[C: CoefficientRing, D: NumericDomain](source: C, value: D.Scalar, proof: _Witness, mut budget: Budget) raises PolyError -> Bool:
    return _check_parts[C, D](source, value, proof.real, proof.imag, budget)


def _enclose[C: CoefficientRing, D: NumericDomain](value: C, context: D.Context, mut budget: Budget) raises PolyError -> Tuple[D.Scalar, _Witness]:
    comptime assert D == BallDomain or D == ComplexBallDomain, "algebraic enclosure needs a ball domain"
    var c = rebind[BallContext](context)
    if not c.precision():
        raise PolyError(PolyErrorKind.DOMAIN, "choose an explicit ball precision")
    var bits = c.precision().value()
    var proof = _parts(value, budget)
    comptime if D == BallDomain:
        if not proof.imag.rational or not proof.imag.low.is_zero():
            raise PolyError(PolyErrorKind.DOMAIN, "nonreal algebraic coefficient in a real domain")
    var re = _localize(_copy_witness(proof.real, budget), budget)
    var im = _localize(_copy_witness(proof.imag, budget), budget)
    var real_width = _relative_width(re, bits, budget)
    var imag_width = _relative_width(im, bits, budget)
    proof.real = _refine(re^, real_width, budget)
    proof.imag = _refine(im^, imag_width, budget)
    var result = numeric_value[D](_box(proof, budget), context, budget)
    if not _check[C, D](value, result, proof, budget):
        raise PolyError(PolyErrorKind.DOMAIN, "algebraic conversion failed enclosure verification")
    return result^, proof^



def _coarse[C: CoefficientRing, D: NumericDomain](value: C, context: D.Context, mut budget: Budget) raises PolyError -> Tuple[D.Scalar, _Witness]:
    # Unrefined placeholders validate original inputs during stream preparation.
    # Refinement events replace these values before the arithmetic graph runs.
    var proof = _parts(value, budget)
    var result = numeric_value[D](_box(proof, budget), context, budget)
    if not _check[C, D](value, result, proof, budget):
        raise PolyError(PolyErrorKind.DOMAIN, "algebraic stream input failed verification")
    return result^, proof^


def _round(value: Real, context: ArithmeticContext, mut budget: Budget) raises PolyError -> Float:
    var witness = _localize(_initial(value, budget), budget)
    if witness.rational:
        return numeric_value[FloatDomain](NumericBox.point(witness.low, QQ(Rational(0)), budget), context, budget)
    var p = context.format().precision()
    var ceiling = context.max_precision().value() if context.max_precision() else max(_product(8, p), _sum(p, 4096))
    var bits = min(_sum(p, 8), ceiling)
    while True:
        budget.charge(1)
        var width = _relative_width(witness, bits, budget)
        witness = _refine(witness^, width, budget)
        if witness.rational:
            return numeric_value[FloatDomain](NumericBox.point(witness.low, QQ(Rational(0)), budget), context, budget)
        var ball_context: BallContext
        try:
            ball_context = BallContext(precision=bits)
        except e:
            raise PolyError(PolyErrorKind.DOMAIN, String(e))
        var ball = numeric_value[BallDomain](NumericBox.real_interval(witness.low, witness.high, budget), ball_context, budget)
        _charge_size(_sum(bits, p), budget)
        var certain: Optional[Float]
        try:
            certain = to_float_if_certain(ball, context=context)
        except e:
            raise PolyError(PolyErrorKind.DOMAIN, String(e))
        if certain:
            var result = certain.take()
            if not result.is_finite():
                raise PolyError(PolyErrorKind.DOMAIN, "algebraic value does not round to a finite Float")
            return result^
        if bits == ceiling:
            raise PolyError(PolyErrorKind.DOMAIN, "algebraic rounding exceeded context max_precision")
        bits = ceiling if bits > ceiling // 2 else bits * 2


def _value[D: NumericDomain, C: CoefficientRing](value: C, context: D.Context, mut budget: Budget) raises PolyError -> D.Scalar:
    comptime if D == BallDomain or D == ComplexBallDomain:
        var pair = _enclose[C, D](value, context, budget)
        # A list provides an owned slot without assuming a default scalar value.
        _allocation[D.Scalar](1, budget)
        _ = numeric_bounds[D](pair[0], budget)
        return pair[0].copy()
    elif D == FloatDomain:
        comptime if C == Real:
            return rebind_var[D.Scalar](_round(rebind[Real](value), rebind[ArithmeticContext](context), budget))
        else:
            ref z = rebind[Complex](value)
            var imaginary = _initial(z.imag(), budget)
            if not imaginary.rational or not imaginary.low.is_zero():
                raise PolyError(PolyErrorKind.DOMAIN, "nonreal algebraic coefficient in a real domain")
            return rebind_var[D.Scalar](_round(z.real(), rebind[ArithmeticContext](context), budget))
    else:
        var c = rebind[ComplexContext](context)
        var real: Float
        var imag: Float
        comptime if C == Real:
            real = _round(rebind[Real](value), c.real(), budget)
            imag = numeric_value[FloatDomain](NumericBox.point(QQ(Rational(0)), QQ(Rational(0)), budget), c.imag(), budget)
        else:
            real = _round(rebind[Complex](value).real(), c.real(), budget)
            imag = _round(rebind[Complex](value).imag(), c.imag(), budget)
        try:
            return rebind_var[D.Scalar](RoundedComplex(real, imag, context=c))
        except e:
            raise PolyError(PolyErrorKind.DOMAIN, String(e))


def _dense[C: CoefficientRing, D: NumericDomain](p: DensePoly[C], context: D.Context, mut budget: Budget, tighten: Bool = True) raises PolyError -> Tuple[NumericDensePoly[D], List[_Witness]]:
    p.check_invariants(budget)
    _allocation[D.Scalar](len(p.coefficients), budget)
    _allocation[_Witness](len(p.coefficients), budget)
    var values = List[D.Scalar](capacity=len(p.coefficients))
    var proofs = List[_Witness](capacity=len(p.coefficients))
    for coefficient_index in range(len(p.coefficients)):
        ref c = p.coefficients[coefficient_index]
        var converted = _enclose[C, D](c, context, budget) if tighten else _coarse[C, D](c, context, budget)
        _ = numeric_bounds[D](converted[0], budget)
        values.append(converted[0].copy())
        proofs.append(_Witness(_copy_witness(converted[1].real, budget), _copy_witness(converted[1].imag, budget)))
    return NumericDensePoly[D].from_coefficients(p.variable, values^, budget), proofs^


def _sparse[C: CoefficientRing, D: NumericDomain](p: SparsePoly[C], context: D.Context, mut budget: Budget, tighten: Bool = True) raises PolyError -> Tuple[NumericSparsePoly[D], List[_Witness]]:
    p.check_invariants(budget)
    _allocation[NumericTerm[D]](len(p.terms), budget)
    _allocation[_Witness](len(p.terms), budget)
    var terms = List[NumericTerm[D]](capacity=len(p.terms))
    var proofs = List[_Witness](capacity=len(p.terms))
    for t in p.terms:
        var converted = _enclose[C, D](t.coefficient, context, budget) if tighten else _coarse[C, D](t.coefficient, context, budget)
        _charge_size(len(t.monomial.pairs), budget)
        _ = numeric_bounds[D](converted[0], budget)
        terms.append(NumericTerm[D](t.monomial.copy(), converted[0].copy()))
        proofs.append(_Witness(_copy_witness(converted[1].real, budget), _copy_witness(converted[1].imag, budget)))
    return NumericSparsePoly[D].from_terms(terms^, p.order, context, budget), proofs^

comptime _Numeric[D: NumericDomain, sparse: Bool] = NumericSparsePoly[D] if sparse else NumericDensePoly[D]
comptime _Exact[C: CoefficientRing, sparse: Bool] = SparsePoly[C] if sparse else DensePoly[C]


@fieldwise_init
struct _Polynomial[D: NumericDomain, sparse: Bool](Copyable, Movable):
    var value: _Numeric[Self.D, Self.sparse]
    var witnesses: List[_Witness]


def _convert[C: CoefficientRing, D: NumericDomain, sparse: Bool](p: _Exact[C, sparse], context: D.Context, mut budget: Budget, tighten: Bool = True) raises PolyError -> _Polynomial[D, sparse]:
    var proofs = List[_Witness]()
    comptime if sparse:
        var data = _sparse[C, D](rebind[SparsePoly[C]](p), context, budget, tighten)
        var value = NumericSparsePoly[D](MonomialOrder())
        swap(value, data[0])
        swap(proofs, data[1])
        return _Polynomial[D, sparse](rebind_var[_Numeric[D, sparse]](value^), proofs^)
    else:
        var data = _dense[C, D](rebind[DensePoly[C]](p), context, budget, tighten)
        var value = NumericDensePoly[D]()
        swap(value, data[0])
        swap(proofs, data[1])
        return _Polynomial[D, sparse](rebind_var[_Numeric[D, sparse]](value^), proofs^)


def _check_polynomial[C: CoefficientRing, D: NumericDomain, sparse: Bool](p: _Exact[C, sparse], result: _Polynomial[D, sparse], mut budget: Budget) raises PolyError -> Bool:
    comptime if sparse:
        ref source = rebind[SparsePoly[C]](p)
        source.check_invariants(budget)
        ref value = rebind[NumericSparsePoly[D]](result.value)
        _charge_size(_sum(len(source.order.variables), len(value.order.variables)), budget)
        if source.order != value.order or len(source.terms) != len(value.terms) or len(result.witnesses) != len(source.terms):
            return False
        for i in range(len(source.terms)):
            budget.charge(1)
            value.terms[i].monomial.check_invariants(budget)
            if source.order.compare(source.terms[i].monomial, value.terms[i].monomial, budget) != 0 or not _check[C, D](source.terms[i].coefficient, value.terms[i].coefficient, result.witnesses[i], budget):
                return False
    else:
        ref source = rebind[DensePoly[C]](p)
        source.check_invariants(budget)
        ref value = rebind[NumericDensePoly[D]](result.value)
        if source.variable != value.variable or len(source.coefficients) != len(value.coefficients) or len(result.witnesses) != len(source.coefficients):
            return False
        for i in range(len(source.coefficients)):
            budget.charge(1)
            if not _check[C, D](source.coefficients[i], value.coefficients[i], result.witnesses[i], budget):
                return False
    return True
