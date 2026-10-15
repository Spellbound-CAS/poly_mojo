"""Mapped representation copies preserve stored values, formats, and metadata."""
from poly import series, basis, rational, sparse, dense, modint
from apn_mojo import Rational, ExactComplex, Float, ArithmeticContext, FloatFormat, BallContext
from polynomials import Budget, PolyErrorKind, NumericSeries, FloatDomain
from std.testing import assert_equal, assert_true


def main() raises:
    var fc = ArithmeticContext(format=FloatFormat(precision=37))
    var lower = ArithmeticContext(format=FloatFormat(precision=11))
    var a = series.asarray([Rational(1, 3), Rational(2)], order=5, context=fc)
    var core = series.to_core(a, variable=-19)
    assert_equal(core.polynomial.variable, Int64(-19))
    var restored = series.from_core(core, context=lower)
    assert_equal(restored.order(), 5)
    assert_equal(restored[0].to_json(), a[0].to_json())
    assert_equal(restored[4].to_rational_exact(), Rational(0))
    assert_equal(series.to_core(restored).polynomial.variable, Int64(0))
    core.polynomial.coefficients[0] = Float(9, context=fc)
    assert_equal(restored[0].to_json(), a[0].to_json())
    var exact_s = series.asarray([1, 2], order=5)
    var exact_core = series.to_core(exact_s, variable=88)
    assert_equal(exact_core.polynomial.variable, Int64(88))
    assert_equal(series.coefficients(series.from_core(exact_core)).to_list(), series.coefficients(exact_s).to_list())
    var b = basis.asarray([Rational(1, 3), Rational(2), Rational(0)], basis="newton", nodes=[Rational(1, 7), Rational(1, 5)], context=fc)
    var bcore = basis.to_core(b, variable=42)
    assert_equal(bcore.variable, Int64(42))
    var bcopy = basis.from_core(bcore, context=lower)
    assert_equal(len(bcopy), 3)
    assert_equal(basis.to_core(bcopy).variable, Int64(0))
    assert_equal(bcopy[0].to_json(), b[0].to_json())
    assert_equal(basis.nodes(bcopy)[0].to_json(), basis.nodes(b)[0].to_json())
    bcore.nodes[0] = Float(9, context=fc)
    assert_equal(basis.nodes(bcopy)[0].to_json(), basis.nodes(b)[0].to_json())
    var exact_b = basis.asarray([1, 2, 0], basis="bernstein")
    var exact_bc = basis.to_core(exact_b, variable=-22)
    assert_equal(exact_bc.variable, Int64(-22))
    assert_equal(basis.coefficients(basis.from_core(exact_bc)).to_list(), basis.coefficients(exact_b).to_list())
    var ids: List[Int64] = [11, -8, 400]
    var f = rational.from_terms([((1, 0, 0), 2), ((0, 1, 0), 3)], [((0, 0, 0), 1), ((0, 1, 0), 1)])
    var fcore = rational.to_core(f, variables=ids)
    assert_equal(fcore.numerator.order.variables, ids)
    assert_equal(rational.numerator(rational.from_core[3](fcore, variables=ids)).to_list(), rational.numerator(f).to_list())
    var zf = rational.asarray([ExactComplex(Rational(1), Rational(2)), ExactComplex(Rational(3))], [ExactComplex(Rational(1))])
    var single_id: List[Int64] = [-33]
    var zcore = rational.to_core(zf, variables=single_id)
    assert_equal(zcore[0].order.variables, single_id)
    assert_equal(rational.numerator(rational.from_core[1](zcore, variables=single_id)).to_list(), rational.numerator(zf).to_list())
    var mf = rational.asarray([modint(2, 7), modint(3, 7)], [modint(1, 7)])
    var mcore = rational.to_core(mf, variables=single_id)
    assert_equal(mcore[0].ring.value, UInt64(7))
    assert_equal(rational.numerator(rational.from_core[1](mcore, variables=single_id)).to_list(), rational.numerator(mf).to_list())
    var nf = rational.asarray(f, context=fc)
    var ncore = rational.to_core(nf, variables=ids)
    assert_equal(ncore[0].order.variables, ids)
    assert_equal(ncore[1].order.variables, ids)
    var ncopy = rational.from_core[3](ncore, variables=ids)
    assert_equal(rational.numerator(ncopy).to_list(), rational.numerator(nf).to_list())
    assert_equal(rational.denominator(ncopy).to_list(), rational.denominator(nf).to_list())
    ncore[0].terms[0].coefficient = Float(9, context=fc)
    assert_equal(rational.numerator(ncopy).to_list(), rational.numerator(nf).to_list())
    var bc = BallContext(precision=47)
    var bf = rational.asarray(f, context=bc)
    var bc_core = rational.to_core(bf, variables=ids)
    assert_equal(rational.numerator(rational.from_core[3](bc_core, variables=ids))[0][1].to_json(), rational.numerator(bf)[0][1].to_json())
    for mode in range(6):
        var rejected = False
        try:
            if mode == 0:
                var bad = series.to_core(a)
                bad.order = -1
                _ = series.from_core(bad, context=fc)
            elif mode == 1:
                var bad = series.to_core(a)
                bad.order = 1
                _ = series.from_core(bad, context=fc)
            elif mode == 2:
                var bad = basis.to_core(b)
                bad.kind = 255
                _ = basis.from_core(bad, context=fc)
            elif mode == 3: _ = rational.from_core[3](ncore)
            elif mode == 4: _ = rational.to_core(nf, variables=List[Int64]([Int64(1), Int64(1), Int64(2)]))
            else:
                var bad = rational.to_core(nf)
                bad[1].terms.clear()
                _ = rational.from_core[3](bad)
        except e: rejected = e.kind == PolyErrorKind.DOMAIN
        assert_true(rejected)
    var meter = Budget()
    _ = series.from_core(core, meter, context=lower)
    var used = UInt64.MAX - meter.remaining()
    var exact = Budget(used)
    _ = series.from_core(core, exact, context=lower)
    assert_equal(exact.remaining(), UInt64(0))
    var short = Budget(used - 1)
    var limited = False
    try: _ = series.from_core(core, short, context=lower)
    except e: limited = e.kind == PolyErrorKind.LIMIT
    assert_true(limited)
    print("surface representation interoperability checks passed")
