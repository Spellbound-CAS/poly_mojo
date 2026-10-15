"""Certified McCallum and Lazard projection with primitive coprime bases.

Source: SMT-RAT fc28c1d, projectionoperator/{McCallum,Lazard}.h.
Using Res(f,f') instead of disc(f) adds only the leading coefficient,
already present in both projection operators.
"""
# Copyright (c) 2015-2022 Florian Corzilius, Gereon Kremer, Jasper Nalbach,
# Sebastian Junges, Stefan Schupp, Erika Abraham.
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
# THE SOFTWARE.
from .errors import Budget, PolyError, PolyErrorKind, _allocation, _sum, _charge_size
from .coefficients import QQ, NoContext
from .sparse import SparsePoly
from .monomial import MonomialOrder
from .evaluate import derivative
from ._views import _split, _coefficient, _lc
from .certificates import CoprimalityCertificate, check_coprimality, _equal
from .gcd import coprimality_certificate
from .multivariate_factor import MultivariateFactorization, factor_multivariate, check_multivariate_factorization, _copy_factorization
from .determinant import ProjectionResultantCertificate, projection_resultant, check_projection_resultant


@fieldwise_init
struct CADProjectionCertificate(Copyable, Movable, Writable):
    """Complete projection closure with determinant and basis normalization proofs."""
    var variables: List[Int64]
    """Variables in base-to-top cylindrical order."""
    var lazard: Bool
    """Whether to use leading/trailing coefficients rather than every coefficient."""
    var levels: List[List[SparsePoly[QQ]]]
    """Primitive squarefree pairwise coprime bases at each main-variable level."""
    var factorizations: List[MultivariateFactorization[QQ]]
    """Product identities for inputs and projection polynomials in schedule order."""
    var resultants: List[ProjectionResultantCertificate]
    """Discriminant and pairwise resultant determinant transcripts."""
    var coprimality: List[CoprimalityCertificate[QQ]]
    """Coprimality of each basis member with its derivative and with earlier members."""


@fieldwise_init
struct _State(Movable):
    var result: CADProjectionCertificate
    var factors: Int
    var resultants: Int
    var coprimality: Int
    var inputs: List[SparsePoly[QQ]]
    var witnesses: List[Int]


def _empty(variables: List[Int64], lazard: Bool, mut budget: Budget) raises PolyError -> CADProjectionCertificate:
    _allocation[Int64](len(variables), budget)
    _allocation[List[SparsePoly[QQ]]](len(variables), budget)
    for i in range(len(variables)):
        budget.charge(1)
        if variables[i] < 0:
            raise PolyError(PolyErrorKind.DOMAIN, "negative CAD variable")
        for j in range(i):
            budget.charge(1)
            if variables[i] == variables[j]:
                raise PolyError(PolyErrorKind.DOMAIN, "duplicate CAD variable")
    return CADProjectionCertificate(variables.copy(), lazard, List[List[SparsePoly[QQ]]](length=len(variables), fill=List[SparsePoly[QQ]]()), List[MultivariateFactorization[QQ]](), List[ProjectionResultantCertificate](), List[CoprimalityCertificate[QQ]]())


def _level(f: SparsePoly[QQ], variables: List[Int64], mut budget: Budget) raises PolyError -> Int:
    var level = -1
    for variable in f.variables(budget):
        budget.charge(1)
        var found = False
        for i in range(len(variables)):
            budget.charge(1)
            if variables[i] == variable:
                level = max(level, i)
                found = True
                break
        if not found:
            raise PolyError(PolyErrorKind.DOMAIN, "CAD variable order does not cover the polynomial")
    return level


def _collect(proof: MultivariateFactorization[QQ], maximum: Int, mut state: _State, mut budget: Budget) raises PolyError:
    if proof.unit.is_zero() and len(proof.factors):
        raise PolyError(PolyErrorKind.DOMAIN, "zero projection input must have no factors")
    for entry in proof.factors:
        budget.charge(1)
        ref factor = entry[0]
        factor.check_invariants(budget)
        var level = _level(factor, state.result.variables, budget)
        if factor.order != MonomialOrder() or level < 0 or level >= maximum or entry[1] <= 0:
            raise PolyError(PolyErrorKind.DOMAIN, "invalid projection factor")
        _charge_size(factor.terms[len(factor.terms) - 1].coefficient.size(), budget)
        if factor.terms[len(factor.terms) - 1].coefficient != QQ.one(NoContext()):
            raise PolyError(PolyErrorKind.DOMAIN, "projection factors must be normalized")
        var found = False
        for known in state.result.levels[level]:
            budget.charge(1)
            if _equal(factor, known, budget):
                found = True
                break
        if not found:
            _allocation[SparsePoly[QQ]](1, budget)
            state.result.levels[level].append(factor._copy(budget))


def _factor[checking: Bool](f: SparsePoly[QQ], maximum: Int, given: CADProjectionCertificate, mut state: _State, mut budget: Budget) raises PolyError:
    comptime if checking:
        if state.factors >= len(given.factorizations) or not check_multivariate_factorization(f, given.factorizations[state.factors], budget):
            raise PolyError(PolyErrorKind.DOMAIN, "invalid projection factorization identity")
        _collect(given.factorizations[state.factors], maximum, state, budget)
        state.factors = _sum(state.factors, 1)
    else:
        var cached = -1
        for i in range(len(state.inputs)):
            budget.charge(1)
            if _equal(f, state.inputs[i], budget):
                cached = state.witnesses[i]
                break
        var proof: MultivariateFactorization[QQ]
        if cached >= 0:
            proof = _copy_factorization(state.result.factorizations[cached], budget)
        else:
            proof = factor_multivariate(f, budget)
            var key = f._copy(budget)
            _allocation[SparsePoly[QQ]](1, budget)
            _allocation[Int](1, budget)
            state.inputs.append(key^)
            state.witnesses.append(len(state.result.factorizations))
        _collect(proof, maximum, state, budget)
        _allocation[MultivariateFactorization[QQ]](1, budget)
        state.result.factorizations.append(proof^)


def _coprime[checking: Bool](f: SparsePoly[QQ], g: SparsePoly[QQ], given: CADProjectionCertificate, mut state: _State, mut budget: Budget) raises PolyError:
    comptime if checking:
        if state.coprimality >= len(given.coprimality) or not check_coprimality(f, g, given.coprimality[state.coprimality], budget):
            raise PolyError(PolyErrorKind.DOMAIN, "projection basis is not primitive squarefree and coprime")
        state.coprimality = _sum(state.coprimality, 1)
    else:
        var proof = coprimality_certificate(f, g, budget)
        if not proof:
            raise PolyError(PolyErrorKind.DOMAIN, "factorization did not yield a projection basis")
        _allocation[CoprimalityCertificate[QQ]](1, budget)
        state.result.coprimality.append(proof.take())


def _resultant[checking: Bool](f: SparsePoly[QQ], g: SparsePoly[QQ], variable: Int64, maximum: Int, given: CADProjectionCertificate, mut state: _State, mut budget: Budget) raises PolyError:
    comptime if checking:
        if state.resultants >= len(given.resultants):
            raise PolyError(PolyErrorKind.DOMAIN, "missing projection determinant")
        ref proof = given.resultants[state.resultants]
        if proof.variable != variable or not check_projection_resultant(f, g, proof, budget):
            raise PolyError(PolyErrorKind.DOMAIN, "incorrect projection determinant")
        _factor[checking](proof.resultant, maximum, given, state, budget)
        state.resultants = _sum(state.resultants, 1)
    else:
        var proof = projection_resultant(f, g, variable, budget)
        _factor[checking](proof.resultant, maximum, given, state, budget)
        _allocation[ProjectionResultantCertificate](1, budget)
        state.result.resultants.append(proof^)


def _project[checking: Bool](inputs: List[SparsePoly[QQ]], variables: List[Int64], lazard: Bool, given: CADProjectionCertificate, mut budget: Budget) raises PolyError -> CADProjectionCertificate:
    var state = _State(_empty(variables, lazard, budget), 0, 0, 0, List[SparsePoly[QQ]](), List[Int]())
    for f in inputs:
        budget.charge(1)
        _ = _level(f, variables, budget)
        _factor[checking](f.reorder(MonomialOrder(), budget), len(variables), given, state, budget)
    for level in range(len(variables) - 1, -1, -1):
        budget.charge(1)
        var variable = variables[level]
        var size = len(state.result.levels[level])
        for i in range(size):
            budget.charge(1)
            # Only lower levels grow during projection; owned copies avoid
            # retaining a list borrow across that mutation.
            var f = state.result.levels[level][i]._copy(budget)
            var differentiated = derivative(f, variable, budget)
            _coprime[checking](f, differentiated, given, state, budget)
            if level:
                if lazard:
                    _factor[checking](_coefficient(f, variable, 0, budget), level, given, state, budget)
                    _factor[checking](_lc(f, variable, budget), level, given, state, budget)
                else:
                    var coefficients = _split(f, variable, budget)
                    for coefficient in coefficients:
                        budget.charge(1)
                        _factor[checking](coefficient, level, given, state, budget)
                _resultant[checking](f, differentiated, variable, level, given, state, budget)
            for j in range(i):
                budget.charge(1)
                var g = state.result.levels[level][j]._copy(budget)
                _coprime[checking](f, g, given, state, budget)
                if level:
                    _resultant[checking](f, g, variable, level, given, state, budget)
    comptime if checking:
        if state.factors != len(given.factorizations) or state.resultants != len(given.resultants) or state.coprimality != len(given.coprimality) or len(given.levels) != len(variables):
            raise PolyError(PolyErrorKind.DOMAIN, "incomplete or extraneous projection witnesses")
        for level in range(len(variables)):
            budget.charge(1)
            if len(given.levels[level]) != len(state.result.levels[level]):
                raise PolyError(PolyErrorKind.DOMAIN, "incorrect projection basis size")
            for i in range(len(given.levels[level])):
                budget.charge(1)
                if not _equal(given.levels[level][i], state.result.levels[level][i], budget):
                    raise PolyError(PolyErrorKind.DOMAIN, "incorrect projection basis member")
    var result = _empty(List[Int64](), False, budget)
    swap(result, state.result)
    return result^


def cad_projection(inputs: List[SparsePoly[QQ]], variables: List[Int64], mut budget: Budget, lazard: Bool = False) raises PolyError -> CADProjectionCertificate:
    """Construct the complete McCallum or Lazard projection closure.

    Args:
        inputs: The rational polynomial family, including repeated factors and zero inputs.
        variables: Distinct variables in base-to-top cylindrical order, covering every input.
        budget: The shared work allowance.
        lazard: Whether to use Lazard leading/trailing coefficient projection.

    Returns:
        Projection bases with product, coprimality and exact determinant certificates.

    Raises:
        PolyError: DOMAIN for malformed inputs or variable order, LIMIT or CAPACITY on work limits.
    """
    var empty = _empty(List[Int64](), False, budget)
    return _project[False](inputs, variables, lazard, empty, budget)


def check_cad_projection(inputs: List[SparsePoly[QQ]], certificate: CADProjectionCertificate, mut budget: Budget) raises PolyError -> Bool:
    """Check projection closure without factorization, gcd search or determinant division.

    Args:
        inputs: The original rational polynomial family.
        certificate: The proposed projection closure and witnesses.
        budget: The shared work allowance.

    Returns:
        Whether every projection obligation and basis property is certified.

    Raises:
        PolyError: LIMIT or CAPACITY on work limits.
    """
    try:
        _ = _project[True](inputs, certificate.variables, certificate.lazard, certificate, budget)
        return True
    except e:
        if e.kind == PolyErrorKind.DOMAIN:
            return False
        raise e^
