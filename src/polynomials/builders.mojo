"""Named polynomial rings with a restricted iterative expression grammar."""
from apn_mojo import ExactComplex
from .coefficients import CoefficientRing
from .exact_complex import QQi
from .sparse import SparsePoly, Term
from .monomial import MonomialOrder
from .errors import Budget, PolyError, PolyErrorKind, _allocation, _charge_size, _sum
from ._coefficient_ops import _exact
from ._text_ops import _literal, _scalar_text, _unsigned


def _letter(byte: UInt8) -> Bool:
    return (byte >= 65 and byte <= 90) or (byte >= 97 and byte <= 122) or byte == 95


def _digit(byte: UInt8) -> Bool:
    return byte >= 48 and byte <= 57


def _space(byte: UInt8) -> Bool:
    return byte == 32 or byte == 9 or byte == 10 or byte == 13


def _precedence(op: UInt8) -> Int:
    return 3 if op >= 125 else 2 if op == 42 or op == 47 else 1


def _apply[C: CoefficientRing](mut values: List[SparsePoly[C]], op: UInt8, mut budget: Budget) raises PolyError:
    budget.charge(1)
    if len(values) < (1 if op >= 125 else 2):
        raise PolyError(PolyErrorKind.ENCODING, "missing polynomial operand")
    var right = values.pop()
    if op >= 125:
        values.append(right.neg(budget) if op == 126 else right^)
        return
    var left = values.pop()
    if op == 43:
        values.append(left.add(right, budget))
    elif op == 45:
        values.append(left.sub(right, budget))
    elif op == 42:
        values.append(left.mul(right, budget))
    else:
        if len(right.terms) > 1 or (len(right.terms) == 1 and len(right.terms[0].monomial.pairs)):
            raise PolyError(PolyErrorKind.DOMAIN, "polynomial grammar divides only by scalar constants")
        if right.is_zero():
            raise PolyError(PolyErrorKind.DOMAIN, "zero scalar divisor")
        _allocation[Term[C]](len(left.terms), budget)
        var terms = List[Term[C]](capacity=len(left.terms))
        for term in left.terms:
            _charge_size(len(term.monomial.pairs), budget)
            terms.append(Term[C](term.monomial.copy(), _exact(term.coefficient, right.terms[0].coefficient, budget)))
        values.append(SparsePoly[C].from_terms(terms^, left.order, left.ring, budget))


@fieldwise_init
struct PolynomialRing[C: CoefficientRing](Copyable, Movable):
    """Bind ASCII variable names to stable consecutive IDs and an exact ring.

    Parameters:
        C: A built-in exact coefficient domain for parsing and formatting.
    """
    var _names: List[String]
    var _ring: Self.C.Context
    var _order: MonomialOrder

    @staticmethod
    def create(names: List[String], ring: Self.C.Context, mut budget: Budget, *, kind: UInt8 = MonomialOrder.LEX) raises PolyError -> Self:
        """Create a named ring with IDs equal to name positions.

        Args:
            names: Distinct ASCII identifiers; i is reserved over QQi.
            ring: The exact coefficient context.
            budget: The shared work allowance.
            kind: The admissible monomial order kind.

        Returns:
            The prepared named-variable builder.

        Raises:
            PolyError: DOMAIN for invalid names or orders, LIMIT or CAPACITY on work limits.
        """
        _allocation[String](len(names), budget)
        _allocation[Int64](len(names), budget)
        var copied = List[String](capacity=len(names))
        var variables = List[Int64](capacity=len(names))
        for i in range(len(names)):
            _charge_size(names[i].byte_length(), budget)
            var bytes = names[i].as_bytes()
            if not len(bytes) or not _letter(bytes[0]):
                raise PolyError(PolyErrorKind.DOMAIN, "invalid polynomial variable name")
            for byte in bytes:
                if not _letter(byte) and not _digit(byte):
                    raise PolyError(PolyErrorKind.DOMAIN, "invalid polynomial variable name")
            comptime if Self.C == QQi:
                if names[i] == "i":
                    raise PolyError(PolyErrorKind.DOMAIN, "i is reserved for the imaginary unit")
            for j in range(i):
                _charge_size(_sum(names[i].byte_length(), names[j].byte_length()), budget)
                if names[i] == names[j]:
                    raise PolyError(PolyErrorKind.DOMAIN, "duplicate polynomial variable name")
            copied.append(names[i].copy())
            if kind != MonomialOrder.KERNEL:
                variables.append(Int64(i))
        return Self(copied^, ring.copy(), MonomialOrder(kind, variables^, budget))

    def variable(self, name: String, mut budget: Budget) raises PolyError -> SparsePoly[Self.C]:
        """Construct a variable by its configured name.

        Args:
            name: A name supplied when creating the ring.
            budget: The shared work allowance.

        Returns:
            The corresponding polynomial generator.

        Raises:
            PolyError: DOMAIN for an unknown name, LIMIT or CAPACITY on work limits.
        """
        for i in range(len(self._names)):
            _charge_size(_sum(name.byte_length(), self._names[i].byte_length()), budget)
            if name == self._names[i]:
                return SparsePoly[Self.C].variable(Int64(i), self._order, self._ring, budget)
        raise PolyError(PolyErrorKind.DOMAIN, "unknown polynomial variable")

    def parse(self, text: String, mut budget: Budget) raises PolyError -> SparsePoly[Self.C]:
        """Parse a restricted polynomial expression without executing code.

        Args:
            text: Integers, names, parentheses, unary signs, +, -, *, scalar / and literal nonnegative powers.
            budget: The shared work allowance.

        Returns:
            The canonical expanded polynomial; exponent chains require parentheses.

        Raises:
            PolyError: ENCODING for malformed grammar, DOMAIN for unsupported division or names, LIMIT or CAPACITY on work limits.
        """
        _charge_size(text.byte_length(), budget)
        var bytes = text.as_bytes()
        var values = List[SparsePoly[Self.C]]()
        var operators = List[UInt8]()
        var position = 0
        var operand = True
        var powered = False
        while position < len(bytes):
            budget.charge(1)
            var byte = bytes[position]
            if _space(byte):
                position += 1
                continue
            if operand:
                if byte == 43 or byte == 45 or byte == 40:
                    _allocation[UInt8](1, budget)
                    operators.append(UInt8(125 if byte == 43 else 126 if byte == 45 else 40))
                    position += 1
                    continue
                var start = position
                if _digit(byte):
                    while position < len(bytes) and _digit(bytes[position]):
                        budget.charge(1)
                        position += 1
                    var scalar = _literal[Self.C](String(text[byte=start:position]), self._ring, budget)
                    _allocation[SparsePoly[Self.C]](1, budget)
                    values.append(SparsePoly[Self.C].constant(scalar, self._order, budget))
                elif _letter(byte):
                    while position < len(bytes) and (_letter(bytes[position]) or _digit(bytes[position])):
                        budget.charge(1)
                        position += 1
                    var name = String(text[byte=start:position])
                    _allocation[SparsePoly[Self.C]](1, budget)
                    comptime if Self.C == QQi:
                        if name == "i":
                            values.append(SparsePoly[Self.C].constant(rebind_var[Self.C](QQi(ExactComplex(0, 1))), self._order, budget))
                        else:
                            values.append(self.variable(name, budget))
                    else:
                        values.append(self.variable(name, budget))
                else:
                    raise PolyError(PolyErrorKind.ENCODING, "expected polynomial operand")
                operand = False
                powered = False
            elif byte == 94:
                if powered:
                    raise PolyError(PolyErrorKind.ENCODING, "parenthesize chained polynomial powers")
                position += 1
                while position < len(bytes) and _space(bytes[position]):
                    budget.charge(1)
                    position += 1
                var start = position
                while position < len(bytes) and _digit(bytes[position]):
                    budget.charge(1)
                    position += 1
                var exponent = _unsigned(String(text[byte=start:position]), UInt64(UInt32.MAX), budget)
                var base = values.pop()
                values.append(base.pow(exponent, budget))
                powered = True
            elif byte == 41:
                while len(operators) and operators[len(operators) - 1] != 40:
                    _apply(values, operators.pop(), budget)
                if not len(operators):
                    raise PolyError(PolyErrorKind.ENCODING, "unmatched closing parenthesis")
                _ = operators.pop()
                position += 1
                powered = False
            elif byte == 43 or byte == 45 or byte == 42 or byte == 47:
                while len(operators) and operators[len(operators) - 1] != 40 and _precedence(operators[len(operators) - 1]) >= _precedence(byte):
                    _apply(values, operators.pop(), budget)
                _allocation[UInt8](1, budget)
                operators.append(byte)
                operand = True
                position += 1
            else:
                raise PolyError(PolyErrorKind.ENCODING, "expected polynomial operator")
        if operand:
            raise PolyError(PolyErrorKind.ENCODING, "incomplete polynomial expression")
        while len(operators):
            var op = operators.pop()
            if op == 40:
                raise PolyError(PolyErrorKind.ENCODING, "unmatched opening parenthesis")
            _apply(values, op, budget)
        if len(values) != 1:
            raise PolyError(PolyErrorKind.ENCODING, "invalid polynomial expression")
        return values.pop()

    def format(self, poly: SparsePoly[Self.C], mut budget: Budget) raises PolyError -> String:
        """Format a polynomial as a parsable named-variable expression.

        Args:
            poly: A polynomial in the configured coefficient ring and order.
            budget: The shared work allowance.

        Returns:
            Exact coefficient text in decreasing term order.

        Raises:
            PolyError: DOMAIN for incompatible variables or contexts, LIMIT or CAPACITY on work limits.
        """
        poly.check_invariants(budget)
        _charge_size(_sum(len(poly.order.variables), len(self._order.variables)), budget)
        if poly.ring != self._ring or poly.order != self._order:
            raise PolyError(PolyErrorKind.DOMAIN, "polynomial and named ring disagree")
        if poly.is_zero():
            return "0"
        var pieces = List[String]()
        for i in range(len(poly.terms) - 1, -1, -1):
            ref term = poly.terms[i]
            var scalar = _scalar_text(term.coefficient, budget)
            comptime if Self.C == QQi:
                var components = scalar.split(",")
                var formatted = String(components[0], " + (", components[1], ")*i")
                scalar = formatted^
            _allocation[String](1, budget)
            var piece = String("(", scalar, ")")
            for pair in term.monomial.pairs:
                budget.charge(1)
                if pair[0] < 0 or pair[0] >= Int64(len(self._names)):
                    raise PolyError(PolyErrorKind.DOMAIN, "unconfigured polynomial variable ID")
                _charge_size(_sum(piece.byte_length(), self._names[Int(pair[0])].byte_length()), budget)
                piece += String("*", self._names[Int(pair[0])])
                if pair[1] != 1:
                    piece += String("^", pair[1])
            pieces.append(piece^)
        var length = 0
        for piece in pieces:
            length = _sum(length, _sum(piece.byte_length(), 3))
        _charge_size(length, budget)
        return String(" + ").join(pieces)
