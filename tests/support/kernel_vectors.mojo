"""Parse independent polynomial operation vectors without computing expected answers."""
from apn_mojo import Rational
from polynomials import Budget, QQ, NoContext, Monomial, MonomialOrder, Term, SparsePoly

def parse_poly(text: String, mut work: Budget) raises -> SparsePoly[QQ]:
    var terms = List[Term[QQ]]()
    if text == "0":
        return SparsePoly[QQ].zero(MonomialOrder(), NoContext(), work)
    for term_view in text.split(";"):
        var term = String(term_view)
        var parts = term.split("[")
        if len(parts) != 2 or not term.endswith("]"):
            raise Error("malformed kernel vector polynomial")
        var coeff = QQ(Rational(String(parts[0])))
        var tail = String(parts[1])
        var body = String(tail[byte=0:tail.byte_length() - 1])
        var pairs = List[Tuple[Int64, UInt32]]()
        if body.byte_length():
            for pair in body.split(","):
                var fields = pair.split(":")
                pairs.append((Int64(Int(fields[0])), UInt32(Int(fields[1]))))
        terms.append(Term[QQ](Monomial(pairs^, work), coeff^))
    return SparsePoly[QQ].from_terms(terms^, MonomialOrder(), NoContext(), work)

def read_lines(path: String) raises -> List[String]:
    var text: String
    with open(path, "r") as f:
        text = f.read()
    var lines = List[String]()
    for line in text.split("\n"):
        if line.byte_length() and not line.startswith("#"):
            lines.append(String(line))
    return lines^

def parse_bytes(text: String) raises -> List[UInt8]:
    if text == "-":
        return List[UInt8]()
    var digits = text.as_bytes()
    if len(digits) % 2:
        raise Error("odd hex length")
    var result = List[UInt8]()
    var value = UInt8(0)
    for i in range(len(digits)):
        var digit = Int(digits[i])
        if digit >= ord("0") and digit <= ord("9"):
            digit -= ord("0")
        elif digit >= ord("a") and digit <= ord("f"):
            digit -= ord("a") - 10
        else:
            raise Error("invalid hex digit")
        value = (value << 4) | UInt8(digit)
        if i % 2:
            result.append(value)
            value = 0
    return result^
