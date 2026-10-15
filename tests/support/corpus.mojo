"""Read the independent corpus while preserving repeated fields and record order."""
import polynomials
from apn_mojo import Rational, Integer, div_rem_euclid
from polynomials import Budget, Monomial, MonomialOrder

@fieldwise_init
struct Case(Copyable, Movable):
    var name: String
    var lines: List[String]

    def field(self, prefix: String) raises -> String:
        for line in self.lines:
            if line.startswith(prefix):
                return String(line[byte=prefix.byte_length():])
        raise Error(String("missing corpus field ", prefix, " in ", self.name))

def read_cases(path: String = "tests/reference/corpus.txt") raises -> List[Case]:
    var text: String
    with open(path, "r") as f:
        text = f.read()
    var result = List[Case]()
    var current = Case("", List[String]())
    var active = False
    for view in text.split("\n"):
        var line = String(view)
        if line.byte_length() == 0 or line.startswith("#"):
            continue
        if line.startswith("case "):
            if active:
                raise Error("nested corpus case")
            current = Case(String(line[byte=5:]), List[String]())
            active = True
        elif line == "end":
            if not active:
                raise Error("end outside corpus case")
            result.append(current.copy())
            active = False
        elif active and (line.startswith("op ") or line.startswith("ring ") or line.startswith("var ") or line.startswith("order ") or line.startswith("in ") or line.startswith("out ")):
            current.lines.append(line)
        else:
            raise Error(String("invalid corpus line: ", line))
    if active:
        raise Error("unterminated corpus case")
    return result^

def parse_order(text: String, mut budget: Budget) raises -> MonomialOrder:
    var parts = text.split(" ")
    var kind = MonomialOrder.LEX
    if String(parts[0]) == "grlex":
        kind = MonomialOrder.GRLEX
    elif String(parts[0]) == "grevlex":
        kind = MonomialOrder.GREVLEX
    elif String(parts[0]) != "lex":
        raise Error("unknown corpus order")
    var variables = List[Int64]()
    for i in range(1, len(parts)):
        variables.append(Int64(Int(parts[i])))
    return MonomialOrder(kind, variables^, budget)

def parse_monomial(text: String, mut budget: Budget) raises -> Monomial:
    var pairs = List[Tuple[Int64, UInt32]]()
    if text == "1":
        return Monomial()
    for factor in text.split("*"):
        var parts = factor.split("^")
        var variable = String(parts[0])
        if not variable.startswith("x"):
            raise Error("expected monomial variable")
        var index = Int64(Int(variable[byte=1:]))
        var exponent = UInt32(1)
        if len(parts) == 2:
            exponent = UInt32(Int(parts[1]))
        pairs.append((index, exponent))
    # Corpus spelling is independent of storage order.
    for i in range(1, len(pairs)):
        var j = i
        while j > 0 and pairs[j][0] < pairs[j - 1][0]:
            var old = pairs[j - 1]
            pairs[j - 1] = pairs[j]
            pairs[j] = old
            j -= 1
    return Monomial(pairs^, budget)


def parse_polynomial(text: String, order: MonomialOrder, mut budget: Budget) raises -> polynomials.SparsePoly[polynomials.QQ]:
    """Parse the corpus polynomial grammar independently of production arithmetic."""
    var terms = List[polynomials.Term[polynomials.QQ]]()
    if text == "0":
        return polynomials.SparsePoly[polynomials.QQ].zero(order, polynomials.NoContext(), budget)
    var separated = text.replace(" - ", " + -")
    for term_view in separated.split(" + "):
        var term = String(term_view)
        var sign = 1
        if term.startswith("-"):
            sign = -1
            var unsigned_term = String(term[byte=1:])
            term = unsigned_term^
        var coefficient_text = String(sign)
        var monomial_text = String("1")
        var pieces = term.split("*")
        if term.startswith("x"):
            monomial_text = term
        else:
            coefficient_text = String("-" if sign < 0 else "", pieces[0])
            if len(pieces) > 1:
                monomial_text = String(pieces[1])
                for i in range(2, len(pieces)):
                    monomial_text += String("*", pieces[i])
        terms.append(polynomials.Term[polynomials.QQ](parse_monomial(monomial_text, budget), polynomials.QQ(Rational(coefficient_text))))
    return polynomials.SparsePoly[polynomials.QQ].from_terms(terms^, order, polynomials.NoContext(), budget)


def as_zz(p: polynomials.SparsePoly[polynomials.QQ], mut work: Budget) raises -> polynomials.SparsePoly[polynomials.ZZ]:
    """Convert an independently parsed integer polynomial to ZZ."""
    var terms = List[polynomials.Term[polynomials.ZZ]]()
    for t in p.terms:
        terms.append(polynomials.Term[polynomials.ZZ](t.monomial.copy(), polynomials.ZZ(t.coefficient.value.to_integer_exact())))
    return polynomials.SparsePoly[polynomials.ZZ].from_terms(terms^, p.order, polynomials.NoContext(), work)


def as_mod(p: polynomials.SparsePoly[polynomials.QQ], ctx: polynomials.Modulus, mut work: Budget) raises -> polynomials.SparsePoly[polynomials.ModInt]:
    """Convert integer corpus coefficients to residues using independent apn division."""
    var terms = List[polynomials.Term[polynomials.ModInt]]()
    for t in p.terms:
        var qr = div_rem_euclid(t.coefficient.value.to_integer_exact(), Integer(Int(ctx.value)))
        var residue = UInt64(Int(qr[1].to_string()))
        terms.append(polynomials.Term[polynomials.ModInt](t.monomial.copy(), polynomials.ModInt(residue, ctx)))
    return polynomials.SparsePoly[polynomials.ModInt].from_terms(terms^, p.order, ctx, work)
