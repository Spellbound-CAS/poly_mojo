"""Every kernel codec vector, malformed payloads, and generated round trips."""
from apn_mojo import Integer, Rational
from polynomials import Budget, PolyErrorKind, QQ, NoContext, MonomialOrder, SparsePoly, to_kernel_bytes, from_kernel_bytes, MAX_KERNEL_BYTES
from support.check import Checker
from support.kernel_vectors import read_lines, parse_poly, parse_bytes
from support.gen import polynomial
from support.rng import SplitMix64

def main() raises:
    var c = Checker("codec")
    var work = Budget()
    var encoded_count = 0
    var bad_count = 0
    for line in read_lines("tests/reference/kernel/poly_ops.txt"):
        var parts = line.split(" ")
        var op = String(parts[0])
        if op == "enc":
            var p = parse_poly(String(parts[1]), work)
            var bytes = parse_bytes(String(parts[2]))
            c.check(to_kernel_bytes(p, work) == bytes, "independent encoding vector")
            var decoded = from_kernel_bytes(Span(bytes), work)
            decoded.check_invariants(work)
            c.check(decoded == p, "independent decoding vector")
            encoded_count += 1
        elif op == "bad":
            var bytes = parse_bytes(String(parts[1]))
            try:
                _ = from_kernel_bytes(Span(bytes), work)
                c.check(False, "malformed vector accepted")
            except e:
                c.check(e.kind == PolyErrorKind.ENCODING, "malformed vector ENCODING")
            bad_count += 1
    c.check(encoded_count == 127 and bad_count == 11, "all codec vectors executed")
    var table_polys: List[String] = ["2[]", "1[0:1]", "-1/2[];1[-1:2]"]
    var table_index = 0
    for line in read_lines("tests/reference/kernel/poly_vectors.txt"):
        var hex_bytes = String("")
        var started = False
        for part in line.split(" "):
            if part.byte_length() == 2:
                var text = String(part)
                if text == "01" or text == "02":
                    started = True
            if started and part.byte_length():
                hex_bytes += String(part)
        var bytes = parse_bytes(hex_bytes)
        c.check(to_kernel_bytes(parse_poly(table_polys[table_index], work), work) == bytes, "kernel table")
        table_index += 1
    c.check(table_index == 3, "all three table rows")
    var rng = SplitMix64(123)
    for _ in range(100):
        var order = MonomialOrder(MonomialOrder.GREVLEX, [2, 1, 0], work)
        var p = polynomial[QQ](NoContext(), order, rng, work, terms=8)
        if p.is_zero():
            continue
        var bytes = to_kernel_bytes(p, work)
        c.check(from_kernel_bytes(Span(bytes), work) == p.reorder(MonomialOrder(), work), "generated reordered roundtrip")
    var x = parse_poly("1[0:1]", work)
    var valid = to_kernel_bytes(x, work)
    for length in range(len(valid)):
        try:
            _ = from_kernel_bytes(Span(valid)[0:length], work)
            c.check(False, "truncated payload accepted")
        except e:
            c.check(e.kind == PolyErrorKind.ENCODING, "truncation ENCODING")
    var trailing = valid.copy()
    trailing.append(0)
    try:
        _ = from_kernel_bytes(Span(trailing), work)
        c.check(False, "trailing byte accepted")
    except e:
        c.check(e.kind == PolyErrorKind.ENCODING, "trailing byte ENCODING")
    for n in range(2):
        var constant = SparsePoly[QQ].constant(QQ(Rational(n)), MonomialOrder(), work)
        try:
            _ = to_kernel_bytes(constant, work)
            c.check(False, "zero/one literal accepted")
        except e:
            c.check(e.kind == PolyErrorKind.ENCODING, "zero/one ENCODING")
    var large = QQ(Rational(-(Integer(1) << 8193), Integer(3)))
    var huge = SparsePoly[QQ].constant(large, MonomialOrder(), work)
    var huge_bytes = to_kernel_bytes(huge, work)
    c.check(from_kernel_bytes(Span(huge_bytes), work) == huge, "large signed rational binary splitting")
    var oversized = List[UInt8](length=MAX_KERNEL_BYTES + 1, fill=0)
    try:
        _ = from_kernel_bytes(Span(oversized), work)
        c.check(False, "oversized payload accepted")
    except e:
        c.check(e.kind == PolyErrorKind.CAPACITY, "kernel ceiling CAPACITY")
    for steps in range(70):
        var allowance = Budget(UInt64(steps))
        try:
            _ = from_kernel_bytes(Span(valid), allowance)
        except e:
            c.check(e.kind == PolyErrorKind.LIMIT, "decoder exhaustion")
        var encode_allowance = Budget(UInt64(steps))
        try:
            _ = to_kernel_bytes(x, encode_allowance)
        except e:
            c.check(e.kind == PolyErrorKind.LIMIT, "encoder exhaustion")
    print("codec vectors:", encoded_count, "valid,", bad_count, "invalid")
    c.finish()
