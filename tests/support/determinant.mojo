"""Independent fraction-free determinants and Sylvester minors for P2 tests."""
from apn_mojo import Integer, div_rem_trunc
from polynomials import ZZ, DensePoly


def determinant(var matrix: List[List[Integer]]) raises -> Integer:
    var size = len(matrix)
    if not size:
        return Integer(1)
    var previous = Integer(1)
    var sign = 1
    for k in range(size - 1):
        var pivot = k
        while pivot < size and matrix[pivot][k] == Integer(0):
            pivot += 1
        if pivot == size:
            return Integer(0)
        if pivot != k:
            var temporary = List[Integer]()
            swap(temporary, matrix[k])
            swap(temporary, matrix[pivot])
            swap(temporary, matrix[k])
            sign = -sign
        for i in range(k + 1, size):
            for j in range(k + 1, size):
                var numerator = matrix[i][j] * matrix[k][k] - matrix[i][k] * matrix[k][j]
                var qr = div_rem_trunc(numerator, previous)
                if qr[1] != Integer(0):
                    raise Error("nonexact independent Bareiss division")
                matrix[i][j] = qr[0].copy()
        previous = matrix[k][k].copy()
        for i in range(k + 1, size):
            matrix[i][k] = Integer(0)
    return Integer(sign) * matrix[size - 1][size - 1]


def subresultant_coefficients(f: DensePoly[ZZ], g: DensePoly[ZZ], j: Int) raises -> List[ZZ]:
    # The first r-1 columns and one of the last j+1 columns of the
    # rectangular Sylvester matrix define each coefficient of S_j.
    var m = f.degree()
    var n = g.degree()
    var rows = m + n - 2 * j
    var columns = m + n - j
    var rectangle = List[List[Integer]]()
    for r in range(rows):
        var row = List[Integer](length=columns, fill=Integer(0))
        if r < n - j:
            for i in range(m + 1):
                row[r + i] = f.coefficients[m - i].value.copy()
        else:
            for i in range(n + 1):
                row[r - (n - j) + i] = g.coefficients[n - i].value.copy()
        rectangle.append(row^)
    var result = List[ZZ]()
    for degree in range(j + 1):
        var square = List[List[Integer]]()
        for r in range(rows):
            var row = List[Integer]()
            for column in range(rows - 1):
                row.append(rectangle[r][column].copy())
            row.append(rectangle[r][columns - 1 - degree].copy())
            square.append(row^)
        result.append(ZZ(determinant(square^)))
    return result^
