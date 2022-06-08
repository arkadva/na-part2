import math

# assuming n by n+1 matrix
def solveMatrix(matrix, size):
    for i in range(size):
        # preprocess the matrix
        currentPivot = abs(matrix[i][i])
        maxRow = i
        row = i + 1
        while row < size:
            if abs(matrix[row][i]) > currentPivot:
                currentPivot = abs(matrix[row][i])
                maxRow = row
            row += 1
        matrix = swapRows(matrix, i, maxRow, size)

        matrix = setPivotToOne(matrix, i, i, size)
        row = i + 1
        while row < size:
            matrix = nullify(matrix, row, i, size, matrix[i][i])
            row += 1
    for i in range(1, size):
        row = i - 1
        while row >= 0:
            matrix = nullify(matrix, row, i, size, matrix[i][i])
            row -= 1

    return matrix

def identityMatrix(size):
    matrix = []
    for i in range(size):
        matrix.append([])
        for j in range(size):
            matrix[i].append(0)

    for i in range(size):
        matrix[i][i] = 1
    return matrix

def swapRows(matrix, row1, row2, size):
    identity = identityMatrix(size)
    identity[row1], identity[row2] = identity[row2], identity[row1]
    return multiplyMatrices(identity, matrix, size)

def editMatrix(matrix, x, y, val):
    matrix[x][y] = val
    return matrix

# assuming nxn+1 matrix
def multiplyMatrices(matrix1, matrix2, size):
    mat = []
    for i in range(size+1):
        mat.append([])

    for i in range(size):
        for j in range(size + 1):
            sum = 0
            for k in range(size):
                sum += matrix1[i][k] * matrix2[k][j]
            mat[i].append(sum)

    printMultiplication(matrix1, matrix2, mat, size, size)
    return mat

def nullify(matrix, x, y, size, pivot):
    identity = identityMatrix(size)
    return multiplyMatrices(editMatrix(identity, x, y, -1*matrix[x][y]/pivot), matrix, size)

def setPivotToOne(matrix, x, y, size):
    identity = identityMatrix(size)
    return multiplyMatrices(editMatrix(identity, x, y, 1/matrix[x][y]), matrix, size)

# assuming nxn+1
def printMultiplication(A, B, res, size, size2):
    for i in range(size):
        print('[', end =' ')
        for j in range(size2):
            print('{:6.4f}'.format(A[i][j]), end=' ')
        print('][', end =' ')
        for j in range(size2 + 1):
            print('{:6.4f}'.format(B[i][j]), end=' ')
        print(']\t[', end=' ')
        for j in range(size2 + 1):
            print('{:6.4f}'.format(res[i][j]), end=' ')
        print(']', end='')
        print()
    print()

def printMatrix(A, size, size2):
    for i in range(size):
        for j in range(size2):
            print(A[i][j], end=' ')
        print()

class Point:
    def __init__(self, x, y):
        self.x = x
        self.y = y

def cubic_spline_interpolation(table, target):
    gamma = []
    mu = []
    d = []
    h = []

    for i in range(0, len(table) - 1):
        print(f'h_{i} = {table[i + 1].x - table[i].x}')
        h.append(table[i + 1].x - table[i].x)

    print()
    for i in range(1, len(table) - 1):
        print(f'g_{i} = {h[i]/(h[i] + h[i - 1])}')
        gamma.append(h[i]/(h[i] + h[i - 1]))
        print(f'mu_{i} = {1 - h[i]/(h[i] + h[i - 1])}')
        mu.append(1 - h[i]/(h[i] + h[i - 1]))
        print(f'd_{i} = {(6/(h[i] + h[i - 1])*((table[i + 1].y - table[i].y)/h[i] - (table[i].y - table[i - 1].y)/h[i - 1]))}')
        d.append((6/(h[i] + h[i - 1])*((table[i + 1].y - table[i].y)/h[i] - (table[i].y - table[i - 1].y)/h[i - 1])))
        print()

    # build matrix
    mat = identityMatrix(len(d))
    for i in range(len(d)):
        mat[i][i] = 2
        if i != 0:
            mat[i][i - 1] = mu[i]
        if i != len(d) - 1:
            mat[i][i + 1] = gamma[i]
        mat[i].append(d[i])

    print("Target matrix:")
    printMatrix(mat, len(d), len(d) + 1)
    # extract result
    m = [0]
    res = solveMatrix(mat, len(d))
    for x in range(len(res) - 1):
        m.append(res[x][-1])
    m.append(0)

    print("Matrix result:", m)
    for y in range(len(table) - 1):
        if target > table[y].x:
            if target < table[y + 1].x:
                res = ((table[y + 1].x - target)**3 * m[y] + (target - table[y].x) ** 3 * m[y + 1])/(6*h[y]) \
                        + ((table[y + 1].x - target) * table[y].y + (target - table[y].x) * table[y + 1].y)/(h[y]) \
                        - h[y]*((table[y + 1].x - target) * m[y] + (target - table[y].x) * m[y + 1])/6
                print(f'Final result: (({table[y + 1].x} - {target})**3 * {m[y]} + ({target} - {table[y].x}) ** 3 * {m[y + 1]})/{(6*h[y])} + (({table[y + 1].x} - {target}) * {table[y].y} + ({target} - {table[y].x}) * {table[y + 1].y})/({h[y]}) - {h[y]}*(({table[y + 1].x} - {target}) * {m[y]} + ({target} - {table[y].x}) * {m[y + 1]})/6 = {res}')
                return res
    print("Target out of bounds")

def lagrange(points, xr):
    n = len(points)
    l_mult = 1  # calculating multiplier in L calculation
    l_sum = 0  # calculating sum in P calculation
    print("P = 0")
    for i in range(n):
        print("L = 1")
        l_mult = 1
        for j in range(n):
            if i != j:
                print(f'L = {l_mult} * ({xr} - {points[j].x}) / ({points[i].x} - {points[j].x}) = {l_mult * ((xr - points[j].x) / (points[i].x - points[j].x))}')
                l_mult *= ((xr - points[j].x) / (points[i].x - points[j].x))
        print(f'P = {l_sum} + {l_mult} * {points[i].y} = {l_sum + l_mult * points[i].y}\n')
        l_sum += l_mult * points[i].y

    return l_sum


point_table = [Point(0, 0), Point(2, 0), Point(2.25, 0.112463), Point(2.3, 0.167996), Point(2.7, 0.222709)]
target_point = 2.4

print("Lagrange start:")
lag_res = lagrange(point_table[1:], target_point)
print("Cubic spline interpolation start:")
spline_res = cubic_spline_interpolation(point_table, target_point)
print("Lagrange interpolation:", lag_res)
print("Cubic spline interpolation:", spline_res)