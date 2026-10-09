# Матричный калькулятор, реализующий 15 функций над матрицами:
# 1. Транспонирование матрицы (Stepskelet)
# 2. Нахождение следа матрицы (Bebrick322 +)
# 3. Умножение матрицы на число (GroM +)
# 4. Сложение матриц (Bebrick322 +)
# 5. Вычитание матриц (Stepskelet)
# 6. Генератор матриц(единичная, нулевая, случайная, случайная-диагональная, случайная-симметричная) (GroM +)
# 7. Экспорт в LaTeX (GroM +)
# 8. Перемножение матриц (GroM +)
# 9. Вычисление определителя (GroM +)
# 10. Решение СЛАУ методом Гаусса (Stepskelet)
# 11. Нахождение степеней матрицы (Stepskelet)
# 12. Решение СЛАУ методом Крамера (Bebrick322 +)
# 13. Нахождение обратной матрицы (Stepskelet)
# 14. Решение СЛАУ матричным методом (Bebrick322 +)
# 15. LUP разложение матрицы (Bebrick322 +)

import random


def matrix_x_number(matrix: list[list], value: int) -> list[list]:
    return [[cell * value for cell in row] for row in matrix]


def generate_identity_matrix(n: int) -> list[list[int]]:
    return [[1 if i == j else 0 for j in range(n)] for i in range(n)]


def generate_zero_matrix(rows: int, cols: int) -> list[list[int]]:
    return [[0 for _ in range(cols)] for _ in range(rows)]


def generate_random_matrix(rows: int, cols: int, min_val: int = -100, max_val: int = 100) -> list[list[int]]:
    return [[random.randint(min_val, max_val) for _ in range(cols)] for _ in range(rows)]


def generate_random_diagonal_matrix(n: int, min_val: int = -100, max_val: int = 100) -> list[list[int]]:
    matrix = [[0] * n for _ in range(n)]
    for i in range(n):
        matrix[i][i] = random.randint(min_val, max_val)
    return matrix


def generate_random_symmetric_matrix(n: int, min_val: int = -100, max_val: int = 100) -> list[list[int]]:
    matrix = [[0] * n for _ in range(n)]
    for i in range(n):
        for j in range(i, n):
            val = random.randint(min_val, max_val)
            matrix[i][j] = val
            matrix[j][i] = val
    return matrix


def matrix_to_latex(matrix: list[list], env: str = "bmatrix") -> str:
    rows = []
    for row in matrix:
        formatted_row = []
        for val in row:
            if isinstance(val, float):
                val_str = f"{val:.4f}".rstrip('0').rstrip('.')
            else:
                val_str = str(val)
            formatted_row.append(val_str)
        rows.append(" & ".join(formatted_row))

    body = " \\\\\n  ".join(rows)
    return f"\\begin{{{env}}}\n  {body}\n\\end{{{env}}}"


def matrix_x_matrix(matrix_a: list[list], matrix_b: list[list]) -> list[list]:
    rows_a = len(matrix_a)
    cols_a = len(matrix_a[0])
    cols_b = len(matrix_b[0])

    result = [[0] * cols_b for _ in range(rows_a)]

    for i in range(rows_a):
        for j in range(cols_b):
            result[i][j] = sum(matrix_a[i][k] * matrix_b[k][j]
                               for k in range(cols_a))
    return result


def calculate_determinant(matrix: list[list]) -> float | int:
    n = len(matrix)

    if n == 1:
        return matrix[0][0]
    if n == 2:
        return matrix[0][0] * matrix[1][1] - matrix[0][1] * matrix[1][0]

    det = 0
    for col in range(n):
        minor = [row[:col] + row[col+1:] for row in matrix[1:]]
        sign = 1 if col % 2 == 0 else -1
        det += sign * matrix[0][col] * calculate_determinant(minor)

    return det


EPS = 1e-12


def trace(matrix: list[list]) -> float | int:
    """Сумма элементов на главной диагонали. Матрица должна быть квадратной"""
    n = len(matrix)
    if n == 0 or any(len(row) != n for row in matrix):
        raise ValueError("Матрица должна быть квадратной и непустой")
    return sum(matrix[i][i] for i in range(n))


def add_matrices(matrix_a: list[list], matrix_b: list[list]) -> list[list]:
    """Поэлементное сложение двух матриц одинакового размера"""
    if len(matrix_a) != len(matrix_b):
        raise ValueError("Число строк не совпадает")
    if not matrix_a:
        return []
    cols = len(matrix_a[0])
    if any(len(row) != cols for row in matrix_a) or any(len(row) != cols for row in matrix_b):
        raise ValueError("Размеры матриц не совпадают")
    return [[matrix_a[i][j] + matrix_b[i][j] for j in range(cols)]
            for i in range(len(matrix_a))]


def lup_decomposition(matrix: list[list]) -> tuple[list[list], list[list], list[list]]:
    n = len(matrix)
    if n == 0 or any(len(row) != n for row in matrix):
        raise ValueError("Матрица должна быть квадратной и непустой")

    U = [[float(x) for x in row] for row in matrix]
    L = [[0.0] * n for _ in range(n)]
    P = [[1.0 if i == j else 0.0 for j in range(n)] for i in range(n)]

    for i in range(n):
        L[i][i] = 1.0

    for k in range(n):
        pivot = max(range(k, n), key=lambda i: abs(U[i][k]))
        if abs(U[pivot][k]) < EPS:
            raise ValueError("Матрица вырождена: нулевой ведущий элемент")

        if pivot != k:
            U[k], U[pivot] = U[pivot], U[k]
            P[k], P[pivot] = P[pivot], P[k]
            L[k][:k], L[pivot][:k] = L[pivot][:k], L[k][:k]

        for i in range(k + 1, n):
            factor = U[i][k] / U[k][k]
            L[i][k] = factor
            for j in range(k, n):
                U[i][j] -= factor * U[k][j]

    return L, U, P


def _det_permutation(P: list[list]) -> float:
    """Знак матрицы перестановок"""
    n = len(P)
    perm = []
    for i in range(n):
        for j in range(n):
            if abs(P[i][j] - 1.0) < EPS:
                perm.append(j)
                break
        else:
            raise ValueError("Некорректная матрица перестановок")
    inversions = 0
    for i in range(n):
        for j in range(i + 1, n):
            if perm[i] > perm[j]:
                inversions += 1
    return -1.0 if inversions % 2 else 1.0


def _solve_lup(L: list[list], U: list[list], P: list[list], b: list) -> list:
    """Решение СЛАУ по готовому LUP-разложению"""
    n = len(L)
    if len(b) != n:
        raise ValueError("Несовпадение размеров")
    pb = [sum(P[i][j] * b[j] for j in range(n)) for i in range(n)]

    y = [0.0] * n
    for i in range(n):
        s = pb[i]
        for j in range(i):
            s -= L[i][j] * y[j]
        y[i] = s / L[i][i]

    x = [0.0] * n
    for i in range(n - 1, -1, -1):
        s = y[i]
        for j in range(i + 1, n):
            s -= U[i][j] * x[j]
        x[i] = s / U[i][i]
    return x


def solve_cramer(matrix_a: list[list], vector_b: list) -> list[float]:
    """Решение СЛАУ A*x = b методом Крамера"""
    n = len(matrix_a)
    if n == 0 or any(len(row) != n for row in matrix_a):
        raise ValueError("Матрица должна быть квадратной и непустой")
    if len(vector_b) != n:
        raise ValueError("Вектор b должен иметь длину n")

    L, U, P = lup_decomposition(matrix_a)
    det_u = 1.0
    for i in range(n):
        det_u *= U[i][i]
    det_a = det_u * _det_permutation(P)
    if abs(det_a) < EPS:
        raise ValueError("Определитель равен нулю — метод Крамера неприменим")

    x = []
    for i in range(n):
        a_i = [row[:] for row in matrix_a]
        for r in range(n):
            a_i[r][i] = vector_b[r]
        Li, Ui, Pi = lup_decomposition(a_i)
        det_u_i = 1.0
        for k in range(n):
            det_u_i *= Ui[k][k]
        det_i = det_u_i * _det_permutation(Pi)
        x.append(det_i / det_a)
    return x


def _inverse_matrix(matrix_a: list[list]) -> list[list]:
    n = len(matrix_a)
    if n == 0 or any(len(row) != n for row in matrix_a):
        raise ValueError("Матрица должна быть квадратной и непустой")

    L, U, P = lup_decomposition(matrix_a)
    inv = [[0.0] * n for _ in range(n)]

    for j in range(n):
        e = [0.0] * n
        e[j] = 1.0
        col = _solve_lup(L, U, P, e)
        for i in range(n):
            inv[i][j] = col[i]
    return inv


def solve_matrix_method(matrix_a: list[list], vector_b: list) -> list[float]:
    """Решение СЛАУ A*x = b матричным методом: x = A^{-1}*b"""
    n = len(matrix_a)
    if n == 0 or any(len(row) != n for row in matrix_a):
        raise ValueError("Матрица должна быть квадратной и непустой")
    if len(vector_b) != n:
        raise ValueError("Вектор b должен иметь длину n")

    inv = _inverse_matrix(matrix_a)
    return [sum(inv[i][j] * vector_b[j] for j in range(n)) for i in range(n)]
