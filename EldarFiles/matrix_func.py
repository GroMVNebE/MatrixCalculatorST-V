# Функции Bebrick322: пункты 2, 4, 12, 14, 15

EPS = 1e-12


# ---------- Пункт 2: след матрицы ----------
def trace(matrix: list[list]) -> float | int:
    """Сумма элементов на главной диагонали. Матрица должна быть квадратной."""
    n = len(matrix)
    if n == 0 or any(len(row) != n for row in matrix):
        raise ValueError("Матрица должна быть квадратной и непустой")
    return sum(matrix[i][i] for i in range(n))


# ---------- Пункт 4: сложение матриц ----------
def add_matrices(matrix_a: list[list], matrix_b: list[list]) -> list[list]:
    """Поэлементное сложение двух матриц одинакового размера."""
    if len(matrix_a) != len(matrix_b):
        raise ValueError("Число строк не совпадает")
    if not matrix_a:
        return []
    cols = len(matrix_a[0])
    if any(len(row) != cols for row in matrix_a) or any(len(row) != cols for row in matrix_b):
        raise ValueError("Размеры матриц не совпадают")
    return [[matrix_a[i][j] + matrix_b[i][j] for j in range(cols)]
            for i in range(len(matrix_a))]


# ---------- Пункт 15: LUP-разложение ----------
def lup_decomposition(matrix: list[list]) -> tuple[list[list], list[list], list[list]]:
    """P·A = L·U. L — нижнетреугольная с 1 на диагонали,
    U — верхнетреугольная, P — матрица перестановок."""
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
    """Знак матрицы перестановок."""
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
    """Решение СЛАУ по готовому LUP-разложению."""
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


# ---------- Пункт 12: метод Крамера ----------
def solve_cramer(matrix_a: list[list], vector_b: list) -> list[float]:
    """Решение СЛАУ A·x = b методом Крамера."""
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


# ---------- Пункт 13 (вспомогательная для п.14): обратная матрица ----------
def inverse_matrix(matrix_a: list[list]) -> list[list]:
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


# ---------- Пункт 14: матричный метод ----------
def solve_matrix_method(matrix_a: list[list], vector_b: list) -> list[float]:
    """Решение СЛАУ A·x = b матричным методом: x = A^{-1}·b."""
    n = len(matrix_a)
    if n == 0 or any(len(row) != n for row in matrix_a):
        raise ValueError("Матрица должна быть квадратной и непустой")
    if len(vector_b) != n:
        raise ValueError("Вектор b должен иметь длину n")

    inv = inverse_matrix(matrix_a)
    return [sum(inv[i][j] * vector_b[j] for j in range(n)) for i in range(n)]