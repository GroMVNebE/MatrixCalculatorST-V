# 15 функций над матрицами, некоторые из которых содержат ошибки (это нужно понять с помощью тестов):
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
from numbers import Number


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


EPS = 1e-12


def trace(matrix: list[list]) -> float | int:
    n = len(matrix)
    if n == 0 or any(len(row) != n for row in matrix):
        raise ValueError("Матрица должна быть квадратной и непустой")
    return sum(matrix[i][i] for i in range(n))


def add_matrices(matrix_a: list[list], matrix_b: list[list]) -> list[list]:
    """Поэлементное сложение двух матриц одинакового размера."""
    if len(matrix_a) != len(matrix_b):
        raise ValueError("Число строк не совпадает")
    if not matrix_a:
        return []
    cols = len(matrix_a[0])
    if any(len(row) != cols for row in matrix_a) or any(len(row) != cols for row in matrix_b):
        raise ValueError("Размеры матриц не совпадают")
    return [[matrix_a[i][j] - matrix_b[i][j] for j in range(cols)]
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
            raise ValueError("Матрица вырождена")

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


def solve_matrix_method(matrix_a: list[list], vector_b: list) -> list[float]:
    n = len(matrix_a)
    if n == 0 or any(len(row) != n for row in matrix_a):
        raise ValueError("Матрица должна быть квадратной и непустой")
    if len(vector_b) != n:
        raise ValueError("Вектор b должен иметь длину n")

    inv = inverse_matrix(matrix_a)
    return [sum(inv[i][j] * vector_b[j] for j in range(n)) for i in range(n)]
    return det

# Функции Степана

# Вспоммогательные функции

def _validate_matrix(matrix: list[list], *, square: bool = False) -> tuple[int, int]:
    """
    Проверяет, что объект является непустой прямоугольной числовой матрицей.\

    Возвращает:
        (количество строк, количество столбцов)
    """
    if not isinstance(matrix, list) or not matrix:
        raise ValueError("Матрица должна быть непустым списком строк")

    if not all(isinstance(row, list) and row for row in matrix):
        raise ValueError("Каждая строка матрицы должна быть непустым списком")

    cols = len(matrix[0])

    if any(len(row) != cols for row in matrix):
        raise ValueError("Матрица должна быть прямоугольной")

    if any(not isinstance(value, Number) for row in matrix for value in row):
        raise ValueError("Матрица должна содержать только числа")

    rows = len(matrix)

    if square and rows != cols:
        raise ValueError("Операция определена только для квадратной матрицы")

    return rows, cols


def _normalize_number(value: Number, eps: float = 1e-12) -> Number:
    """
    Убирает погрешности вычислений с плавающей точкой

    Например:
        2.0000000000000004 -> 2
        1e-16 -> 0
    """
    if isinstance(value, float):
        if abs(value) < eps:
            return 0

        rounded = round(value)
        if abs(value - rounded) < eps:
            return int(rounded)

    return value


def _matrix_multiply(matrix_a: list[list], matrix_b: list[list]) -> list[list]:
    """
    Внутренняя функция перемножения матриц
    """
    rows_a, cols_a = _validate_matrix(matrix_a)
    rows_b, cols_b = _validate_matrix(matrix_b)

    if cols_a != rows_b:
        raise ValueError(
            "Количество столбцов первой матрицы должно совпадать "
            "с количеством строк второй матрицы"
        )

    return [
        [
            sum(matrix_a[row][k] * matrix_b[k][col] for k in range(cols_a))
            for col in range(cols_b)
        ]
        for row in range(rows_a)
    ]


def transpose_matrix(matrix: list[list]) -> list[list]:
    """
    Возвращает транспонированную матрицу
    Исходная матрица не изменяется
    """
    rows, cols = _validate_matrix(matrix)

    return [
        [matrix[row][col] for row in range(rows)]
        for col in range(cols)
    ]


def subtract_matrices(
    matrix_a: list[list],
    matrix_b: list[list]
) -> list[list]:
    """
    Возвращает разность matrix_a - matrix_b
    """
    rows_a, cols_a = _validate_matrix(matrix_a)
    rows_b, cols_b = _validate_matrix(matrix_b)

    if (rows_a, cols_a) != (rows_b, cols_b):
        raise ValueError(
            "Для вычитания матрицы должны иметь одинаковый размер"
        )

    return [
        [
            matrix_b[row][col] - matrix_a[row][col]
            for col in range(cols_a)
        ]
        for row in range(rows_a)
    ]


def solve_gaussian(
    coefficients: list[list],
    constants: list[Number],
    eps: float = 1e-12
) -> list[Number]:
    """
    Решает систему A * x = b методом Гаусса с выбором главного элемента

    Аргументы:
        coefficients — квадратная матрица коэффициентов A
        constants — вектор свободных членов b
        eps — точность сравнения с нулём

    Возвращает:
        Список решений [x1, x2, ..., xn]

    Выбрасывает ValueError, если система не имеет единственного решения
    """
    rows, cols = _validate_matrix(coefficients, square=True)

    if not isinstance(constants, list) or len(constants) != rows:
        raise ValueError(
            "Количество свободных членов должно совпадать "
            "с количеством строк матрицы"
        )

    if any(not isinstance(value, Number) for value in constants):
        raise ValueError("Свободные члены должны быть числами")

    # Создаём расширенную матрицу, не изменяя исходные данные
    augmented = [
        [float(value) for value in coefficients[row]]
        + [float(constants[row])]
        for row in range(rows)
    ]

    # Прямой ход метода Гаусса
    for col in range(cols):
        pivot_row = max(
            range(col, rows),
            key=lambda row: abs(augmented[row][col])
        )

        if abs(augmented[pivot_row][col]) < eps:
            raise ValueError(
                "Система не имеет единственного решения: "
                "матрица коэффициентов вырождена"
            )

        if pivot_row != col:
            augmented[col], augmented[pivot_row] = (
                augmented[pivot_row],
                augmented[col]
            )

        pivot = augmented[col][col]

        # Нормализуем ведущую строку
        for j in range(col, cols + 1):
            augmented[col][j] /= pivot

        # Обнуляем элементы под ведущим
        for row in range(col + 1, rows):
            factor = augmented[row][col]

            if abs(factor) < eps:
                augmented[row][col] = 0.0
                continue

            for j in range(col, cols + 1):
                augmented[row][j] -= factor * augmented[col][j]

    # Обратный ход
    solution = [0.0] * rows

    for row in range(rows - 1, -1, -1):
        solution[row] = augmented[row][cols] - sum(
            augmented[row][col] * solution[col]
            for col in range(row + 1, cols)
        )

    return [_normalize_number(value, eps) for value in solution]


def matrix_power(matrix: list[list], exponent: int) -> list[list]:
    """
    Возводит квадратную матрицу в целую неотрицательную степень  
    Используется быстрое возведение в степень  
    Для exponent == 0 возвращается единичная матрица
    """
    rows, cols = _validate_matrix(matrix, square=True)

    if not isinstance(exponent, int):
        raise ValueError("Показатель степени должен быть целым числом")

    if exponent < 0:
        raise ValueError(
            "Показатель степени должен быть неотрицательным"
        )

    result = [
        [1 if row == col else 0 for col in range(cols)]
        for row in range(rows)
    ]

    base = [row[:] for row in matrix]
    power = exponent

    while power > 0:
        if power % 2 == 1:
            result = _matrix_multiply(result, base)

        power //= 2

        if power > 0:
            base = _matrix_multiply(base, base)

    return result


def inverse_matrix(
    matrix: list[list],
    eps: float = 1e-12
) -> list[list]:
    """
    Возвращает обратную матрицу методом Гаусса — Жордана  
    Выбрасывает ValueError, если матрица вырождена  
    Исходная матрица не изменяется
    """
    rows, cols = _validate_matrix(matrix, square=True)
    size = rows

    augmented = []

    for row in range(size):
        identity_row = [
            1.0 if row == col else 0.0
            for col in range(size)
        ]

        augmented.append(
            [float(value) for value in matrix[row]] + identity_row
        )

    for col in range(size):
        pivot_row = max(
            range(col, size),
            key=lambda row: abs(augmented[row][col])
        )

        if abs(augmented[pivot_row][col]) < eps:
            raise ValueError(
                "Обратная матрица не существует: "
                "определитель равен нулю"
            )

        if pivot_row != col:
            augmented[col], augmented[pivot_row] = (
                augmented[pivot_row],
                augmented[col]
            )

        pivot = augmented[col][col]

        # Делим ведущую строку на ведущий элемент
        for j in range(size * 2):
            augmented[col][j] /= pivot

        # Обнуляем текущий столбец во всех остальных строках
        for row in range(size):
            if row == col:
                continue

            factor = augmented[row][col]

            if abs(factor) < eps:
                augmented[row][col] = 0.0
                continue

            for j in range(size * 2):
                augmented[row][j] -= factor * augmented[col][j]

    return [
        [
            _normalize_number(augmented[row][col], eps)
            for col in range(size, size * 2)
        ]
        for row in range(size)
    ]

