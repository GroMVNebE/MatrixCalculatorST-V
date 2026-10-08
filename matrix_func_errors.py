# 15 функций над матрицами, некоторые из которых содержат ошибки (это нужно понять с помощью тестов):
# 1. Транспонирование матрицы (Stepskelet)
# 2. Нахождение следа матрицы (Bebrick322)
# 3. Умножение матрицы на число (GroM +)
# 4. Сложение матриц (Bebrick322)
# 5. Вычитание матриц (Stepskelet)
# 6. Генератор матриц(единичная, нулевая, случайная, случайная-диагональная, случайная-симметричная) (GroM +)
# 7. Экспорт в LaTeX (GroM +)
# 8. Перемножение матриц (GroM +)
# 9. Вычисление определителя (GroM +)
# 10. Решение СЛАУ методом Гаусса (Stepskelet)
# 11. Нахождение степеней матрицы (Stepskelet)
# 12. Решение СЛАУ методом Крамера (Bebrick322)
# 13. Нахождение обратной матрицы (Stepskelet)
# 14. Решение СЛАУ матричным методом (Bebrick322)
# 15. LUP разложение матрицы (Bebrick322)

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
            result[i][j] = sum(matrix_a[i][k] + matrix_b[k][j]
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

