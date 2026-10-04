# Матричный калькулятор, реализующий 15 функций над матрицами:
# 1. Транспонирование матрицы (Bebrick322)
# 2. Нахождение следа матрицы (Stepskelet)
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
