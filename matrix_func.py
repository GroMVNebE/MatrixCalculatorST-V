# 1. Транспонирование матрицы (не назначено)
# 2. Нахождение следа матрицы (не назначено)
# 3. Умножение матрицы на число (GroM)
# 4. Сложение матриц (не назначено)
# 5. Вычитание матриц (не назначено)
# 6. Генератор матриц(единичная, нулевая, случайная, диагональная, симметричная) (GroM)
# 7. Экспорт в LaTeX (GroM)
# 8. Перемножение матриц (GroM)
# 9. Вычисление определителя (GroM)
# 10. Решение СЛАУ методом Гаусса (не назначено)
# 11. Нахождение степеней матрицы (не назначено)
# 12. Решение СЛАУ методом Крамера (не назначено)
# 13. Нахождение обратной матрицы (не назначено)
# 14. Решение СЛАУ матричным методом (не назначено)
# 15. LUP разложение матрицы (не назначено)

import random


def matrix_x_number(matrix: list[list], value: int) -> list[list]:
    return [[cell * value for cell in row] for row in matrix]


def generate_identity_matrix(n: int) -> list[list[int]]:
    return [[1 if i == j else 0 for j in range(n)] for i in range(n)]


def generate_zero_matrix(rows: int, cols: int) -> list[list[int]]:
    return [[0 for _ in range(cols)] for _ in range(rows)]


def generate_random_matrix(rows: int, cols: int, min_val: int = -100, max_val: int = -100) -> list[list[int]]:
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
