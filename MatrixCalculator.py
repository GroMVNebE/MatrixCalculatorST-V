# Матричный калькулятор, реализующий 15 функций над матрицами:
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

def exmpl_func(x: int) -> int:
    """Функция-пример, реализующая возведение в квадрат
    и содержащая намеренную ошибку при попытке возведения в квадрат числа 5"""
    if x == 5:
        return 228
    return x**2


# ===================
# Модуль ввода-вывода
# ===================

class IO_Unit:
    def __init__(self) -> None:
        self.create_matrix()

    def _prompt_int(self, prompt_text: str, min_val: int = None) -> int:
        while True:
            try:
                val = int(input(prompt_text).strip())
                if min_val is not None and val < min_val:
                    print(f'Число должно быть не меньше {min_val}!')
                    continue
                return val
            except ValueError:
                print('Некорректный ввод! Введите целое число.')

    def create_matrix(self):
        # Ввод размера матрицы
        print('Введите Матрицу, над которой хотите проводить операции')
        print('Начнём с размера: введите через пробел кол-во строк и столбцов в Матрице (пр. "3 4")')
        rows, cols = 0, 0
        while True:
            try:
                parts = input().strip().split()
                if len(parts) != 2:
                    raise ValueError
                rows, cols = int(parts[0]), int(parts[1])
                if rows <= 0 or cols <= 0:
                    raise ValueError
                break
            except ValueError:
                print(
                    'Некорректный ввод! Введите два целых положительных числа через пробел:')
        self.rows, self.cols = rows, cols

        # Создание пустой матрицы
        self.matrix = [[None] * self.cols for _ in range(self.rows)]
        # Заполнение матрицы значениями
        for r in range(self.rows):
            for c in range(self.cols):
                self.draw_matrix()
                val = self._prompt_int(
                    f'Введите значение ячейки [{r + 1}, {c + 1}]: ')
                self.matrix[r][c] = val
        print('Итоговая матрица:')
        self.draw_matrix()

    def draw_matrix(self):
        first_none = True
        formatted_grid = []
        for row in self.matrix:
            formatted_row = []
            for val in row:
                if val is not None:
                    formatted_row.append(str(val))
                else:
                    formatted_row.append('X' if first_none else 'O')
                    first_none = False
            formatted_grid.append(formatted_row)

        cell_width = max(
            len(cell) for row in formatted_grid for cell in row
        )
        cell_width = max(cell_width, 1) + 2

        border = '+' + '+'.join(['-' * cell_width] * self.cols) + '+'

        print(border)
        for row in formatted_grid:
            row_str = '|' + '|'.join(cell.center(cell_width)
                                     for cell in row) + '|'
            print(row_str)
            print(border)


unit = IO_Unit()
