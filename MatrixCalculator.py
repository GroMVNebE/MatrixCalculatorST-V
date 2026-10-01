# Матричный калькулятор, реализующий 15 функций над матрицами:
# 1. Транспонирование матрицы (не назначено)
# 2. Нахождение следа матрицы (не назначено)
# 3. Умножение матрицы на число (GroM +)
# 4. Сложение матриц (не назначено)
# 5. Вычитание матриц (не назначено)
# 6. Генератор матриц(единичная, нулевая, случайная, случайная-диагональная, случайная-симметричная) (GroM +)
# 7. Экспорт в LaTeX (GroM +)
# 8. Перемножение матриц (GroM)
# 9. Вычисление определителя (GroM)
# 10. Решение СЛАУ методом Гаусса (не назначено)
# 11. Нахождение степеней матрицы (не назначено)
# 12. Решение СЛАУ методом Крамера (не назначено)
# 13. Нахождение обратной матрицы (не назначено)
# 14. Решение СЛАУ матричным методом (не назначено)
# 15. LUP разложение матрицы (не назначено)

from matrix_func import *


class IO_Unit:
    def __init__(self) -> None:
        self.commands = [
            {'name': 'Ввести матрицу', 'command': self.create_matrix_cmd},
            {'name': 'Вывести матрицу', 'command': self.draw_matrix_cmd,
                'matrix_required': True},
            {'name': 'Умножить матрицу на число', 'command': self.matrix_x_number_cmd,
                'matrix_required': True},
            {'name': 'Сгенерировать новую матрицу',
                'command': self.generate_matrix_cmd},
            {'name': 'Экспорт матрицы в LaTeX',
                'command': self.export_to_latex_cmd, 'matrix_required': True}
        ]
        self.matrix = None
        self.rows, self.cols = None, None
        self.run()

    def _prompt_int(self, prompt_text: str, min_val: int = None, max_val: int = None) -> int:
        while True:
            try:
                val = int(input(prompt_text).strip())
                if min_val is not None and val < min_val or max_val is not None and val > max_val:
                    print(
                        f'Число должно быть {f"не меньше {min_val}" if min_val is not None else ""} {"и" if all(val_lim is not None for val_lim in [min_val, max_val]) else ""} {f"не больше {max_val}" if max_val is not None else ""}')
                    continue
                return val
            except ValueError:
                print('Некорректный ввод! Введите целое число')

    def _validate_command(self, command: dict) -> bool:
        if command.get('matrix_required') is True and self.matrix is None:
            return False
        return True

    def run(self):
        while True:
            self.request_command()

    def request_command(self):
        idx = 0
        ava_commands = []
        print('Доступные команды:')
        for command in self.commands:
            if self._validate_command(command):
                print(f'{idx+1}. {command.get("name")}')
                ava_commands += [command]
                idx += 1
        cmd_idx = self._prompt_int(
            "Введите номер выбранной команды: ", 1, idx) - 1
        ava_commands[cmd_idx].get('command')()

    def create_matrix_cmd(self):
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
                self.draw_matrix_cmd()
                val = self._prompt_int(
                    f'Введите значение ячейки [{r + 1}, {c + 1}]: ')
                self.matrix[r][c] = val
        print('Итоговая матрица:')
        self.draw_matrix_cmd()

    def draw_matrix_cmd(self, matrix: list[list] = None):
        if matrix is None:
            matrix = self.matrix
        first_none = True
        # Сбор содержимого матрицы как строк
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
        # Вычисление ширины ячейки
        cell_width = max(
            len(cell) for row in formatted_grid for cell in row
        )
        cell_width = max(cell_width, 1) + 2
        # Создание горизонтальных границ
        border = '+' + '+'.join(['-' * cell_width] * len(matrix[0])) + '+'
        # Отрисовка матрицы в виде таблицы
        print(border)
        for row in formatted_grid:
            row_str = '|' + '|'.join(cell.center(cell_width)
                                     for cell in row) + '|'
            print(row_str)
            print(border)

    def matrix_x_number_cmd(self):
        print('Исходная матрица:')
        self.draw_matrix_cmd()

        value = self._prompt_int(
            'Введите целое число, на которое будет умножена матрица: ')
        try:
            self.matrix = matrix_x_number(self.matrix, value)
            print('Умножение матрицы на число успешно выполнено! Результат:')
            self.draw_matrix_cmd()
        except Exception as e:
            print(f'В ходе умножения матрицы на число произошла ошибка: {e}')

    def generate_matrix_cmd(self):
        print("1. Единичная матрица (квадратная)")
        print("2. Нулевая матрица")
        print("3. Случайная матрица")
        print("4. Случайная диагональная матрица")
        print("5. Случайная симметричная матрица")

        choice = self._prompt_int("Выберите тип матрицы (1-5): ", 1, 5)
        min_rnd, max_rnd = -100, 100
        if choice >= 3:
            min_rnd = self._prompt_int(
                'Введите минимальное значение элемента: ')
            max_rnd = self._prompt_int(
                'Введите максимальное значение элемента: ')
            min_rnd, max_rnd = min(min_rnd, max_rnd), max(min_rnd, max_rnd)

        match(choice):
            case 1:
                n = self._prompt_int(
                    "Введите размерность квадратной матрицы n: ", 1)
                self.rows, self.cols = n, n
                self.matrix = generate_identity_matrix(n)
            case 2:
                rows = self._prompt_int("Введите количество строк: ", 1)
                cols = self._prompt_int("Введите количество столбцов: ", 1)
                self.rows, self.cols = rows, cols
                self.matrix = generate_zero_matrix(rows, cols)
            case 3:
                rows = self._prompt_int("Введите количество строк: ", 1)
                cols = self._prompt_int("Введите количество столбцов: ", 1)
                self.rows, self.cols = rows, cols
                self.matrix = generate_random_matrix(
                    rows, cols, min_rnd, max_rnd)
            case 4:
                n = self._prompt_int(
                    "Введите размерность квадратной матрицы n: ", 1)
                self.rows, self.cols = n, n
                self.matrix = generate_random_diagonal_matrix(
                    n, min_rnd, max_rnd)
            case 5:
                n = self._prompt_int(
                    "Введите размерность квадратной матрицы n: ", 1)
                self.rows, self.cols = n, n
                self.matrix = generate_random_symmetric_matrix(
                    n, min_rnd, max_rnd)

        print("Сгенерированная матрица:")
        self.draw_matrix_cmd()

    def export_to_latex_cmd(self):
        print("Выберите стиль скобок:")
        print("1. Квадратные скобки [ ] (bmatrix)")
        print("2. Круглые скобки ( ) (pmatrix)")
        print("3. Определитель | | (vmatrix)")
        print("4. Без скобок (matrix)")

        choice = self._prompt_int("Выберите вариант (1-4): ", 1, 4)

        match choice:
            case 1: env = "bmatrix"
            case 2: env = "pmatrix"
            case 3: env = "vmatrix"
            case 4: env = "matrix"
            case _: env = "bmatrix"

        latex_code = matrix_to_latex(self.matrix, env)

        print("Сгенерированный LaTeX-код:")
        print(latex_code)
        print("Скопируйте этот код для вставки в любой LaTeX-редактор")


unit = IO_Unit()
