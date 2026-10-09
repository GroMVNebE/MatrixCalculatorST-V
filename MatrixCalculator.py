# Матричный калькулятор, реализующий 15 функций над матрицами:
# 1. Транспонирование матрицы (Stepskelet +)
# 2. Нахождение следа матрицы (Bebrick322 +)
# 3. Умножение матрицы на число (GroM +)
# 4. Сложение матриц (Bebrick322 +)
# 5. Вычитание матриц (Stepskelet +)
# 6. Генератор матриц(единичная, нулевая, случайная, случайная-диагональная, случайная-симметричная) (GroM +)
# 7. Экспорт в LaTeX (GroM +)
# 8. Перемножение матриц (GroM +)
# 9. Вычисление определителя (GroM +)
# 10. Решение СЛАУ методом Гаусса (Stepskelet +)
# 11. Нахождение степеней матрицы (Stepskelet +)
# 12. Решение СЛАУ методом Крамера (Bebrick322 +)
# 13. Нахождение обратной матрицы (Stepskelet +)
# 14. Решение СЛАУ матричным методом (Bebrick322 +)
# 15. LUP разложение матрицы (Bebrick322 +)

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
                'command': self.export_to_latex_cmd, 'matrix_required': True},
            {'name': 'Перемножить матрицы', 'command': self.matrix_x_matrix_cmd,
                'matrix_required': True},
            {'name': 'Вычислить определитель', 'command': self.calculate_determinant_cmd,
                'matrix_required': True, 'is_square_matrix': True},
            {'name': 'Найти след матрицы', 'command': self.trace_cmd,
                'matrix_required': True, 'is_square_matrix': True},
            {'name': 'Сложить матрицы', 'command': self.add_matrices_cmd,
                'matrix_required': True},
            {'name': 'Решить СЛАУ методом Крамера', 'command': self.cramer_cmd,
                'matrix_required': True, 'is_square_matrix': True},
            {'name': 'Решить СЛАУ матричным методом', 'command': self.matrix_method_cmd,
                'matrix_required': True, 'is_square_matrix': True},
            {'name': 'LUP-разложение матрицы', 'command': self.lup_cmd,
                'matrix_required': True, 'is_square_matrix': True},

            # Команды функций Stepskelet
            {'name': 'Транспонировать матрицу', 'command': self.transpose_matrix_cmd,
             'matrix_required': True},
            {'name': 'Вычесть матрицу', 'command': self.subtract_matrix_cmd,
             'matrix_required': True},
            {'name': 'Решить СЛАУ методом Гаусса', 'command': self.solve_gaussian_cmd,
             'matrix_required': True, 'is_square_matrix': True},
            {'name': 'Возвести матрицу в степень', 'command': self.matrix_power_cmd,
             'matrix_required': True, 'is_square_matrix': True},
            {'name': 'Найти обратную матрицу', 'command': self.inverse_matrix_cmd,
             'matrix_required': True, 'is_square_matrix': True},
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

    def _prompt_vector(self, n: int, name: str = "b") -> list:
        print(
            f'Введите вектор {name} длины {n} (каждое значение — с новой строки):')
        vector = []
        for i in range(n):
            val = self._prompt_int(f'{name}[{i + 1}] = ')
            vector.append(val)
        return vector

    def _validate_command(self, command: dict) -> bool:
        if command.get('matrix_required') is True and self.matrix is None:
            return False
        if command.get('is_square_matrix') is True:
            if self.matrix is None or self.rows != self.cols:
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
        print('Введите новую матрицу, которая будет использоваться в операциях')
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

        self.matrix = [[None] * self.cols for _ in range(self.rows)]
        for r in range(self.rows):
            for c in range(self.cols):
                self.draw_matrix_cmd()
                val = self._prompt_int(
                    f'Введите значение ячейки [{r + 1}, {c + 1}]: ')
                self.matrix[r][c] = val
        print('Введённая матрица:')
        self.draw_matrix_cmd()

    def draw_matrix_cmd(self, matrix: list[list] = None):
        if matrix is None:
            matrix = self.matrix
        first_none = True
        formatted_grid = []
        for row in matrix:
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
        border = '+' + '+'.join(['-' * cell_width] * len(matrix[0])) + '+'
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

    def matrix_x_matrix_cmd(self):
        matrix_a = [row[:] for row in self.matrix]
        cols_a = self.cols

        self.create_matrix_cmd()
        matrix_b = self.matrix

        if cols_a != self.rows:
            print(f"Ошибка: Перемножение невозможно!")
            print(
                f"Число столбцов первой матрицы ({cols_a}) должно совпадать с числом строк второй матрицы ({self.rows})")
            self.matrix = matrix_a
            return

        try:
            result = matrix_x_matrix(matrix_a, matrix_b)
            self.matrix = result
            self.rows = len(result)
            self.cols = len(result[0])
            print("Результат перемножения матриц:")
            self.draw_matrix_cmd()
        except Exception as e:
            print(f"В ходе перемножения матриц произошла ошибка: {e}")
            self.matrix = matrix_a

    def calculate_determinant_cmd(self):
        print('Исходная матрица:')
        self.draw_matrix_cmd()
        try:
            det = calculate_determinant(self.matrix)
            print(f"Определитель матрицы det(A) = {det}")
        except Exception as e:
            print(f"При вычислении определителя произошла ошибка: {e}")

    def trace_cmd(self):
        print('Исходная матрица:')
        self.draw_matrix_cmd()
        try:
            tr = trace(self.matrix)
            print(f"След матрицы tr(A) = {tr}")
        except Exception as e:
            print(f"При вычислении следа произошла ошибка: {e}")

    def add_matrices_cmd(self):
        matrix_a = [row[:] for row in self.matrix]
        rows_a, cols_a = self.rows, self.cols

        print('Введите вторую матрицу для сложения:')
        self.create_matrix_cmd()
        matrix_b = self.matrix

        if rows_a != self.rows or cols_a != self.cols:
            print("Ошибка: размеры матриц не совпадают!")
            print(
                f"Первая: {rows_a}x{cols_a}, вторая: {self.rows}x{self.cols}")
            self.matrix = matrix_a
            self.rows, self.cols = rows_a, cols_a
            return

        try:
            result = add_matrices(matrix_a, matrix_b)
            self.matrix = result
            self.rows = len(result)
            self.cols = len(result[0])
            print("Результат сложения матриц:")
            self.draw_matrix_cmd()
        except Exception as e:
            print(f"В ходе сложения матриц произошла ошибка: {e}")

    def cramer_cmd(self):
        print('Матрица системы:')
        self.draw_matrix_cmd()
        b = self._prompt_vector(self.rows, "b")
        try:
            x = solve_cramer(self.matrix, b)
            print("Решение СЛАУ методом Крамера:")
            for i, val in enumerate(x):
                print(f"  x{i + 1} = {val}")
        except Exception as e:
            print(f"При решении СЛАУ методом Крамера произошла ошибка: {e}")

    def matrix_method_cmd(self):
        print('Матрица системы:')
        self.draw_matrix_cmd()
        b = self._prompt_vector(self.rows, "b")
        try:
            x = solve_matrix_method(self.matrix, b)
            print("Решение СЛАУ матричным методом:")
            for i, val in enumerate(x):
                print(f"  x{i + 1} = {val}")
        except Exception as e:
            print(f"При решении СЛАУ матричным методом произошла ошибка: {e}")

    def lup_cmd(self):
        print('Исходная матрица:')
        self.draw_matrix_cmd()
        try:
            L, U, P = lup_decomposition(self.matrix)
            print("Матрица L (нижнетреугольная):")
            self.draw_matrix_cmd(L)
            print("Матрица U (верхнетреугольная):")
            self.draw_matrix_cmd(U)
            print("Матрица перестановок P:")
            self.draw_matrix_cmd(P)
        except Exception as e:
            print(f"При LUP-разложении произошла ошибка: {e}")
    # Интерфейс функций Stepskelet

    def _input_matrix_of_size(self, rows: int, cols: int, title: str) -> list[list]:
        """
        Вспомогательный интерфейсный метод

        Считывает дополнительную матрицу, но не заменяет self.matrix  
        Это важно для операций с двумя матрицами
        """
        print(title)
        matrix = [[None] * cols for _ in range(rows)]

        for row in range(rows):
            for col in range(cols):
                self.draw_matrix_cmd(matrix)
                matrix[row][col] = self._prompt_int(
                    f'Введите значение ячейки [{row + 1}, {col + 1}]: '
                )

        return matrix

    def transpose_matrix_cmd(self):
        """Интерфейс транспонирования матрицы"""
        print('Исходная матрица:')
        self.draw_matrix_cmd()

        try:
            result = transpose_matrix(self.matrix)

            self.matrix = result
            self.rows, self.cols = self.cols, self.rows

            print('Транспонирование успешно выполнено. Результат:')
            self.draw_matrix_cmd()
        except (TypeError, ValueError) as error:
            print(f'Не удалось транспонировать матрицу: {error}')

    def subtract_matrix_cmd(self):
        """Интерфейс вычитания из текущей матрицы другой матрицы"""
        print('Первая матрица A:')
        self.draw_matrix_cmd()

        matrix_b = self._input_matrix_of_size(
            self.rows,
            self.cols,
            'Введите матрицу B того же размера. Будет вычислено A - B.'
        )

        print('Вторая матрица B:')
        self.draw_matrix_cmd(matrix_b)

        try:
            result = subtract_matrices(self.matrix, matrix_b)

            self.matrix = result

            print('Вычитание успешно выполнено. Результат A - B:')
            self.draw_matrix_cmd()
        except (TypeError, ValueError) as error:
            print(f'Не удалось выполнить вычитание: {error}')

    def solve_gaussian_cmd(self):
        """
        Интерфейс решения СЛАУ

        Текущая матрица используется как матрица коэффициентов A  
        Решение не заменяет текущую матрицу
        """
        print('Матрица коэффициентов системы A:')
        self.draw_matrix_cmd()

        constants = []

        print('Введите столбец свободных членов b:')

        for row in range(self.rows):
            value = self._prompt_int(
                f'b[{row + 1}] = '
            )
            constants.append(value)

        try:
            solution = solve_gaussian(self.matrix, constants)

            print('Решение системы:')
            for index, value in enumerate(solution, start=1):
                print(f'x{index} = {value}')

        except (TypeError, ValueError) as error:
            print(f'Не удалось решить систему методом Гаусса: {error}')

    def matrix_power_cmd(self):
        """Интерфейс возведения матрицы в степень"""
        print('Исходная матрица:')
        self.draw_matrix_cmd()

        exponent = self._prompt_int(
            'Введите целую неотрицательную степень: ',
            min_val=0
        )

        try:
            result = matrix_power(self.matrix, exponent)

            self.matrix = result
            self.rows = len(result)
            self.cols = len(result[0])

            print(
                f'Матрица успешно возведена в степень {exponent}. Результат:')
            self.draw_matrix_cmd()
        except (TypeError, ValueError) as error:
            print(f'Не удалось возвести матрицу в степень: {error}')

    def inverse_matrix_cmd(self):
        """Интерфейс нахождения обратной матрицы"""
        print('Исходная матрица:')
        self.draw_matrix_cmd()

        try:
            result = inverse_matrix(self.matrix)

            self.matrix = result
            self.rows = len(result)
            self.cols = len(result[0])

            print('Обратная матрица успешно найдена:')
            self.draw_matrix_cmd()
        except (TypeError, ValueError) as error:
            print(f'Не удалось найти обратную матрицу: {error}')


unit = IO_Unit()
