import unittest

from matrix_func_errors import (
    inverse_matrix,
    matrix_power,
    solve_gaussian,
    subtract_matrices,
    transpose_matrix,
)


class TestMatrixFunctions(unittest.TestCase):

    def assertListAlmostEqual(self, actual, expected, places=9):
        self.assertEqual(len(actual), len(expected))
        for x, y in zip(actual, expected):
            self.assertAlmostEqual(x, y, places=places)

    def test_transpose_matrix(self):
        original = [[1, 2, 3], [4, 5, 6]]

        self.assertEqual(transpose_matrix(original), [[1, 4], [2, 5], [3, 6]])
        self.assertEqual(original, [[1, 2, 3], [4, 5, 6]])
        self.assertEqual(transpose_matrix(
            transpose_matrix(original)), original)

    def test_subtract_matrices(self):
        a = [[5, 7, 9], [4, 6, 8]]
        b = [[1, 2, 3], [6, 5, 4]]

        self.assertEqual(subtract_matrices(a, b), [[4, 5, 6], [-2, 1, 4]])
        self.assertEqual(subtract_matrices(b, a), [[-4, -5, -6], [2, -1, -4]])
        with self.assertRaises(ValueError):
            subtract_matrices([[1, 2]], [[1], [2]])

    def test_solve_gaussian(self):
        # 2x + y - z = 8;  -3x - y + 2z = -11;  -2x + y + 2z = -3
        a = [[2, 1, -1], [-3, -1, 2], [-2, 1, 2]]
        self.assertListAlmostEqual(solve_gaussian(a, [8, -11, -3]), [2, 3, -1])

        self.assertListAlmostEqual(solve_gaussian(
            [[0, 1], [1, 0]], [2, 3]), [3, 2])

        with self.assertRaises(ValueError):
            solve_gaussian([[1, 2], [2, 4]], [1, 2])

    def test_matrix_power(self):
        a = [[1, 2], [3, 4]]

        self.assertEqual(matrix_power(a, 0), [[1, 0], [0, 1]])
        self.assertEqual(matrix_power(a, 1), a)
        self.assertEqual(matrix_power(a, 3), [[37, 54], [81, 118]])
        with self.assertRaises(ValueError):
            matrix_power(a, -1)

    def test_inverse_matrix(self):
        a = [[4, 7], [2, 6]]  # det = 10
        inv = inverse_matrix(a)

        self.assertListAlmostEqual(inv[0], [0.6, -0.7])
        self.assertListAlmostEqual(inv[1], [-0.2, 0.4])

        product = [
            [sum(a[i][k] * inv[k][j] for k in range(2)) for j in range(2)]
            for i in range(2)
        ]
        self.assertListAlmostEqual(product[0], [1, 0])
        self.assertListAlmostEqual(product[1], [0, 1])

        with self.assertRaises(ValueError):
            inverse_matrix([[1, 2], [2, 4]])
