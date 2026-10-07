import pytest
from matrix_func import (
    trace, add_matrices,
    lup_decomposition, solve_cramer,
    solve_matrix_method, inverse_matrix,
)


# ---------- Пункт 2: след ----------
class TestTrace:
    def test_trace_2x2(self):
        assert trace([[1, 2], [3, 4]]) == 5

    def test_trace_identity_3x3(self):
        assert trace([[1, 0, 0], [0, 1, 0], [0, 0, 1]]) == 3

    def test_trace_zero_matrix(self):
        assert trace([[0, 0], [0, 0]]) == 0

    def test_trace_non_square_raises(self):
        with pytest.raises(ValueError):
            trace([[1, 2, 3], [4, 5, 6]])

    def test_trace_empty_raises(self):
        with pytest.raises(ValueError):
            trace([])


# ---------- Пункт 4: сложение ----------
class TestAddMatrices:
    def test_add_simple(self):
        a = [[1, 2], [3, 4]]
        b = [[5, 6], [7, 8]]
        assert add_matrices(a, b) == [[6, 8], [10, 12]]

    def test_add_with_zero(self):
        a = [[1, 2], [3, 4]]
        z = [[0, 0], [0, 0]]
        assert add_matrices(a, z) == a

    def test_add_negative(self):
        a = [[1, -2], [3, 4]]
        b = [[-1, 2], [-3, -4]]
        assert add_matrices(a, b) == [[0, 0], [0, 0]]

    def test_add_row_count_mismatch_raises(self):
        with pytest.raises(ValueError):
            add_matrices([[1, 2], [3, 4]], [[1, 2]])

    def test_add_shape_mismatch_raises(self):
        with pytest.raises(ValueError):
            add_matrices([[1, 2]], [[1], [2]])


# ---------- Пункт 12: метод Крамера ----------
class TestCramer:
    def test_cramer_2x2(self):
        a = [[2, 1], [1, 3]]
        b = [5, 7]
        x = solve_cramer(a, b)
        assert x[0] == pytest.approx(1.6)
        assert x[1] == pytest.approx(1.8)

    def test_cramer_3x3(self):
        a = [[2, -1, 1], [3, 3, 9], [3, 3, 5]]
        b = [2, -1, 4]
        x = solve_cramer(a, b)
        assert x[0] == pytest.approx(-1.0)
        assert x[1] == pytest.approx(2.0)
        assert x[2] == pytest.approx(2.0)

    def test_cramer_singular_raises(self):
        with pytest.raises(ValueError):
            solve_cramer([[1, 2], [2, 4]], [1, 2])

    def test_cramer_wrong_b_length_raises(self):
        with pytest.raises(ValueError):
            solve_cramer([[1, 2], [3, 4]], [1, 2, 3])

    def test_cramer_non_square_raises(self):
        with pytest.raises(ValueError):
            solve_cramer([[1, 2, 3], [4, 5, 6]], [1, 2])


# ---------- Пункт 14: матричный метод ----------
class TestMatrixMethod:
    def test_matrix_method_2x2(self):
        a = [[2, 1], [1, 3]]
        b = [5, 7]
        x = solve_matrix_method(a, b)
        assert x[0] == pytest.approx(1.6)
        assert x[1] == pytest.approx(1.8)

    def test_matrix_method_identity(self):
        a = [[1, 0, 0], [0, 1, 0], [0, 0, 1]]
        b = [1, 2, 3]
        x = solve_matrix_method(a, b)
        assert x == pytest.approx([1.0, 2.0, 3.0])

    def test_matrix_method_singular_raises(self):
        with pytest.raises(ValueError):
            solve_matrix_method([[1, 2], [2, 4]], [1, 2])

    def test_matrix_method_wrong_b_length_raises(self):
        with pytest.raises(ValueError):
            solve_matrix_method([[1, 2], [3, 4]], [1, 2, 3])


# ---------- Пункт 15: LUP ----------
class TestLUP:
    def test_lup_reconstruct(self):
        a = [[4, 3], [6, 3]]
        L, U, P = lup_decomposition(a)
        n = len(a)
        pa = [[sum(P[i][k] * a[k][j] for k in range(n)) for j in range(n)] for i in range(n)]
        lu = [[sum(L[i][k] * U[k][j] for k in range(n)) for j in range(n)] for i in range(n)]
        for i in range(n):
            for j in range(n):
                assert pa[i][j] == pytest.approx(lu[i][j])

    def test_lup_lower_has_ones_on_diagonal(self):
        a = [[1, 2, 3], [4, 5, 6], [7, 8, 10]]
        L, U, P = lup_decomposition(a)
        for i in range(len(a)):
            assert L[i][i] == pytest.approx(1.0)

    def test_lup_lower_is_lower_triangular(self):
        a = [[1, 2, 3], [4, 5, 6], [7, 8, 10]]
        L, U, P = lup_decomposition(a)
        n = len(a)
        for i in range(n):
            for j in range(i + 1, n):
                assert L[i][j] == pytest.approx(0.0)

    def test_lup_upper_is_upper_triangular(self):
        a = [[1, 2, 3], [4, 5, 6], [7, 8, 10]]
        L, U, P = lup_decomposition(a)
        n = len(a)
        for i in range(n):
            for j in range(i):
                assert U[i][j] == pytest.approx(0.0)

    def test_lup_singular_raises(self):
        with pytest.raises(ValueError):
            lup_decomposition([[1, 2], [2, 4]])

    def test_inverse_matrix(self):
        a = [[2, 1], [1, 3]]
        inv = inverse_matrix(a)
        assert inv[0][0] == pytest.approx(0.6)
        assert inv[0][1] == pytest.approx(-0.2)
        assert inv[1][0] == pytest.approx(-0.2)
        assert inv[1][1] == pytest.approx(0.4)