import pytest
from matrix_func_errors import (
    trace, add_matrices,
    lup_decomposition, solve_cramer,
    solve_matrix_method, inverse_matrix,
)


class TestTrace:
    def test_trace_2x2(self):
        assert trace([[1, 2], [3, 4]]) == 5

    def test_trace_identity_3x3(self):
        assert trace([[1, 0, 0], [0, 1, 0], [0, 0, 1]]) == 3


class TestAddMatrices:
    # ЭТОТ ТЕСТ ПАДАЕТ на версии с ошибкой
    def test_add_simple(self):
        a = [[1, 2], [3, 4]]
        b = [[5, 6], [7, 8]]
        assert add_matrices(a, b) == [[6, 8], [10, 12]]

    def test_add_with_zero(self):
        a = [[1, 2], [3, 4]]
        z = [[0, 0], [0, 0]]
        assert add_matrices(a, z) == a


class TestCramer:
    def test_cramer_2x2(self):
        a = [[2, 1], [1, 3]]
        b = [5, 7]
        x = solve_cramer(a, b)
        assert x[0] == pytest.approx(1.6)
        assert x[1] == pytest.approx(1.8)


class TestMatrixMethod:
    def test_matrix_method_2x2(self):
        a = [[2, 1], [1, 3]]
        b = [5, 7]
        x = solve_matrix_method(a, b)
        assert x[0] == pytest.approx(1.6)
        assert x[1] == pytest.approx(1.8)


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