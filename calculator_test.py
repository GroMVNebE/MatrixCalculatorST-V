# Тесты над матричным калькулятором
# Здесь проводится тестирование MatrixCalculator.py

import pytest
from MatrixCalculator import *

def test_exmpl_func():
    """Тест для функции :func:`exmpl_func`"""
    assert exmpl_func(4) == 16
    assert exmpl_func(5) == 25
