import ytypes
from ytypes import Vector, check_hints, Matrix, Scalar, List, Either, Any, Tuple, TypeHintError
import unittest
import numpy as np


ytypes.enabled = True


class TestClass:
    ...


class TestNoHints(unittest.TestCase):
    def test_args_only(self):
        @check_hints
        def one_arg(a):
            return

        @check_hints
        def two_args(a, b):
            return

        self.assertIsNone(one_arg(1))
        self.assertIsNone(one_arg(1.1))
        self.assertIsNone(one_arg('1'))

        self.assertIsNone(two_args(1, 1))
        self.assertIsNone(two_args(1.1, TestClass()))
        self.assertIsNone(two_args('1', ['9', 91]))

    def test_kwargs_only(self):
        @check_hints
        def one_kwarg(*, a):
            return

        self.assertIsNone(one_kwarg(a=1))
        self.assertIsNone(one_kwarg(a=1.1))
        self.assertIsNone(one_kwarg(a='1'))


class TestArgs(unittest.TestCase):
    def test_scalar(self):
        @check_hints
        def one_hint(scalar1: Scalar, scalar2):
            return

        self.assertIsNone(one_hint(1, 2))
        self.assertIsNone(one_hint(1.1, 2))
        self.assertIsNone(one_hint(1, 2.1))
        self.assertIsNone(one_hint(1.1, 2.1))
        with self.assertRaises(TypeHintError):
            one_hint([1, 2], 1.1)
        with self.assertRaises(TypeHintError):
            one_hint('1', 1.1)

    def test_pytype_str(self):
        @check_hints
        def _str(a: str):
            return
        self.assertIsNone(_str('abc'))
        with self.assertRaises(TypeHintError):
            _str(1)

    def test_pytype_float(self):
        @check_hints
        def _float(a: float):
            return
        self.assertIsNone(_float(1.2))
        self.assertIsNone(_float(1e-1))
        self.assertIsNone(_float(1e10))
        with self.assertRaises(TypeHintError):
            _float(1)
        with self.assertRaises(TypeHintError):
            _float('abc')

    def test_pytype_int(self):
        @check_hints
        def _int(a: int):
            return
        self.assertIsNone(_int(-1))
        self.assertIsNone(_int(10))
        with self.assertRaises(TypeHintError):
            _int(1.31)
        with self.assertRaises(TypeHintError):
            _int(-1.3e-23)
        with self.assertRaises(TypeHintError):
            _int('abc')

    def test_pytype_list(self):
        @check_hints
        def _list(a: list):
            return
        self.assertIsNone(_list([1, 2, 3]))
        self.assertIsNone(_list([]))
        with self.assertRaises(TypeHintError):
            _list(())
        with self.assertRaises(TypeHintError):
            _list((1, [0, 1]))

    def test_either(self):
        @check_hints
        def _either_str_float(a: Either(str, float)):
            return
        self.assertIsNone(_either_str_float(1.1))
        self.assertIsNone(_either_str_float('abc'))
        with self.assertRaises(TypeHintError):
            _either_str_float(())
        with self.assertRaises(TypeHintError):
            _either_str_float(1)

    def test_List(self):
        @check_hints
        def _List_int(a: List(int)):
            return
        self.assertIsNone(_List_int([1, 2, 3, 4]))
        self.assertIsNone(_List_int(np.array([1, 2, 3, 4])))
        self.assertIsNone(_List_int((1, 2, 3, 4)))
        self.assertIsNone(_List_int({1, 2, 3, 4}))
        with self.assertRaises(TypeHintError):
            _List_int([1.1, 2, 5])
        with self.assertRaises(TypeHintError):
            _List_int(1)

        @check_hints
        def _List_int(a: List(Any)):
            return
        self.assertIsNone(_List_int([1, 2.1, '3', TestClass()]))


unittest.main()
