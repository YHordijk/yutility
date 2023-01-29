from collections.abc import Iterable
from typing import get_type_hints
import numpy as np
import inspect
import typing
from functools import wraps
import warnings


enabled: bool = False


def check_hints(func):
    '''Decorator used to enforce type-hints. Apply to a 
    function and at run-time it will check if the values 
    given to the function are of the correct type. Running 
    python with the ``-o`` flag will disable type-checking 
    as will setting ``ychem.ytypes.enabled`` to ``False``.'''
    if enabled and __debug__:
        @wraps(func)
        def inner(*args, **kwargs):
            arg_names = inspect.getfullargspec(func).args
            kwargs_ = kwargs.copy()
            kwargs_.update({arg_name: arg for arg_name,
                           arg in zip(arg_names, args)})
            _check_types(func, **kwargs_)
            ret = func(*args, **kwargs)
            _check_types(func, **{'return': ret})
            return ret
        return inner
    else:
        return func


class TypeHintError(Exception):
    def __init__(self, cls, var_name, val, msg=None):
        self.cls = cls
        self.var_name = var_name
        self.val = val
        self.msg = msg

    def __str__(self):
        if self.msg is None:
            return f"Value {self.var_name} = {self.val} with type '{type(self.val).__name__}' is not of type '{self.cls}'"
        return self.msg


class ClassPrinter(type):
    def __repr__(cls):
        return cls.__name__.split('.')[-1]


class GenericType(metaclass=ClassPrinter):
    def __repr__(self):
        return type(self).__name__.split('.')[-1]


class Vector(GenericType):
    def __init__(self, *args):
        self.shape = args
        self.check = self.check_init

    @property
    def __name__(self):
        return repr(self)

    @classmethod
    def __class_name__(cls):
        return 'Vector'

    def __repr__(self):
        return f'{type(self).__name__.split(".")[-1]}({",".join(str(s) for s in self.shape)})'

    def check_init(self, var_name, val):
        if isinstance(val, np.ndarray):
            if val.ndim > 1:
                raise TypeHintError(type(self), var_name, val, f'Vector {var_name} with shape {val.shape} has more than one dimension')
            val = val.tolist()
        if not isinstance(val, Iterable):
            raise TypeHintError(type(self), var_name, val)
        if not all([type(x) in [int, float] for x in val]):
            raise TypeHintError(type(self), var_name, val)

        if np.array(val).size != self.shape[0] and self.shape[0] is not ...:
            raise TypeHintError(type(self), val, f'Size of Vector {var_name} is {np.array(val).size} instead of {self.shape[0]}')
        return True

    @classmethod
    def check(cls, var_name, val):
        if isinstance(val, np.ndarray):
            if val.ndim > 1:
                raise TypeHintError(cls, var_name, val, f'Vector {var_name} with shape {val.shape} has more than one dimension')
            val = val.tolist()
        if not isinstance(val, Iterable):
            raise TypeHintError(cls, var_name, val)
        if not all([type(x) in [int, float] for x in val]):
            raise TypeHintError(cls, var_name, val)
        return True


class Matrix(GenericType):
    def __init__(self, *args):
        self.shape = args
        self.check = self.check_init

    def check_init(self, var_name, val):
        if isinstance(val, np.ndarray):
            if val.ndim > 2:
                raise TypeHintError(type(self), var_name, val, f'Matrix {var_name} with shape {val.shape} has more than two dimensions')
            val = val.tolist()
        if not (isinstance(val, Iterable) and isinstance(val[0], Iterable)):
            raise TypeHintError(type(self), var_name, val)
        if not all([[type(x) in [int, float] for x in row] for row in val]):
            raise TypeHintError(type(self), var_name, val)

        mval = np.array(val, ndmin=2)
        if mval.shape[0] != self.shape[0] and self.shape[0] is not ...:
            raise TypeHintError(type(self), var_name, val, f'Dimesion 0 of Matrix {var_name} has size {mval.shape[0]} instead of {self.shape[0]}')
        if mval.shape[1] != self.shape[1] and self.shape[1] is not ...:
            raise TypeHintError(type(self), var_name, val, f'Dimesion 1 of Matrix {var_name} has size {mval.shape[1]} instead of {self.shape[1]}')
        return True

    @classmethod
    def check(cls, var_name, val):
        if isinstance(val, np.ndarray):
            if val.ndim > 2:
                raise TypeHintError(cls, var_name, val, f'Matrix {var_name} with shape {val.shape} has more than two dimensions')
            val = val.tolist()
        if not (isinstance(val, Iterable) and isinstance(val[0], Iterable)):
            raise TypeHintError(cls, var_name, val)
        if not all([[type(x) in [int, float] for x in row] for row in val]):
            raise TypeHintError(cls, var_name, val)
        return True



# class Scalar(GenericType):
    # @classmethod
    # def check(cls, var_name, val):
    #     if isinstance(val, np.ndarray):
    #         if type(val[0].item()) not in [int, float]:
    #             raise TypeHintError(cls, var_name, val)
    #     else:
    #         if type(val) not in [int, float]:
    #             raise TypeHintError(cls, var_name, val)
    #     return True



class List(GenericType):
    def __init__(self, *args):
        self.types = args

    def __repr__(self):
        return f'{type(self)}({", ".join(typ.__name__ for typ in self.types)})'

    @property
    def __name__(self):
        return repr(self)

    def check(self, var_name, val):
        if isinstance(val, np.ndarray):
            val = val.tolist()
        if not isinstance(val, (list, tuple, set)):
            raise TypeHintError(repr(self), var_name, val)
        if not all([self.check_val(var_name, x, i) for i, x in enumerate(val)]):
            raise TypeHintError(repr(self), var_name, val)
        return True

    def check_val(self, var_name, val, index=None):
        for typ in self.types:
            try:
                _check(typ, var_name, val)
                return True
            except TypeHintError:
                pass
        if index is None:
            raise TypeHintError(repr(self), var_name, val)
        raise TypeHintError(repr(self), var_name + f'[{index}]', val)


class Tuple(GenericType):
    def __init__(self, *args):
        self.types = args

    def __repr__(self):
        return f'{type(self)}({", ".join(typ.__name__ for typ in self.types)})'

    @property
    def __name__(self):
        return repr(self)

    def check(self, var_name, val):
        if not isinstance(val, tuple):
            raise TypeHintError(type(self), var_name, val)
        if not all([_check(typ, var_name, x) for x, typ in zip(val, self.types)]):
            raise TypeHintError(type(self), var_name, val)
        return True


class Dict(GenericType):
    def __init__(self, *args):
        self.kv_pairs = args

    def __repr__(self):
        return f'{type(self)}({", ".join(typ[0].__name__ + ": " + typ[0].__name__ for typ in self.kv_pairs)})'

    @property
    def __name__(self):
        return repr(self)

    def check(self, var_name, val):
        if not isinstance(val, dict):
            raise TypeHintError(type(self), var_name, val)
        for k, v in val.items():
            if not any(_check(kallow, var_name, k) and _check(vallow, var_name, v) for kallow, vallow in self.kv_pairs):
                raise TypeHintError(type(self), var_name, val)
        return True


class String(GenericType):
    @classmethod
    def check(cls, var_name, val):
        if not isinstance(val, str):
            raise TypeHintError(cls, var_name, val)
        return True


class Either(GenericType):
    def __init__(self, *args):
        self.types = args

    def check(self, var_name, val):
        for typ in self.types:
            try:
                _check(typ, var_name, val)
                return True
            except TypeHintError as err:
                pass
        raise TypeHintError(typ, var_name, val)

    def __repr__(self):
        return f'Either({", ".join([repr(typ) for typ in self.types])})'


class Any(GenericType):
    @classmethod
    def check(cls, var_name, val):
        return True


Scalar = Either(float, int)


def _check(typ, key, val):
    # check if typ is actually a class or object
    if not isinstance(typ, type):  # True if object
        typ_ = type(typ)
    else:
        typ_ = typ

    # None is special
    if typ is None and val is None:
        return True

    if typ_ in GenericType.__subclasses__():
        typ.check(key, val)
    else:
        if not isinstance(val, typ):
            raise TypeHintError(typ, key, val)
    return True


def _check_types(obj, **kwargs):
    # thint = typing.get_type_hints(obj, include_extras=True)
    thint = obj.__annotations__
    for key, val in kwargs.items():
        if val is None:
            return
        if key not in thint:
            continue
        typ = thint[key]
        # check if typ is optional
        if hasattr(typ, '__origin__'):
            if typ.__origin__ is typing.Union:
                if val is None:
                    return
                allowed = typ.__args__

                for allowtyp in allowed:
                    try:
                        _check(allowtyp, key, val)
                        return
                    except BaseException:
                        raise
                raise TypeHintError(typ, key, val)

        _check(thint[key], key, val)


def is_ytype(obj, typs):
    if not isinstance(typs, (tuple, list)):
        typs = [typs]
    for typ in typs:
        try:
            _check(typ, '', obj)
            return True
        except:
            pass
    return False


def is_iterable(obj):
    return hasattr(obj, '__iter__') and not isinstance(obj, str)


if __name__ == '__main__':
    print(is_ytype(np.array([0, 1, 2, 3]), List(int)))
