from time import perf_counter
import numpy as np
from yutility import units, log
import sys


def add_to_func(**kwargs):
    def decorator(func):
        for key, val in kwargs.items():
            setattr(func, key, val)
        return func
    return decorator


def timer(*args, **kwargs):
    if len(args) == 1 and len(kwargs) == 0 and callable(args[0]):
        func = args[0]

        def inner(*args, **kwargs):
            start = perf_counter()
            res = func(*args, **kwargs)
            time = perf_counter() - start
            prefix, factor = units.get_si_prefix(time)
            log.log(f't({func.__name__}) = {time*factor:.2f} {prefix}s')
            return res

        return inner
    else:
        repeats = args[0]

        def decorator(func):
            def inner(*args, **kwargs):
                times = []
                for _ in range(repeats):
                    start = perf_counter()
                    res = func(*args, **kwargs)
                    times.append(perf_counter() - start)
                time_mean = np.mean(times)
                prefix, factor = units.get_si_prefix(time_mean)
                log.log(
                    f't({func.__name__}, {repeats}) = {np.mean(times)*factor:.2f} +- {np.std(times)*factor:.2f} {prefix}s')
                return res
            return inner
        return decorator


def return_time(func):
    def inner(*args, **kwargs):
        unit = units.Time('s')
        start = perf_counter()
        func(*args, **kwargs)
        return unit.string(perf_counter() - start, 2)
    return inner


def raw(func):
    func.is_raw = True
    return func


def data(func):
    func.is_data = True
    return func


def unit(unit):
    def decorator(func):
        func.unit = unit
        return func
    return decorator


def plot_label(plot_label):
    def decorator(func):
        func.plot_label = plot_label
        return func
    return decorator


def identifier(func):
    func.is_identifier = True
    return func


def info(func):
    func.is_info = True
    return func


if __name__ == '__main__':
    @timer(10)
    def test_func(n, m=1):
        s = 0
        for x in range(n):
            s += x
            if x % 2 == 3:
                s *= x // 3
        return s

    log.logfile = sys.stdout
    test_func(100000)
