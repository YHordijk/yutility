from time import perf_counter
from functools import wraps
import numpy as np
from yutility import log
'''
Implements a decorator that records and stores the number of times a function
has been called and also the amount of time taken in total/per call
'''
times = {}
exec_start = perf_counter()

enabled: bool = True



class Timer:
    def __init__(self, name=None):
        self.__qualname__ = name
        times.setdefault(self.__qualname__, {'calls': 0, 'timings': []})

    def __enter__(self):
        self.start = perf_counter()
        return self

    def __exit__(self, *args, **kwargs):
        times[self.__qualname__]['calls'] += 1
        times[self.__qualname__]['timings'].append(perf_counter() - self.start)


def time(func):
    times[func] = {'calls': 0, 'timings': []}

    @wraps(func)
    def inner(*args, **kwargs):
        if enabled and __debug__:
            start = perf_counter()
            ret = func(*args, **kwargs)
            times[func]['calls'] += 1
            times[func]['timings'].append(perf_counter() - start)
            return ret
        else:
            return func(*args, **kwargs)
    return inner


def print_timings():
    funcs = [func if not hasattr(func, '__qualname__') else func.__qualname__ for func in times.keys()]
    ncalls = [str(time['calls']) for time in times.values()]
    means = ['None' if time['calls'] == 0 else f"{np.mean(time['timings']): .3f}" for time in times.values()]
    stds = ['None' if time['calls'] == 0 else f"{np.std(time['timings']): .3f}" for time in times.values()]

    total_time = [sum(time['timings']) for time in times.values()]
    total_total_time = perf_counter() - exec_start
    rel_time = [f'{tt: .3f} ({tt/total_total_time:.1%})' for tt in total_time]

    header = ['Function', 'Calls', 'Mean (s)', 'Std. dev. (s)', 'Time Spent (s)']
    data = [x for x in zip(funcs, ncalls, means, stds, rel_time)]
    log.print_list(data, header=header)


if __name__ == '__main__':
    from time import sleep

    @time
    def test():
        sleep(.141 + np.random.rand()/7)

    @time
    def test2():
        sleep(.231 + np.random.rand()/7)

    [test() for _ in range(3)]
    [test2() for _ in range(5)]


    for i in range(100):
        with Timer('inner'):
            x = np.random.randn(1000000)

    print_timings()
