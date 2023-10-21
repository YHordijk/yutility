from time import perf_counter
from functools import wraps
import numpy as np
from yutility import log, listfunc
'''
Implements a decorator that records and stores the number of times a function
has been called and also the amount of time taken in total/per call
'''
times = {}
exec_start = perf_counter()

enabled: bool = True


class Timer:
    def __init__(self, name=None):
        self.name = name
        times.setdefault(self.name, {'calls': 0, 'timings': []})

    def __enter__(self):
        self.start = perf_counter()
        return self

    def __exit__(self, *args, **kwargs):
        times[self.name]['calls'] += 1
        times[self.name]['timings'].append(perf_counter() - self.start)


def Time(func):
    times.setdefault(func.__qualname__, {'calls': 0, 'timings': []})

    @wraps(func)
    def inner(*args, **kwargs):
        if enabled and __debug__:
            start = perf_counter()
            ret = func(*args, **kwargs)
            times[func.__qualname__]['calls'] += 1
            times[func.__qualname__]['timings'].append(perf_counter() - start)
            return ret
        else:
            return func(*args, **kwargs)
    return inner


def print_timings():
    funcs = [func if not hasattr(func, '__qualname__') else func.__qualname__ for func in times.keys()] + ['TOTAL']
    ncalls = [str(time['calls']) for time in times.values()]
    ncalls = ncalls + [str(sum(int(n) for n in ncalls))]
    means = ['None' if time['calls'] == 0 else f"{np.mean(time['timings']): .3f}" for time in times.values()] + ['']
    stds = ['None' if time['calls'] == 0 else f"{np.std(time['timings']): .3f}" for time in times.values()] + ['']

    total_time = [sum(time['timings']) for time in times.values()]
    total_total_time = perf_counter() - exec_start
    rel_time = [f'{tt: .3f} ({tt/total_total_time:.1%})' for tt in total_time] + [f'{total_total_time: .3f} (100%)']

    header = ['Function', 'Calls', 'Mean (s)', 'Std. dev. (s)', 'Time Spent (s)']
    data = [x for x in zip(funcs, ncalls, means, stds, rel_time)]
    log.print_list(data, header=header, hline=[-1])


def print_timings2():
    names = list(times.keys())
    names = sorted(names)
    names_splits = [name.split('.') for name in names]
    ps = []
    for i, name in enumerate(names_splits):
        parents = []
        for name_ in names_splits[:i][::-1]:
            parents.append(".".join(listfunc.common_list(name, name_)))
        parents = [parent for parent in parents if parent]
        if parents != []:
            longest_parent = sorted(parents, key=lambda x: len(x))[-1]
            if longest_parent not in names:
                names.append(longest_parent)
                ps.append(longest_parent)
    names = sorted(names)
    names_splits = [name.split('.') for name in names]
    # get the parents:
    parents = []
    parent_levels = []
    for i, (name, splits) in enumerate(zip(names, names_splits)):
        parent = 'TOTAL'
        parent_level = 0
        for name_, splits_ in zip(names[:i], names_splits[:i]):
            if listfunc.startswith(splits, splits_):
                parent = name_
                parent_level += 1
        parents.append(parent)
        parent_levels.append(parent_level)

    parent_times = {parent: {'timings': [sum([sum(times.get(name, {'timings': [0]})['timings']) for name in names if name.startswith(parent)])]} for parent in parents}
    for parent in parents:
        if parent in times:
            parent_times[parent] = {'timings': times[parent]['timings'], 'calls': times[parent]['calls']}
    parent_times['TOTAL'] = {'timings': [perf_counter() - exec_start], 'calls': 1}

    header = ['Function', 'Calls', 'Mean (s)', 'Time Spent (s)', 'Rel. Time']
    times_ = parent_times.copy()
    times_.update(times)
    lines = []
    for name, parent, level in zip(names, parents, parent_levels):
        line = []
        line.append(name.replace(parent, ' > ' * level))  # function name
        mean = np.mean(times_[name]["timings"])
        total = np.sum(times_[name]["timings"])
        calls = times_[name].get("calls", 0)
        if calls == 0:
            line.append('')
            line.append('')
        else:
            line.append(str(calls))
            line.append(f'{mean:.3f}')

        line.append(f'{total:.3f}')

        rel = sum(times_[name]["timings"])/sum(parent_times[parent]["timings"])*100
        rel_total = sum(times_[name]["timings"])/sum(parent_times["TOTAL"]["timings"])*100
        if parent != 'TOTAL':
            line.append(f'{rel_total: >3.0f}% ({rel: >3.0f}%)')
        else:
            line.append(f'{rel_total: >3.0f}%')

        lines.append(line)
    lines.append(['TOTAL', '', '', f'{sum(times_["TOTAL"]["timings"]):.3f}', '100%'])
    
    log.print_list(lines, header=header, hline=[-2])


if __name__ == '__main__':
    from time import sleep

    @Time
    def test():
        sleep(.141 + np.random.rand()/7)

    @Time
    def test2():
        sleep(.231 + np.random.rand()/7)

    [test() for _ in range(3)]
    [test2() for _ in range(5)]

    for i in range(100):
        with Timer('inner'):
            x = np.random.randn(1000000)

    print_timings()
    print()
    print_timings2()
