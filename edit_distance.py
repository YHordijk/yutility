import numpy as np


def lev(a, b, max_dist=None):
    if len(b) == 0:
        return len(a)
    if len(a) == 0:
        return len(b)

    def tail(s):
        return s[1:]

    if a[0] == b[0]:
        return lev(tail(a), tail(b))

    d = 1 + min(lev(tail(a), b), lev(a, tail(b)), lev(tail(a), tail(b)))
    return d


def WagnerFischer(a, b):
    d = np.zeros((len(a)+1, len(b)+1)).astype(int)

    for i in range(len(a)):
        d[i+1, 0] = i+1
    for i in range(len(b)):
        d[0, i+1] = i+1

    for i in range(1, len(a)+1):
        for j in range(1, len(b)+1):
            if a[i-1] == b[j-1]:
                cost = 0
            else:
                cost = 1

            d[i, j] = min(d[i-1, j] + 1, d[i, j-1] + 1, d[i-1, j-1] + cost)

    return d[-1, -1]


def get_closest(a, others, compare_func=WagnerFischer, ignore_case=False, ignore_chars='', maximum_distance=None):
    if ignore_case:
        a = a.lower()
    a = a.replace(ignore_chars, '')

    dists = []
    for other in others:
        if ignore_case:
            other = other.lower()
        other = other.replace(ignore_chars, '')
        dists.append(compare_func(a, other))

    lowest_strs = [other for dist, other in zip(dists, others) if dist <= max(maximum_distance, min(dists))]
    return lowest_strs


if __name__ == '__main__':
    from tcutility.data import functionals
    from yutility import timer
    from tcutility import log

    def get(functional_name: str):
        '''
        Return information about a given functional.

        Args:
            functional_name: the name of the functional. It should exist in the get_available_functionals keys.

        Return:
            A :class:`results.Result` object containing information about the functional if it exists. Else it will return ``None``.
        '''
        funcs = functionals.get_available_functionals()
        ret = funcs[functional_name]
        if not ret:
            funcs.prune()
            closest = get_closest(functional_name, funcs.keys(), ignore_case=True, ignore_chars='-', maximum_distance=2)
            if len(closest) > 1:
                closest_string = " or ".join([", ".join(closest[:-1]), closest[-1]])
            else:
                closest_string = closest[0]
            log.warn(f'Could not find functional "{functional_name}". Did you mean {closest_string}?')

        return ret

    get('lyp')
    get('blyp-d3(bj)')

    for _ in range(200):
        with timer.Timer('Wagner-Fischer'):
            WagnerFischer('sitting', 'kitten')
    
    for _ in range(200):
        with timer.Timer('Naïve'):
            lev('sitting', 'kitten')

    timer.print_timings2()
