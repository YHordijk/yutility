ensure_list = lambda x: [x] if not isinstance(x, (list, tuple, set)) else x
squeeze_list = lambda x: x[0] if len(x) == 1 else x


def argsort(lst, key=None):
    if key is None:
        key = lambda x: x

    enumerated_key = lambda x: key(x[1])
    return [x[0] for x in sorted(enumerate(lst), key=enumerated_key)]


def remove(lst, x):
    x = ensure_list(x)
    for x_ in x:
        try:
            lst.remove(x_)
        except:
            pass


def move(lst, x, newidx):
    remove(lst, x)
    lst.insert(newidx, x)


def get_first_truthy(lst):
    # returns the first element from left that evaluates to True
    # it therefore ignores None, False, '', [], {}, etc.
    if not any(lst):
        return

    for item in lst:
        if item:
            return item


def split(lst, val):
    # splits a list based on value, similar to str.split
    idx = lst.index(val)
    return lst[:idx], lst[idx+1:]


def indices(lst, val):
    idxs = []
    for i, item in enumerate(lst):
        if item == val:
            idxs.append(i)
    return idxs


def replace(lst, old, new):
    idx = lst.index(old)
    lst[idx] = new
    return lst



if __name__ == '__main__':
    print(split([0, 1, 2, 3, 4, 5, 6], 3))

    print(indices([0, 1, 2, 1, 1, 4, 5], 1))
