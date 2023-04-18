ensure_list = lambda x: [x] if not isinstance(x, (list, tuple, set)) else x


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
