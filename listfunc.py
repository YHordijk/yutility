ensure_list = lambda x: [x] if not isinstance(x, (list, tuple, set)) else x


def get_first_truthy(lst):
    if not any(lst):
        return

    for item in lst:
        if item:
            return item
