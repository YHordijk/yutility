def as_correct_type(s :str):
    return parse_str(s)


def parse_str(s: str):
    # checks if string should be an int, float, bool or string
    if s == '':
        return None
    if not isinstance(s, str):
        return s
    try:
        return int(s)
    except ValueError:
        pass
    try:
        return float(s)
    except:
        pass
    if s in ['True', 'False']:
        return bool(s)
    return s


def format_string(s: str):
    return eval('f' + '"' + s + '"')


if __name__ == '__main__':
    class Test:
        def __init__(self):
            self.x = 1.423

        def __repr__(self):
            return format_string('Test({self.x})')

    t = Test()
    print(repr(t))
