from yutility.ytypes import check_hints, List, Either, Any


Key = Either(str, List(str))


@check_hints
def get_subkeys(headkey: Key, all_keys: List(Key)):
    if headkey == '*':
        return all_keys

    headkey = parse_key(headkey)
    subkeys = []
    for key in all_keys:
        key = parse_key(key)
        for keypart, headkeypart in zip(key, headkey):
            if keypart != headkeypart:
                break
        else:
            subkeys.append(key_to_string(key))
    return subkeys


@check_hints
def parse_key(key: Key) -> List(str):
    '''
    Parses a key from specific format.

    Keys are given as key1:key2:...:keyn or as an iterable [key1, key2, ..., keyn].
    This function ensures it is given in list format.
    '''
    if isinstance(key, str):
        return key.split(':')
    else:
        return key


@check_hints
def key_to_string(key: Key) -> str:
    '''
    Returns `key` in string format with key parts separated by `:`.
    '''
    key = parse_key(key)
    return ':'.join(key)    


@check_hints
def keys_to_dict(keys: List(Key), values=None) -> dict:
    if values is None:
        values = [None] * len(keys)
    else:
        assert len(values) == len(keys)

    d = {}
    for key, val in zip(keys, values):
        key = parse_key(key)
        d_ = d
        for keypart in key[:-1]:
            d_ = d_.setdefault(keypart, {})
        d_[key[-1]] = val
    return d


@check_hints
def add_to_dict(key, dictionary, value):
    d_ = dictionary
    for keypart in parse_key(key):
        key = parse_key(key)
        for keypart in key[:-1]:
            d_ = d_.setdefault(keypart, {})
        d_[key[-1]] = value


@check_hints
def get_from_dict(key: Key, dictionary: dict) -> Any:
    for keypart in parse_key(key):
        if keypart not in dictionary:
            raise KeyError(f'Could not find key {key}')
        dictionary = dictionary.get(keypart)
    return dictionary


if __name__ == '__main__':
    from pprint import pprint
    keys = ['test:test1:test2', 'test:a:b', 'test:a:c']
    values = [1, 2, 3]
    d = keys_to_dict(keys, values)
    add_to_dict('test:a:a', d, 'HLLO')
    val = get_from_dict('test:a:a', d)
    print(val)
