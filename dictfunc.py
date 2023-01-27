def dict_match(a, b):
    matched_keys = [key for key in a.keys() if key in b.keys()]
    match = True
    for key in matched_keys:
        aval, bval = a[key], b[key]
        if isinstance(aval, dict) and isinstance(bval, dict):
            match = dict_match(aval, bval)
        else:
            if aval != bval:
                match = False
    return match


def dict_match_strict(a, b):
    matched_keys = [key for key in a.keys() if key in b.keys()]
    match = True
    for key in matched_keys:
        aval, bval = a[key], b[key]
        if isinstance(aval, dict) and isinstance(bval, dict):
            match = dict_match(aval, bval)
        else:
            if aval != bval:
                match = False
    return match


def dict_to_list(a: dict):
    l = []
    if not isinstance(a, dict):
        return [[a]]
    for k, v in a.items():
        if isinstance(v, dict):
            if v == {}:
                l.append([k, {}])
            else:
                [l.append([k, *x]) for x in dict_to_list(v)]
        else:
            l.append([k, v])
    return l


def list_to_dict(a: list):
    d = {}
    for lst in a:
        d_ = d
        for i, key in enumerate(lst):
            if i == len(lst) - 2:
                d_[lst[-2]] = lst[-1]
                break
            else:
                d_ = d_.setdefault(key, {})
    return d


def common_dict(dicts):
    if len(dicts) == 0:
        return {}
    _d = {}
    _d.update(dicts[0])
    for d in dicts:
        for k1, v1 in d.items():
            if k1 not in _d.keys():
                _d[k1] = v1
            elif isinstance(v1, dict):
                _d[k1].update(v1)
    return _d


def get_inverse(key, dic):
    keys, values = list(dic.keys()), list(dic.values())
    idx = values.index(key)
    return keys[idx]


def invert(dic):
    return {val: key for key, val in dic.items()}


def key_from_value(value, dic):
    for key, val in dic.items():
        if val == value:
            return key


def remove_false_keys(dic):
    to_remove = []
    for key, value in dic.items():
        if not value:
            to_remove.append(key)
    for key in to_remove:
        del(dic[key])
    return dic


def merged_dict(a, b):
    a.update(b)
    return a


if __name__ == '__main__':
    d = {'name': 'transitionstate',
 'reactants': {'catalyst': 'catalyst',
               'radical': 'methyl',
               'substrate': 'acrolein_E'},
 'reaction': 'Radical Addition',
 'reaction_specific': {},
 'settings_preset': 'BLYP-D3(BJ)/TZ2P/Good',
 'solvent': 'vacuum',
 'substituents': {'catalyst': {'Rcat': 'I2'},
                  'radical': {},
                  'substrate': {'R1': 'pyrrolidine',
                                'R2': 'm-HOPh',
                                'R3': 'H',
                                'R4': 'H',
                                'Rch': 'O'}}}
    lst = dict_to_list(d)
    print(d)
    print(list_to_dict(lst))