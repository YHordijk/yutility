import os
from yutility import log
import re

j = os.path.join


def check_paths(pathlist, make=True):
    if not all(os.path.exists(p) for p in pathlist if len(p.split('.')) == 1):
        log.log('Warning, not all paths exist:')
        log.tab_level += 1
        for p in pathlist:
            if len(p.split('.')) == 1:
                if not os.path.exists(p):
                    log.log(p)
                    os.makedirs(p, exist_ok=True)
        log.tab_level -= 1
    else:
        print('Paths OK ✔️')


def get_size(start_path):
    # return size of path in Bytes
    if os.path.isfile(start_path):
        return os.path.getsize(start_path)
    
    total_size = 0
    for dirpath, dirnames, filenames in os.walk(start_path):
        for f in filenames:
            fp = os.path.join(dirpath, f)
            # skip if it is symbolic link
            if not os.path.islink(fp):
                total_size += os.path.getsize(fp)

    return total_size


def get_path_versions(path):
    '''Returns all versions of a path
    e.g. in a tree such as 
    root
    |- dir_or_file
    |- dir_or_file.002
    |- dir_or_file.003
    |- ...
    |- dir_or_file.n

    returns all directories starting with root/dir_or_file'''

    root = os.path.dirname(path)
    ret = []
    for d in os.listdir(root):
        if j(root, d).startswith(path):
            ret.append(j(root, d))
    return ret


def next_path_version(path):
    i = 1
    while os.path.exists(path + [f'.{str(i).zfill(3)}', ''][i == 1]):
        i += 1
    return path + [f'.{str(i).zfill(3)}', ''][i == 1]


def print_paths(pathlist):
    from yutility import log, units
    header = [log.emojis['empty'], 'Path', 'Size']
    ls = []
    for p in pathlist:
        l = []
        if os.path.exists(p):
            l.append(log.emojis['good'])
        else:
            l.append(log.emojis['fail'])
        l.append(p)
        v, u = units.Binary().convert(get_size(p), use_si=True)
        l.append(f'{v:.1f} {u}')
        ls.append(l)
    log.print_list(ls, header=header)


def split_all(path):
    parts = []
    while True:
        a, b = os.path.split(path)
        if not a or not b:
            return parts[::-1]
        parts.append(b)
        path = a


def get_subdirectories(root, include_intermediates=False):
    dirs = [root]
    subdirs = set()

    while len(dirs) > 0:
        _dirs = []
        for cdir in dirs:
            csubdirs = [j(cdir, d) for d in os.listdir(cdir) if os.path.isdir(j(cdir, d))]
            if len(csubdirs) == 0:
                subdirs.add(cdir)
            else:
                if include_intermediates:
                    subdirs.add(cdir)
                _dirs.extend(csubdirs)

        dirs = _dirs

    return subdirs


class PathMatches:
    def __init__(self):
        self.directories = []
        self.data = {}
        self._keys = set()

    def add(self, directory, **data):
        self.directories.append(directory)
        self.data[directory] = data
        [self._keys.add(k) for k in data.keys()]

    def __iter__(self):
        return iter(self.data.keys())

    def __getitem__(self, key):
        if key in self._keys:
            return {d: self.data[d][key] for d in self}
        if isinstance(key, int):
            d = self.directories[key]
            ret = self.data[d].copy()
            ret['directory'] = d
            return ret

    def __getattr__(self, key):
        return self[key]

    def keys(self):
        return ['directory'] + list(self._keys)

    def items(self):
        return self.data.items()

    def rows(self):
        keys = self.keys()[1:]
        ret = []
        for directory, data in self.items():
            ret.append([directory] + [data[key] for key in keys])
        return ret

    def print(self):
        log.table(self.rows(), self.keys())

    def __len__(self):
        return len(self.directories)


def match(root, pattern):
    '''
    Find and return PathMatches object of root that match the given pattern.
    '''

    # get the number and names of substitutions in the given pattern
    substitutions = re.findall(r'{(\w+)}', pattern)
    # the pattern should resolve to words and may contain - and _
    # replace them here
    for sub in substitutions:
        pattern = pattern.replace('{' + sub + '}', '([a-zA-Z0-9_-]+)')

    ret = PathMatches()
    root_length = len(split_all(root))
    subdirs = get_subdirectories(root, include_intermediates=True)
    subdirs = [j(*split_all(subdir)[root_length:]) for subdir in subdirs if len(split_all(subdir)[root_length:]) > 0]
    for subdir in subdirs:
        match = re.fullmatch(pattern, subdir)
        if match:
            p = j(root, subdir)
            ret.add(p, **{substitutions[i]: match.group(i+1) for i in range(len(substitutions))})

    return ret


if __name__ == '__main__':
    from tcutility import log

    path_matches = match('pathfunc_test', '{system}/{functional}_{basis_set}')
    path_matches.print()

    # for p in path_matches:
    #     print(p, path_matches.basis_set[p])

    print(path_matches.directories)

    for i in range(len(path_matches)):
        print(path_matches[i])