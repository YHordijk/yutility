import os
from yutility import log
import re
from typing import Union
from tcutility import results

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


def match(root, pattern: Union[str, list[str]]):
    '''
    Find and return PathMatches object of root that match the given pattern.
    '''

    # get the number and names of substitutions in the given pattern
    substitutions = re.findall(r'{(\w+[+*?]?)}', pattern)
    # the pattern should resolve to words and may contain - and _
    # replace them here
    for sub in substitutions:
        quantifier = sub[-1] if sub[-1] in '+*?' else '+'
        pattern = pattern.replace('{' + sub + '}', f'([a-zA-Z0-9_-]{quantifier})')

    ret = results.Result()
    # root dir can be any level deep. We should count how many directories are in root
    root_length = len(split_all(root))
    # get all subdirectories first, we can loop through them later
    subdirs = get_subdirectories(root, include_intermediates=True)
    # remove the root from the subdirectories. We cannot use str.removeprefix because it was added in python 3.9
    subdirs = [j(*split_all(subdir)[root_length:]) for subdir in subdirs if len(split_all(subdir)[root_length:]) > 0]
    for subdir in subdirs:
        # check if we get a match with our pattern
        match = re.fullmatch(pattern, subdir)
        if not match:
            continue

        p = j(root, subdir)
        # get the group data and add it to the return dictionary. We skip the first group because it is the full directory path
        ret[p] = results.Result(directory=p, **{substitutions[i]: match.group(i+1) for i in range(len(substitutions))})

    return ret


if __name__ == '__main__':
    from tcutility import log

    systems = match('tmp', '{system}/EDA')
    systems.print()
    print(systems.systems)

    path_matches = match('tmp', '{system}/EDA/frag_{fragment}')
    path_matches.print()

    path_matches = match('tmp', '{system}/EDA/complex{suffix*}')
    path_matches.print()
