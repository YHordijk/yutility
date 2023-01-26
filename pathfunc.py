import os
from yutility import log


def get_server_path(path):
    server = j('/scistor', 'tc', 'yhk800', 'PhD', 'ychem')
    rel = os.path.relpath(path, base)
    return j(server, rel).replace('\\', r'/')


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
