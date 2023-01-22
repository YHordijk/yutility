import os

j = os.path.join
dn = os.path.dirname

base = dn(dn(__file__))
test_molecules = j(base, 'data', 'molecules')
tmp = j(base, 'tmp')
manager_start = j(base, 'reaction_generation', 'start.sh')
SG_molecules = j(base, 'reaction_generation', 'reactants')
SG_substituents = j(base, 'reaction_generation', 'substituents')
SG_input_mols = j(base, 'reaction_generation', 'input_mols')
SG_log = j(base, 'reaction_generation', 'RG_log.txt')
calculations = j(base, 'calculations')
calculations_test = j(base, 'calculations.bak')
viewer_data = j(base, 'viewer', 'data')
configs = j(base, 'configs')
manager_logs = j(base, 'reaction_generation', 'manager_logs')
reactions = j(base, 'results', 'reactions.txt')
gui_resources = j(base, 'gui', 'resources')
reservations = j(base, 'reaction_generation', 'reservations')
connect = j(base, 'connect')
plot = j(base, 'plot')
sp_database = j(base, 'results', 'sp_database.csv')
rxn_database = j(base, 'results', 'rxn_database.csv')
gui_active_manager_info = j(base, 'gui', 'resources', 'managers.json')
database = j(base, 'results', 'database')
rxns = j(base, 'results', 'rxns')
raw = j(base, 'results', 'raw')
sps_zip = j(base, 'results', 'sps.tar.gz')
test = j(base, 'tests')


__all__ = [base,
           test_molecules,
           tmp,
           SG_molecules,
           SG_substituents,
           SG_input_mols,
           SG_log,
           calculations,
           calculations_test,
           viewer_data,
           gui_resources,
           reservations,
           configs,
           manager_logs,
           reactions,
           connect,
           plot,
           sp_database,
           rxn_database,
           gui_active_manager_info,
           database,
           raw,
           sps_zip,
           test,
           ]


def get_server_path(path):
    server = j('/scistor', 'tc', 'yhk800', 'PhD', 'ychem')
    rel = os.path.relpath(path, base)
    return j(server, rel).replace('\\', r'/')


def check_paths(make=True):
    if not all(map(os.path.exists, [
               p for p in __all__ if len(p.split('.')) == 1])):
        print('[paths.py]: Warning, not all paths exist.')
        for p in __all__:
            if len(p.split('.')) == 1:
                if not os.path.exists(p):
                    print(p)
                    os.makedirs(p, exist_ok=True)
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


def print_paths():
    from ychem.utility import log, units
    header = [log.emojis['empty'], 'Path', 'Size']
    ls = []
    unit = units.Binary()
    for p in __all__:
        l = []
        if os.path.exists(p):
            l.append(log.emojis['good'])
        else:
            l.append(log.emojis['fail'])
        l.append(p)
        v, u = unit.convert(get_size(p), use_si=True)
        l.append(f'{v:.1f} {u}')
        ls.append(l)
    log.print_list(ls, header=header)


if __name__ == '__main__':
    # check_paths()
    print_paths()

    # print(get_path_versions(r"D:\Users\Yuman\Desktop\PhD\ychem\calculations\04129edec013cd795578878661f10004e9ccc71d7350ec3e5602a4c8a6c70670\test"))
