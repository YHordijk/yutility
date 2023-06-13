import os
import sys
import json
import psutil

j = os.path.join

RESERVATION_DIR = j(os.path.split(__file__)[0], 'reservations')

USING_SLURM = 'SLURM_JOB_ID' in os.environ


def get_slurmids():
    '''
    Return a list of SLURM-IDs owned by the current user running on this server
    '''
    squeue_lines = [line for line in os.popen('squeue --me').read().split('\n') if line]
    return [int(line.split()[0]) for line in squeue_lines[1:]]


def get_pids():
    '''
    Returns all process-IDs for processes that are running Python
    '''
    return [proc.pid for proc in psutil.process_iter() if 'python' in proc.name().lower()]


def get_ids():
    if USING_SLURM:
        return get_slurmids()
    return get_pids()


def get_slurmid():
    return int(os.environ['SLURM_JOB_ID'])


def get_pid():
    '''
    Return process-ID of the current process
    '''
    return os.getpid()


def get_id():
    '''
    Get an ID for the current process. This will be a SLURM-ID if slurm is active, else it will be the PID
    '''
    if USING_SLURM:
        return get_slurmid()
    return get_pid()


def get_reservation_files():
    ''' 
    Returns a list of files in RESERVATION_DIR if they are currently in used by a process (if filename equals an active PID)
    '''
    return [j(RESERVATION_DIR, file) for file in os.listdir(RESERVATION_DIR) if file != '.DS_Store' and int(file) in get_ids()]

def get_reservations():
    ''' 
    Returns a list of files in RESERVATION_DIR if they are currently in used by a process (if filename equals an active PID)
    '''
    ret = []
    for file in get_reservation_files():
        with open(file) as r:
            content = json.loads(r.read())
        ret.append(content)
    return ret


def is_reserved(obj):
    ''' 
    Returns whether a json-dumpable object `obj` is reserved by an active ID
    '''
    return obj in get_reservations()


def make_reservation(obj):
    ''' 
    Create a file with the name of the ID of this process containing json-dumpable obj
    '''
    file = j(RESERVATION_DIR, str(get_id()))
    os.makedirs(RESERVATION_DIR, exist_ok=True)
    with open(file, 'w+') as f:
        f.write(json.dumps(obj, indent=4, sort_keys=True))


def clean_reservations():
    ''' 
    Remove all files in paths.reservations that are not used by an active process
    '''
    active_files = get_reservation_files()
    for file in os.listdir(RESERVATION_DIR):
        if j(RESERVATION_DIR, file) in active_files:
            continue
        os.remove(j(RESERVATION_DIR, file))


if __name__ == '__main__':
    # print(USING_SLURM)
    # print(RESERVATION_DIR)
    # print(os.getpid())
    # # print(get_pids())
    make_reservation({'test': 123})
    # print(get_reservation_files())
    print(is_reserved({'test': 124}))
    # # print(is_reserved_from_path(r"D:\Users\Yuman\Desktop\PhD\ychem\calculations\388662676b6f5a5e7b3c11ba90a16454d4c0672281b3bd175ddc76297fedb09e"))
    # clean_reservations()
