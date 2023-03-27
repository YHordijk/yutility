from yutility import log

ensure_list = lambda x: [x] if not isinstance(x, (list, tuple, set)) else x


def download_from_github(file, out, repo=None, token=None):
    '''
    Download a file in path ``file`` in repository with name ``repo`` into path ``out``
    '''
    from github import Github
    import requests

    if repo is None:
        repo, file = file.split(':')

    git = Github(token)
    user = git.get_user()
    repo = user.get_repo(repo)
    content = repo.get_contents(file)
    url = content.download_url
    r = requests.get(url, allow_redirects=True)
    with open(out, 'wb') as outf:
        outf.write(r.content)


def dowhile(main_func, condition_func):
    main_func()
    while condition_func():
        main_func()


def print_kf(path, print_variables=False):
    from yutility import units
    from scm import plams
    import sys
    kf = plams.KFFile(path)

    lines = []
    unit = units.Binary('B')

    header = ['Section', 'Size', 'Unit']
    if print_variables:
        header = ['Section/Variable', 'Size', 'Unit']

    for section in kf.sections():
        total, un = unit.convert(sum(sys.getsizeof(x) for x in kf.read_section(section).values()), 'B', use_si=True)
        lines.append((section, round(total, 1), un))
        if print_variables:
            for var, val in kf.read_section(section).items():
                size, un = unit.convert(sys.getsizeof(val), 'B', use_si=True)
                lines.append(('    ' + var, round(size, 1), un))
    log.print_list(lines, header=header)


if __name__ == '__main__':
    x = 10


    def f():
        global x
        x *= 2


    def cond():
        global x
        return x*2 <= 20480*2*2


    dowhile(f, cond)
    print(x)
