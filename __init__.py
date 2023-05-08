from yutility import log, pathfunc, listfunc

ensure_list = listfunc.ensure_list
squeeze_list = listfunc.squeeze_list


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


def make_gif(outf, frames, fps=14):
    import moviepy.editor as mvp
    clip = mvp.ImageSequenceClip(frames, fps=fps)
    clip.write_gif(outf, fps=fps)


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

    header = ['Section', 'Size', 'Unit', 'Usage']
    if print_variables:
        header = ['Section/Variable', 'Size', 'Unit', 'Usage']

    total_size = pathfunc.get_size(path)
    total_size_sum = 0
    for section in kf.sections():
        size = sum(sys.getsizeof(x) for x in kf.read_section(section).values())
        total_size_sum += size
        total, un = unit.convert(size, 'B', use_si=True)
        lines.append((section, round(total, 1), un, f'{size/total_size*100:.2f}%'))
        if print_variables:
            for var, val in kf.read_section(section).items():
                varsize = sys.getsizeof(val)
                size_, un = unit.convert(varsize, 'B', use_si=True)
                lines.append(('    ' + var, round(size_, 1), un, f'{varsize/size*100:.2f}%'))
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
