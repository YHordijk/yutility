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
