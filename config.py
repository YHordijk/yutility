from yutility import paths

_config = {}
with open(paths.config) as conf:
    lines = [line.strip() for line in conf.readlines()]

for line in lines:
    _config[line.split('=')[0]] = line.split('=')[1]


def get(key):
    return _config.get(key)
