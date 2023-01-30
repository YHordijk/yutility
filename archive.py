import tarfile
import os
from ychem.utility import paths

j = os.path.join


def archive_dir(path, arcname=None, ftype='gz'):
    assert ftype in ('gz', 'bz2', 'xz')
    with tarfile.open(arcname, f'w|{ftype}') as tar:
        tar.add(path, os.path.split(path)[1])


if __name__ == '__main__':
    archive_dir(paths.database, arcname=paths.database + '.tar.gz')
