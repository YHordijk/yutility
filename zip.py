import tarfile
from ychem.utility import paths
import os

j = os.path.join


def zip_sps():
    with tarfile.open(j(paths.sps_zip), 'w:gz') as tar:
        tar.add(paths.sps, arcname=os.path.basename(paths.sps), recursive=True)


def unzip_sps():
    with tarfile.open(paths.sps_zip, 'r:gz') as tar:
        print(os.path.dirname(paths.sps), paths.sps)
        tar.extractall(path=os.path.dirname(paths.sps))


# zip_sps()
unzip_sps()
