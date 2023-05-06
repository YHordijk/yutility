from yutility.orbitals import sfo, mo, info
from scm import plams


class Orbitals:
    def __init__(self, path):
        self.reader = plams.KFReader(path)
        self.mos = mo.MOs(reader=self.reader)
        self.sfos = sfo.SFOs(reader=self.reader)


def load(path):
    return Orbitals(path)


if __name__ == '__main__':
    p = '../test/orbitals/rkf/methyl.rkf'
    orbs = Orbitals(p)
    print(orbs.mos.mos)
    print(orbs.sfos.sfos)
