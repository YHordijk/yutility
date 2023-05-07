from yutility.orbitals import sfo, mo, info
from scm import plams
import numpy as np


class Orbitals:
    def __init__(self, path, moleculename=None):
        self.reader = plams.KFReader(path)
        self.mos = mo.MOs(reader=self.reader, moleculename=moleculename)
        self.sfos = sfo.SFOs(reader=self.reader)
        self.rename_fragments = self.sfos.rename_fragments

    def mulliken_contribution(self, sfo, mo):
        if sfo.symmetry != mo.symmetry:
            return 0

        if sfo.spin != mo.spin:
            return 0

        n_i = sfo.occupation
        c_ie = np.array(mo.coeffs)
        c_iu = mo @ sfo
        S_ue = np.array(sfo.overlaps)

        return n_i * np.sum(c_ie*c_iu*S_ue)


def load(path):
    return Orbitals(path)


if __name__ == '__main__':
    from yutility import ensure_list

    p = '../test/orbitals/rkf/methyl_unrestricted_frags_symm.rkf'
    orbs = Orbitals(p)
    mos = ensure_list(orbs.mos['2AA1_B'])
    sfos = orbs.sfos.sfos

    for mo_ in mos:
        for sfo_ in sfos:
            if orbs.mulliken_contribution(sfo_, mo_) == 0:
                continue
            print(mo_, sfo_, f'{orbs.mulliken_contribution(sfo_, mo_):.3%}')
