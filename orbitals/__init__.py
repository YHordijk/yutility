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
        r'''
        Calculate the Mulliken contribution of a selected SFO to a selected MO.
        The Mulliken contribution originates from the Mulliken population analysis,
        but before scaling the result by the number of electrons occupying the MO.

        The Mulliken contribution of SFO $\mu$ to MO $i$ is given as

            $$\hat{C}_{i\mu} = c_{i\mu}^2 + \sum_{\nu \neq \mu} c_{i\mu} c_{i\nu} S_{\mu\nu}
                             = \sum_\nu c_{i\mu} c_{i\nu} S_{\mu\nu}$$

            Where the index $\nu$ denotes all SFOs, $c_{i\nu}$ is the coefficient 
            of SFO $\nu$ in MO $i$ and $S_{\mu\nu}$ is the overlap between SFOs
            $\mu$ and $\nu$.

        The marginals of the resulting Mulliken contribution matrix $\hat{C}_{i\mu}$ should all equal one.
        '''
        # coefficient of all SFOs contributing to the selected MO
        c_iv = np.array(mo @ self.sfos.sfos)
        # coefficient of the selected SFO
        c_iu = mo @ sfo
        # overlaps between selected SFO and all other SFOs
        S_uv = np.array(sfo @ self.sfos.sfos)
        # calculate the mulliken contribution
        return np.sum(c_iv*c_iu*S_uv)

    def mulliken_population(self, sfo, mo):
        r'''
        Calculate the Mulliken population of a selected SFO to a selected MO.
        The Mulliken population is the product of the occupation of the MO with the 
        Mulliken contribution of the SFO to the MO.

        The Mulliken population of SFO $\mu$ to MO $i$ is then
            
            $$P_{i\mu} = n_i\hat{C}_{i\mu}$$

            Where $n_i$ is the occupation of MO $i$ and $\hat{C}_{i\mu}$ is the 
            Mulliken contribution of SFO $\mu$ to MO $i$.

        The marginals of the Mulliken population matrix $P_{i\mu}$ should equal
        the occupation numbers of the MOs.
        '''
        return mo.occupation * self.mulliken_contribution(sfo, mo)


def load(path):
    return Orbitals(path)


if __name__ == '__main__':
    import matplotlib.pyplot as plt
    from yutility import ensure_list

    p = '../test/orbitals/rkf/methyl_unrestricted_frags_symm.rkf'
    orbs = Orbitals(p)

    contributions = []
    for mo_ in orbs.mos:
        contributions.append([])
        for sfo_ in orbs.sfos:
            contributions[-1].append(orbs.mulliken_contribution(sfo_, mo_))

    print(np.sum(contributions, axis=1))
    plt.imshow(contributions)
    plt.show()
