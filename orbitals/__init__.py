from yutility.orbitals import sfo, mo, info, indirect
from scm import plams
import numpy as np
import matplotlib.pyplot as plt
from yutility import plot, ensure_list


class Orbitals:
    def __init__(self, path, moleculename=None):
        if isinstance(path, (plams.KFReader, plams.KFFile)):
            self.reader = path
        else:
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
        the occupation numbers of the MOs on axis 1 and the 
        '''
        return mo.occupation * self.mulliken_contribution(sfo, mo)

    @property
    def fragments(self):
        return self.sfos.fragments

    @property
    def spins(self):
        return self.sfos.spins


def load(path):
    return Orbitals(path)


def sort_orb_pairs(orbs1, orbs2, prop=None):
    '''
    Sort pairs from sfos1 and sfos2 based on the values prop(sfos1, sfos2)
    args:
        sfos1, sfos2: lists of SFO objects
        prop:         function taking SFO objects or lists of SFO objects
    return:
        list of tuples containing (sfo1, sfo2, prop(sfo1, sfo2)) sorted by prop(sfo1, sfo2)
        here sfo1 and sfo2 are taken from sfos1 and sfos2
    '''
    M = prop(orbs1, orbs2)
    ret = []
    for i, orb1 in enumerate(ensure_list(orbs1)):
        for j, orb2 in enumerate(ensure_list(orbs2)):
            if np.isnan(M[i, j]):
                continue
            ret.append((orb1, orb2, M[i, j]))

    ret = sorted(ret, key=lambda pair: pair[-1])
    return ret


def plot_property(orbs1, orbs2, prop=None, cmap='Greens', title=None, unit=None, use_relname=False, use_indexname=False, scale=None):
    if cmap is None:
        cmap = 'Greens'
        if hasattr(prop, 'cmap'):
            cmap = prop.cmap

    if scale is None:
        scale = 1
        if hasattr(prop, 'scale'):
            scale = prop.scale

    prop_name = ''
    if hasattr(prop, '__name__'):
        prop_name = prop.__name__

    if title is None:
        title = prop_name
        if hasattr(prop, 'title'):
            title = prop.title

    if unit is None:
        if hasattr(prop, 'unit'):
            unit = prop.unit

    if unit is not None:
        unit = '(' + unit + ')'
    else:
        unit = ''
 
    if callable(prop):
        M = prop(orbs1, orbs2)
    else:
        M = prop

    plotname = orbs1[0].spin + ' ' + orbs1[0].kfpath
    plt.figure(figsize=(10, 8), label=f'{prop_name} {plotname}')
    occ_virt_border1 = [i for i in range(1, len(orbs1)) if orbs1[i-1].occupation != orbs1[i].occupation]
    occ_virt_border1 = 0 if len(occ_virt_border1) == 0 else occ_virt_border1[0]
    occ_virt_border2 = [i for i in range(1, len(orbs2)) if orbs2[i-1].occupation != orbs2[i].occupation]
    occ_virt_border2 = 0 if len(occ_virt_border2) == 0 else occ_virt_border2[0]
    plt.imshow(M, origin='lower', cmap=cmap)
    # gridlines
    plt.hlines(y=np.arange(0, len(orbs1))+0.5, xmin=np.full(len(orbs1), -0.5), xmax=np.full(len(orbs1), len(orbs2)-0.5), color="w", linewidth=1.5)
    plt.vlines(x=np.arange(0, len(orbs2))+0.5, ymin=np.full(len(orbs2), -0.5), ymax=np.full(len(orbs2), len(orbs1)-0.5), color="w", linewidth=1.5)
    # occ_virt border lines
    plt.vlines(occ_virt_border2-.5, -.5, len(orbs1)-.5, colors='k', linewidth=2)
    plt.hlines(occ_virt_border1-.5, -.5, len(orbs2)-.5, colors='k', linewidth=2)
    # text inside cells
    for i in range(len(orbs1)):
        for k in range(len(orbs2)):
            val = M[i, k]
            if np.isnan(val):
                continue
            color = 'w' if val > np.nanmax(M) / 2 else 'k'
            plt.gca().text(k, i, f'{val*scale:.2f}', ha="center", va="center", color=color, fontsize=8)

    try:
        psi1 = r'\phi_{' + orbs1[0].fragment_unique_name + r'}'
        psi2 = r'\phi_{' + orbs2[0].fragment_unique_name + r'}'
    except:
        try:
            psi1 = r'\phi_{' + orbs1[0].moleculename + r'}'
            psi2 = r'\phi_{' + orbs2[0].moleculename + r'}'
        except:
            psi1 = r'\phi_{' + 'Frag1' + r'}'
            psi2 = r'\phi_{' + 'Frag2' + r'}'

    plt.xlabel('$'+psi2+'$', fontsize=16)
    plt.ylabel('$'+psi1+'$', fontsize=16)
    yticks = range(len(orbs1))
    xticks = range(len(orbs2))
    if use_relname:
        plt.xticks(xticks, [orb.relative_name for orb in orbs2], rotation=90)
        plt.yticks(yticks, [orb.relative_name for orb in orbs1], rotation=0)
    elif use_indexname:
        plt.xticks(xticks, [orb.index_name for orb in orbs2], rotation=90)
        plt.yticks(yticks, [orb.index_name for orb in orbs1], rotation=0)
    else:
        plt.xticks(xticks, [repr(orb) for orb in orbs2], rotation=90)
        plt.yticks(yticks, [repr(orb) for orb in orbs1], rotation=0)
    plt.title(title + r'$(' + psi1 + r', ' + psi2 + r')$ ' + unit, fontsize=16)
    plt.tight_layout()

    return plot.ShowCaller()


# @decorators.add_to_func(title=r'$\Delta E_{oi}$', scale=1e3, unit=r'10$^3$ eV$^{-1}$')
# def orbint(sfos1: list[SFO] or SFO, sfos2: list[SFO] or SFO, use_mask: bool = True) -> float or np.ndarray:
#     S = overlap(sfos1, sfos2)
#     dE = energy_gap(sfos1, sfos2)
#     oi = np.array(S**2/dE)
#     if use_mask:
#         mask = occ_virt_mask(sfos1, sfos2)
#         oi[np.logical_not(mask)] = None
#     return oi.squeeze()


# def highest_contr(sfos: list[sfo.SFO] or sfo.SFO, mos: list[mo.MO] or mo.MO) -> float or np.ndarray:



if __name__ == '__main__':
    p = '../test/orbitals/rkf/BH3NH3.rkf'
    orbs = Orbitals(p)
    print(orbs.fragments)

    sfos1 = orbs.sfos[:'Donor(LUMO+4)']
    sfos2 = orbs.sfos[:'Acceptor(LUMO+4)']

    plot_property(sfos1, sfos2, sfo.orbint, use_relname=True).show()
    best_pair = sort_orb_pairs(sfos1, sfos2, sfo.orbint)[-1]
    best_pair[1].generate_orbital().show()


    p = '../test/orbitals/rkf/substrate_cat_complex.rkf'
    orbs = Orbitals(p)
    print(orbs.fragments)

    sfos = orbs.sfos['C:1(1P)']
    for sfo in sfos:
        sfo.generate_orbital().show()
    mos = orbs.mos['HOMO-10':'LUMO+10']

    plot_property(sfos1, sfos2, sfo.orbint, use_relname=True).show()
    pairs = sort_orb_pairs(sfos1, sfos2, sfo.orbint)
    print(pairs[-1])

