from scm import plams
import numpy as np
import os
import matplotlib.pyplot as plt
from yutility import log, units, ensure_list, squeeze_list, print_kf, run, plot, symmetry, decorators, orbitals
from yutility.ytypes import Either

j = os.path.join


class SFOs:
    def __init__(self, kfpath=None, reader=None):
        assert reader or kfpath, 'Please provide a KFReader or path to a kf-file'
        self.reader = reader or plams.KFReader(kfpath)
        self.kfpath = kfpath
        if not self.kfpath:
            self.kfpath = self.reader.path
        self.get_calc_info()
        self.get_sfos()

    def _decode_key(self, key):
        '''
        Keys are given in the following format:

            {fragname}[:{fragidx}]({orbname})[_{spin}]

        Where [:fragidx] is optional
        '''

        # get spin from the key
        if '_' in key:
            orbname, spin = key.split('_')
        else:
            spin = None
            orbname = key

        # split key into fragment name and orbname 
        if '(' in orbname:
            fragname, orbname = orbname.split('(')
            orbname = orbname.strip(')')
        else:
            fragname = None

        # extract fragment index from fragment name if present
        if fragname is not None and ':' in fragname:
            fragname, fragidx = fragname.split(':')
            fragidx = int(fragidx)
        else:
            fragidx = None

        return {'fragname': fragname, 'fragidx': fragidx, 'orbname': orbname, 'spin': spin}

    def get_sfo(self, orbname=None, fragname=None, fragidx=None, spin=None, index=None):
        ret = []
        for sfo in self.sfos:
            if orbname is not None and orbname not in [sfo.name, sfo.relname]:
                if ':' in sfo.name and sfo.name.startswith(orbname):
                    pass
                else:
                    continue

            if fragname is not None and fragname != sfo.fragment:
                continue

            if spin is not None and spin != sfo.spin:
                continue

            if fragidx is not None and fragidx != sfo.fragment_index:
                continue

            if index is not None and index != sfo.index:
                continue

            ret.append(sfo)

        return ret

    def __getitem__(self, key):
        if isinstance(key, int):
            return self.sfos[key - 1]

        if isinstance(key, (tuple, list)):
            ret = []
            for key_ in key:
                ret.extend(ensure_list(self.__getitem__(key_)))    
            return squeeze_list(ret)

        if isinstance(key, str):
            if key in self.fragments:
                return squeeze_list(self.get_fragment_sfos(key))

            decoded = self._decode_key(key)
            ret = self.get_sfo(**decoded)
            return squeeze_list(ret)

        if isinstance(key, slice):
            start_sfo = ensure_list(self.__getitem__(key.start or 1))
            stop_sfo = ensure_list(self.__getitem__(key.stop or 0))

            start_index = min(sfo.index for sfo in start_sfo)
            stop_index = max(sfo.index for sfo in stop_sfo)

            ret = []
            for index in range(start_index, stop_index + 1):
                ret.extend(self.get_sfo(index=index))
            return squeeze_list(ret)

    def get_fragment_sfos(self, fragname):
        ret = [sfo for sfo in self.sfos if sfo.fragment == fragname]
        if len(ret) == 0:
            ret = [sfo for sfo in self.sfos if sfo.fragment.startswith(fragname)]
        return ret

    @property
    def fragments(self):
        return set([sfo.fragment for sfo in self.sfos])

    def get_calc_info(self):
        calc_info = orbitals.info.get_calc_info(self.reader)
        self.is_unrestricted = calc_info['unrestricted_sfos']
        self.uses_molecular_fragments = calc_info['used_regions']
        self.uses_atomic_fragments = not calc_info['used_regions']
        self.is_relativistic = calc_info['relativistic']

        self.spins = ['A', 'B'] if self.is_unrestricted else ['AB']
        self.symlabels = calc_info['symlabels']

    def get_sfos(self):
        data = orbitals.info.read_SFO_data(self.reader)

        sfo_data = []
        for idx in range(data['nsfo']['total']):
            for spin in self.spins:
                symlabel = data['symlabel_by_sfo'][idx]
                isfo = data['isfo'][idx] - 1
                subspecies = data['subspecies'][idx]
                ifo = data['ifo'][idx]
                relindex = data['relindices'][spin][idx]
                if relindex > 0:
                    if relindex == 1:
                        relname = 'LUMO'
                    else:
                        relname = f'LUMO+{relindex-1}'
                else:
                    if relindex == 0:
                        relname = 'HOMO'
                    else:
                        relname = f'HOMO-{abs(relindex)}'

                sfo_data.append({
                    'index':                idx + 1,
                    'relindex':             relindex,
                    'name':                 f'{ifo}{subspecies}',
                    'relname':              relname,
                    'fragment_index':       data['fragidx'][idx],
                    'fragment':             data['fragtypes'][idx],
                    'fragment_unique_name': data['fraguniquenames'][idx],
                    'fragment_orb_index':   data['fragorb'][idx],
                    'symmetry_type_index':  isfo,
                    'symmetry':             symlabel,
                    'energy':               data['energy'][spin][idx] * 27.21139664,
                    'spin':                 spin,
                    'reader':               self.reader,
                    'kfpath':               self.kfpath,
                    'overlaps':             data['overlaps'][symlabel][spin][isfo],
                    'occupation':           data['occupations'][spin][idx],
                    'atomic_fragments':     self.uses_atomic_fragments,
                })
        self.sfos = [SFO(**sfo_datum) for sfo_datum in sfo_data]


class SFO:
    '''
    lightweight class holding SFO data
    '''

    def __init__(self, **kwargs):
        for key, value in kwargs.items():
            setattr(self, key, value)

    def __repr__(self):
        return self.full_name

    def __matmul__(self, *args, **kwargs):
        return self.get_overlap(*args, **kwargs)

    def __rmatmul__(self, *args, **kwargs):
        return self.__matmul__(*args, **kwargs)

    def __sub__(self, other):
        if isinstance(other, list):
            return [self.__sub__(sfo) for sfo in other]
        return abs(self.energy - other.energy)

    def __rsub__(self, *args, **kwargs):
        return self.__sub__(*args, **kwargs)

    def get_overlap(self, other: 'SFO') -> float:
        '''
        Return overlap between two SFO objects
        SFOs must have same symmetry and spin, else overlap is 0
        '''
        if isinstance(other, list):
            return [self.get_overlap(sfo) for sfo in other]
        else:
            if self.symmetry != other.symmetry:
                return 0
            if self.spin != other.spin:
                return 0

            return self.overlaps[other.symmetry_type_index]

    @property
    def full_name(self):
        spin_part = ''
        if self.spin != 'AB':
            spin_part = f'_{self.spin}'

        return f'{self.fragment_unique_name}({self.name}){spin_part}'

    @property
    def index_name(self):
        spin_part = ''
        if self.spin != 'AB':
            spin_part = f'_{self.spin}'

        return f'{self.index}{self.symmetry}{spin_part}' 

    @property
    def relative_name(self):
        spin_part = ''
        if self.spin != 'AB':
            spin_part = f'_{self.spin}'

        return f'{self.fragment_unique_name}({self.relname}){spin_part}'

    @property
    def occupied(self):
        return self.occupation > 0

    def generate_orbital(self):
        return run.orbital_cub(self.rkf_path, self.index, orbtype='SFO', symlabel=self.symmetry)


def occ_virt_mask(sfos1: list[SFO] or SFO, sfos2: list[SFO] or SFO) -> float or np.ndarray:
    ret = []
    for sfo1 in ensure_list(sfos1):
        ret.append([])
        for sfo2 in ensure_list(sfos2):
            ret[-1].append(sfo1.occupied != sfo2.occupied)
    return np.array(ret).squeeze()


def occ_occ_mask(sfos1: list[SFO] or SFO, sfos2: list[SFO] or SFO) -> float or np.ndarray:
    ret = []
    for sfo1 in ensure_list(sfos1):
        ret.append([])
        for sfo2 in ensure_list(sfos2):
            ret[-1].append(sfo1.occupied and sfo2.occupied)
    return np.array(ret).squeeze()


@decorators.add_to_func(title=r'$S$')
def overlap(sfos1: list[SFO] or SFO, sfos2: list[SFO] or SFO) -> float or np.ndarray:
    ret = []
    for sfo1 in ensure_list(sfos1):
        ret.append([])
        for sfo2 in ensure_list(sfos2):
            ret[-1].append(abs(sfo1 @ sfo2))
    return np.array(ret).squeeze()


@decorators.add_to_func(title=r'$\epsilon$')
def energy_gap(sfos1: list[SFO] or SFO, sfos2: list[SFO] or SFO) -> float or np.ndarray:
    ret = []
    for sfo1 in ensure_list(sfos1):
        ret.append([])
        for sfo2 in ensure_list(sfos2):
            ret[-1].append(abs(sfo1.energy - sfo2.energy))
    return np.array(ret).squeeze()


@decorators.add_to_func(title=r'$\Delta E_{oi}$')
def orbint(sfos1: list[SFO] or SFO, sfos2: list[SFO] or SFO, use_mask: bool = True) -> float or np.ndarray:
    S = overlap(sfos1, sfos2)
    dE = energy_gap(sfos1, sfos2)
    oi = np.array(S**2/dE)
    if use_mask:
        mask = occ_virt_mask(sfos1, sfos2)
        oi[np.logical_not(mask)] = None
    return oi.squeeze()


@decorators.add_to_func(title=r'$\Delta E_{Pauli}$', cmap='Reds')
def pauli(sfos1: list[SFO] or SFO, sfos2: list[SFO] or SFO, use_mask: bool = True) -> float or np.ndarray:
    S = overlap(sfos1, sfos2)
    pauli = np.array(S**2)
    if use_mask:
        mask = occ_occ_mask(sfos1, sfos2)
        pauli[np.logical_not(mask)] = None
    return pauli.squeeze()


def sort_sfo_pairs(sfos1, sfos2, prop=orbint):
    '''
    Sort pairs from sfos1 and sfos2 based on the values prop(sfos1, sfos2)
    args:
        sfos1, sfos2: lists of SFO objects
        prop:         function taking SFO objects or lists of SFO objects
    return:
        list of tuples containing (sfo1, sfo2, prop(sfo1, sfo2)) sorted by prop(sfo1, sfo2)
        here sfo1 and sfo2 are taken from sfos1 and sfos2
    '''
    M = prop(sfos1, sfos2)
    ret = []
    for i, sfo1 in enumerate(ensure_list(sfos1)):
        for j, sfo2 in enumerate(ensure_list(sfos2)):
            if np.isnan(M[i, j]):
                continue
            ret.append((sfo1, sfo2, M[i, j]))

    ret = sorted(ret, key=lambda pair: pair[-1])
    return ret


def plot_sfos_prop(sfos1, sfos2, prop=orbint, cmap='Greens', title=None, use_relname=False, use_indexname=False):
    if cmap is None:
        cmap = 'Greens'
        if hasattr(prop, 'cmap'):
            cmap = prop.cmap

    if title is None:
        title = prop.__name__
        if hasattr(prop, 'title'):
            title = prop.title

    if callable(prop):
        M = prop(sfos1, sfos2)
    else:
        M = prop

    plotname = sfos1[0].spin + ' ' + sfos1[0].kfpath
    plt.figure(figsize=(10, 8), num=f'{prop.__name__} {plotname}')
    occ_virt_border1 = [i for i in range(1, len(sfos1)) if sfos1[i-1].occupation != sfos1[i].occupation]
    occ_virt_border1 = 0 if len(occ_virt_border1) == 0 else occ_virt_border1[0]
    occ_virt_border2 = [i for i in range(1, len(sfos2)) if sfos2[i-1].occupation != sfos2[i].occupation]
    occ_virt_border2 = 0 if len(occ_virt_border2) == 0 else occ_virt_border2[0]
    plt.imshow(M, origin='lower', cmap=cmap)
    # gridlines
    plt.hlines(y=np.arange(0, len(sfos1))+0.5, xmin=np.full(len(sfos1), -0.5), xmax=np.full(len(sfos1), len(sfos2)-0.5), color="w", linewidth=1.5)
    plt.vlines(x=np.arange(0, len(sfos2))+0.5, ymin=np.full(len(sfos2), -0.5), ymax=np.full(len(sfos2), len(sfos1)-0.5), color="w", linewidth=1.5)
    # occ_virt border lines
    plt.vlines(occ_virt_border2-.5, -.5, len(sfos1)-.5, colors='k', linewidth=2)
    plt.hlines(occ_virt_border1-.5, -.5, len(sfos2)-.5, colors='k', linewidth=2)
    # text inside cells
    for i in range(len(sfos1)):
        for k in range(len(sfos2)):
            val = M[i, k]
            if np.isnan(val):
                continue
            color = 'w' if val > np.nanmax(M) / 2 else 'k'
            plt.gca().text(k, i, f'{val:.2f}', ha="center", va="center", color=color, fontsize=8)

    psi1 = r'\phi_{' + sfos1[0].fragment_unique_name + r'}'
    psi2 = r'\phi_{' + sfos2[0].fragment_unique_name + r'}'
    plt.xlabel('$'+psi2+'$', fontsize=16)
    plt.ylabel('$'+psi1+'$', fontsize=16)
    yticks = range(len(sfos1))
    xticks = range(len(sfos2))
    if use_relname:
        plt.xticks(xticks, [orb.relname for orb in sfos2], rotation=45)
        plt.yticks(yticks, [orb.relname for orb in sfos1], rotation=0)
    elif use_indexname:
        plt.xticks(xticks, [orb.index_name for orb in sfos2], rotation=45)
        plt.yticks(yticks, [orb.index_name for orb in sfos1], rotation=0)
    else:
        plt.xticks(xticks, [repr(orb) for orb in sfos2], rotation=45)
        plt.yticks(yticks, [repr(orb) for orb in sfos1], rotation=0)
    plt.title(title + r'$(' + psi1 + r', ' + psi2 + r')$', fontsize=16)
    plt.tight_layout()

    return plot.ShowCaller()


if __name__ == '__main__':
    p = '../test/orbitals/rkf/BH3NH3.rkf'
    sfos = SFOs(kfpath=p)

    sfos_donor = sfos[:'Donor(LUMO+2)']
    sfos_acceptor = sfos['Acceptor(1A)':'Acceptor(LUMO+2)']
    sfo_donor_best, sfo_acceptor_best, oi = sort_sfo_pairs(sfos_donor, sfos_acceptor, orbint)[-1]
    print(sfo_donor_best, sfo_acceptor_best, oi)
    plot_sfos_prop(sfos_donor, sfos_acceptor, orbint, use_relname=True).hold()

    p = '../test/orbitals/rkf/methyl.rkf'
    sfos = SFOs(kfpath=p)

    sfos_c = sfos['C(1S)', 'C(2S)', 'C(1P)']
    sfos_h = sfos.get_fragment_sfos('H')
    sfos_c_best, sfos_h_best, oi = sort_sfo_pairs(sfos_c, sfos_h, orbint)[-1]
    plot_sfos_prop(sfos_c, sfos_h, overlap, use_relname=False).hold()