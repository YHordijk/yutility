from scm import plams
import numpy as np
import os
import matplotlib.pyplot as plt
from yutility import log, units, ensure_list, print_kf, run, plot, symmetry, decorators, volume
from yutility.ytypes import Either

j = os.path.join


def get_calc_info(reader):
    '''
    Function to read useful info about orbitals from kf reader
    '''
    ret = {}

    # determine if calculation used relativistic corrections
    # if it did, variable 'escale' will be present in 'SFOs'
    # if it didnt, only variable 'energy' will be present
    ret['relativistic'] = ('SFOs', 'escale') in reader

    # determine if SFOs are unrestricted or not
    ret['unrestricted_sfos'] = ('SFOs', 'energy_B') in reader

    # get the symmetry labels
    ret['symlabels'] = reader.read('Symmetry', 'symlab').strip().split()

    # determine if MOs are unrestricted or not
    ret['unrestricted_mos'] = (ret['symlabels'][0], 'eps_B') in reader

    # determine if the calculation used regions or not
    natoms = reader.read('Geometry', 'nr of atoms')
    nfrags = reader.read('Geometry', 'nr of fragments')
    ret['used_regions'] = natoms != nfrags

    # get symmetry labels

    return ret


def read_MO_data(reader):  # noqa: N802
    calc_info = get_calc_info(reader)

    ret = {}
    # symlabels
    ret['symlabels'] = reader.read('Symmetry', 'symlab').strip().split()

    # number of MOs
    ret['nmo'] = {}
    for symlabel in ret['symlabels']:
        if calc_info['unrestricted_mos']:
            ret['nmo'][symlabel] = {
                'A': reader.read(symlabel, 'nmo_A'),
                'B': reader.read(symlabel, 'nmo_B')
            }
        else:
            ret['nmo']['symlabel'] = reader.read(symlabel, 'nmo_A')

    # MO energies
    ret['energy'] = {}
    energyprefix = 'escale' if calc_info['relativistic'] else 'eps'
    for symlabel in ret['symlabels']:
        if calc_info['unrestricted_mos']:
            ret['energy'][symlabel] = {
                'A': reader.read(symlabel, f'{energyprefix}_A'),
                'B': reader.read(symlabel, f'{energyprefix}_B')
            }
        else:
            ret['energy']['symlabel'] = reader.read(symlabel, f'{energyprefix}_A')

    # MO occupations
    ret['occs'] = {}
    for symlabel in ret['symlabels']:
        if calc_info['unrestricted_mos']:
            ret['occs'][symlabel] = {
                'A': reader.read(symlabel, 'froc_A'),
                'B': reader.read(symlabel, 'froc_B')
            }
        else:
            ret['occs']['symlabel'] = reader.read(symlabel, 'froc_A')
 
    # MO coefficients
    ret['coeffs'] = {}
    for symlabel in ret['symlabels']:
        if calc_info['unrestricted_mos']:
            ret['coeffs'][symlabel] = {
                'A': reader.read(symlabel, 'Eig-CoreSFO_A'),
                'B': reader.read(symlabel, 'Eig-CoreSFO_B')
            }
        else:
            ret['coeffs']['symlabel'] = reader.read(symlabel, 'Eig-CoreSFO_A')

    return ret


def read_SFO_data(reader):  # noqa: N802
    def square_overlaps(S, nmo):
        Srows = []
        for i in range(nmo):
            # start index will be the number of elements before this row
            minidx1 = i * (i+1) // 2
            # stop index will be the number of elements of the next row
            maxidx1 = (i+1) * (i+2) // 2
            Srows.append(S[minidx1:maxidx1])
        # then we go through rows again and add the remaining terms
        Srowsfixed = []
        for i, row in enumerate(Srows):
            Srowsfixed.append(row + [row2[i] for row2 in Srows[i+1:]])
        return Srowsfixed

    calc_info = get_calc_info(reader)

    ret = {}
    # symlabels
    ret['symlabels'] = reader.read('Symmetry', 'symlab').strip().split()
    ret['symlabel_by_sfo'] = [] 

    # number of SFOs per symlabel
    ret['nsfo'] = {}
    norb = ensure_list(reader.read('Symmetry', 'norb'))
    for i, symlabel in enumerate(ret['symlabels']):
        ret['nsfo'][symlabel] = norb[i]
        ret['symlabel_by_sfo'].extend(norb[i] * [symlabel])
    ret['nsfo']['total'] = reader.read('SFOs', 'number')

    # isfo, index of sfo in symlabel
    ret['isfo'] = reader.read('SFOs', 'isfo')

    # fragment, index of fragment sfo belongs to
    ret['fragidx'] = reader.read('SFOs', 'fragment')

    # fragmenttype, name of fragment sfo belongs to
    ret['fragtypes'] = reader.read('SFOs', 'fragtype').strip().split()

    # fragorb, index of sfo in fragment
    ret['fragorb'] = reader.read('SFOs', 'fragorb')

    # subspecies sfo belongs to
    ret['subspecies'] = reader.read('SFOs', 'subspecies').strip().split()

    # ifo, order of sfo
    ret['ifo'] = reader.read('SFOs', 'ifo')

    # if the calculation did not use regions we have to
    # fix the ordering of the atoms (i.e. fragment indices)
    if not calc_info['used_regions']:
        atom_order_index = reader.read('Geometry', 'atom order index')
        atom_order_index = atom_order_index[:len(atom_order_index)//2]
        ret['fragidx'] = [atom_order_index.index(i) + 1 for i in ret['fragidx']]

    # SFO energies
    # and get index if energies are sorted
    energyprefix = 'escale' if calc_info['relativistic'] else 'energy'
    if calc_info['unrestricted_sfos']:
        ret['energy'] = {
            'A': reader.read('SFOs', f'{energyprefix}'),
            'B': reader.read('SFOs', f'{energyprefix}_B')
        }
        ret['energyidx'] = {
            'A': np.argsort(ret['energy']['A']),
            'B': np.argsort(ret['energy']['B']),
        }
    else:
        ret['energy'] = {'AB': reader.read('SFOs', f'{energyprefix}')}
        ret['energyidx'] = {'AB': np.argsort(ret['energy']['AB'])}

    # SFO occupations
    if calc_info['unrestricted_sfos']:
        ret['occupations'] = {
            'A': np.array(reader.read('SFOs', 'occupation')),
            'B': np.array(reader.read('SFOs', 'occupation_B'))
        }
    else:
        ret['occupations'] = {'AB': np.array(reader.read('SFOs', 'occupation'))}

    # get indices relative to HOMO for each fragment
    def get_rel_indices(spin):
        relindices = []
        newidx = 0
        for i in range(len(ret['ifo'])-1):
            occ = ret['occupations'][spin][i]
            occ2 = ret['occupations'][spin][i+1]
            if occ > 0 and occ2 == 0:
                homoidx = i
            if (occ == 0 and occ2 > 0) or i == len(ret['ifo'])-2:
                relindices.extend(np.arange(newidx - homoidx, i - homoidx + 2))
                newidx = i + 1
        return np.array(relindices)

    if calc_info['unrestricted_sfos']:
        ret['relindices'] = {
            'A': get_rel_indices('A'),
            'B': get_rel_indices('B')
        }
    else:
        ret['relindices'] = {'AB': get_rel_indices('AB')}

    # overlaps
    ret['overlaps'] = {}
    for symlabel in ret['symlabels']:
        nsfo = ret['nsfo'][symlabel]
        if calc_info['unrestricted_sfos']:
            S_A = ensure_list(reader.read(symlabel, 'S-CoreSFO'))
            S_B = ensure_list(reader.read(symlabel, 'S-CoreSFO_B'))
            ret['overlaps'][symlabel] = {
                'A': square_overlaps(S_A, nsfo),
                'B': square_overlaps(S_B, nsfo)
            }
        else:
            S = ensure_list(reader.read(symlabel, 'S-CoreSFO'))
            ret['overlaps'][symlabel] = {'AB': square_overlaps(S, nsfo)}

    return ret


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
            if orbname is not None and orbname not in [sfo.name]:
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

        if isinstance(key, str):
            decoded = self._decode_key(key)
            ret = self.get_sfo(**decoded)
            return ret[0] if len(ret) == 1 else ret

        if isinstance(key, slice):
            start_decoded = self._decode_key(key.start)
            stop_decoded = self._decode_key(key.stop)

            start_sfo = ensure_list(self.get_sfo(**start_decoded))
            stop_sfo = ensure_list(self.get_sfo(**stop_decoded))

            start_index = min(sfo.index for sfo in start_sfo)
            stop_index = min(sfo.index for sfo in stop_sfo)

            ret = []
            for index in range(start_index, stop_index + 1):
                ret.extend(self.get_sfo(index=index))
            return ret

    def get_calc_info(self):
        calc_info = get_calc_info(self.reader)
        self.is_unrestricted = calc_info['unrestricted_sfos']
        self.uses_molecular_fragments = calc_info['used_regions']
        self.uses_atomic_fragments = not calc_info['used_regions']
        self.is_relativistic = calc_info['relativistic']

        self.spins = ['A', 'B'] if self.is_unrestricted else ['AB']
        self.symlabels = calc_info['symlabels']

    def get_sfos(self):
        data = read_SFO_data(self.reader)

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
        return self.get_overlap(*args, **kwargs)

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

        fragidx_part = ''
        if self.atomic_fragments:
            fragidx_part = f':{self.fragment_index}'

        return f'{self.fragment}{fragidx_part}({self.name}){spin_part}'

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

        fragidx_part = ''
        if self.atomic_fragments:
            fragidx_part = f':{self.fragment_index}'

        return f'{self.fragment}{fragidx_part}({self.relname}){spin_part}'

    @property
    def occupied(self):
        return self.occupation > 0

    def generate_orbital(self):
        return run.orbital_cub(self.kfpath, self.index, orbtype='SFO', symlabel=self.symmetry)

    def show(self):
        self.generate_orbital().show()


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
def orbint(sfos1: list[SFO] or SFO, sfos2: list[SFO] or SFO) -> float or np.ndarray:
    S = overlap(sfos1, sfos2)
    dE = energy_gap(sfos1, sfos2)
    oi = np.array(S**2/dE)
    mask = occ_virt_mask(sfos1, sfos2)
    oi[np.logical_not(mask)] = None
    return oi.squeeze()


@decorators.add_to_func(title=r'$\Delta E_{Pauli}$', cmap='Reds')
def pauli(sfos1: list[SFO] or SFO, sfos2: list[SFO] or SFO) -> float or np.ndarray:
    S = overlap(sfos1, sfos2)
    pauli = np.array(S**2)
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


def plot_sfos_prop(sfos1, sfos2, prop=orbint, cmap=None, use_relname=False, use_indexname=False):
    cmap = 'Greens'
    if hasattr(prop, 'cmap'):
        cmap = prop.cmap

    title = prop.__name__
    if hasattr(prop, 'title'):
        title = prop.title

    M = prop(sfos1, sfos2)

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

    psi1 = r'\Psi_{' + sfos2[0].fragment + r'}'
    psi2 = r'\Psi_{' + sfos1[0].fragment + r'}'
    plt.xlabel('$'+psi1+'$')
    plt.ylabel('$'+psi2+'$')
    yticks = range(len(sfos1))
    xticks = range(len(sfos2))
    if use_relname:
        plt.xticks(xticks, [orb.relname for orb in sfos2], rotation=90)
        plt.yticks(yticks, [orb.relname for orb in sfos1], rotation=0)
    elif use_indexname:
        plt.xticks(xticks, [orb.index_name for orb in sfos2], rotation=90)
        plt.yticks(yticks, [orb.index_name for orb in sfos1], rotation=0)
    else:
        plt.xticks(xticks, [repr(orb) for orb in sfos2], rotation=90)
        plt.yticks(yticks, [repr(orb) for orb in sfos1], rotation=0)
    plt.title(title + r'$(' + psi1 + r', ' + psi2 + r')$')
    plt.tight_layout()

    return plot.ShowCaller()


if __name__ == '__main__':
    p = 'test/orbitals/rkf/methyl_unrestricted_frags_symm.rkf'
    p = 'test/orbitals/rkf/methyl.rkf'
    p = 'test/orbitals/rkf/BH3NH3.rkf'
    sfos = SFOs(kfpath=p)

    # for sfo in sfos.sfos:
    #     print(sfo)

    # sfo1 = sfos[61]
    # print(sfos['C(2S)'].index_name)
    # print([sfo.index_name for sfo in sfos['H(1S)']])
    # print(sfos['H(1S)'] @ sfos['C(2S)'])

    # print(overlap(sfos['H(1S)'], sfos['C(2S)']))

    sfos_donor = sfos['Donor(1A)':'Donor(10A)']
    sfos_acceptor = sfos['Acceptor(1A)':'Acceptor(10A)']
    plot_sfos_prop(sfos_donor, sfos_acceptor, orbint, use_relname=True).show()
    sfo_donor_best, sfo_acceptor_best, oi = sort_sfo_pairs(sfos_donor, sfos_acceptor, orbint)[-1]
    print(sfo_donor_best, sfo_acceptor_best)
    volume.show_multiple([sfo_donor_best.generate_orbital(), sfo_acceptor_best.generate_orbital()])
    # print(sfo2.full_name)

    # reader = plams.KFReader(p)
    # # print(read_MO_data(reader))
    # log.log(read_SFO_data(reader))

    # p = '/Users/yumanhordijk/PhD/yutility/test/orbitals/rkf/methyl_unrestricted_frags_symm.rkf'
    # reader = plams.KFReader(p)
    # # print(read_MO_data(reader))
    # print(read_SFO_data(reader))   
