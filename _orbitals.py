from scm import plams
import numpy as np
import os
import matplotlib.pyplot as plt
from yutility import log, units, ensure_list, print_kf, run, plot, symmetry
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

    # determine if MOs are unrestricted or not
    symlabels = reader.read('Symmetry', 'symlab').strip().split()
    ret['unrestricted_mos'] = (symlabels[0], 'eps_B') in reader

    return ret


class Orbitals:
    def __init__(self, kfpath, moleculename=None):
        self.kfpath = kfpath
        self.moleculename = moleculename
        self.fmos = FMOs(kfpath, moleculename=moleculename)
        self.sfos = SFOs(kfpath)

    def rename_fragments(self, old, new):
        self.sfos.rename_fragments(old, new)

    @property
    def is_unrestricted(self):
        return self.sfos.sfos[0].spin != 'AB'

    @property
    def spins(self):
        return ['A', 'B'] if self.is_unrestricted else ['AB']


class FMOs:
    def __init__(self, kfpath, moleculename=None):
        self.fmos = _get_all_FMOs(kfpath)
        self.moleculename = moleculename
        for fmo in self.fmos:
            fmo.fragment = moleculename

    def __iter__(self):
        return iter(self.fmos)

    def __getitem__(self, key):
        if isinstance(key, slice):
            start = key.start or 0
            stop = key.stop or -1
            step = key.step or 1

            startidx = min([fmo.index for fmo in ensure_list(self[start])])
            stopidx = max([fmo.index for fmo in ensure_list(self[stop])]) + 1

            ret = []
            for i in range(startidx, stopidx, step):
                ret.extend(ensure_list(self[i]))

            if len(ret) == 1:
                return ret[0]
            return ret


        if isinstance(key, int):
            ret = []
            for fmo in self.fmos:
                if fmo.index == key:
                    ret.append(fmo)
            if len(ret) == 1:
                return ret[0]
            return ret


        if isinstance(key, str):
            # check if key corresponds to name of fmo or AMSlevels name
            # also check if spin is given
            key_parts = key.split()
            if len(key_parts) > 1:
                key, spin = key_parts
            else:
                spin = 'AB'

            ret = []
            for fmo in self.fmos:
                if fmo.name == key and fmo.spin in spin:
                    ret.append(fmo)
                if fmo.AMSlevels_name == key:
                    ret.append(fmo)
            if len(ret) == 1:
                return ret[0]
            return ret

    @property
    def is_unrestricted(self):
        return self.fmos[0].spin != 'AB'

    def mulliken_contributions(self, sfos):
        ...


class FMO:
    def __init__(self, index, offset, energy, name, coeffs, occupation, rkf, rkf_path=None, symmetry=None, spin=None, overlaps=None, fragment=None):
        self.index = index
        self.offset = offset
        self.energy = energy
        self.name = name
        self.coeffs = coeffs
        self.occupation = occupation
        self.fragment = fragment
        self.rkf = rkf
        self.rkf_path = rkf_path
        self.symmetry = symmetry
        self.spin = spin
        self.overlaps = overlaps

    @property
    def coeff_norm(self):
        return 1 / sum(c**2 for c in self.coeffs)

    @property
    def occupied(self):
        return self.occupation > 0

    @property
    def AMSlevels_name(self):
        if self.spin != 'AB':
            return f'{self.index}{self.symmetry}_{self.spin}'
        return f'{self.index}{self.symmetry}'

    @property
    def relname(self):
        return self.__repr__()

    @property
    def AMSlevels_relname(self):
        return self.relname

    def __repr__(self):
        if self.spin != 'AB':
            return self.name + ' ' + self.spin
        return self.name

    def __getitem__(self, key):
        if isinstance(key, SFO):
            return self.coeffs[key.index]

    def get_coeff(self, sfo):
        assert self.rkf_path == sfo.rkf_path
        return self.coeffs[sfo.index]

    def generate_orbital(self):
        return run.orbital_cub(self.rkf_path, self.index, orbtype='SCF', symlabel=self.symmetry)


def _get_all_FMOs(kfpath):
    reader = plams.KFReader(kfpath)
    def get_FMOs_of_spin_and_symmetry(symmetry, spin):
        nmo = len(reader.read(symlabel, f'eps_{spin}'))
        eps = reader.read(symlabel, f'eps_{spin}')
        occs = np.array(reader.read(symlabel, f'froc_{spin}'))
        coeffs = reader.read(symlabel, f'Eig-CoreSFO_{spin}')

        # make coeffs square
        coeffs = np.array(coeffs).reshape(nmo, nmo)
        fmos = []
        for index in range(nmo):
            # set the name and offset of the orbital
            offset = index - homo_idx
            if index == somo_idx:
                orbname = 'SOMO'
            elif offset > 0:
                if somo_idx:
                    offset -= 1
                orbname = f'LUMO+{offset-1}' if offset > 1 else 'LUMO'
            else:
                orbname = f'HOMO-{abs(offset)}' if offset < 0 else 'HOMO'

            fmos.append(FMO(index + 1, offset, eps[index], orbname, coeffs[index], occs[index], reader, rkf_path=kfpath, symmetry=symlabel, spin=spin))
        return fmos

    if ('Symmetry', 'symlab') in reader:
        symmlabels = reader.read('Symmetry', 'symlab').strip().split()
    else:
        symmlabels = symmetry.labels[reader.read('Geometry', 'grouplabel').strip()]

    fmos = []
    for symlabel in symmlabels:
        is_unrestricted = (symlabel, 'eps_B') in reader
        somo_idx = None
        if is_unrestricted:
            occs = np.array(reader.read(symlabel, 'froc_A'))
            homo_idx_A = [i for i, occ in enumerate(occs) if occ > 0][-1]

            occs = np.array(reader.read(symlabel, 'froc_B'))
            homo_idx_B = [i for i, occ in enumerate(occs) if occ > 0][-1]
            homo_idx = min(homo_idx_A, homo_idx_B)
            if homo_idx_B != homo_idx_A:
                somo_idx = max(homo_idx_A, homo_idx_B)

            fmos.extend(get_FMOs_of_spin_and_symmetry(symlabel, 'A'))
            fmos.extend(get_FMOs_of_spin_and_symmetry(symlabel, 'B'))
        else:
            occs = np.array(reader.read(symlabel, 'froc_A'))
            occupied = [i for i, occ in enumerate(occs) if occ > 0]
            homo_idx = 0 if len(occupied) == 0 else occupied[-1]
            fmos_ = get_FMOs_of_spin_and_symmetry(symlabel, 'A')
            for fmo in fmos_:
                fmo.spin = 'AB'
            fmos.extend(fmos_)
            
    return fmos


class SFOs:
    def __init__(self, kfpath):
        self.sfos = _get_all_SFOs(kfpath)

    def __iter__(self):
        return iter(self.sfos)

    def get_orbital(self, index, fragidx=None, spin=None):
        ret = []
        for sfo in self.sfos:
            # print(index, fragidx, sfo.index, sfo.fragmentindex)
            if sfo.index != index:
                continue
            if fragidx is not None and sfo.fragmentindex != fragidx:
                continue
            if spin is not None and sfo.spin != spin:
                continue
            ret.append(sfo)
            # print(index, fragidx, sfo.index, sfo.fragmentindex)
        return ret

    def _decode_key(self, key):
        '''
        Keys are given in the following format:

            {frag}[:{fragidx}]({orbname}{symmlabel}_{spin})

        Where :fragidx is optional
        '''
        if isinstance(key, int):
            key = str(key)

        key_splits = key.split('(')
        if len(key_splits) > 1:
            frag, orbname = key_splits
            orbname = orbname.strip(')')
        else:
            frag = None
            orbname = key_splits[0]

        fragidx = None
        if frag:
            frag_splits = frag.split(':')
            if len(frag_splits) > 1:
                frag, fragidx = frag_splits
                fragidx = int(fragidx)

        orbname_splits = orbname.split('_')

        if len(orbname_splits) > 1:
            orbname, spin = orbname_splits
        else:
            spin = None

        ''' 
        get symmetry label from orbname
            orbname: 'HOMOA',     symmlabel: 'A'
            orbname: 'LUMO+4  B1.g', symmlabel: 'B1.g'
            orbname: '29AA',       symmlabel: 'AA'
        '''
        symmlabel = None
        on = orbname.replace('HOMO', '').replace('LUMO', '').replace('+', '').replace('-', '')
        for i, char in enumerate(on):
            if not char.isnumeric():
                symmlabel = on[i:]
                if ':' in symmlabel:
                    symmlabel = None
                break

        # print(frag, fragidx, orbname, symmlabel, spin)
        # assert frag is None or frag in self.fragments, f'Fragment {frag} is not present, must be one of [{" ".join(self.fragments)}]'
        # assert spin is None or spin in self.spins, f'Spin state {spin} is not present for {"un"*self.is_unrestricted}restricted orbitals'
        # assert symmlabel is None or symmlabel in self.symmetry_labels, f'Symmetry species {symmlabel} is not present, must be one of [{" ".join(self.symmetry_labels)}]'
        return frag, fragidx, orbname, symmlabel, spin

    def get_sfo(self, frag=None, fragidx=None, orbname=None, symmlabel=None, spin=None, index=None):
        ret = []
        for sfo in self.sfos:
            # print(sfo)
            if frag is not None and sfo.fragment != frag:
                continue

            if fragidx is not None and sfo.fragmentindex != fragidx:
                continue

            if spin is not None and sfo.spin != spin:
                continue

            if symmlabel is not None and sfo.symmetry != symmlabel:
                continue

            if index is not None and sfo.index == index:
                ret.append(sfo)
                continue

            if orbname is not None and orbname == sfo.name:
                ret.append(sfo)
                continue

            if orbname is not None and orbname == sfo._relname:
                ret.append(sfo)
                continue

            if orbname is not None and orbname == sfo.AMSlevels_indexname:
                ret.append(sfo)
                continue

            if orbname is not None and sfo.name.startswith(orbname):
                ret.append(sfo)
                continue

        if len(ret) == 1:
            return ret[0]
        return ret

    def __getitem__(self, key):
        if isinstance(key, str):
            args = self._decode_key(key)
            return self.get_sfo(*args)
        elif isinstance(key, slice):
            startargs = self._decode_key(key.start)
            stopargs = self._decode_key(key.stop)

            start_sfo = ensure_list(self.get_sfo(*startargs))
            stop_sfo = ensure_list(self.get_sfo(*stopargs))

            start_idx = min([sfo.index for sfo in start_sfo])
            stop_idx = max([sfo.index for sfo in stop_sfo]) + 1

            frag, fragidx, _, symmlabel, spin = startargs

            ret = []
            for idx in range(start_idx, stop_idx):
                ret.extend(ensure_list(self.get_sfo(frag, fragidx, symmlabel=symmlabel, spin=spin, index=idx)))
            return ret

    @property
    def fragments(self):
        return set([sfo.fragment for sfo in self.sfos])

    @property
    def spins(self):
        return set([sfo.spin for sfo in self.sfos])

    @property
    def is_unrestricted(self):
        return 'AB' not in self.spins

    @property
    def symmetry_labels(self):
        return list(set([sfo.symmetry for sfo in self.sfos]))

    @property
    def from_fragment_analysis(self):
        return self.sfos[0].from_fragment_analysis

    def rename_fragments(self, old, new):
        for sfo in self:
            if sfo.fragment not in old:
                continue
            idx = old.index(sfo.fragment)
            sfo.fragment = new[idx]


class SFO:
    def __init__(self, 
                 index, 
                 fragorbindex, 
                 name, 
                 relname=None,
                 coeffs=None, 
                 fragment=None, 
                 fragmentindex=None, 
                 rkf=None, 
                 rkf_path=None, 
                 energy=None, 
                 occupation=None, 
                 overlaps=None,
                 spin=None,
                 symmetry=None,
                 from_fragment_analysis=None):
        self.index = index
        self.fragorbindex = fragorbindex
        self.name = name
        self._relname = relname
        self.coeffs = coeffs
        self.fragment = fragment
        self.fragmentindex = fragmentindex
        self.rkf = rkf
        self.rkf_path = rkf_path
        self.energy = energy
        self.occupation = occupation
        self.overlaps = overlaps
        self.spin = spin
        self.symmetry = symmetry
        self.from_fragment_analysis = from_fragment_analysis

    def __repr__(self):
        if self.spin != 'AB':
            return f'{self.fragment}({self.name}_{self.spin})'
        return f'{self.fragment}({self.name})'

    def __getitem__(self, key):
        if isinstance(key, FMO):
            return self.coeffs[key.index-1]

    @property
    def AMSlevels_name(self):
        if self.spin != 'AB':
            return f'{self.name}_{self.spin}'
        return f'{self.name}'

    @property
    def AMSlevels_indexname(self):
        return f'{self.index}{self.symmetry}'

    @property
    def AMSlevels_relname(self):
        return self._relname

    @property
    def index_name(self):
        if self.spin != 'AB':
            return f'{self.index}{self.symmetry}_{self.spin}'
        return f'{self.index}{self.symmetry}'

    @property
    def relname(self):
        if self.spin != 'AB':
            return f'{self.fragment}({self._relname}_{self.spin})'
        return f'{self.fragment}({self._relname})'

    def coeff(self, fmo):
        assert self.rkf_path == fmo.rkf_path
        return self.coeffs[fmo.index-1]

    def overlap(self, sfo):
        assert self.rkf_path == sfo.rkf_path
        if self.symmetry != sfo.symmetry:
            return 0
        if self.spin != sfo.spin:
            return 0
        return self.overlaps[sfo.index-1]

    def __matmul__(self, other):
        if isinstance(other, FMO):
            return self.coeff(other)
        if isinstance(other, SFO):
            return self.overlap(other)

    def __rmatmul__(self, other):
        return self.__matmul__(other)

    def generate_orbital(self):
        return run.orbital_cub(self.rkf_path, self.index, orbtype='SFO', symlabel=self.symmetry)


def _get_all_SFOs(kfpath):
    reader = plams.KFReader(kfpath)
    nsfos       = reader.read('SFOs', 'number')
    fragtypes   = reader.read('SFOs', 'fragtype').split()
    fragidx     = np.array(reader.read('SFOs', 'fragment'))

    from_fragment_analysis = True
    natoms = reader.read('Geometry', 'nr of atoms')
    if max(fragidx) == natoms:
        atom_order_index = reader.read('Geometry', 'atom order index')
        atom_order_index = atom_order_index[:len(atom_order_index)//2]
        fragidx     = [atom_order_index.index(i) + 1 for i in fragidx]
        from_fragment_analysis = False

    fragnames    = [name for name, idx in zip(fragtypes, fragidx)]
    if ('SFOs', 'escale') in reader:
        energies = reader.read('SFOs', 'escale')
    else:
        energies = reader.read('SFOs', 'energy')
        if ('SFOs', 'energy_B') in reader:
            energies_B = reader.read('SFOs', 'energy_B')

    occupations  = reader.read('SFOs', 'occupation')
    ifo          = reader.read('SFOs', 'ifo')
    subspecies   = reader.read('SFOs', 'subspecies').split()
    sfonames     = [f'{ifo_}{subsp}' for ifo_, subsp in zip(ifo, subspecies)]
    indices      = reader.read('SFOs', 'fragorb')

    if ('Symmetry', 'symlab') in reader:
        symmlabels = reader.read('Symmetry', 'symlab').strip().split()
        symmnorb = ensure_list(reader.read('Symmetry', 'norb'))
    else:
        symmlabels = symmetry.labels[reader.read('Geometry', 'grouplabel').strip()]
        symmnorb = [reader.read(symmlabel, 'nmo') for symmlabel in symmlabels]

    sfo_symmlabel = [symmlabel for i, symmlabel in enumerate(symmlabels) for _ in range(symmnorb[i])]
    coeffs = {}
    S = {}
    for symmlabel in symmlabels:
        is_unrestricted = (symmlabel, 'eps_B') in reader
        coeffs[symmlabel] = {}
        for spin in ['A', 'B'] if is_unrestricted else ['A']:
            nmo = len(reader.read(symmlabel, f'eps_{spin}'))
            coeffs_ = reader.read(symmlabel, f'Eig-CoreSFO_{spin}')
            # make coeffs square
            coeffs[symmlabel][spin if is_unrestricted else "AB"] = np.array(coeffs_).reshape(nmo, nmo)

        '''
        read overlap, overlap is saved as follows:
            start       S       Compressed
             idx  
              0    | 1        |    | 1 |
              1    | 2 3      |    | 2 |
              3    | 4 5 6    | -> | 3 |
              6    | 7 8 9 10 |    | 4 |
                                   | 5 |
                                   |...|
        '''
        S[symmlabel] = {}
        for spin in ['A', 'B'] if is_unrestricted else ['A']:
            if spin == 'A':
                overlaps = reader.read(symmlabel, 'S-CoreSFO')
            else:
                try:
                    overlaps = reader.read(symmlabel, 'S-CoreSFO_B')
                except:
                    overlaps = reader.read(symmlabel, 'S-CoreSFO')
            # read each row
            ov = []
            for i in range(nmo):
                # start index will be the number of elements before this row
                minidx1 = i * (i+1) // 2
                # stop index will be the number of elements of the next row
                maxidx1 = (i+1) * (i+2) // 2
                ov.append(overlaps[minidx1:maxidx1])
            # then we go through rows again and add the remaining terms
            ov2 = []
            for i, row in enumerate(ov):
                ov2.append(row + [row2[i] for row2 in ov[i+1:]])
            S[symmlabel][spin if is_unrestricted else "AB"] = np.array(ov2)

    sfos = []
    for i in range(nsfos):
        for spin in ['A', 'B'] if is_unrestricted else ['AB']:
            symm = sfo_symmlabel[i]
            symmidx = indices[i]

            sfo = SFO(i+1, 
                      indices[i], 
                      sfonames[i],
                      coeffs=coeffs[symm][spin][:, i],
                      fragment=fragnames[i],
                      fragmentindex=fragidx[i],
                      rkf=reader,
                      rkf_path=kfpath,
                      energy=energies[i],
                      occupation=occupations[i], 
                      symmetry=symm,
                      spin=spin,
                      overlaps=S[symm][spin][:, i],
                      from_fragment_analysis=from_fragment_analysis)
            sfos.append(sfo)

    hofo_indices = {idx: [sfo.index for sfo in sfos if sfo.fragmentindex == idx and sfo.occupation > 0][-1] for idx in set(fragidx)}
    for sfo in sfos:
        rel = sfo.index - hofo_indices[sfo.fragmentindex]
        if rel > 1:
            sfo._relname = f'LUMO+{rel-1}'
        elif rel == 1:
            sfo._relname = 'LUMO'
        elif rel == 0:
            sfo._relname = 'HOMO'
        elif rel < 0:
            sfo._relname = f'HOMO{rel}'

    return sfos



def S(sfos1, sfos2):
    Sm = np.zeros((len(sfos1), len(sfos2)))
    for i, sfo1 in enumerate(sfos1):
        for k, sfo2 in enumerate(sfos2):
            Sm[i, k] = abs(sfo1 @ sfo2)
    return Sm


def dE(sfos1, sfos2):
    dEm = np.zeros((len(sfos1), len(sfos2)))
    for i, sfo1 in enumerate(sfos1):
        for k, sfo2 in enumerate(sfos2):
            dEm[i, k] = abs(sfo1.energy - sfo2.energy) * 27.2114079527
    return dEm


def occ_and(sfos1, sfos2):
    occm = np.zeros((len(sfos1), len(sfos2)))
    for i, sfo1 in enumerate(sfos1):
        for k, sfo2 in enumerate(sfos2):
            occm[i, k] = sfo1.occupation and sfo2.occupation
    return occm.astype(bool)

def occ_or(sfos1, sfos2):
    occm = np.zeros((len(sfos1), len(sfos2)))
    for i, sfo1 in enumerate(sfos1):
        for k, sfo2 in enumerate(sfos2):
            occm[i, k] = sfo1.occupation or sfo2.occupation
    return occm.astype(bool)


def orbint(sfos1, sfos2):
    orbintm = S(sfos1, sfos2)**2 / abs(dE(sfos1, sfos2))
    mask = occ_or(sfos1, sfos2) ^ occ_and(sfos1, sfos2)
    orbintm[mask == 0] = None
    return orbintm


def pauli(sfos1, sfos2):
    paulim = S(sfos1, sfos2)**2
    mask = occ_and(sfos1, sfos2)
    paulim[mask == 0] = None
    return paulim


def orbint_interactions(sfos1, sfos2):
    orbintm = orbint(sfos1, sfos2)
    ret = []
    for i, sfo1 in enumerate(sfos1):
        for k, sfo2 in enumerate(sfos2):
            if np.isnan(orbintm[i, k]):
                continue
            ret.append((sfo1, sfo2, orbintm[i, k]))
    # set highest interaction first
    ret = sorted(ret, key=lambda x: -x[2])
    return ret


def plot_sfos_prop(sfos1, sfos2, prop=S, cmap=None, use_relname=False, use_indexname=False):
    # assert prop in [orbint, pauli, S, dE]

    if cmap is None:
        if prop is pauli:
            cmap = 'Reds'
        else:
            cmap = 'Greens'

    M = prop(sfos1, sfos2)

    plotname = sfos1[0].spin + ' ' + sfos1[0].rkf_path
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
        plt.xticks(xticks, [orb.AMSlevels_relname for orb in sfos2], rotation=90)
        plt.yticks(yticks, [orb.AMSlevels_relname for orb in sfos1], rotation=0)
    if use_indexname:
        plt.xticks(xticks, [orb.index_name for orb in sfos2], rotation=90)
        plt.yticks(yticks, [orb.index_name for orb in sfos1], rotation=0)
    else:
        plt.xticks(xticks, [repr(orb) for orb in sfos2], rotation=90)
        plt.yticks(yticks, [repr(orb) for orb in sfos1], rotation=0)
    plt.title(prop.__name__ + r'$(' + psi1 + r', ' + psi2 + r')$')
    plt.tight_layout()

    return plot.ShowCaller()


def plot_pauli(sfos1, sfos2):
    pauli = pauli(sfos1, sfos2) * 1000

    plotname = sfos1[0].rkf_path
    plt.figure(figsize=(7, 7), num='Pauli ' + plotname)
    occ_virt_border1 = [i for i in range(1, len(sfos1)) if sfos1[i-1].occupation != sfos1[i].occupation]
    occ_virt_border1 = 0 if len(occ_virt_border1) == 0 else occ_virt_border1[0]
    occ_virt_border2 = [i for i in range(1, len(sfos2)) if sfos2[i-1].occupation != sfos2[i].occupation]
    occ_virt_border2 = 0 if len(occ_virt_border2) == 0 else occ_virt_border2[0]
    plt.imshow(pauli, origin='lower', cmap='Reds')
    # gridlines
    plt.hlines(y=np.arange(0, len(sfos1))+0.5, xmin=np.full(len(sfos1), -0.5), xmax=np.full(len(sfos1), len(sfos2)-0.5), color="w", linewidth=1.5)
    plt.vlines(x=np.arange(0, len(sfos2))+0.5, ymin=np.full(len(sfos2), -0.5), ymax=np.full(len(sfos2), len(sfos1)-0.5), color="w", linewidth=1.5)
    # occ_virt border lines
    plt.vlines(occ_virt_border2-.5, -.5, len(sfos1)-.5, colors='k', linewidth=2)
    plt.hlines(occ_virt_border1-.5, -.5, len(sfos2)-.5, colors='k', linewidth=2)
    # text inside cells
    for i in range(len(sfos1)):
        for k in range(len(sfos2)):
            val = pauli[i, k]
            if np.isnan(val):
                continue
            color = 'w' if val > np.nanmax(pauli) / 2 else 'k'
            plt.gca().text(k, i, f'{val:.1f}', ha="center", va="center", color=color, fontsize=8)

    psi1 = r'\Psi_{' + sfos2[0].fragment + r'}'
    psi2 = r'\Psi_{' + sfos1[0].fragment + r'}'
    plt.xlabel('$'+psi1+'$')
    plt.ylabel('$'+psi2+'$')
    yticks = range(len(sfos1))
    xticks = range(len(sfos2))
    plt.xticks(xticks, [orb.name for orb in sfos2], rotation=90)
    plt.yticks(yticks, [orb.name for orb in sfos1], rotation=0)
    plt.title(r'Pauli-Repulsions  $\propto$  $\langle' + psi1 + ' | ' + psi2 + r' \rangle^2$   $(10^3)$')
    plt.tight_layout()

    return plot.ShowCaller()


def plot_orbint(sfos1, sfos2):
    orbintM = orbint(sfos1, sfos2) * 1000

    plotname = sfos1[0].rkf_path
    plt.figure(figsize=(7, 7), num='OrbInt ' + plotname)
    occ_virt_border1 = [i for i in range(1, len(sfos1)) if sfos1[i-1].occupation != sfos1[i].occupation]
    occ_virt_border1 = 0 if len(occ_virt_border1) == 0 else occ_virt_border1[0]
    occ_virt_border2 = [i for i in range(1, len(sfos2)) if sfos2[i-1].occupation != sfos2[i].occupation]
    occ_virt_border2 = 0 if len(occ_virt_border2) == 0 else occ_virt_border2[0]
    plt.imshow(orbintM, origin='lower', cmap='Greens')
    # gridlines
    plt.hlines(y=np.arange(0, len(sfos1))+0.5, xmin=np.full(len(sfos1), -0.5), xmax=np.full(len(sfos1), len(sfos2)-0.5), color="w", linewidth=1.5)
    plt.vlines(x=np.arange(0, len(sfos2))+0.5, ymin=np.full(len(sfos2), -0.5), ymax=np.full(len(sfos2), len(sfos1)-0.5), color="w", linewidth=1.5)
    # occ_virt border lines
    plt.vlines(occ_virt_border2-.5, -.5, len(sfos1)-.5, colors='k', linewidth=2)
    plt.hlines(occ_virt_border1-.5, -.5, len(sfos2)-.5, colors='k', linewidth=2)
    # text inside cells
    for i in range(len(sfos1)):
        for k in range(len(sfos2)):
            val = orbintM[i, k]
            if np.isnan(val):
                continue
            color = 'w' if val > np.nanmax(orbintM) / 2 else 'k'
            plt.gca().text(k, i, f'{val:.1f}', ha="center", va="center", color=color, fontsize=8)

    psi1 = r'\Psi_{' + sfos2[0].fragment + r'}'
    psi2 = r'\Psi_{' + sfos1[0].fragment + r'}'
    eps1 = r'\epsilon_{' + sfos2[0].fragment + r'}'
    eps2 = r'\epsilon_{' + sfos1[0].fragment + r'}'
    plt.xlabel('$'+psi1+'$')
    plt.ylabel('$'+psi2+'$')
    yticks = range(len(sfos1))
    xticks = range(len(sfos2))
    plt.xticks(xticks, [orb.name for orb in sfos2], rotation=90)
    plt.yticks(yticks, [orb.name for orb in sfos1], rotation=0)
    plt.title(r'Orbital-Interactions  $\propto$  $\frac{\langle' + psi1 + ' | ' + psi2 + r' \rangle^2}{|' + eps1 + ' - ' + eps2 + r'|}$   $(10^3/eV)$')
    plt.tight_layout()

    return plot.ShowCaller()


if __name__ == '__main__':
    p = '/Users/yumanhordijk/PhD/LewisAcid_coordination/calculations_final/EDA_vacuum/I2_N_pi/full/adf.rkf'
    orbs = Orbitals(p)
    orbs.rename_fragments(['subsystem1', 'subsystem2'], ['LA', 'Substrate'])
    print('PI')
    print(orbs.sfos['LA(54A)'].relname, orbs.sfos['Substrate(14A)'].relname)
    print((orbs.sfos['LA(54A)'].energy - orbs.sfos['Substrate(14A)'].energy) * 27.2114079527)
    print(orbs.sfos['LA(54A)'] @ orbs.sfos['Substrate(14A)'])
    print()