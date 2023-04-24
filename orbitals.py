from scm import plams
import numpy as np
import os
import matplotlib.pyplot as plt
from yutility import log, units, ensure_list, print_kf, run, plot, symmetry
from yutility.ytypes import Either

j = os.path.join


class Orbitals:
    def __init__(self, kfpath):
        self.kfpath = kfpath
        self.fmos = FMOs(kfpath)
        self.sfos = SFOs(kfpath)

    def rename_fragments(self, old, new):
        self.sfos.rename_fragments(old, new)

    @property
    def is_unrestricted(self):
        return self.sfos.sfos[0].spin != 'AB'


class FMOs:
    def __init__(self, kfpath):
        self.fmos = _get_all_FMOs(kfpath)

    def __iter__(self):
        return iter(self.fmos)

    def __getitem__(self, key):
        if isinstance(key, slice):
            start = key.start or 0
            stop = key.stop or -1
            step = key.step or 1

            startidx = min([fmo.index for fmo in ensure_list(self[start])])
            stopidx = max([fmo.index for fmo in ensure_list(self[stop])])

            ret = []
            for i in range(startidx, stopidx, step):
                ret.extend(self[i])

            return ret

        if isinstance(key, int):
            ret = []
            for fmo in self.fmos:
                if fmo.index == key:
                    ret.append(fmo)
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
    def __init__(self, index, offset, energy, name, coeffs, occupation, rkf, rkf_path=None, symmetry=None, spin=None, overlaps=None):
        self.index = index
        self.offset = offset
        self.energy = energy
        self.name = name
        self.coeffs = coeffs
        self.occupation = occupation
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

            fmos.append(FMO(index, offset, eps[index], orbname, coeffs[index], occs[index], reader, rkf_path=kfpath, symmetry=symlabel, spin=spin))
        return fmos

    if ('Symmetry', 'symlab') in reader:
        symmlabels = reader.read('Symmetry', 'symlab').strip()
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
            homo_idx = [i for i, occ in enumerate(occs) if occ > 0][-1]
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

    def __decode_key(self, key):
        frag, orbname = key.split('(')
        frag_splits = frag.split(':')
        if len(frag_splits) > 1:
            frag, fragidx = frag_splits
            fragidx = int(fragidx)
        else:
            fragidx = None

        orbname = orbname.strip(')')
        orbname_splits = orbname.split('_')

        if len(orbname_splits) > 1:
            orbname, spin = orbname_splits
        else:
            spin = None

        return frag, fragidx, orbname, spin

    def get_sfo(self, frag, fragidx, orbname, spin):
        ret = []
        for sfo in self.sfos:
            if sfo.fragment != frag:
                continue
            if fragidx is not None and sfo.fragmentindex != fragidx:
                continue
            if spin is not None and sfo.spin != spin:
                continue

            if orbname == sfo.name:
                ret.append(sfo)

            if orbname == sfo._relname:
                ret.append(sfo)

        return ret

    def __getitem__(self, key):
        '''
        
        '''
        if isinstance(key, str):
            args = self.__decode_key(key)
            return self.get_sfo(*args)
        elif isinstance(key, slice):
            startargs = self.__decode_key(key.start)
            stopargs = self.__decode_key(key.stop)
            



        # if isinstance(key, slice):
        #     start = key.start or 0
        #     stop = key.stop or -1
        #     step = key.step or 1

        #     start_sfo = self[start]
        #     stop_sfo = self[stop]

        #     fragidx = None
        #     for fragment in self.fragments:
        #         if start.startswith(fragment):
        #             fragidx = start_sfo.fragmentindex
        #             assert stop.startswith(fragment), 'Slice start and stop keys must have same fragment'
        #             break
        #     spin = start_sfo.spin
        #     assert stop_sfo.spin == spin



        #     startidx = min([fmo.index for fmo in ensure_list(start_sfo)])
        #     stopidx = max([fmo.index for fmo in ensure_list(stop_sfo)])
        #     ret = []
        #     for idx in range(startidx, stopidx+1):
        #         ret.extend(self.get_orbital(idx, fragidx, spin))

        #     return ret

        # if isinstance(key, int):
        #     ret = []
        #     for sfo in self.sfos:
        #         if sfo.index == key:
        #             ret.append(sfo)

        # elif isinstance(key, str):
        #     # check if key corresponds to name of sfo or AMSlevels name
        #     ret = []
        #     for sfo in self.sfos:
        #         if repr(sfo) == key:
        #             ret.append(sfo)
        #         if sfo.AMSlevels_name == key:
        #             ret.append(sfo)
        #         if sfo.relname == key:
        #             ret.append(sfo)
        #         if sfo.AMSlevels_relname == key:
        #             ret.append(sfo)
        # if len(ret) == 1:
        #     return ret[0]
        # return ret

    @property
    def fragments(self):
        return set([sfo.fragment for sfo in self.sfos])

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
                 symmetry=None):
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
    def AMSlevels_relname(self):
        return self._relname

    @property
    def relname(self):
        if self.spin != 'AB':
            return f'{self.fragment}({self._relname}_{self.spin})'
        return f'{self.fragment}({self._relname})'

    def coeff(self, fmo):
        assert self.rkf_path == fmo.rkf_path
        return self.coeffs[fmo.spin][fmo.index-1]

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
    fragidx     = reader.read('SFOs', 'fragment')
    fragnames   = [name for name, idx in zip(fragtypes, fragidx)]
    energies    = reader.read('SFOs', 'escale')
    occupations = reader.read('SFOs', 'occupation')
    ifo         = reader.read('SFOs', 'ifo')
    subspecies  = reader.read('SFOs', 'subspecies').split()
    sfonames    = [f'{ifo_}{subsp}' for ifo_, subsp in zip(ifo, subspecies)]
    indices     = reader.read('SFOs', 'fragorb')
    
    if ('Symmetry', 'symlab') in reader:
        symmlabels = reader.read('Symmetry', 'symlab').strip().split()
        symmnorb = ensure_list(reader.read('Symmetry', 'norb'))
    else:
        symmlabels = symmetry.labels[reader.read('Geometry', 'grouplabel').strip()]
        symmnorb = [reader.read(symlabel, 'nmo') for symlabel in symmlabels]

    sfo_symmlabel = [symlabel for i, symlabel in enumerate(symmlabels) for _ in range(symmnorb[i])]
    coeffs = {}
    S = {}
    for symlabel in symmlabels:
        is_unrestricted = (symlabel, 'eps_B') in reader
        coeffs[symlabel] = {}
        for spin in ['A', 'B'] if is_unrestricted else ['A']:
            nmo = len(reader.read(symlabel, f'eps_{spin}'))
            coeffs_ = reader.read(symlabel, f'Eig-CoreSFO_{spin}')
            # make coeffs square
            coeffs[symlabel][spin if is_unrestricted else "AB"] = np.array(coeffs_).reshape(nmo, nmo)

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
        S[symlabel] = {}
        for spin in ['A', 'B'] if is_unrestricted else ['A']:
            if spin == 'A':
                overlaps = reader.read(symlabel, 'S-CoreSFO')
            else:
                overlaps = reader.read(symlabel, 'S-CoreSFO_B')
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
            S[symlabel][spin if is_unrestricted else "AB"] = np.array(ov2)

        # plt.imshow(np.abs(S[symlabel]))
        # plt.show()
        # plt.imshow(S[symlabel])
        # plt.show()

    sfos = []
    for i in range(nsfos):
        for spin in ['A', 'B'] if is_unrestricted else ['AB']:
            symm = sfo_symmlabel[i]
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
                      overlaps=S[symm][spin][:, i])
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
            dEm[i, k] = abs(sfo1.energy - sfo2.energy)
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
    orbintm = S(sfos1, sfos2)**2 / (abs(dE(sfos1, sfos2)) * 27.2114079527)
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


def plot_sfos_prop(sfos1, sfos2, prop=S, cmap=None, use_relname=False):
    assert prop in [orbint, pauli, S, dE]

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
    else:
        plt.xticks(xticks, [orb.AMSlevels_name for orb in sfos2], rotation=90)
        plt.yticks(yticks, [orb.AMSlevels_name for orb in sfos1], rotation=0)
    plt.title(r'$' + prop.__name__ + r'(' + psi1 + r', ' + psi2 + r')$')
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

    p = '/Users/yumanhordijk/PhD/LewisAcid_coordination/calculations_final/EDA_vacuum/I2_N_sigma/full/adf.rkf'
    orbs = Orbitals(p)
    orbs.rename_fragments(['subsystem1', 'subsystem2'], ['LA', 'Substrate'])
    print('SIGMA')
    print(orbs.sfos['LA(54A)'].relname, orbs.sfos['Substrate(15A)'].relname)
    print((orbs.sfos['LA(54A)'].energy - orbs.sfos['Substrate(15A)'].energy) * 27.2114079527)
    print(orbs.sfos['LA(54A)'] @ orbs.sfos['Substrate(15A)'])


    print('\nBF3')
    p = '/Users/yumanhordijk/PhD/LewisAcid_coordination/calculations_final/EDA_vacuum/BF3_N_pi/full/adf.rkf'
    orbs = Orbitals(p)
    orbs.rename_fragments(['subsystem1', 'subsystem2'], ['LA', 'Substrate'])
    print('PI')
    print(orbs.sfos['LA(17A)'].relname, orbs.sfos['Substrate(15A)'].relname)
    print((orbs.sfos['LA(17A)'].energy - orbs.sfos['Substrate(15A)'].energy) * 27.2114079527)
    print(orbs.sfos['LA(17A)'] @ orbs.sfos['Substrate(15A)'])
    print()
    
    p = '/Users/yumanhordijk/PhD/LewisAcid_coordination/calculations_final/EDA_vacuum/BF3_N_sigma/full/adf.rkf'
    orbs = Orbitals(p)
    orbs.rename_fragments(['subsystem1', 'subsystem2'], ['LA', 'Substrate'])
    print('SIGMA')
    print(orbs.sfos['LA(17A)'].relname, orbs.sfos['Substrate(15A)'].relname)
    print((orbs.sfos['LA(17A)'].energy - orbs.sfos['Substrate(15A)'].energy)* 27.2114079527)
    print(orbs.sfos['LA(17A)'] @ orbs.sfos['Substrate(15A)'])
    