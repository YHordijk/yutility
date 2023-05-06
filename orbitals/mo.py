from scm import plams
import numpy as np
import os
import matplotlib.pyplot as plt
from yutility import log, units, ensure_list, squeeze_list, print_kf, run, plot, symmetry, decorators, orbitals
from yutility.ytypes import Either

j = os.path.join


class MOs:
    def __init__(self, kfpath=None, reader=None):
        assert reader or kfpath, 'Please provide a KFReader or path to a kf-file'
        self.reader = reader or plams.KFReader(kfpath)
        self.kfpath = kfpath
        if not self.kfpath:
            self.kfpath = self.reader.path
        self.get_calc_info()
        self.get_mos()

    def __iter__(self):
        return iter(self.mos)

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
        calc_info = orbitals.info.get_calc_info(self.reader)
        self.is_unrestricted = calc_info['unrestricted_mos']
        self.uses_molecular_fragments = calc_info['used_regions']
        self.uses_atomic_fragments = not calc_info['used_regions']
        self.is_relativistic = calc_info['relativistic']

        self.spins = ['A', 'B'] if self.is_unrestricted else ['AB']
        self.symlabels = calc_info['symlabels']

    def get_mos(self):
        data = orbitals.info.read_MO_data(self.reader)

        mo_data = []
        for symlabel in data['symlabels']:
            for spin in self.spins:
                for idx in range(data['nmo'][symlabel][spin]):
                    energy = data['energy'][symlabel][spin][idx]
                    mo_index = data['energy_idx'][symlabel][spin][idx]
                    relindex = mo_index + 1 - data['noccs']
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
                    mo_data.append({
                        'index':                mo_index + 1,
                        'relindex':             relindex,
                        'index_in_symlabel':    idx,
                        'name':                 f'{idx + 1}{symlabel}',
                        'relname':              relname,
                        # 'fragment_index':       data['fragidx'][idx],
                        # 'fragment':             data['fragtypes'][idx],
                        # 'fragment_orb_index':   data['fragorb'][idx],
                        # 'symmetry_type_index':  isfo,
                        'symmetry':             symlabel,
                        'energy':               energy * 27.21139664,
                        'spin':                 spin,
                        'reader':               self.reader,
                        'kfpath':               self.kfpath,
                        # 'overlaps':             data['overlaps'][symlabel][spin][isfo],
                        'occupation':           data['occs'][symlabel][spin][idx],
                        'atomic_fragments':     self.uses_atomic_fragments,
                    })
        self.mos = [MO(**mo_datum) for mo_datum in mo_data]


class MO:
    '''
    lightweight class holding MO data
    '''

    def __init__(self, **kwargs):
        for key, value in kwargs.items():
            setattr(self, key, value)

    def __repr__(self):
        return self.full_name

    @property
    def full_name(self):
        spin_part = ''
        if self.spin != 'AB':
            spin_part = f'_{self.spin}'

        return f'{self.name}{spin_part}'

    @property
    def relative_name(self):
        spin_part = ''
        if self.spin != 'AB':
            spin_part = f'_{self.spin}'

        return f'{self.relname}{spin_part}'

    @property
    def index_name(self):
        spin_part = ''
        if self.spin != 'AB':
            spin_part = f'_{self.spin}'

        return f'{self.index}{self.symmetry}{spin_part}' 

    @property
    def occupied(self):
        return self.occupation > 0

    def generate_orbital(self):
        return run.orbital_cub(self.rkf_path, self.index, orbtype='SCF', symlabel=self.symmetry)


if __name__ == '__main__':
    p = '../test/orbitals/rkf/methyl_unrestricted_frags_symm.rkf'
    mos = MOs(kfpath=p)
    for mo in sorted(mos, key=lambda mo: mo.energy):
        print(mo, mo.occupied, mo.relative_name)
