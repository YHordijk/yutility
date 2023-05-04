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

    # determine if the calculation used regions or not
    natoms = reader.read('Geometry', 'nr of atoms')
    nfrags = reader.read('Geometry', 'nr of fragments')
    ret['used_regions'] = natoms != nfrags

    return ret


def read_MO_data(reader):
    calc_info = get_calc_info(reader)

    ret = {}
    ## symlabels
    ret['symlabels'] = reader.read('Symmetry', 'symlab').strip().split()

    ## number of MOs
    ret['nmo'] = {}
    for symlabel in ret['symlabels']:
        if calc_info['unrestricted_mos']:
            ret['nmo'][symlabel] = {
                'A': reader.read(symlabel, 'nmo_A'),
                'B': reader.read(symlabel, 'nmo_B')
            }
        else:
            ret['nmo']['symlabel'] = reader.read(symmlabel, 'nmo_A')

    ## MO energies
    ret['energy'] = {}
    energyprefix = 'escale' if calc_info['relativistic'] else 'eps'
    for symlabel in ret['symlabels']:
        if calc_info['unrestricted_mos']:
            ret['energy'][symlabel] = {
                'A': reader.read(symlabel, f'{energyprefix}_A'),
                'B': reader.read(symlabel, f'{energyprefix}_B')
            }
        else:
            ret['energy']['symlabel'] = reader.read(symmlabel, f'{energyprefix}_A')

    ## MO occupations
    ret['occs'] = {}
    for symlabel in ret['symlabels']:
        if calc_info['unrestricted_mos']:
            ret['occs'][symlabel] = {
                'A': reader.read(symlabel, 'froc_A'),
                'B': reader.read(symlabel, 'froc_B')
            }
        else:
            ret['occs']['symlabel'] = reader.read(symmlabel, 'froc_A')
 
    ## MO coefficients
    ret['coeffs'] = {}
    for symlabel in ret['symlabels']:
        if calc_info['unrestricted_mos']:
            ret['coeffs'][symlabel] = {
                'A': reader.read(symlabel, 'Eig-CoreSFO_A'),
                'B': reader.read(symlabel, 'Eig-CoreSFO_B')
            }
        else:
            ret['coeffs']['symlabel'] = reader.read(symmlabel, 'Eig-CoreSFO_A')

    
    return ret


def read_SFO_data(reader):
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
        return  Srowsfixed


    calc_info = get_calc_info(reader)

    ret = {}
    ## symlabels
    ret['symlabels'] = reader.read('Symmetry', 'symlab').strip().split()

    ## number of SFOs per symlabel
    ret['nsfo'] = {}
    norb = ensure_list(reader.read('Symmetry', 'norb'))
    for i, symlabel in enumerate(ret['symlabels']):
        ret['nsfo'][symlabel] = norb[i]

    ## isfo, index of sfo in symlabel
    ret['isfo'] = reader.read('SFOs', 'isfo')

    ## fragment, index of fragment sfo belongs to
    ret['fragidx'] = reader.read('SFOs', 'fragment')

    ## fragmenttype, name of fragment sfo belongs to
    ret['fragtypes'] = reader.read('SFOs', 'fragtype').strip().split()

    ## fragorb, index of sfo in fragment
    ret['fragorb'] = reader.read('SFOs', 'fragorb')

    ## subspecies sfo belongs to
    ret['subspecies'] = reader.read('SFOs', 'subspecies').strip().split()

    ## ifo, order of sfo
    ret['ifo'] = reader.read('SFOs', 'ifo')

    # if the calculation did not use regions we have to
    # fix the ordering of the atoms (i.e. fragment indices)
    if not calc_info['used_regions']:
        atom_order_index = reader.read('Geometry', 'atom order index')
        atom_order_index = atom_order_index[:len(atom_order_index)//2]
        ret['fragidx'] = [atom_order_index.index(i) + 1 for i in ret['fragidx']]

    ## SFO energies
    energyprefix = 'escale' if calc_info['relativistic'] else 'energy'
    if calc_info['unrestricted_sfos']:
        ret['energy'] = {
            'A': reader.read('SFOs', f'{energyprefix}'),
            'B': reader.read('SFOs', f'{energyprefix}_B')
        }
    else:
        ret['energy'] = reader.read('SFOs', f'{energyprefix}'),

    ## overlaps
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
            ret['overlaps']['symlabel'] = square_overlaps(S, nsfo)

    return ret



if __name__ == '__main__':
    p = '/Users/yumanhordijk/PhD/yutility/test/orbitals/rkf/methyl.rkf'
    reader = plams.KFReader(p)
    # print(read_MO_data(reader))
    log.log(read_SFO_data(reader))

    # p = '/Users/yumanhordijk/PhD/yutility/test/orbitals/rkf/methyl_unrestricted_frags_symm.rkf'
    # reader = plams.KFReader(p)
    # # print(read_MO_data(reader))
    # print(read_SFO_data(reader))   
