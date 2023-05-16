import numpy as np
from yutility import ensure_list, symmetry


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
    if ('Symmetry', 'symlab') in reader:
        ret['symlabels'] = reader.read('Symmetry', 'symlab').strip().split()
    else:
        ret['symlabels'] = symmetry.labels[reader.read('Geometry', 'grouplabel').strip()]

    # determine if MOs are unrestricted or not
    ret['unrestricted_mos'] = (ret['symlabels'][0], 'eps_B') in reader

    # determine if the calculation used regions or not
    natoms = reader.read('Geometry', 'nr of atoms')
    nfrags = reader.read('Geometry', 'nr of fragments')
    ret['used_regions'] = natoms != nfrags

    # get symmetry labels

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
    if ('Symmetry', 'symlab') in reader:
        ret['symlabels'] = reader.read('Symmetry', 'symlab').strip().split()
    else:
        ret['symlabels'] = symmetry.labels[reader.read('Geometry', 'grouplabel').strip()]

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

    # get fragment + index unique pairs to construct unique fragment names
    ret['fraguniquenames'] = []
    unqs = set([(name, i) for name, i in zip(ret['fragtypes'], ret['fragidx'])])
    for name, i in zip(ret['fragtypes'], ret['fragidx']):
        # check if fragname is unique
        if len([pair for pair in unqs if pair[0] == name]) == 1:
            ret['fraguniquenames'].append(name)
        else:
            ret['fraguniquenames'].append(f'{name}:{i}')
            
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
            if (occ == 0 and occ2 > 0):
                relindices.extend(np.arange(newidx - homoidx, i - homoidx + 1))
                newidx = i + 1
            if i == len(ret['ifo'])-2:
                relindices.extend(np.arange(newidx - homoidx, i - homoidx + 2))
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


def read_MO_data(reader):  # noqa: N802
    calc_info = get_calc_info(reader)

    ret = {}
    # symlabels
    if ('Symmetry', 'symlab') in reader:
        ret['symlabels'] = reader.read('Symmetry', 'symlab').strip().split()
    else:
        ret['symlabels'] = symmetry.labels[reader.read('Geometry', 'grouplabel').strip()]

    # number of MOs
    ret['nmo'] = {}
    for symlabel in ret['symlabels']:
        if calc_info['unrestricted_mos']:
            ret['nmo'][symlabel] = {
                'A': reader.read(symlabel, 'nmo_A'),
                'B': reader.read(symlabel, 'nmo_B')
            }
        else:
            ret['nmo'][symlabel] = {'AB': reader.read(symlabel, 'nmo_A')}
    ret['nmo']['total'] = sum(ensure_list(reader.read('Symmetry', 'norb')))

    # MO energies
    ret['energy'] = {}
    ret['energy']['sorted'] = {'A': [], 'B': [], 'AB': []}
    energyprefix = 'escale' if calc_info['relativistic'] else 'eps'
    for symlabel in ret['symlabels']:
        if calc_info['unrestricted_mos']:
            ret['energy'][symlabel] = {
                'A': ensure_list(reader.read(symlabel, f'{energyprefix}_A')),
                'B': ensure_list(reader.read(symlabel, f'{energyprefix}_B'))
            }
            ret['energy']['sorted']['A'].extend(ensure_list(reader.read(symlabel, f'{energyprefix}_A')))
            ret['energy']['sorted']['B'].extend(ensure_list(reader.read(symlabel, f'{energyprefix}_B')))
        else:
            ret['energy'][symlabel] = {'AB': ensure_list(reader.read(symlabel, f'{energyprefix}_A'))}
            ret['energy']['sorted']['AB'].extend(ensure_list(reader.read(symlabel, f'{energyprefix}_A')))
    ret['energy']['sorted']['A'] = sorted(ret['energy']['sorted']['A'])
    ret['energy']['sorted']['B'] = sorted(ret['energy']['sorted']['B'])
    ret['energy']['sorted']['AB'] = sorted(ret['energy']['sorted']['AB'])

    ret['energy_idx'] = {}
    all_idx = {'A': [], 'B': [], 'AB': []}
    for symlabel in ret['symlabels']:
        if calc_info['unrestricted_mos']:
            energies_A = ensure_list(reader.read(symlabel, f'{energyprefix}_A'))
            idx_A = []
            for energy_A in energies_A:
                idx = ret['energy']['sorted']['A'].index(energy_A)
                while idx in all_idx['A']:
                    idx += 1
                idx_A.append(idx) 
                all_idx['A'].append(idx)

            energies_B = ensure_list(reader.read(symlabel, f'{energyprefix}_B'))
            idx_B = []
            for energy_B in energies_B:
                idx = ret['energy']['sorted']['B'].index(energy_B)
                while idx in all_idx['B']:
                    idx += 1
                idx_B.append(idx)
                all_idx['B'].append(idx)

            ret['energy_idx'][symlabel] = {
                'A': idx_A,
                'B': idx_B
            }
        else:
            energies = ensure_list(reader.read(symlabel, f'{energyprefix}_A'))
            idxs = []
            for energy_B in energies:
                idx = ret['energy']['sorted']['AB'].index(energy_B)
                while idx in all_idx['AB']:
                    idx += 1
                idxs.append(idx)
                all_idx['AB'].append(idx)
            ret['energy_idx'][symlabel] = {'AB': idxs}

    # MO occupations
    ret['occs'] = {}
    ret['noccs'] = {}
    for symlabel in ret['symlabels']:
        if calc_info['unrestricted_mos']:
            ret['occs'][symlabel] = {
                'A': ensure_list(reader.read(symlabel, 'froc_A')),
                'B': ensure_list(reader.read(symlabel, 'froc_B'))
            }

            ret['noccs'][symlabel] = {
                'A': len([occ for occ in ensure_list(reader.read(symlabel, 'froc_A')) if occ > 0]),
                'B': len([occ for occ in ensure_list(reader.read(symlabel, 'froc_B')) if occ > 0])
            }
        else:
            ret['occs'][symlabel] = {'AB': ensure_list(reader.read(symlabel, 'froc_A'))}
            ret['noccs'][symlabel] = {'AB': len([occ for occ in ensure_list(reader.read(symlabel, 'froc_A')) if occ > 0]),}
 
    # MO coefficients
    ret['coeffs'] = {}
    for symlabel in ret['symlabels']:
        if calc_info['unrestricted_mos']:
            c_A = ensure_list(reader.read(symlabel, 'Eig-CoreSFO_A'))
            c_B = ensure_list(reader.read(symlabel, 'Eig-CoreSFO_B'))
            ret['coeffs'][symlabel] = {
                'A': np.array(c_A).reshape(ret['nmo'][symlabel]['A'], ret['nmo'][symlabel]['A']),
                'B': np.array(c_B).reshape(ret['nmo'][symlabel]['B'], ret['nmo'][symlabel]['B'])
            }
        else:
            c = ensure_list(reader.read(symlabel, 'Eig-CoreSFO_A'))
            ret['coeffs'][symlabel] = {'AB': np.array(c).reshape(ret['nmo'][symlabel]['AB'], ret['nmo'][symlabel]['AB'])}

    # get index of MO in symmetry label
    ret['symmidx'] = {}


    return ret
