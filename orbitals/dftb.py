from yutility.orbitals import info
import numpy as np
from yutility import ensure_list
from scm import plams



def read_SFO_data(reader):  # noqa: N802
    calc_info = info.get_calc_info(reader)

    assert calc_info['engine'] == 'dftb'

    ret = {}

    # number of SFOs and number of fragments
    nsfo = reader.read('FragmentOrbitals', 'nOrbitals')
    ret['nsfo'] = {'A': nsfo, 'total': nsfo}
    ret['nfrags'] = reader.read('FragmentOrbitals', 'nFragments')

    ret['symlabel_by_sfo'] = ['A'] * nsfo
    ret['isfo'] = list(range(1, nsfo+1))

    # read the fragment index the orbital belongs to
    ret['fragidx'] = reader.read('FragmentOrbitals', 'Fragment')

    # read names of the fragments (use element names if we did not do a fragment analysis)
    if calc_info['used_regions']:
        ret['fragtypes'] = [reader.read('FragmentOrbitals', 'BaseNameFragFile').strip().split()[i-1] for i in ret['fragidx']]
        # fragment names are already uniques if regions were used
        ret['fraguniquenames'] = ret['fragtypes']
    else:
        ret['fragtypes'] = [reader.read('FragmentOrbitals', 'FragmentSymbols').strip().split()[i-1] for i in ret['fragidx']]
        # make fragment names unique as {element}:{index}
        ret['fraguniquenames'] = [f'{name}:{i}' for i, name in zip(ret['fragidx'], ret['fragtypes'])]

    # fragorb, index of sfo in the fragment its belongs to
    ret['fragorb'] = reader.read('FragmentOrbitals', 'iFO')
    ret['ifo'] = reader.read('FragmentOrbitals', 'iFO')

    # subspecies sfo belongs to
    ret['subspecies'] = reader.read('FragmentOrbitals', 'SubSpecies').strip().split()

    # SFO energies
    # and get index if energies are sorted
    ret['energy'] = {'AB': reader.read('FragmentOrbitals', 'Energies(1)')}
    ret['energyidx'] = np.argsort(ret['energy'])

    # SFO occupations
    ret['occupations'] = {'AB': np.array(reader.read('FragmentOrbitals', 'Occupations(1)'))}

    # get the indices of the orbitals relative to the HOMO for each fragment
    # loop through each fragment
    ret['relindices'] = {'AB': []}
    for fragidx in range(1, ret['nfrags']+1):
        # first make a list of bools signifying if the fragorb is occupied
        is_occupied = []
        for i, occ in enumerate(ret['occupations']['AB']):
            if ret['fragidx'][i] != fragidx:
                continue
            is_occupied.append(occ > 0)
        # build the indices for each fragorb. The index is determined by the number of occupied orbitals higher than this orb
        # if there are only unoccupied orbitals higher than this orb the index is determined by the number of fragments remaining
        relindices = [-sum(is_occupied[i:]) + 1 if sum(is_occupied[i:]) > 0 else i - sum(is_occupied) + 1 for i in range(len(is_occupied))]
        ret['relindices']['AB'].extend(relindices)

    # overlaps
    ret['overlaps'] = {'A': {'AB': np.array(reader.read('FragmentOrbitals', 'Overlaps(1)')).reshape(nsfo, nsfo)}}

    return ret


def read_MO_data(reader):  # noqa: N802
    calc_info = info.get_calc_info(reader)

    ret = {}
    # symlabels
    if ('Symmetry', 'symlab') in reader:
        ret['symlabels'] = reader.read('Symmetry', 'symlab').strip().split()
    elif ('Geometry', 'grouplabel') in reader:
        ret['symlabels'] = symmetry.labels[reader.read('Geometry', 'grouplabel').strip()]
    else:
        ret['symlabels'] = symmetry.labels['NOSYM']

    # number of MOs
    ret['nmo'] = {}
    for symlabel in ret['symlabels']:
        if calc_info['unrestricted_mos']:
            ret['nmo'][symlabel] = {
                'A': int(np.sqrt(len(reader.read(symlabel, 'Eig-CoreSFO_A')))),
                'B': int(np.sqrt(len(reader.read(symlabel, 'Eig-CoreSFO_B'))))
            }
        else:
            ret['nmo'][symlabel] = {'AB': int(np.sqrt(len(reader.read(symlabel, 'Eig-CoreSFO_A'))))}

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

if __name__ == '__main__':
    reader = plams.KFReader('../test/orbitals/rkf/substrate_dftb.rkf')
    reader = plams.KFReader('../test/orbitals/rkf/substrate_cat_complex_dftb.rkf')
    print(read_SFO_data(reader))
