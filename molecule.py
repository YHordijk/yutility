import scm.plams as plams
import json
import numpy as np
from ReactionRunner.utility import parse_string


def el2num(elements):
    pd = plams.PeriodicTable
    return [pd.get_atomic_number(el) for el in elements]


def get_xyz(mol):
    s = f'{len(mol.atoms)}\n\n'
    for e, p in zip([atom.symbol for atom in mol.atoms], [atom.coords for atom in mol.atoms]):
        s += f'{e:2}\t{p[0]:>10.6f}\t{p[1]:>10.6f}\t{p[2]:>10.6f}\n'
    return s


def load(path):
    data = {}
    with open(path) as f:
        mol = plams.Molecule()
        lines = [line.strip() for line in f.readlines()]
        natoms = int(lines[0])
        comment = lines[1]
        atom_lines = lines[2:natoms+2]
        for line in atom_lines:
            symbol, x, y, z, *args = line.split()
            atom_identifier_data = {}
            atom_identifier_data['labels'] = []
            for arg in args:
                if '=' in arg:
                    key, value = arg.split('=')
                    atom_identifier_data[key.strip()] = value.strip()
                else:
                    atom_identifier_data['labels'].append(arg.strip())

            atom = plams.Atom(symbol=symbol, coords=(float(x), float(y), float(z)))
            atom.identifier_data = atom_identifier_data
            mol.add_atom(atom)

        flag_lines = lines[natoms+2:]
        flags = {}
        for line in flag_lines:
            line = line.strip()
            if line == '':
                continue
            # all flags are in format key=value1_value2_value3
            key, val = line.split('=')
            val = val.replace("'", '"')
            # split value1_value2_value3 into [value1, value2, value3] and cast to int
            val = [i if not i.isnumeric() else float(i) for i in val.strip().split('_')]
            flags[key.strip()] = val

        mol_identifier_data = {'charge': [0.], 'spinpol': [0.]}
        for flag, val in flags.items():
            # substituents have either 2 or 3 elements so special case
            if flag.startswith('R'):
                for i in val[:2]:
                    mol.atoms[int(i)-1].identifier_data['labels'].append(flag)
            # simple flags
            elif flag in ['center', 'conn', 'delete']:
                for i in val:
                    mol.atoms[int(i)-1].identifier_data['labels'].append(flag)

            # flags that are ordered
            elif flag in ['align', 'plane', 'active_atoms']:
                for j, i in enumerate(val):
                    mol.atoms[int(i)-1].identifier_data['labels'].append(f'{flag}_{j}')

            # elif flag in ['charge', 'spinpol']:
            #     mol_identifier_data[flag] = float(val[0])
            # elif flag in ['task', 'pre_rotation']:
            #     mol_identifier_data[flag] = val[0]

            # elif flag.startswith('default_R'):
            #     mol_identifier_data[flag] = val[0]
            elif flag  == 'reaction_specific':
                mol_identifier_data[flag] = json.loads(val[0].strip())
                # print(flag, json.loads(val[0].strip()), type(json.loads(val[0].strip())))
            else:
                mol_identifier_data[flag] = val

        for flag, val in mol_identifier_data.items():
            if len(val) == 1:
                mol_identifier_data[flag] = parse_string.parse_str(val[0])
        mol.identifier_data = mol_identifier_data

        data['natoms'] = natoms
        data['comment'] = comment
        data['molecule'] = mol
        flags_ = flags.copy()
        flags_.update(mol_identifier_data)
        data['flags'] = flags_

    return data


def copy_data(mol1, mol2):
    mol1['flags'] = mol2['flags']
    mol1['molecule'].identifier_data = mol2['molecule'].identifier_data
    for a1, a2 in zip(mol1['molecule'].atoms, mol2['molecule'].atoms):
        a1.identifier_data = a2.identifier_data


def get_labeled_atoms(molecule, label, origin=None, return_idx=False):
    atoms = []
    for atom in molecule.atoms:
        if origin is None or atom.identifier_data['origin'] == origin:
            if label in atom.identifier_data['labels']:
                atoms.append(atom)

    if len(atoms) == 0:
        for atom in molecule.atoms:
            if origin is None or atom.identifier_data['origin'] == origin:
                for flag in atom.identifier_data['labels']:
                    if flag.startswith(label):
                        i, flag_ = flag[::-1].split('_', 1)
                        flag_ = flag_[::-1]
                        i = i[::-1]
                        atoms.append((int(i), atom))
        atoms = [atom[1] for atom in sorted([atom for atom in atoms], key=lambda atom: atom[0])]
    if return_idx:
        return [molecule.atoms.index(a) for a in atoms]
    return atoms


def TSRC(molecule):
    atoms = get_labeled_atoms(molecule, 'active_atoms')
    tsats = {label.split('_')[-1] for atom in atoms for label in atom.identifier_data['labels'] if label.startswith('active_atoms')}
    dists = {}
    for i in tsats:
        atoms_ = get_labeled_atoms(molecule, f'active_atoms_{i}')
        if len(atoms_) != 2:
            continue
        idxs = get_labeled_atoms(molecule, f'active_atoms_{i}', return_idx=True)
        dists[i] = idxs, np.linalg.norm(np.array(atoms_[0].coords) - np.array(atoms_[1].coords))

    return dists


def TSRC_idx(molecule):
    atoms = get_labeled_atoms(molecule, 'active_atoms')
    tsats = {label.split('_')[-1] for atom in atoms for label in atom.identifier_data['labels'] if label.startswith('active_atoms')}
    for i in tsats:
        idxs = get_labeled_atoms(molecule, f'active_atoms_{i}', return_idx=True)
        if len(idxs) != 2:
            continue
        yield idxs



if __name__ == '__main__':
    mol = load(r"D:\Users\Yuman\Desktop\PhD\ychem\reaction_generation\input_mols\3d9e03fb7103af580fd95a3fa8a820856bcd70bea2d2363709438e1ad5ba10df.xyz")
    print(mol)
import scm.plams as plams
import json
import numpy as np
from ReactionRunner.utility import parse_string


def el2num(elements):
    pd = plams.PeriodicTable
    return [pd.get_atomic_number(el) for el in elements]


def get_xyz(mol):
    s = f'{len(mol.atoms)}\n\n'
    for e, p in zip([atom.symbol for atom in mol.atoms], [atom.coords for atom in mol.atoms]):
        s += f'{e:2}\t{p[0]:>10.6f}\t{p[1]:>10.6f}\t{p[2]:>10.6f}\n'
    return s


def load(path):
    data = {}
    with open(path) as f:
        mol = plams.Molecule()
        lines = [line.strip() for line in f.readlines()]
        natoms = int(lines[0])
        comment = lines[1]
        atom_lines = lines[2:natoms+2]
        for line in atom_lines:
            symbol, x, y, z, *args = line.split()
            atom_identifier_data = {}
            atom_identifier_data['labels'] = []
            for arg in args:
                if '=' in arg:
                    key, value = arg.split('=')
                    atom_identifier_data[key.strip()] = value.strip()
                else:
                    atom_identifier_data['labels'].append(arg.strip())

            atom = plams.Atom(symbol=symbol, coords=(float(x), float(y), float(z)))
            atom.identifier_data = atom_identifier_data
            mol.add_atom(atom)

        flag_lines = lines[natoms+2:]
        flags = {}
        for line in flag_lines:
            line = line.strip()
            if line == '':
                continue
            # all flags are in format key=value1_value2_value3
            key, val = line.split('=')
            val = val.replace("'", '"')
            # split value1_value2_value3 into [value1, value2, value3] and cast to int
            val = [i if not i.isnumeric() else float(i) for i in val.strip().split('_')]
            flags[key.strip()] = val

        mol_identifier_data = {'charge': [0.], 'spinpol': [0.]}
        for flag, val in flags.items():
            # substituents have either 2 or 3 elements so special case
            if flag.startswith('R'):
                for i in val[:2]:
                    mol.atoms[int(i)-1].identifier_data['labels'].append(flag)
            # simple flags
            elif flag in ['center', 'conn', 'delete']:
                for i in val:
                    mol.atoms[int(i)-1].identifier_data['labels'].append(flag)

            # flags that are ordered
            elif flag in ['align', 'plane', 'active_atoms']:
                for j, i in enumerate(val):
                    mol.atoms[int(i)-1].identifier_data['labels'].append(f'{flag}_{j}')

            # elif flag in ['charge', 'spinpol']:
            #     mol_identifier_data[flag] = float(val[0])
            # elif flag in ['task', 'pre_rotation']:
            #     mol_identifier_data[flag] = val[0]

            # elif flag.startswith('default_R'):
            #     mol_identifier_data[flag] = val[0]
            elif flag  == 'reaction_specific':
                mol_identifier_data[flag] = json.loads(val[0].strip())
                # print(flag, json.loads(val[0].strip()), type(json.loads(val[0].strip())))
            else:
                mol_identifier_data[flag] = val

        for flag, val in mol_identifier_data.items():
            if len(val) == 1:
                mol_identifier_data[flag] = parse_string.parse_str(val[0])
        mol.identifier_data = mol_identifier_data

        data['natoms'] = natoms
        data['comment'] = comment
        data['molecule'] = mol
        flags_ = flags.copy()
        flags_.update(mol_identifier_data)
        data['flags'] = flags_

    return data


def copy_data(mol1, mol2):
    mol1['flags'] = mol2['flags']
    mol1['molecule'].identifier_data = mol2['molecule'].identifier_data
    for a1, a2 in zip(mol1['molecule'].atoms, mol2['molecule'].atoms):
        a1.identifier_data = a2.identifier_data


def get_labeled_atoms(molecule, label, origin=None, return_idx=False):
    atoms = []
    for atom in molecule.atoms:
        if origin is None or atom.identifier_data['origin'] == origin:
            if label in atom.identifier_data['labels']:
                atoms.append(atom)

    if len(atoms) == 0:
        for atom in molecule.atoms:
            if origin is None or atom.identifier_data['origin'] == origin:
                for flag in atom.identifier_data['labels']:
                    if flag.startswith(label):
                        i, flag_ = flag[::-1].split('_', 1)
                        flag_ = flag_[::-1]
                        i = i[::-1]
                        atoms.append((int(i), atom))
        atoms = [atom[1] for atom in sorted([atom for atom in atoms], key=lambda atom: atom[0])]
    if return_idx:
        return [molecule.atoms.index(a) for a in atoms]
    return atoms


def TSRC(molecule):
    atoms = get_labeled_atoms(molecule, 'active_atoms')
    tsats = {label.split('_')[-1] for atom in atoms for label in atom.identifier_data['labels'] if label.startswith('active_atoms')}
    dists = {}
    for i in tsats:
        atoms_ = get_labeled_atoms(molecule, f'active_atoms_{i}')
        if len(atoms_) != 2:
            continue
        idxs = get_labeled_atoms(molecule, f'active_atoms_{i}', return_idx=True)
        dists[i] = idxs, np.linalg.norm(np.array(atoms_[0].coords) - np.array(atoms_[1].coords))

    return dists


def TSRC_idx(molecule):
    atoms = get_labeled_atoms(molecule, 'active_atoms')
    tsats = {label.split('_')[-1] for atom in atoms for label in atom.identifier_data['labels'] if label.startswith('active_atoms')}
    for i in tsats:
        idxs = get_labeled_atoms(molecule, f'active_atoms_{i}', return_idx=True)
        if len(idxs) != 2:
            continue
        yield idxs



if __name__ == '__main__':
    mol = load(r"D:\Users\Yuman\Desktop\PhD\ychem\reaction_generation\input_mols\3d9e03fb7103af580fd95a3fa8a820856bcd70bea2d2363709438e1ad5ba10df.xyz")
    print(mol)
