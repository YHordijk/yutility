import numpy as np
import scm.plams as plams
from ychem.utility.ytypes import check_hints, List, Either, Tuple, Vector, Matrix, Scalar
from ychem.utility import geometry
import networkx as nx
import itertools


@check_hints
def mol2graph(mol: plams.Molecule) -> nx.Graph:
    G = nx.Graph()
    G.add_nodes_from(mol.atoms)
    nx.set_node_attributes(G, {at: {'atnum': at.atnum} for at in mol})
    G.add_edges_from([(bond.atom1, bond.atom2) for bond in mol.bonds])
    return G


def atom_index(mol, atom):
    return mol.atoms.index(atom)


def atom_neighbour_tree(mol, atom):
    succs = dict(
        nx.bfs_successors(
            mol2graph(mol),
            atom,
            sort_neighbors=lambda lst: sorted(
                lst,
                key=lambda at: -
                at.atnum)))
    print(succs)


def get_priorities(mol, atom):
    G = mol2graph(mol)
    succs = dict(
        nx.bfs_successors(
            G,
            atom,
            sort_neighbors=lambda lst: sorted(
                lst,
                key=lambda at: -
                at.atnum)))

    def get_at_depth(depth):
        # returns the atoms in the groups at a certain depth
        groups = [[at] for at in succs[atom]]
        for i in range(depth):
            for j, group in enumerate(groups):
                group_ = []
                for at in group:
                    group_ = group_ + succs.get(at, [])
                groups[j] = group_
        return groups

    def get_next(depth, level):
        # returns the atoms at certain depth and level
        # the depth is the distance from the central atom
        # the level is the index of zip_longest of the groups
        # at that depth
        groups = get_at_depth(depth)
        return [group[level] if len(group) >= (
            level + 1) else None for group in groups]

    def max_levels(depth):
        groups = get_at_depth(depth)
        return max(len(group) for group in groups)

    def max_depth():
        i = 0
        groups = get_at_depth(i)
        while max(len(group) for group in groups) > 0:
            i += 1
            groups = get_at_depth(i)
        return i

    def atnums(ats):
        return [None if at is None else at.atnum for at in ats]

    # # check if there are 4 bonds to this atom, if not it cannot be chiral
    # if len(succs[atom]) <= 3:
    #     return None
    # ranks will be stored here
    # will contain information on which group it originates from
    # which atom was selected as endpoint and the depth and level it occured
    ranking = []
    # track if done to escape from double loop
    done = False
    # look through all depths which is the bond distance from the source
    for d in range(max_depth()):
        # and levels, which is the max of the number of atoms in each group at
        # the current depth
        for l in range(max_levels(d)):
            # get the atoms at this depth and level but only if the atom number is from a group that is not finished
            # if we dont do this there is a chance that finished groups interfere with non-finished groups
            # (by providing atnums that would otherwise be unique)
            ats = [
                at if not any(
                    i == rank[0] for rank in ranking) else None for i, at in enumerate(
                    get_next(
                        d, l))]
            atns = atnums(ats)
            # check if there are unique atnums here, if there are unique ones
            # we can stop checking that group
            unqs = [
                at for at in ats if at is not None and atns.count(
                    at.atnum) == 1]
            unqidxs = [ats.index(unq) for unq in unqs]
            # for each unique atom we append the ranking
            [ranking.append((i, succs[atom][i], d, l, ats[i]))
             for i in unqidxs if not any(rank[0] == i for rank in ranking)]
            # stopping condition
            if len(ranking) == len(succs[atom]):
                done = True
                break
        if done:
            break

    # the groups are still not ordered
    # get the paths between root and the decisionmaking atom
    paths = [nx.shortest_path(G, atom, rank[4]) for rank in ranking]
    # convert to atomic numbers
    patnums = [[at.atnum for at in path] for path in paths]
    # array used to store the order of the groups
    curr = np.arange(len(paths))
    # go through the paths level by level
    for patnum in list(itertools.zip_longest(*patnums))[1:]:
        # if list is too short there will be a None in the list
        # if so, the group is taken care of and we can mark it as done (False
        # in undone)
        undone = np.array([curr.tolist().index(i)
                          for i, atn in enumerate(patnum) if atn is not None])
        undone_idxs = np.array(
            [i for i, atn in enumerate(patnum) if atn is not None])
        sidxs = np.argsort(-np.array(patnum)[undone_idxs])
        curr[undone] = undone_idxs[sidxs]

    # get the group roots
    ranking = [ranking[c][1] for c in curr]
    return ranking


def Nneighbours(mol, atom):
    return len([bond for bond in mol.bonds if atom in bond])


@check_hints
def enantiomer(mol: plams.Molecule, atom: plams.Atom) -> Either(str, None):
    if Nneighbours(mol, atom) != 4:
        return None
    ranking = get_priorities(mol, atom)
    if len(ranking) != Nneighbours(mol, atom):
        return None
    # now to check which way the groups are turning we first
    # gather the bond vectors originating from the central atom
    vecs = [np.array(at.coords) - np.array(atom.coords) for at in ranking]
    # the align the vectors such that groups 3 and 4 are in the xy plane
    # using the right-hand rule we know that when group 1 is in the positive
    # direction of the normal vector of the plane formed by 3 and 4 it must be R
    # if it is opposite it must be S
    vecs = geometry.align_vecs_to_plane(vecs, vecs[-1], vecs[-2])
    if vecs[0][2] > 0:
        return 'R'
    return 'S'


@check_hints
def get_enantiomers(mol: plams.Molecule) -> List(Tuple(plams.Atom, str)):
    ret = []
    for atom in mol.atoms:
        enant = enantiomer(mol, atom)
        if enant is not None:
            ret.append((atom, enant))
    return ret


@check_hints
def print_chirality(mol: plams.Molecule):
    chirs = get_enantiomers(mol)
    ats = [(mol.atoms[i], chir)
           for i, chir in enumerate(chirs) if chir is not None]
    if len(ats) == 0:
        print('This molecule is not chiral')
    elif len(ats) == 1:
        print(
            f'This molecule is chiral with {len(ats)} chiral centrum [{ats[0][0].symbol + str(mol.atoms.index(ats[0][0])+1) + "=" + chirs[mol.atoms.index(ats[0][0])]}]')
    else:
        print(
            f'This molecule is chiral with {len(ats)} chiral centra [{", ".join(at.symbol + str(mol.atoms.index(at)+1) + "=" + chir for at, chir in ats)}]')


@check_hints
def EZ_isomer(mol: plams.Molecule, bond: plams.Bond) -> Either(str, None):
    if bond.order != 2:
        return None

    rank1 = [at for at in get_priorities(mol, bond.atom1) if at != bond.atom2]
    rank2 = [at for at in get_priorities(mol, bond.atom2) if at != bond.atom1]

    if len(rank1) == 0 or len(rank2) == 0:
        return None
    high1, high2 = atom_index(mol, rank1[0]), atom_index(mol, rank2[0])

    geometry.align_to_plane(mol, [rank1[0], bond.atom1, bond.atom2])
    geometry.center(mol, bond.atom1)
    geometry.align_to_axis(mol, [bond.atom1, bond.atom2])
    if mol.atoms[high1].coords[1] / mol.atoms[high2].coords[1] < 0:
        return 'E'
    return 'Z'


@check_hints
def compare_atoms(mol1: plams.Molecule, at1: plams.Atom,
                  mol2: plams.Molecule, at2: plams.Atom) -> bool:
    if at1.atnum != at2.atnum:
        return False
    if Nneighbours(mol1, at1) != Nneighbours(mol2, at2):
        return False

    print(atom_neighbour_tree(mol1, at1))
    prio1 = [at.atnum for at in get_priorities(mol1, at1)]
    prio2 = [at.atnum for at in get_priorities(mol2, at2)]
    return prio1 == prio2


@check_hints
def compare_molecules(mol1: plams.Molecule, mol2: plams.Molecule,
                      rmsd_thresh: Scalar = 1e-3) -> bool:
    G1, G2 = mol2graph(mol1), mol2graph(mol2)
    if not nx.is_isomorphic(G1, G2, node_match=lambda a1,
                            a2: a1['atnum'] == a2['atnum']):
        return False

    if plams.Molecule.rmsd(mol1, mol2) < rmsd_thresh:
        return True

    enant1, enant2 = get_enantiomers(mol1), get_enantiomers(mol2)
    if len(enant1) != len(enant2):
        return False
    for e1, e2 in zip(enant1, enant2):
        if not compare_atoms(mol1, e1[0], mol2, e2[0]):
            return False
    return True


@check_hints
def get_ez_isomers(mol: plams.Molecule) -> List(Tuple(plams.Bond, str)):
    ret = []
    for bond in mol.bonds:
        ez = EZ_isomer(mol, bond)
        if ez is not None:
            ret.append((bond, ez))
    return ret


@check_hints
def print_ez_isomerism(mol: plams.Molecule):
    ezs = get_ez_isomers(mol)

    bonds = ezs
    if len(bonds) == 0:
        print('This molecule does not have E/Z isomerism')
        return
    bondstrs = [
        f'{bond[0].atom1.symbol}({atom_index(mol, bond[0].atom1)+1})={bond[0].atom2.symbol}({atom_index(mol, bond[0].atom2)+1})' for bond in bonds]
    print(
        f'This molecule has {len(bonds)} E/Z isomerisms with [{", ".join(bondstrs[i] + ":" + bonds[i][1] for i in range(len(bonds)))}]')


mol1 = plams.Molecule('../tests/data/compare_mols/testa_1.xyz')
mol1.guess_bonds()
# print_chirality(mol1)
mol2 = plams.Molecule('../tests/data/compare_mols/testa_2.xyz')
mol2.guess_bonds()
# print_chirality(mol2)

print(compare_molecules(mol1, mol2))
