import numpy as np
from scm import plams
from yviewer import viewer


class CubeFile:
    def __init__(self, path):
        self.path = path
        self.load(path)

    def load(self, path):
        with open(path) as cub:
            lines = [line.strip() for line in cub.readlines()]

        natoms, origin = int(lines[2].split()[0]), np.array(lines[2].split()[1:]).astype(float) * 0.52918
        xvec, yvec, zvec = np.array([line.split()[1:] for line in lines[3:6]]).astype(float) * 0.52918
        xn, yn, zn = np.array([line.split()[0] for line in lines[3:6]]).astype(int)
        atomcoords = []
        atnums = []
        for line in lines[6:6+natoms]:
            atomcoords.append(np.array(line.split()[2:]).astype(float) * 0.52918)
            atnums.append(int(line.split()[0]))

        atomcoords = np.array(atomcoords)
        points = []
        for line in lines[6+natoms:]:
            points.extend(line.split())
        points = np.array(points).astype(float)

        pos = []
        for x in range(xn):
            for y in range(yn):
                for z in range(zn):
                    pos.append(origin + x*xvec + y*yvec + z*zvec - np.mean(atomcoords, axis=0))
        pos = np.array(pos)
        pos[:, 1] *= -1

        self.positions = pos
        self.values = points
        self.atomnumbers = atnums
        self.atomcoords = atomcoords - np.mean(atomcoords, axis=0)
        self.xvec = xvec
        self.yvec = yvec
        self.zvec = zvec
        self.xn = xn
        self.yn = yn
        self.zn = zn

    def get_molecule(self):
        mol = plams.Molecule()
        for atnum, atcoord in zip(self.atomnumbers, self.atomcoords):
            mol.add_atom(plams.Atom(atnum=atnum, coords=atcoord))
        mol.guess_bonds()
        return mol

    def set_cutoff(self, cutoff=.03):
        self.cutoffindices = np.arange(len(self.values))[np.abs(self.values) >= cutoff]
        return self.cutoffindices


    def set_colors(self, colors=[(255, 0, 0), (0, 0, 255)]):
        self.colors = []
        for value in self.values:
            if value >= 0:
                self.colors.append(colors[0])
            elif value < 0:
                self.colors.append(colors[1])
        self.colors = np.array(self.colors)
        return self.colors


    def show(self, **kwargs):
        self.set_cutoff(kwargs.get('cutoff', .03))
        self.set_colors(**kwargs)
        idx = self.cutoffindices
        viewer.show(self.get_molecule(), molinfo=[{'cub': [self.positions[idx], self.colors[idx]]}])


def show_multiple(cubs, molinfo=None, **kwargs):
    mols = []
    molinfo = molinfo or [{} for _ in cubs]
    for i, cub in enumerate(cubs):
        idx = cub.set_cutoff(kwargs.get('cutoff', .03))
        cub.set_colors(**kwargs)
        molinfo[i]['cub'] = [cub.positions[idx], cub.colors[idx]]
        mols.append(cub.get_molecule())
    viewer.show(mols, molinfo=molinfo)
