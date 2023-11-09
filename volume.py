import numpy as np
from scm import plams
from yviewer import viewer
import scipy.spatial
import scipy.cluster
import math
from tcintegral import grid


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
        self.origin = origin - np.mean(atomcoords, axis=0)
        self.values = points
        self.atomnumbers = atnums
        self.atomcoords = atomcoords - np.mean(atomcoords, axis=0)
        self.xvec = xvec
        self.yvec = yvec
        self.zvec = zvec
        self.xn = xn
        self.yn = yn
        self.zn = zn
        self.cell_volume = np.linalg.norm(np.cross(self.xvec*self.xn, self.yvec*self.yn) @ (self.zvec*self.zn))  # volume of a trapezoid

        self.set_cutoff(0)  # by default include all points

    @property
    def volume(self):
        return self.cutoffindices.size / self.positions.size * self.cell_volume

    def get_molecule(self):
        mol = plams.Molecule()
        for atnum, atcoord in zip(self.atomnumbers, self.atomcoords):
            mol.add_atom(plams.Atom(atnum=atnum, coords=atcoord))
        mol.guess_bonds()
        return mol

    def set_cutoff(self, cutoff=.03):
        self.cutoff = cutoff
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
        self.set_cutoff(kwargs.pop('cutoff', self.cutoff))
        self.set_colors(**kwargs)
        idx = self.cutoffindices
        viewer.show(self.get_molecule(), molinfo=[{'cub': [self.positions[idx], self.colors[idx]]}])

    def to_grid(self):
        gridd = grid.Grid(sum([self.xvec, self.yvec, self.zvec]).tolist())
        gridd += grid.Cube(self.origin.tolist(), sum([self.xvec*self.xn, self.yvec*self.yn, self.zvec*self.zn]).tolist())
        gridd.values = self.values
        return gridd


    # def set_clusters(self, maxdist=0.4):
    #     def select_cluster(idx, cluster_list):
    #         for cluster in cluster_list:
    #             dists = D[:, idx][cluster]
    #             if any(dists < maxdist):
    #                 return cluster

    #     D = scipy.spatial.distance_matrix(self.positions, self.positions)

    #     pos_clusters = []
    #     neg_clusters = []

    #     for i, (pos, val) in enumerate(zip(self.positions, self.values)):
    #         print(i)
    #         if val > 0:
    #             c = select_cluster(i, pos_clusters)
    #             if not c:
    #                 pos_clusters.append([i])
    #             else:
    #                 c.append(i)
    #         else:
    #             c = select_cluster(i, neg_clusters)
    #             if not c:
    #                 neg_clusters.append([i])
    #             else:
    #                 c.append(i)
        
    #     print(pos_clusters)
    #     print(len(pos_clusters))



    def distance_matrix(self):
        return scipy.spation.distance_matrix(self.positions, self.positions)


    def get_hull(self):
        for cluster in self.clusters:
            points = self.positions[cluster]
            hull = scipy.spatial.ConvexHull(points)


# def group_points(points, dlimit):
#     groups = []
#     while points:
#         far_points = []
#         ref = points.pop()
#         groups.append([ref])
#         for point in points:
#             d = sum([(x1 - x2)**2 for x1, x2 in zip(point, ref)])
#             if d < dlimit:
#                 groups[-1].append(point)
#             else:
#                 far_points.append(point)

#         points = far_points

#     # perform average operation on each group
#     return [list(np.mean(x, axis=1).astype(int)) for x in groups]


# def get_distance(ref, point):
#     # print('ref: {} , point: {}'.format(ref, point))
#     x1, y1, z1 = ref
#     x2, y2, z2 = point
#     return math.hypot(x2 - x1, y2 - y1)



def show_multiple(cubs, molinfo=None, **kwargs):
    mols = []
    molinfo = molinfo or [{} for _ in cubs]
    for i, cub in enumerate(cubs):
        idx = cub.set_cutoff(kwargs.get('cutoff', .03))
        cub.set_colors(**kwargs)
        # cub.set_clusters()
        molinfo[i]['cub'] = [cub.positions[idx], cub.colors[idx]]
        mols.append(cub.get_molecule())
    viewer.show(mols, molinfo=molinfo)
