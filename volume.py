import numpy as np


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
        for line in lines[6:6+natoms]:
            atomcoords.append(np.array(line.split()[2:]).astype(float) * 0.52918)

        atomcoords = np.array(atomcoords)
        points = []
        for line in lines[6+natoms:]:
            points.extend(line.split())
        points = np.array(points).astype(float)

        limit = .03
        keep_idx = []
        colors = []
        for i, point in enumerate(points):
            if point >= limit:
                colors.append((0, 219, 100))
                keep_idx.append(i)
            elif point < -limit:
                colors.append((255, 165, 0))
                keep_idx.append(i)
            else:
                colors.append((0, 0, 0))

        pos = []
        for x in range(xn):
            for y in range(yn):
                for z in range(zn):
                    pos.append(origin + x*xvec + y*yvec + z*zvec - np.mean(atomcoords, axis=0))
        pos = np.array(pos)
        pos[:, 1] *= -1

        colors = np.array(colors)
        pos = pos[keep_idx]
        colors = colors[keep_idx]
        
        self.positions = pos
        self.values = points
        self.xvec = xvec
        self.yvec = yvec
        self.zvec = zvec
        self.xn = xn
        self.yn = yn
        self.zn = zn
