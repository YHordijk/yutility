import numpy as np
import scm.plams as plams
from math import sin, cos
from yutility import molecule as yu_molecule
from yutility.ytypes import Matrix, Either, Vector


normal_to_plane = {'xy': np.array([0., 0., 1.]),
                   'yx': np.array([0., 0., 1.]),
                   'xz': np.array([0., 1., 0.]),
                   'zx': np.array([0., 1., 0.]),
                   'yz': np.array([1., 0., 0.]),
                   'zy': np.array([1., 0., 0.])}

axis_directions = {'x': np.array([1.0, 0.0, 0.0]),
                   'y': np.array([0.0, 1.0, 0.0]),
                   'z': np.array([0.0, 0.0, 1.0])}


class Transform:
    def __init__(self, 
                 R: Matrix(3, 3) = None, 
                 T: Vector(3) = None, 
                 S: Vector(3) = None):
        self.M = self.build_matrix(R, T, S)

    def __str__(self):
        return str(self.M)

    def translate(self, 
                  T: Vector(3) = None, 
                  x: float = None, 
                  y: float = None, 
                  z: float = None):
        '''
        Add a translation component to the transformation matrix.
        Arguments can be given as a container of x, y, z values. They can also be given separately.
        You can also specify x, y and z components separately

        Example usage:
            Transform.translate([2, 3, 0])
            Transform.translate(x=2, y=3)
        '''

        if T is None:
            T = [x or 0, y or 0, z or 0]

        self.M = self.M @ self.build_matrix(T=T)

    def rotate(self, 
               R: Matrix(3, 3) = None, 
               x: float = None, 
               y: float = None, 
               z: float = None):
        r'''
        Add a rotational component to transformation matrix.
        Arguments can be given as a rotation matrix R \in R^3x3 or by specifying the angle to rotate along the x, y or z axes
        
        Example usage:
            Transform.rotate(get_rotation_matrix(x=1, y=-1))
            Transform.rotate(x=1, y=-1)
        '''
        if R is None:
            R = get_rotation_matrix(x=x, y=y, z=z)

        self.M = self.M @ self.build_matrix(R=R)

    def scale(self, 
              S: Vector(3) = None, 
              x: float = None, 
              y: float = None, 
              z: float = None):
        '''
        Add a scaling component to the transformation matrix.
        Arguments can be given as a container of x, y, z values.
        You can also specify x, y and z components separately

        Example usage:
            Transform.scale([0, 0, 3])
            Transform.scale(z=3)
        '''

        if S is None:
            S = [x or 0, y or 0, z or 0]

        self.M = self.M @ self.build_matrix(S=S)

    def build_matrix(self, 
                     R: Matrix(3, 3) = None, 
                     T: Vector(3) = None, 
                     S: Vector(3) = None):
        r'''
        Build and return a transformation matrix. 
        This matrix encodes in one matrix rotations, translations and scaling.

        M = | R   r |, where R \in R^3x3, r \in R^3, 0_3 = [0, 0, 0] \in R^3 and 1 \in R
            | 0_3 1 |

        When applied to a positional vector [x, y, z, 1] it will apply these 
        transformations simultaneously. 

        Note: This matrix is an element of the special SE(3) Lie group.
        '''
        R = R if R is not None else get_rotation_matrix()
        T = T if T is not None else np.array([0, 0, 0])
        S = S if S is not None else np.array([1, 1, 1])

        return np.array([
            [R[0, 0]*S[0], R[0, 1], R[0, 2], T[0]],
            [R[1, 0], R[1, 1]*S[1], R[1, 2], T[1]],
            [R[2, 0], R[2, 1], R[2, 2]*S[2], T[2]],
            [0, 0, 0, 1]
        ])
        
    def __call__(self, *args, **kwargs):
        return self.apply(*args, **kwargs)

    def apply(self, v: Matrix(3, ...)) -> Matrix(3, ...):
        r'''
        Applies the transformation matrix to vector(s) v \in R^3xN
        v should be a series of column vectors.

        Application is a three-step process.
         1) Append row vector of ones to the bottom of v
         2) Apply the transformation matrix M \dot v
         3) Remove the bottom row vector of ones and return the result
        '''
        v = np.asarray(v)
        N = v.shape[1]
        v = np.vstack([v, np.ones(N)])
        vprime = self.M @ v
        return vprime[:3, :]

    def __matmul__(self, other):
        if isinstance(other, Transform):
            return self.combine_transforms(other)

    def combine_transforms(self, other: 'Transform') -> 'Transform':
        '''
        Combine two different transform objects. This involves creating 
        a new Transform object and multiplying the two transform matrices
        and assigning it to the new object.
        '''
        new = Transform()
        new.M = self.M @ other.M
        return new


def align_to_plane(mol, atoms, plane='xy'):
    plane_normal = normal_to_plane[plane]
    plane_coords = [np.array(a.coords) for a in atoms]
    plane_vecs = (plane_coords[1] - plane_coords[0]), (plane_coords[2] - plane_coords[0])
    normal = np.cross(plane_vecs[0], plane_vecs[1])
    R = plams.rotation_matrix(normal, plane_normal)
    mol.rotate(R)


def align_vecs_to_plane(vecs, x, y, plane='xy'):
    # aligns the vectors in vecs to a plane by laying x and y in the plane
    vecs = np.array(vecs)
    plane_norm = normal_to_plane[plane]
    norm = np.cross(x, y)
    R = plams.rotation_matrix(norm, plane_norm)
    return (R @ vecs.T).T


def align_vecs_to_axis(vecs, vec, axis='x', plane='xy'):
    # aligns the vectors in vecs to a plane by laying x and y in the plane
    vecs = np.array(vecs)
    angle = plams.angle(vec, axis_directions[axis])
    R = plams.axis_rotation_matrix(normal_to_plane[plane], -angle)
    return (R @ vecs.T).T


def align_to_axis(mol, atoms, axis='x', plane='xy'):
    align_coords = [np.array(a.coords) for a in atoms]
    vec = align_coords[1] - align_coords[0]
    angle = plams.angle(vec, axis_directions[axis])
    print(angle)
    R = plams.axis_rotation_matrix(normal_to_plane[plane], -angle)
    mol.rotate(R)


def center(mol, atom):
    center_coords = np.array(atom.coords)
    mol.translate(-center_coords)


def random_point_on_sphere(radius, dim=3):
    x = np.random.randn(dim)
    x = x/np.linalg.norm(x) * radius
    return x


def random_point_in_sphere(max_radius, min_radius=0, dim=3):
    radius = np.random.rand() * (max_radius - min_radius) + min_radius
    return random_point_on_sphere(radius, dim=dim)


def get_rotation_matrix(x=None, y=None, z=None):
    R = np.eye(3)

    if x is not None:
        c = cos(x)
        s = sin(x)
        R = R @ np.array(([1, 0, 0],
                          [0, c, -s],
                          [0, s, c]))

    if y is not None:
        c = cos(y)
        s = sin(y)
        R = R @ np.array(([c, 0, s],
                          [0, 1, 0],
                          [-s, 0, c]))

    if z is not None:
        c = cos(z)
        s = sin(z)
        R = R @ np.array(([c, -s, 0],
                          [s, c, 0],
                          [0, 0, 1]))

    return R


def apply_rotmat(coords, R):
    return (R @ coords.T).T


def rotate(coords, x=None, y=None, z=None):
    return apply_rotmat(coords, get_rotation_matrix(x, y, z))


def center_mol(coords):
    c = np.mean(coords, axis=0)
    return coords - c


def align_molecule(molecule, origin=None):
    if len(yu_molecule.get_labeled_atoms(
            molecule, 'plane', origin=origin)) > 0:
        align_to_plane(
            molecule,
            yu_molecule.get_labeled_atoms(
                molecule,
                'plane',
                origin=origin))
    if len(yu_molecule.get_labeled_atoms(
            molecule, 'align', origin=origin)) > 0:
        align_to_axis(
            molecule,
            yu_molecule.get_labeled_atoms(
                molecule,
                'align',
                origin=origin))
    if len(yu_molecule.get_labeled_atoms(
            molecule, 'center', origin=origin)) > 0:
        center(
            molecule,
            yu_molecule.get_labeled_atoms(
                molecule,
                'center',
                origin=origin)[0])


if __name__ == '__main__':
    # import matplotlib.pyplot as plt
    # xs = [random_point_on_sphere(5, dim=1)[0] for _ in range(1000)]
    # plt.hist(xs)
    # plt.show()
    # print()

    T = Transform()
    T.translate([0, 0, 3])

    v = np.arange(12).reshape(3, 4)
    print(T(v))
