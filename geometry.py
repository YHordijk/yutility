import numpy as np
import scm.plams as plams
from math import sin, cos
from yutility import molecule as yu_molecule
from yutility.ytypes import Matrix, Either, Vector
import itertools
import scipy


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
            S = [x or 1, y or 1, z or 1]
        elif isinstance(S, (float, int)):
            S = [S, S, S]

        self.M = self.M @ self.build_matrix(S=S)

    def build_matrix(self, 
                     R: Matrix(3, 3) = None, 
                     T: Vector(3) = None, 
                     S: Vector(3) = None) -> Matrix(4, 4):
        r'''
        Build and return a transformation matrix. 
        This 4x4 matrix encodes rotations, translations and scaling.

        M = | R@diag(S)   T |, where R \in R^3x3, r \in R^3, 0_3 = [0, 0, 0] \in R^3 and 1 \in R
            | 0_3         1 |

        When applied to a positional vector [x, y, z, 1] it will apply these 
        transformations simultaneously. 

        Note: This matrix is an element of the special SE(3) Lie group.
        '''

        # set the default matrix and vectors
        R = R if R is not None else get_rotation_matrix()
        T = T if T is not None else np.array([0, 0, 0])
        S = S if S is not None else np.array([1, 1, 1])

        return np.array([
            [R[0, 0]*S[0], R[0, 1],      R[0, 2],      T[0]],
            [R[1, 0],      R[1, 1]*S[1], R[1, 2],      T[1]],
            [R[2, 0],      R[2, 1],      R[2, 2]*S[2], T[2]],
            [0,            0,            0,            1]])
        
    def __call__(self, *args, **kwargs):
        return self.apply(*args, **kwargs)

    def apply(self, v: Matrix(..., 3)) -> Matrix(..., 3):
        r'''
        Applies the transformation matrix to vector(s) v \in R^Nx3
        v should be a series of column vectors.

        Application is a three-step process.
         1) Append row vector of ones to the bottom of v
         2) Apply the transformation matrix M \dot v
         3) Remove the bottom row vector of ones and return the result
        '''
        v = np.atleast_2d(v)
        v = np.asarray(v).T
        N = v.shape[1]
        v = np.vstack([v, np.ones(N)])
        vprime = self.M @ v
        return vprime[:3, :].T

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

    def get_rotation_matrix(self):
        return self.M[:3, :3]


class KabschTransform(Transform):
    def __init__(self, 
                 X: Matrix(..., 3), 
                 Y: Matrix(..., 3)):
        '''
        Use Kabsch-Umeyama algorithm to calculate the optimal rotation matrix and translation
        to superimpose X onto Y. 

        It is numerically stable and works when the covariance matrix is singular.
        Both sets of points must be the same size for this algorithm to work.
        The coordinates are first centered onto their centroids before determining the 
        optimal rotation matrix.

        References
            https://en.wikipedia.org/wiki/Orthogonal_Procrustes_problem
            https://en.wikipedia.org/wiki/Kabsch_algorithm
        '''

        # make sure arrays are 2d and the same size
        X, Y = np.atleast_2d(X), np.atleast_2d(Y)
        assert X.shape == Y.shape, f"Matrices X with shape {X.shape} and Y with shape {Y.shape} are not the same size"

        # center the coordinates
        centroid_x = np.mean(X, axis=0)
        centroid_y = np.mean(Y, axis=0)
        Xc = X - centroid_x
        Yc = Y - centroid_y

        # calculate covariance matrix
        H = Xc.T @ Yc

        # first do single value decomposition on covariance matrix
        # this step ensures that the algorithm is numerically stable
        # and removes problems with singular covariance matrices
        U, _, V = scipy.linalg.svd(H)

        # get the sign of the determinant of V.T @ U.T
        sign = np.sign(np.linalg.det(V.T@U.T))
        # then build the optimal rotation matrix
        d = np.diag([1, 1, sign])
        R = V.T @ d @ U.T

        # build the transformation:
        # for a sequence of transformation operations we have to invert their order
        # We have that Y ~= (R @ (X - centroid_x).T).T + centroid_y
        # the normal order is to first translate X by -centroid_x
        # then rotate with R
        # finally translate by +centroid_y
        self.M = self.build_matrix()
        self.translate(centroid_y)
        self.rotate(R)
        self.translate(-centroid_x)


def RMSD(X: Matrix(..., 3), 
         Y: Matrix(..., 3) = None,
         axis: int = None) -> float:
    '''
    Calculate Root Mean Squared Deviations between two sets of points
    X and Y. if Y is not given, the RMSD of X will be evaluated w.r.t. the origin.
    Optionally the axis can be given to calculate the RMSD along different axes.
    '''
    Y = Y if Y is not None else 0
    return np.sqrt(np.sum((X - Y)**2, axis=axis)/X.shape[0])


def RMSD_kabsch(X: Matrix(..., 3), 
                Y: Matrix(..., 3)) -> float:
    '''
    Calculate the RMSD after Kabsch-transformation of the coordinates.
    This will yield the smallest possible RMSD for the two sets of coordinates
    '''
    Tkabsch = KabschTransform(X, Y)
    return RMSD(Tkabsch(X), Y)


def RMSD_combinatorial(X: Matrix(..., 3), 
                       Y: Matrix(..., 3),
                       use_kabsch: bool = True) -> float:
    '''
    Exhaustively look through permutations of X to minimize the RMSD.
    This can get very expensive very quickly. Optionally we can use the Kabsch
    algorithm to calculate the smallest RMSD.

    Be careful! This method scales O(n!)
    '''

    RMSD_func = RMSD_kabsch if use_kabsch else RMSD
    return min([RMSD_func(X, y) for y in itertools.permutations(Y)])


@plams.add_to_class(plams.Molecule)
def apply_transform(self, transform: Transform):
    for atom in self.atoms:
        atom.coords = transform.apply(atom.coords)[0]


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
    R = plams.axis_rotation_matrix(normal_to_plane[plane], -angle)
    mol.rotate(R)


def center(mol, atom):
    center_coords = np.array(atom.coords)
    mol.translate(-center_coords)


def random_point_on_sphere(radius, N=1, dim=3):
    x = np.random.randn(dim, N)
    x = x/np.linalg.norm(x, axis=0) * radius
    return x.T


def random_point_in_sphere(max_radius, N=1, min_radius=0, dim=3):
    radii = np.random.rand(N) * (max_radius - min_radius) + min_radius
    x = np.random.randn(dim, N)
    x = x/np.linalg.norm(x, axis=0) * radii
    return random_point_on_sphere(radii, dim=dim, N=N)


def get_rotation_matrix(x=None, y=None, z=None):
    # create a rotation matrix based on Tait-Bryant angles
    # in this system, x, y, and z are angles of rotation around
    # the corresponding axes. This function uses right-handed 
    # convention
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


def vector_align_rotmat(a, b):
    # normalize the vecs
    a = np.array(a) / np.linalg.norm(a)
    b = np.array(b) / np.linalg.norm(b)

    v = np.cross(a, b)
    skew = np.array([
        [0, -v[2], v[1]],
        [v[2], 0, -v[0]],
        [-v[1], v[0], 0]
    ])

    c = a @ b

    R = np.eye(3) + skew + skew @ skew / (1 + c)
    return R







# def align_molecule(molecule, origin=None):
#     if len(yu_molecule.get_labeled_atoms(
#             molecule, 'plane', origin=origin)) > 0:
#         align_to_plane(
#             molecule,
#             yu_molecule.get_labeled_atoms(
#                 molecule,
#                 'plane',
#                 origin=origin))
#     if len(yu_molecule.get_labeled_atoms(
#             molecule, 'align', origin=origin)) > 0:
#         align_to_axis(
#             molecule,
#             yu_molecule.get_labeled_atoms(
#                 molecule,
#                 'align',
#                 origin=origin))
#     if len(yu_molecule.get_labeled_atoms(
#             molecule, 'center', origin=origin)) > 0:
#         center(
#             molecule,
#             yu_molecule.get_labeled_atoms(
#                 molecule,
#                 'center',
#                 origin=origin)[0])


def cartesian_to_spherical(coords):
    x, y, z = coords[:, 0], coords[:, 1], coords[:, 2]
    r = np.linalg.norm(coords)
    theta = np.arccos(z/r)
    phi = y/np.abs(y) * np.arccos(x/np.sqrt(x**2 + y**2))
    return np.vstack([x, y, z]).T

if __name__ == '__main__':
    import matplotlib.pyplot as plt
    
    # X = np.array([[-1, 0, 0], [0, 0, 0], [0, 1, -1]])
    X = plams.Molecule('/Users/yumanhordijk/PhD/ychem/calculations2/c1d4ca95a3911eb1f79bf4ef91cc7a88b479d7dc8357860bfdb3e577747ebc3a/catalyst/input_mol.xyz').as_array()
    T = Transform()
    T.translate(y=-1.5)
    T.rotate(z=.3)
    Y = T(X)
    R = []

    trans = np.linspace(-5, 5, 100)
    rot = np.linspace(-np.pi, np.pi, 100)

    for trans_ in trans:
        R.append([])
        for rot_ in rot:
            T = Transform()
            T.translate(x=trans_)
            T.rotate(z=rot_)
            R[-1].append(RMSD(X, T(Y)))

    plt.title(r'RMSD between two ZnCl$_2$ molecules first Translated, then Rotated')
    im = plt.imshow(R, extent=(trans.min(), trans.max(), rot.min(), rot.max()), origin='lower', cmap='coolwarm', aspect='auto')
    cbar = plt.colorbar(im, label=r'RMSD (Angstrom)')
    plt.xlabel('Translation on X-axis (Angstrom)')
    plt.ylabel('Rotation around Z-axis')
    plt.tight_layout()
    plt.show()



    # demo for Transforms and Kabsch algorithm
    # we have a random Transform which rotates, translates and scales
    T = Transform()
    T.rotate(x=1, y=1, z=-1)
    T.scale(x=20, y=-100)
    T.translate(x=3, y=0, z=-10)
    T.rotate(x=-1, y=-1, z=1)

    # we generate some coordinates X
    # X = np.random.randn(5, 3)
    X = np.arange(5*3).reshape(5, 3)
    # apply the transformation to X to obtain Y
    Y = T(X)

    # get the kabsch optimal transformation
    Tkabsch = KabschTransform(X, Y)

    # print the RMSD between Tkabsch(X) and Y
    # small value indicates that alignment went well
    # this also means that Tkabsch and T are equivalent transformations
    print('RMSD =', RMSD(Tkabsch(X), Y))
    print('Y =', Y)
    print('Tkabsch(X) =', Tkabsch(X))

    mol = plams.Molecule(r"D:\Users\Yuman\Desktop\PhD\ychem\calculations2\1d0673a30f1785b890fc1007d2248a01755aa7b12f69edb22fd52f836fda8dcd.003\radical\input_mol.xyz")
    T = Transform()
    T.translate(x=10)
    print(mol)
    mol.apply_transform(T)
    print(mol)




