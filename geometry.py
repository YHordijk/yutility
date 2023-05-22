import numpy as np
import scm.plams as plams
from math import sin, cos
from yutility import molecule as yu_molecule


normal_to_plane = {'xy': np.array([0., 0., 1.]),
                   'yx': np.array([0., 0., 1.]),
                   'xz': np.array([0., 1., 0.]),
                   'zx': np.array([0., 1., 0.]),
                   'yz': np.array([1., 0., 0.]),
                   'zy': np.array([1., 0., 0.])}

axis_directions = {'x': np.array([1.0, 0.0, 0.0]),
                   'y': np.array([0.0, 1.0, 0.0]),
                   'z': np.array([0.0, 0.0, 1.0])}


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
