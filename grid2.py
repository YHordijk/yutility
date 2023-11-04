import numpy as np
from yutility import ensure_list
from collections.abc import Iterable


class Grid:
    '''
    Class that defines positions and values on a grid.
    '''
    def __init__(self, spacing: Iterable[float] or float = None):
        '''
        Defines a grid with given spacing.
        
        Args:
            spacing: the distance between gridpoints. If a single value, the distance will be the same in all directions. A list/tuple value will define the grid along each direction.

        '''
        self.spacing = ensure_list(spacing)
        self.origin = None
        self.points = None
        self.values = None

    @property
    def ndims(self):
        if self.origin:
            return len(self.origin)

    def __sub__(self, other):
        if isinstance(other, Grid):
            self.remove_grid(other)
        return self

    def __add__(self, other):
        if isinstance(other, Grid):
            self.add_grid(other)
        return self

    def add_grid(self, other):
        other.spacing = self.spacing
        other.set_points()
        if self.points is None:
            self.points = other.points
        else:
            self.points = np.vstack([self.points, other.points])

        if self.values is None:
            self.values = other.values
        else:
            self.values = np.hstack([self.values, other.values])

        self.remove_duplicates()

    def remove_grid(self, other):
        other.spacing = self.spacing
        other.set_points()
        if self.points is None:
            return

        idxs_to_keep = np.unique(np.where(1-np.isin(self.points, other.points))[0])
        self.points = self.points[idxs_to_keep]
        self.values = self.values[idxs_to_keep]

    def remove_duplicates(self):
        _, idxs = np.unique(self.points, return_index=True, axis=0)
        self.points = self.points[idxs]
        self.values = self.values[idxs]



class Cube(Grid):
    def __init__(self, origin: Iterable[float] or float = None, 
                       extent: Iterable[float] or float = None,
                       *args, **kwargs):
        '''
        Build a grid of points in a cube.
        
        Args:
            origin: The origin of the cube to be added.
            extent: The distance the cube goes from the origin. For example, for a 2D box, the extent would be (width, height).
        '''
        super().__init__(*args, **kwargs)
        self.origin = ensure_list(origin)
        self.extent = ensure_list(extent)


    def set_points(self):
        if len(self.extent) == 1 and len(self.extent) != self.ndims:
            self.extent = self.extent * self.ndims

        if len(self.spacing) == 1 and len(self.spacing) != self.ndims:
            self.spacing = self.spacing * self.ndims

        # first build the axes
        # for a cube, the axes are linearly spaced starting in the origin and ending in origin + extent. 
        # The number of points along each axis is given by spacing
        axes = []
        for dim in range(self.ndims):
            axis = np.arange(min(0, self.extent[dim]), max(0, self.extent[dim]), self.spacing[dim]) + self.origin[dim]
            axes.append(axis)

        # after getting all axes we create a meshgrid
        meshed_axes = np.meshgrid(*axes)
        meshed_axes = [axis.flatten() for axis in meshed_axes]
        # flattening and stacking gives our points
        self.points = np.vstack(meshed_axes).T
        # self.incorporate_sub_grid()
        self.values = np.zeros(len(self.points))


class Sphere(Grid):
    def __init__(self, origin: Iterable[float] or float = None, 
                       radius: float = None,
                       *args, **kwargs):
        '''
        Build a grid of points in a cube.
        
        Args:
            origin: The origin of the cube to be added.
            radius: The distance from the origin of the sphere to its edge. Can be tuple or single value.
        '''
        super().__init__(*args, **kwargs)
        self.origin = ensure_list(origin)
        self.radius = radius

    def set_points(self):
        if len(self.spacing) == 1 and len(self.spacing) != self.ndims:
            self.spacing = self.spacing * self.ndims

        # first build the axes for a cube
        # for a cube, the axes are linearly spaced starting in the origin and ending in origin + radius. 
        # The number of points along each axis is given by spacing
        axes = []
        for dim in range(self.ndims):
            axis = np.arange(-self.radius, self.radius, self.spacing[dim]) + self.origin[dim]
            axes.append(axis)

        # after getting all axes we create a meshgrid
        meshed_axes = np.meshgrid(*axes)
        meshed_axes = [axis.flatten() for axis in meshed_axes]
        # flattening and stacking gives our points for the cube
        cube_points = np.vstack(meshed_axes).T

        dists = np.linalg.norm(cube_points - self.origin, axis=1)
        self.points = cube_points[dists <= self.radius]
        self.values = np.zeros(len(self.points))



if __name__ == '__main__':
    import matplotlib.pyplot as plt


    grid = Grid(.5)
    grid = grid + Sphere((0, 0), 5) - Sphere((0, 0), 2)
    # grid.add_cube((-10, -10), 20)
    # grid.remove_cube((-5, -5), 10)
    # grid.add_cube((-4, 0), 4)
    # print(grid.points)

    plt.scatter(grid.points[:,0], grid.points[:,1])
    plt.show()

