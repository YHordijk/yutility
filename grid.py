import numpy as np
from yutility import ensure_list
from collections.abc import Iterable


class Grid:
    '''
    Class that defines positions and values on a grid.
    '''
    def __init__(self, spacing: Iterable[float] or float = None, 
                       origin: Iterable[float] or float = None, 
                       ndims: int = None):
        if ndims is None:
            # check if we can get the number of dimensions from spacing or origin
            if spacing:
                ndims = len(ensure_list(spacing))
            elif origin:
                ndims = len(ensure_list(origin))
            # if not, we default to a 3D grid
            else:
                ndims = 3

        self.ndims = ndims

        self.spacing = np.array(ensure_list(spacing or [1] * ndims))
        self.origin = np.array(ensure_list(origin or [0] * ndims))
        if len(self.spacing) != ndims:
            self.spacing = np.array([self.spacing[0]] * ndims)

        self.sub_grids = []

    # arithmetic and logic
    def __add__(self, other):
        return self.__and__(other)

    def __and__(self, other):
        if isinstance(other, Grid):
            self.sub_grids.append((other, '+'))
        return self

    def __sub__(self, other):
        return self.__or__(other)

    def __or__(self, other):
        if isinstance(other, Grid):
            self.sub_grids.append((other, '-'))
        return self

    # necessary functions
    def set_points(self, parent=None):
        self._points = None
        self.incorporate_sub_grid()

    def incorporate_sub_grid(self):
        for sub_grid, sign in self.sub_grids:
            if self._points is None:
                self._points = sub_grid.points(self)
                continue

            if sign == '+':
                unique_indices = np.invert((sub_grid.points(self) == self._points[:,None]).all(2).any(0))
                self._points = np.append(self._points, sub_grid.points(self)[unique_indices], axis=0)

            if sign == '-':
                unique_indices = np.invert((self._points == sub_grid.points(self)[:,None]).all(2).any(0))
                self._points = self._points[unique_indices]

                duplicate_indices = (sub_grid.points(self) == self._points[:,None]).all(2).any(0)
                self._points = np.append(self._points, sub_grid.points(self)[duplicate_indices], axis=0)

        # if self._points is None:
        #     self._points = np.array([])

        # return self._points


    def points(self, parent=None):
        if not hasattr(self, '_points'):
            self.set_points(parent=parent)
        return self._points


class Cube(Grid):
    ### grid additive methods
    def __init__(self, origin: Iterable[float] or float = None, 
                       extent: Iterable[float] or float = None,
                       spacing: Iterable[float] or float = None, 
                       *args, **kwargs):
        '''
        Add points in a cube to the Grid.
        
        Args:
            origin: The origin of the cube to be added.
            extent: The distance the cube goes from the origin. For example, for a 2D box, the extent would be (width, height).
        '''
        super().__init__(spacing, origin, *args, **kwargs)
        self.extent = ensure_list(extent)

    def set_points(self, parent=None):
        if parent:
            origin = self.origin - parent.origin
            spacing = parent.spacing
        else:
            origin = self.origin
            spacing = self.spacing

        assert len(origin) == self.ndims
        assert len(self.extent) == self.ndims

        # first build the axes
        # for a cube, the axes are linearly spaced starting in the origin and ending in origin + extent. 
        # The number of points along each axis is given by self.spacing
        axes = []
        for dim in range(self.ndims):
            axis = np.arange(0, self.extent[dim], spacing[dim])
            axis = axis + origin[dim]
            axes.append(axis.round(10))

        # after getting all axes we create a meshgrid
        meshed_axes = np.meshgrid(*axes)
        meshed_axes = [axis.flatten() for axis in meshed_axes]
        # flattening and stacking gives our points
        self._points = np.vstack(meshed_axes).T
        # self.incorporate_sub_grid()
        self.values = np.zeros(len(self._points))


if __name__ == '__main__':
    import matplotlib.pyplot as plt
    main_grid = Grid(1, ndims=2)
    cube = Cube((0, 0), (10, 10))
    cube2 = Cube((-10, -10), (30, 30))

    main_grid = main_grid + cube - cube2
    print(main_grid.sub_grids)
    plt.scatter(*main_grid.points().T)
    plt.show()
