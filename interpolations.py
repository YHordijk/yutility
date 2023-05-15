'''
Module used for interpolation of surfaces
'''

import numpy as np


class Interpolator:
	def __init__(self, array, extent=None, interpolation_method=None):
		'''
		Array is a numpy array that contains the values at the array positions.
		The coordinates of the array is given by extent as a container of tuples
		defining the min and max values of the axes.
		'''
		self.array = array
		self.ndim = array.ndim
		self.extent = extent



