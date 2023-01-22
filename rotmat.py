import numpy as np
from math import sin, cos

def get_rotmat(x=0, y=0, z=0):
	Rx = np.array([[1, 0, 0],
				   [0, cos(x), -sin(x)],
				   [0, sin(x), cos(x)]])
	Ry = np.array([[cos(y), 0, sin(y)],
				   [0, 1, 0],
				   [-sin(y), 0, cos(y)]])
	Rz = np.array([[cos(z), -sin(z), 0],
				   [sin(z), cos(z), 0],
				   [0, 0, 1]])
	return Rx @ Ry @ Rz

def unit():
	return get_rotmat(0, 0, 0)

