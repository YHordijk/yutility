import os
import json
import periodictable
from pprint import pprint
import numpy as np
from yutility import geometry, timer
import matplotlib.pyplot as plt
from yviewer import viewer
from scm import plams

j = os.path.join



class BasisSet:
	def __init__(self, bs_name):
		self.bs_name = bs_name
		self.load()

	def load(self):
		file = os.path.abspath(j(__file__, '..', 'basissets', self.bs_name + '.json'))
		with open(file) as inp:
			content = inp.read()
			self.data = json.loads(content)

	def get_ao(self, element, principal, angular, magnetic):
		ao_data = self.data['elements']
		if isinstance(element, int):
			atnum = element
		else:
			atnum = getattr(periodictable.elements, element).number
		ao_data = ao_data[str(atnum)]['electron_shells']
		ao_data = ao_data[angular]
		primitive_type = ao_data['function_type']
		exponents = [float(x) for x in ao_data['exponents']]
		coefficients = [float(x) for x in ao_data['coefficients'][principal]]
		return AtomicOrbital(principal, angular, magnetic, coefficients, exponents, primitive_type)


class AtomicOrbital:
	def __init__(self, principal, angular, magnetic, coefficients, exponents, primitive_type):
		self.principal = principal
		self.angular = angular
		self.magnetic = magnetic
		self.coefficients = coefficients
		self.exponents = exponents
		self.primitive_type = primitive_type
		self.center = np.array([0, 0, 0])

		if angular == 0:
			self.basis_vectors = np.array([[0, 0, 0],
										   [0, 0, 0],
										   [0, 0, 0]], dtype='float')
		if angular == 1:
			self.basis_vectors = np.array([[1, 0, 0],
										   [0, 1, 0],
										   [0, 0, 1]], dtype='float')

		self.rotation_matrix = geometry.get_rotation_matrix()

	def set_center(self, center):
		self.center = center

	# def get_radial(self, p):
	# 	r = np.linalg.norm(p-self.center, axis=1)
	# 	R = np.zeros_like(r)
	# 	for exponent in self.exponents:
	# 		R = np.exp(-r * exponent)
	# 	# wf = r * R * self.coeffs[0]
	# 	return R

	def evaluate(self, p):
		r = np.linalg.norm(p-self.center, axis=1).reshape(-1, 1)

		# wf = np.zeros((p.shape[0], 3), dtype=float)
		# r = np.linalg.norm(p-self.center, axis=1)
		# R = r * np.exp(-r*3)
		basis_vector = (self.basis_vectors @ self.rotation_matrix)[:, self.magnetic]
		R = r**(self.principal - 1)
		R = R.reshape(-1, 1)
		Y = np.zeros((p.shape[0], 1), dtype=float)
		if self.primitive_type == 'gto':
			for exp, coeff in zip(self.exponents, self.coefficients):
				Y += coeff * np.exp(-exp * r**2)
		if self.angular == 0:
			wf = np.sum(Y * R , axis=1)
		if self.angular == 1:
			wf = np.sum(Y * basis_vector * (p - self.center) * R , axis=1)
		return wf / np.sqrt(sum(wf**2))




if __name__ == '__main__':
	bs = BasisSet('cc-pvdz.1')

	H2s = bs.get_ao('C', 1, 1, 0)
	# R = geometry.get_rotation_matrix(y=rot)
	# H2s.rotation_matrix = R
	# X, Y, Z = np.meshgrid(np.linspace(-3, 3, 100), np.linspace(-3, 3, 100), np.array([0]))
	# p = np.vstack([X.flatten(), Y.flatten(), Z.flatten()]).T
	# wf = H2s.evaluate(p)
	# plt.figure()
	# plt.imshow(wf.reshape(100, 100))
	# plt.show()

	X, Y, Z = np.meshgrid(np.linspace(-6, 6, 20), np.linspace(-6, 6, 20), np.linspace(-6, 6, 20))
	p = np.vstack([X.flatten(), Y.flatten(), Z.flatten()]).T
	wf = H2s.evaluate(p)

	print(sum(wf**2))

	wf_abs = abs(wf)
	idx = np.where(wf_abs > (wf_abs.max()/4))

	RED = np.array((255, 0, 0))
	BLUE = np.array((0, 0, 255))
	mol = plams.Molecule()
	mol.add_atom(plams.Atom(symbol='C'))
	viewer.show([mol], molinfo=[{'cub': [p[idx], np.where(wf[idx]>=0, 0, 1).reshape(-1, 1) * RED + np.where(wf[idx]<0, 0, 1).reshape(-1, 1) * BLUE]}])

	timer.print_timings2()
