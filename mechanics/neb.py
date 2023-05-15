from yutility import numdiff
import numpy as np


class NEB:
	'''
	Nudged Elastic Band (NEB) calculator class. 
	NEB gives us the lowest energy path on a surface defined by a certain function. It does so
	by defining a chain of nodes connected by springs on this surface. The nodes are then iteratively
	subject to two forces. The first is the force experienced by the gradient of the function at the
	nodal position. The second force it the spring force due to the positions of the previous and next nodes.
	After some iterations the nodes experience an equal spring and functional force and the band will have converged.
	This method ensures that nodes are roughly equally spaces along the path and the path is therefore sampled nicely.
	It also ensure the energy is minimized along the path, which is usefull in chemical applications
	'''
	
	def __init__(self, start, end, function, n_nodes=10, h=1e-5, spring_constant=15, set_anchors=False, update_strength=1e-2):
		'''
		Start and end are vectors of the initial and final positions
		n_nodes are the number of nodes linking them together
		'''
		self.start = start
		self.end = end
		self.function = function
		self.n_nodes = n_nodes
		self.h = h
		self.spring_constant = spring_constant
		self.set_anchors = set_anchors
		self.update_strength = update_strength
		self.set_nodes()

	def set_nodes(self):
		xs = []
		for axis in range(len(self.start)):
			xs.append(np.linspace(self.start[axis], self.end[axis], self.n_nodes))

		self.node_pos = np.vstack(xs).T
		self.ideal_spring_distance = np.linalg.norm(np.array(start) - np.array(end))/n_nodes

	def get_force(self):
		# first get the functional gradients, we will get these numerically
		F_force = []
		for pos in self.node_pos:
			g = -numdiff.gradient(self.function, pos, h=self.h)
			F_force.append(g)
		F_force = np.vstack(F_force)

		# get spring forces, we will calculate these exactly
		# using hooke's law
		E_force = []
		E_force.append(np.array([0] * len(self.start)))
		for i, b in enumerate(self.node_pos[1:-1], start=1):
			a, c = self.node_pos[i-1], self.node_pos[i+1]
			r_ab = np.linalg.norm(a - b)
			r_bc = np.linalg.norm(b - c)

			n_ab = (a - b)/r_ab
			n_bc = (c - b)/r_bc

			x_ab = abs(r_ab - self.ideal_spring_distance)
			x_bc = abs(r_bc - self.ideal_spring_distance)
			E_force.append(self.spring_constant * (x_ab * n_ab + x_bc * n_bc))

		E_force.append(np.array([0] * len(self.start)))
		E_force = np.vstack(E_force)

		# the total force is simply the sum of the functional and spring forces
		self.F_force = F_force
		self.E_force = E_force
		force = F_force + E_force

		# we can choose to set the start and end positions as anchors
		# this entails setting their gradients to zero so that they do not move
		if self.set_anchors:
			force[0] = np.array([0] * len(self.start))
			force[-1] = np.array([0] * len(self.start))

		return force

	def step(self):
		grads = self.get_force()
		self.node_pos = self.node_pos + grads * self.update_strength
