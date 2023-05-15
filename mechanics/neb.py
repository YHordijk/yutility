from yutility import numdiff


class NEB:
	def __init__(self, start, end, model, n_nodes=10,):
		'''
		Start and end are vectors of the initial and final positions
		n_nodes are the number of nodes linking them together
		'''
		self.start = start
		self.end = end
		self.model = model
		self.n_nodes = n_nodes
		self.set_nodes()
		self.ideal_spring_distance = np.linalg.norm(np.array(start) - np.array(end))/n_nodes

	def set_nodes(self):
		xs = []
		for axis in range(len(self.start)):
			xs.append(np.linspace(self.start[axis], self.end[axis], self.n_nodes))

		self.node_pos = np.vstack(xs).T
		self.node_pos = torch.tensor(self.node_pos, requires_grad=True)
		self.node_pos.retain_grad()

	def get_grad(self, h=1e-5, k=15, set_anchors=False):
		# first get the structural gradients
		self.F = self.model(self.node_pos)
		grads = []
		for axis in range(len(self.start)):
			z = torch.zeros(self.node_pos.shape)
			z[:, axis] = h

			F_disp = self.model(self.node_pos + z)
			grads.append(-(F_disp - self.F)/h)

		grads = torch.hstack(grads)

		# get elastic gradients
		E_grads = []
		E_grads.append(torch.tensor([0] * len(self.start)))
		for i, b in enumerate(self.node_pos[1:-1], start=1):
			a, c = self.node_pos[i-1], self.node_pos[i+1]
			r_ab = torch.linalg.norm(a - b)
			r_bc = torch.linalg.norm(b - c)

			n_ab = (a - b)/r_ab
			n_bc = (c - b)/r_bc

			x_ab = abs(r_ab - self.ideal_spring_distance)
			x_bc = abs(r_bc - self.ideal_spring_distance)
			E_grads.append(k * (x_ab * n_ab + x_bc * n_bc))

		E_grads.append(torch.tensor([0] * len(self.start)))
		E_grads = torch.vstack(E_grads)

		grads = grads + E_grads
		if set_anchors:
			grads[0] = torch.tensor([0] * len(self.start))
			grads[-1] = torch.tensor([0] * len(self.start))
		return grads

	def step(self, h=1e-5, update_strength=1e-2, set_anchors=False):
		grads = self.get_grad(h, set_anchors=set_anchors)
		self.node_pos = self.node_pos + grads * update_strength
