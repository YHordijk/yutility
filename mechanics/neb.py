from yutility import numdiff, peaks
from yutility.mechanics import optimize
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
    
    def __init__(self, start, end, function, n_nodes=10, h=1e-5, spring_constant=None, minimize_anchors=False, update_strength=1e-2):
        '''
        Start and end are vectors of the initial and final positions
        n_nodes are the number of nodes linking them together
        '''
        self.start = start
        self.end = end
        self.function = function
        self.n_nodes = n_nodes
        self.h = h
        self.spring_constant = spring_constant or 2500/n_nodes
        if minimize_anchors:
            self.start = optimize.minimize(start, function)
            self.end = optimize.minimize(end, function)
        self.update_strength = update_strength
        self.set_nodes()

    def set_nodes(self):
        xs = []
        for axis in range(len(self.start)):
            xs.append(np.linspace(self.start[axis], self.end[axis], self.n_nodes))

        self.node_pos = np.vstack(xs).T
        self.ideal_spring_distance = np.linalg.norm(np.array(self.start) - np.array(self.end))/self.n_nodes

    def get_force(self):
        # first get the functional gradients, we will get these numerically
        F_force = []
        self.E = []
        for pos in self.node_pos:
            self.E.append(self.function(pos))
            g = numdiff.gradient(self.function, pos, h=self.h)
            F_force.append(g)
        F_force = np.vstack(F_force)
        F_force = F_force/np.linalg.norm(F_force, axis=1).reshape(-1, 1)

        path_tangents = []
        for i, pos in enumerate(self.node_pos):
            if i == self.n_nodes - 1:
                path_tangents.append(pos - self.node_pos[i-1])
                continue
            if i == 0:
                path_tangents.append(self.node_pos[i+1] - pos)
                continue

            path_tangents.append((self.node_pos[i+1] - self.node_pos[i-1]))

        path_tangents = np.vstack(path_tangents)
        path_tangents = path_tangents/np.linalg.norm(path_tangents, axis=1).reshape(-1, 1)

        F_force_tang = []
        for force, tangent in zip(F_force, path_tangents):
            F_force_tang.append(force - tangent * (tangent @ force))
        F_force = np.vstack(F_force_tang)

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

            x_ab = abs(r_ab - self.ideal_spring_distance)/self.ideal_spring_distance
            x_bc = abs(r_bc - self.ideal_spring_distance)/self.ideal_spring_distance
            E_force.append(self.spring_constant * (x_ab * n_ab + x_bc * n_bc))

        E_force.append(np.array([0] * len(self.start)))
        E_force = np.vstack(E_force)

        # the total force is simply the sum of the functional and spring forces
        self.F_force = F_force
        self.F_force[np.argmax(self.E)] *= -1
        self.E_force = E_force
        self.force = E_force - F_force
        self.path_tangents = path_tangents

        self.force = self.force/np.linalg.norm(self.force, axis=1).reshape(-1, 1)

        # we can choose to set the start and end positions as anchors
        # this entails setting their gradients to zero so that they do not move
        self.force[0] = np.array([0] * len(self.start))
        self.force[-1] = np.array([0] * len(self.start))

        return self.force

    def step(self):
        grads = self.get_force()
        self.node_pos = self.node_pos + grads * self.update_strength
        # for pos, tang in zip(self.node_pos, self.path_tangents):
        #     plt.plot([pos[0], (pos+tang/3)[0]], [pos[1], (pos+tang/3)[1]])

    def optimize(self, max_steps=100):
        for i in range(max_steps):
            self.step()
            print(self.node_pos)


if __name__ == '__main__':
    import matplotlib.pyplot as plt
    cmap = plt.get_cmap('coolwarm')
    # S = lambda P: np.sin(P[0] - P[1]) * P[0]**2 * P[1]**2 + P[1] - P[0]**2
    S = lambda P: P[0]**4 + P[1]**4 - (P[0]**2 + P[1]**2)

    X, Y = np.meshgrid(np.linspace(-1, 1, 100), np.linspace(-1, 1, 100))
    plt.imshow(S((X, Y)), extent=(X.min(), X.max(), Y.min(), Y.max()), origin='lower', aspect='auto')

    neb = NEB(start=(-.5, -.75), end=(.75, -.05), function=S, n_nodes=10, spring_constant=.5, minimize_anchors=True)
    plt.scatter(*neb.start, c='r')
    plt.scatter(*neb.end, c='r')
    # plt.scatter([p[0] for p in neb.node_pos], [p[1] for p in neb.node_pos])
    max_steps = 100
    Es = []
    for i in range(1, max_steps+1):
        neb.step()
        Es.append(neb.E)
        plt.scatter([p[0] for p in neb.node_pos], [p[1] for p in neb.node_pos], color=[cmap(i/max_steps) for _ in range(neb.n_nodes)], alpha=i/max_steps*.5 + .5)
        # plt.plot([p[0] for p in neb.node_pos], [p[1] for p in neb.node_pos], c=(i/max_steps, i/max_steps, i/max_steps), alpha=.1)
    # plt.scatter([p[0] for p in neb.node_pos], [p[1] for p in neb.node_pos])
    plt.xlim(X.min(), X.max())
    plt.ylim(Y.min(), Y.max())
        
    plt.figure()
    for i in range(1, max_steps+1):
        peaks.find_peaks(Es[i-1])
        plt.plot(range(neb.n_nodes), Es[i-1], c=cmap(i/max_steps), alpha=i/max_steps*.5 + .5)
        plt.scatter(np.argmax(Es[i-1]), np.max(Es[i-1]), c='k')
        # plt.plot(range(neb.n_nodes-1), np.diff(Es[i-1]), c=cmap(i/max_steps), alpha=1)

    plt.show()


