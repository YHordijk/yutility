import numpy as np


def gradient(function, point, h=1e-5):
	'''
	Calculate the gradient of a function $f$ at a certain point $p$
	
		$ \nabla_i f|_{x_i=p_i} = \frac{f(p+h_i) - f(p)}{h} $

	The gradient is the vector which points in the direction of greatest increase of the function.
	We calculate it by using the limit-definition of the gradient. The accuracy is controlled by
	the parameter $h$. By default it is 1e-5, however for functions with high kinetic energy it 
	may be set to a smaller number. Set it too small and you will run into precision errors, however.
	'''
	ndims = len(point)
	F = function(point)
	grads = []
	for dim in range(ndims):
		z = np.zeros(ndims)
		z[dim] = h

		F_disp = function(point + z)
		grads.append((F_disp - F)/h)

	return np.hstack(grads)


def hessian(function, point, h=1e-5):
	'''
	Calculate the Hessian matrix H of a function $f$ at a certain point $p$

		$H_{ij}|_{x_i=p_i, x_j=p_j} = \frac{f(p - h_i - h_j) - f(p + h_i - h_j) - f(p - h_i + h_j) + f(p + h_i + h_j)}{4h^2}$

	The Hessian is the matrix built by taking the second derivatives of the function.
	We calculate it similarly to the gradient, but take into account the second order displacements.
	'''
	ndims = len(point)
	H = np.zeros((ndims, ndims))
	for dim1 in range(ndims):
		for dim2 in range(ndims):
			z1 = np.zeros(ndims)
			z1[dim1] = h
			z2 = np.zeros(ndims)
			z2[dim2] = h

			if dim1 == dim2:
				H[dim1, dim2] = (function(point + z1) - 2*function(point) + function(point - z1))/h**2
			else:
				H[dim1, dim2] = (function(point - z1 - z2) - function(point + z1 - z2) - function(point - z1 + z2) + function(point + z1 + z2))/(4*h**2)
	return H
